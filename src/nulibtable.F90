!-*-f90-*-
module nulibtable

  implicit none

  integer :: nulibtable_number_species
  integer :: nulibtable_number_groups
  
  real*8, allocatable,save :: nulibtable_logrho(:)
  real*8, allocatable,save :: nulibtable_logtemp(:)
  real*8, allocatable,save :: nulibtable_ye(:)
  
  real*8, allocatable,save :: nulibtable_energies(:)
  real*8, allocatable,save :: nulibtable_inv_energies(:)
  real*8, allocatable,save :: nulibtable_ewidths(:)
  real*8, allocatable,save :: nulibtable_ebottom(:)
  real*8, allocatable,save :: nulibtable_etop(:)

  real*8, allocatable,save :: nulibtable_emissivities(:,:,:,:)
  real*8, allocatable,save :: nulibtable_absopacity(:,:,:,:)
  real*8, allocatable,save :: nulibtable_scatopacity(:,:,:,:)

  integer :: nulibtable_nrho
  integer :: nulibtable_ntemp
  integer :: nulibtable_nye

  real*8 :: nulibtable_logrho_min
  real*8 :: nulibtable_logrho_max

  real*8 :: nulibtable_logtemp_min
  real*8 :: nulibtable_logtemp_max

  real*8 :: nulibtable_ye_min
  real*8 :: nulibtable_ye_max

  integer :: nulibtable_number_easvariables

end module nulibtable

!this takes rho,temp,ye,species and energy and return eas
subroutine nulibtable_single_species_single_energy(xrho,xtemp,xye,lns,lng,eas,eas_n1)
  
  use nulibtable
  implicit none

  real*8, intent(in) :: xrho, xtemp, xye !inputs
  real*8 :: xlrho, xltemp !log versions
  integer, intent(in) :: lns, lng
  integer, intent(in) :: eas_n1
  real*8, intent(out) :: eas(eas_n1)
  integer :: startindex,endindex

  if (eas_n1.ne.nulibtable_number_easvariables) stop "supplied array dimensions (1) is not commensurate with table"
  
  xlrho = log10(xrho)
  xltemp = log10(xtemp)

  if (xlrho.lt.nulibtable_logrho_min) stop "density below nulib table minimum rho"
  if (xlrho.gt.nulibtable_logrho_max) stop "density above nulib table maximum rho"
  if (xltemp.lt.nulibtable_logtemp_min) stop "temperature below nulib table minimum temp"
  if (xltemp.gt.nulibtable_logtemp_max) stop "temperature above nulib table maximum temp"
  if (xye.lt.nulibtable_ye_min) stop "ye below nulib table minimum ye"
  if (xye.gt.nulibtable_ye_max) stop "ye above nulib table maximum ye"

  startindex = (lns-1)*nulibtable_number_groups+(lng-1)+1
  endindex = startindex

  call intp3d_many_mod(xlrho,xltemp,xye,eas(1), &
       nulibtable_emissivities(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  call intp3d_many_mod(xlrho,xltemp,xye,eas(2), &
       nulibtable_absopacity(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  call intp3d_many_mod(xlrho,xltemp,xye,eas(3), &
       nulibtable_scatopacity(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  eas(:) = 10.0d0**eas(:)

end subroutine nulibtable_single_species_single_energy

!this takes rho,temp,ye,species and return eas over energy range
subroutine nulibtable_single_species_range_energy(xrho,xtemp,xye,lns,eas,eas_n1,eas_n2)
  
  use nulibtable
  implicit none

  real*8, intent(in) :: xrho, xtemp, xye !inputs
  real*8 :: xlrho, xltemp !log versions
  integer, intent(in) :: lns
  integer, intent(in) :: eas_n1,eas_n2
  real*8, intent(out) :: eas(eas_n1,eas_n2)
  integer :: ing
  real*8 :: xeas(eas_n1)
  integer :: startindex,endindex

  if(size(eas,1).ne.nulibtable_number_groups) then
     stop "nulibtable_single_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif
  if(size(eas,2).ne.nulibtable_number_easvariables) then
     stop "nulibtable_single_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif
  xlrho = log10(xrho)
  xltemp = log10(xtemp)

  if (xlrho.lt.nulibtable_logrho_min) stop "density below nulib table minimum rho"
  if (xlrho.gt.nulibtable_logrho_max) stop "density above nulib table maximum rho"
  if (xltemp.lt.nulibtable_logtemp_min) stop "temperature below nulib table minimum temp"
  if (xltemp.gt.nulibtable_logtemp_max) stop "temperature above nulib table maximum temp"
  if (xye.lt.nulibtable_ye_min) stop "ye below nulib table minimum ye"
  if (xye.gt.nulibtable_ye_max) stop "ye above nulib table maximum ye"

  startindex = (lns-1)*nulibtable_number_groups+1
  endindex = startindex + nulibtable_number_groups - 1

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas, &
       nulibtable_emissivities(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)
  
  eas(:,1) = 10.0d0**xeas(:)

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas, &
       nulibtable_absopacity(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)
  
  eas(:,2) = 10.0d0**xeas(:)

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas, &
       nulibtable_scatopacity(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)
  
  eas(:,3) = 10.0d0**xeas(:)

end subroutine nulibtable_single_species_range_energy

!this takes rho,temp,ye,species and return eas over energy and species range
subroutine nulibtable_range_species_range_energy(xrho,xtemp,xye,eas,eas_n1,eas_n2,eas_n3)
  
  use nulibtable
  implicit none

  real*8, intent(in) :: xrho, xtemp, xye !inputs
  real*8 :: xlrho, xltemp !log versions
  integer, intent(in) :: eas_n1,eas_n2,eas_n3
  real*8, intent(out) :: eas(eas_n1,eas_n2,eas_n3)
  integer :: ins,ing
  real*8 :: xeas(eas_n1*eas_n2)
  integer :: index

  if(size(eas,1).ne.nulibtable_number_species) then
     stop "nulibtable_range_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif
  if(size(eas,2).ne.nulibtable_number_groups) then
     stop "nulibtable_range_species_range_energy: supplied array dimensions (3) is not commensurate with table"
  endif  
  if(size(eas,3).ne.nulibtable_number_easvariables) then
     stop "nulibtable_range_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif

  xlrho = log10(xrho)
  xltemp = log10(xtemp)

  if (xlrho.lt.nulibtable_logrho_min) stop "density below nulib table minimum rho"
  if (xlrho.gt.nulibtable_logrho_max) stop "density above nulib table maximum rho"
  if (xltemp.lt.nulibtable_logtemp_min) stop "temperature below nulib table minimum temp"
  if (xltemp.gt.nulibtable_logtemp_max) stop "temperature above nulib table maximum temp"
  if (xye.lt.nulibtable_ye_min) stop "ye below nulib table minimum ye"
  if (xye.gt.nulibtable_ye_max) stop "ye above nulib table maximum ye"

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas,nulibtable_emissivities,nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1*eas_n2,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  do ins=1,nulibtable_number_species
     do ing=1,nulibtable_number_groups
        index = (ins-1)*nulibtable_number_groups + (ing-1) + 1
        eas(ins,ing,1) = 10.0d0**xeas(index)
     enddo
  enddo

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas,nulibtable_absopacity,nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1*eas_n2,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  do ins=1,nulibtable_number_species
     do ing=1,nulibtable_number_groups
        index = (ins-1)*nulibtable_number_groups + (ing-1) + 1
        eas(ins,ing,2) = 10.0d0**xeas(index)
     enddo
  enddo

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas,nulibtable_scatopacity,nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1*eas_n2,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  do ins=1,nulibtable_number_species
     do ing=1,nulibtable_number_groups
        index = (ins-1)*nulibtable_number_groups + (ing-1) + 1
        eas(ins,ing,3) = 10.0d0**xeas(index)
     enddo
  enddo

end subroutine nulibtable_range_species_range_energy
