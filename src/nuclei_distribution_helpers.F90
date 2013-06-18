!-*-f90-*- 

!If you are using Hempel EOS, these routines use his, go download them
!from http://phys-merger.physik.unibas.ch/~hempel/eos.html.
!Currently, this is setup to use the sfho_frdm_comp.zip composition
!data, modify appropriately

subroutine set_up_Hempel

#if NUCLEI_HEMPEL
  use sfho_frdm_composition_module
  use nulib, only : hempel_lookup_table
  implicit none

  !local
  integer :: highestA, highestZ
  integer :: i

  call compdata_readin

  !highest A
  highestA = maxval(az(1:kmax,1))
  !highest Z
  highestZ = maxval(az(1:kmax,2))
  
  allocate(hempel_lookup_table(highestA,highestZ))

  hempel_lookup_table = 0

  do i=1,kmax
     hempel_lookup_table(az(i,1),az(i,2)) = i
  enddo

#endif

end subroutine set_up_Hempel

subroutine nuclei_distribution_Hempel(number_nuclei,nuclei_A,nuclei_Z,mass_fractions,number_densities,eos_variables)

#if NUCLEI_HEMPEL
   use sfho_frdm_composition_module
   use nulib
   implicit none

   !inputs and outputs
   integer, intent(in) :: number_nuclei
   integer, intent(in), dimension(number_nuclei) :: nuclei_A
   integer, intent(in), dimension(number_nuclei) :: nuclei_Z
   real*8, intent(out), dimension(number_nuclei) :: mass_fractions
   real*8, intent(out), dimension(number_nuclei) :: number_densities !number/fm^3
   real*8, intent(in), dimension(total_eos_variables) :: eos_variables
   
   !local
   real*8 :: t,ye,nb
   real*8, dimension(kmax) :: naz, xaz
   real*8 :: xn,xp,nn,np
   integer :: sflag
   integer :: i,inuc
   
   t=eos_variables(tempindex)
   ye=eos_variables(yeindex)
   nb=eos_variables(rhoindex)*1.0d-39/(m_ref*mev_to_gram)

   sflag = 0

   call sub_dist_interpol(t,ye,nb,xaz,xn,xp,naz,nn,np,sflag)

   if (sflag.eq.1) then
      !worked as expected
   else
      stop "nuclei_distribution_Hempel: interpolation failed"
   endif

   do i=1,number_nuclei
      inuc = hempel_lookup_table(nuclei_A(i),nuclei_Z(i))
      number_densities(i) = naz(inuc) !number/fm^3
      mass_fractions(i) = xaz(inuc)
   enddo

#endif
      
end subroutine nuclei_distribution_Hempel

