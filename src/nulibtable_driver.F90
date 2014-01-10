!-*-f90-*-
program driver

  use nulib  
  use nulibtable
  implicit none

  real*8, allocatable :: eas(:)
  real*8, allocatable :: eas_energy(:,:)
  real*8, allocatable :: eas_species_energy(:,:,:)
  real*8 :: scatop,absop,time
  real*8 :: absop_array(2500),scatop_array(2500),emiss_array(2500)

  integer i,j,k,ns,ng,line
  real*8 rho,temp,ye
  logical :: cont,files
  character*1024, dimension(3,3) :: filenames
  character*1024, dimension(3) :: corrections
  real*8, dimension(3,3000,3,12) :: eas_all
  real*8, dimension(3,3000,3,12) :: eas_allbutec
  real*8, dimension(3,3000,3,12) :: eas_ec
  real*8, dimension(3,3000) :: mean_opacity_ec
  real*8, dimension(3,3000) :: rosseland_mean_opacity_ec
  real*8, dimension(3,3000,4) :: state_vars
  real*8, allocatable,dimension(:) :: eos_variables
  real*8, dimension(3,3000) :: total_blackbody
  real*8, dimension(3,3000,3,12) :: blackbody_spectra
  !NuLib parameters file (weak rates and EOS)
  character*200 :: parameters_filename = "/projects/ceclub/gr1dnulib/GitHub/NuLib/parameters"
  !nulib setup
  integer :: keytemp,keyerr
  real*8 :: dxfac,mindx
  real*8 :: matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe
  real*8 :: precision = 1.0d-10
  real*8  :: min_logrho,max_logrho
  real*8  :: min_logtemp,max_logtemp
  real*8  :: min_ye,max_ye


  files=.true.
  filenames(1,1) = "/projects/ceclub/gr1dnulib/GitHub/NuLib/tables/NuLib_Hempel_rho82_temp65_ye51_ng12_ns3_Itemp65_Ieta61_version1.0_20131113_-10.h5"
  filenames(1,2) = "/mnt/simulations/ceclub/sullivan/gr1d/dyedtruns/Data_decfactor10/rho_temp_ye_c.dat"
  filenames(1,3) = "/mnt/simulations/ceclub/sullivan/gr1d/dyedtruns/Data_decfactor10/meanopacities_rho_c.dat"
  filenames(2,1) = "/projects/ceclub/gr1dnulib/GitHub/NuLib/tables/NuLib_Hempel_rho82_temp65_ye51_ng12_ns3_Itemp65_Ieta61_version1.0_20131111_lessng.h5"
  filenames(2,2) = "/mnt/simulations/ceclub/sullivan/gr1d/dyedtruns/Data_all_vel_inel_12_tvd/rho_temp_ye_c.dat"
  filenames(2,3) = "/mnt/simulations/ceclub/sullivan/gr1d/dyedtruns/Data_all_vel_inel_12_tvd/meanopacities_rho_c.dat"
  filenames(3,1) = "/projects/ceclub/gr1dnulib/GitHub/NuLib/tables/NuLib_Hempel_rho82_temp65_ye51_ng12_ns3_Itemp65_Ieta61_version1.0_20131113_m+10.h5"
  filenames(3,2) = "/mnt/simulations/ceclub/sullivan/gr1d/dyedtruns/Data_incfactor10/rho_temp_ye_c.dat"
  filenames(3,3) = "/mnt/simulations/ceclub/sullivan/gr1d/dyedtruns/Data_incfactor10/meanopacities_rho_c.dat"
  corrections(1) = "/projects/ceclub/gr1dnulib/GitHub/NuLib/NuLib_Hempel_rho82_temp65_ye51_ng12_ns3_Itemp65_Ieta61_version1.0_20131118-m-10_noec.h5"
  corrections(2) = "/projects/ceclub/gr1dnulib/GitHub/NuLib/NuLib_Hempel_rho82_temp65_ye51_ng12_ns3_Itemp65_Ieta61_version1.0_20131118-m0_noec.h5"
  corrections(3) = "/projects/ceclub/gr1dnulib/GitHub/NuLib/NuLib_Hempel_rho82_temp65_ye51_ng12_ns3_Itemp65_Ieta61_version1.0_20131118-m+10_noec.h5"
  state_vars = -9999.9d0
  blackbody_spectra = 0.0d0
  total_blackbody = 0.0d0
  rosseland_mean_opacity_ec = 0.0d0
  eas_all = 0.0d0
  eas_allbutec = 0.0d0
  eas_ec = 0.0d0

  call input_parser(parameters_filename)
  call initialize_nulib(1,6,12)
  call readtable(eos_filename)
  m_ref = m_amu !for SFHo_EOS (Hempel)

  allocate(eos_variables(total_eos_variables))
  eos_variables = 0.0d0
  min_ye = 0.035d0
  max_ye = 0.55d0
  min_logrho = 6.0d0
  max_logrho = 15.5d0
  min_logtemp = log10(0.05d0)
  max_logtemp = log10(150.0d0)


  !set up energies bins
  do_integrated_BB_and_emissivity = .false.
  mindx = 2.0d0
  bin_bottom(1) = 0.0d0 !MeV
  bin_bottom(2) = 2.0d0 !MeV
  bin_bottom(3) = bin_bottom(2)+mindx
  bin_bottom(number_groups) = 250.0d0

    call nulib_series2(number_groups-1,bin_bottom(2),bin_bottom(number_groups),mindx,dxfac)
  do i=4,number_groups
     bin_bottom(i) = bin_bottom(i-1)+(bin_bottom(i-1)-bin_bottom(i-2))*dxfac
  enddo
  
  !calculate bin widths & energies from the bottom of the bin & energy at top on bin
  do i=1,number_groups-1
     energies(i) = (bin_bottom(i)+bin_bottom(i+1))/2.0d0
     bin_widths(i) = bin_bottom(i+1)-bin_bottom(i)
     bin_top(i) = bin_bottom(i+1)
  enddo
  energies(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups-1)*dxfac/2.0d0
  bin_widths(number_groups) = 2.0*(energies(number_groups)-bin_bottom(number_groups))
  bin_top(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups)

  do i=1,3
     
     call nulibtable_reader(filenames(i,1))
     
     allocate(eas(nulibtable_number_easvariables))
     allocate(eas_energy(nulibtable_number_groups,nulibtable_number_easvariables))
     allocate(eas_species_energy(nulibtable_number_species, &
          nulibtable_number_groups,nulibtable_number_easvariables))
     
     ns = 1
     line = 0
     
     open(11,file=filenames(i,2),status='old')
     do while(.true.)
        read(11,*,end=10) time,rho,temp,ye
        call nulibtable_single_species_range_energy(rho,temp,ye,ns,eas_energy, &
             nulibtable_number_groups,nulibtable_number_easvariables)     
        line = line + 1
        
        state_vars(i,line,1)=time
        state_vars(i,line,2)=rho
        state_vars(i,line,3)=temp
        state_vars(i,line,4)=ye
        eas_all(i,line,1,:)=eas_energy(:,1)
        eas_all(i,line,2,:)=eas_energy(:,2)
        eas_all(i,line,3,:)=eas_energy(:,3)
        write(*,*) line        
     end do

10   close(11)


!      deallocate(eas)
!      deallocate(eas_energy)
!      deallocate(eas_species_energy)
!      deallocate(nulibtable_energies)
!      deallocate(nulibtable_inv_energies)
!      deallocate(nulibtable_ewidths)
!      deallocate(nulibtable_ebottom)
!      deallocate(nulibtable_etop)     
!      deallocate(nulibtable_logrho)
!      deallocate(nulibtable_logtemp)
!      deallocate(nulibtable_ye)          
!      deallocate(nulibtable_emissivities)
!      deallocate(nulibtable_absopacity)
!      deallocate(nulibtable_scatopacity)
!      line = 0


!      call nulibtable_reader(corrections(i))

!      allocate(eas(nulibtable_number_easvariables))
!      allocate(eas_energy(nulibtable_number_groups,nulibtable_number_easvariables))
!      allocate(eas_species_energy(nulibtable_number_species, &
!           nulibtable_number_groups,nulibtable_number_easvariables))
!      eas_energy = 0.0d0
!      eas = 0.0d0

!      open(11,file=filenames(i,2),status='old')
!      do while(.true.)
!         read(11,*,end=15) time,rho,temp,ye
!         call nulibtable_single_species_range_energy(rho,temp,ye,ns,eas_energy, &
!              nulibtable_number_groups,nulibtable_number_easvariables)     
!         line = line + 1
!         eas_allbutec(i,line,1,:)=eas_energy(:,1)
!         eas_allbutec(i,line,2,:)=eas_energy(:,2)
!         eas_allbutec(i,line,3,:)=eas_energy(:,3)
!         eas_ec(i,line,2,:)=eas_all(i,line,2,:) - eas_energy(:,2)
        
!         !write(*,*) eas_all(i,line,2,12),eas_ec(i,line,2,12),eas_energy(12,2)
!         write(*,*) line
!      end do
! 15   close(11)

     deallocate(eas)
     deallocate(eas_energy)
     deallocate(eas_species_energy)
     deallocate(nulibtable_energies)
     deallocate(nulibtable_inv_energies)
     deallocate(nulibtable_ewidths)
     deallocate(nulibtable_ebottom)
     deallocate(nulibtable_etop)     
     deallocate(nulibtable_logrho)
     deallocate(nulibtable_logtemp)
     deallocate(nulibtable_ye)          
     deallocate(nulibtable_emissivities)
     deallocate(nulibtable_absopacity)
     deallocate(nulibtable_scatopacity)
     absop_array = 0.0d0
     scatop_array = 0.0d0
     line = 0


  end do

  eas_ec = eas_all

  do i=1,3
     open(22,file=filenames(i,3))
     do j=1,3000
        if(state_vars(i,j,1).ne.-9999.9d0)then           
           eos_variables = 0.0d0
           eos_variables(rhoindex) = state_vars(i,j,2)
           eos_variables(tempindex) = state_vars(i,j,3)
           eos_variables(yeindex) = state_vars(i,j,4)

           !! EOS stuff
           keytemp = 1
           keyerr = 0

           call nuc_eos_full(eos_variables(rhoindex),eos_variables(tempindex), &
                eos_variables(yeindex),eos_variables(energyindex),matter_prs, &
                eos_variables(entropyindex),matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe, &
                eos_variables(xaindex),eos_variables(xhindex),eos_variables(xnindex), &
                eos_variables(xpindex),eos_variables(abarindex),eos_variables(zbarindex), &
                eos_variables(mueindex),eos_variables(munindex),eos_variables(mupindex), &
                eos_variables(muhatindex),keytemp,keyerr,precision)
           if (keyerr.ne.0) then
              write(*,*) "rho: ", eos_variables(rhoindex)
              write(*,*) "temperature: ", eos_variables(tempindex)
              write(*,*) "ye: ", eos_variables(yeindex)
              write(*,*) "eos error", keyerr
              stop "set_eos_variables: us eos error"
           endif
           if(eos_variables(xhindex).lt.1.0d-15) then
              eos_variables(xhindex) = 0.0d0
           endif

           call return_blackbody_spectra(blackbody_spectra(i,j,:,:),eos_variables)           
           total_blackbody(i,j) = Sum(blackbody_spectra(i,j,1,:)*bin_widths(:))
!           mean_opacity_ec(i,j) = Sum(blackbody_spectra(i,j,1,:)*bin_widths(:)/eas_ec(i,j,2,:))
           do k=1,number_groups
              rosseland_mean_opacity_ec(i,j) = rosseland_mean_opacity_ec(i,j) + blackbody_spectra(i,j,1,k)*bin_widths(k)/eas_ec(i,j,2,k)
           end do
           rosseland_mean_opacity_ec(i,j) = rosseland_mean_opacity_ec(i,j)/total_blackbody(i,j)
           rosseland_mean_opacity_ec(i,j) = 1/rosseland_mean_opacity_ec(i,j)
           write(*,*) rosseland_mean_opacity_ec(i,j)
           write(22,*) state_vars(i,j,2),rosseland_mean_opacity_ec(i,j)
        end if
     end do
     close(22)
  end do
     

end program driver
