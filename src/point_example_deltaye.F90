!-*-f90-*-
program point_example

  use nulib
  use weak_rates
  implicit none

  !many people use different number of species, this is to denote how they are devided up.
  ! mytable_neutrino_scheme = 1 (three output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon+tau neutrino+antineutrino
  ! neutrino_scheme = 2 (four output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon+tau neutrino             #4 mu and tau antineutrino
  ! neutrino_scheme = 3 (six output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon neutrino                 #4 mu antineutrino
  !         #5: tau neutrino                  #6 tau antineutrino
  integer :: mypoint_neutrino_scheme = 1

  !number of species for nulib to calculate interactions for, must
  !be six currently, average is done via above parameter
  integer :: mypoint_number_species = 6

  !number of species for nulib to output interactions for, must
  !be commensurate with neutrino_scheme above
  integer :: mypoint_number_output_species = 3

  !number of energy groups
  integer :: mypoint_number_groups = 24

  !EOS table
  character*200 :: eos_filename = "/user/sullivan/gr1dnulib/GitHub/NuLib/src/extra_code_and_tables/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"

  !Weak rate data (currently LMP rates only)
  character*200 :: weakrates_filename = "/projects/ceclub/gr1dnulib/GitHub/NuLib/src/extra_code_and_tables/rates-ext.out"

  !local variables
  real*8, allocatable,dimension(:,:) :: local_emissivity
  real*8, allocatable,dimension(:,:) :: local_absopacity
  real*8, allocatable,dimension(:,:) :: local_scatopacity
  real*8, allocatable,dimension(:,:) :: local_Phi0, local_Phi1
  real*8, allocatable,dimension(:,:) :: blackbody_spectra
  real*8, allocatable,dimension(:) :: eos_variables
  real*8 :: matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe
  integer :: keytemp,keyerr
  real*8 :: precision = 1.0d-10
  real*8 :: xrho, xtemp, xye, xtime
  integer i,j
  real*8 dxfac,mindx

  real*8, allocatable,dimension(:) :: mass_fractions
  real*8, allocatable,dimension(:) :: number_densities
  real*8, allocatable,dimension(:) :: probability_dist

  logical :: cont
  character*200 :: time_string
  integer :: nindex
  real*8 :: current_time
  real*8 :: current_energy
  integer :: IO


  !allocate the arrays for the point values
  allocate(local_emissivity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_absopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_scatopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(blackbody_spectra(mypoint_number_output_species,mypoint_number_groups))

  !this sets up many cooefficients and creates the energy grid (one
  !zone + log spacing) see nulib.F90:initialize_nulib
  call initialize_nulib(mypoint_neutrino_scheme,mypoint_number_species,mypoint_number_groups)

  !read in EOS table & set reference mass
  call readtable(eos_filename)
  m_ref = m_amu !for SFHo_EOS (Hempel)
  ! m_ref = m_n !for LS220

  !read in weak rates table and build interpolant functions
  call readrates(weakrates_filename,table_bounds)

  allocate(mass_fractions(nspecies))
  allocate(number_densities(nspecies))
  allocate(probability_dist(number_groups))

  !set up energies bins
  do_integrated_BB_and_emissivity = .false.
  mindx = 1.0d0
  bin_bottom(1) = 0.0d0 !MeV
  bin_bottom(2) = 1.0d0 !MeV
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

  
  allocate(eos_variables(total_eos_variables))


  open(1,file="src/extra_code_and_tables/rho_temp_ye_c.dat",status='old')
  open(2,file="src/extra_code_and_tables/dye.dat")
  open(3,file="src/extra_code_and_tables/M1_nue_enspectra_cen.xg",status='old')


  cont = .true.
  read(3,"(A)") time_string
  nindex = 0

  do while(cont)
     !read eos vars from gr1d evolution

     read(time_string(index(time_string(1:),"= ")+5:index(time_string(1:),"= ")+28),*) current_time
     if(current_time.gt.1.0d0) stop

     do i=1,number_groups
        read(3,*) current_energy,probability_dist(i)
     end do   
     read(3,*)
     read(3,*)
     
     !calculate the blackbody function (dimensionless fermi function in this 
     !case because other terms cancel in the division of E/B (f_avg / f_neq)     
     do i=1,number_groups 
        blackbody_spectrum(i) = fermidirac_dimensionless(energies(i)/eos_variables(tempindex),eos_variables(mueindex)/eos_variables(tempindex))
     end do
     !!insert emmisivity code from microphysical_ec


  integer, intent(in) :: neutrino_species
  real*8, dimension(nspecies,number_groups) :: emissivity
  real*8, dimension(number_groups) :: emissivity_ni56
  real*8, dimension(number_groups) :: blackbody_spectrum

  !Hempel EOS and number of species are set up in readrates
  call nuclei_distribution_Hempel(nspecies,nuclei_A,nuclei_Z,mass_fractions,number_densities,eos_variables)
  emissivity = 0.0d0
  emissivity = 0.0d0
  emissivity_ni56 = 0.0d0
  do i=1,nspecies
     if(i.eq.1)then !if LMP data is not provided for a given nucleus, we will use the rates for 56Ni
         emissivity_ni56(:) = emissivity_from_weak_interaction_rates(56,28,1.0d0,eos_variables,neutrino_species) 
      endif
      if(nucleus_index(nuclei_A(i),nuclei_Z(i)) == 0)then
         emissivity(i,:) = emissivity_ni56(:)*number_densities(i)
      else
        emissivity(i,:) = emissivity_from_weak_interaction_rates(nuclei_A(i),nuclei_Z(i),number_densities(i),&
             eos_variables,neutrino_species)
     end if
  end do



     !!then calculate f_eq
     !!then crunch deltaYe
     

     !read next time and stop if at EOF
     read(3,"(A)",IOSTAT=IO) time_string
     if (IO.lt.0) cont = .false.
     
  enddo

  close(1)
  close(2)
  close(3)

  stop





     eos_variables = 0.0d0
     eos_variables(rhoindex) = xrho
     eos_variables(tempindex) = xtemp
     eos_variables(yeindex) = xye

     !! EOS stuff
     keytemp = 1
     keyerr = 0
     call nuc_eos_full(eos_variables(rhoindex),eos_variables(tempindex), &
          eos_variables(yeindex),eos_variables(energyindex),matter_prs, &
          matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe, &
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

!     write(2,*) eos_variables(abarindex);
     call nuclei_distribution_Hempel(nspecies,nuclei_A,nuclei_Z,mass_fractions,number_densities,eos_variables)
     
     write(2,*) Sum(nuclei_A(:)*number_densities(:))/Sum(number_densities(:))
     




10 close(1)
  close(2)



end program point_example
  
