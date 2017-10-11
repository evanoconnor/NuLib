!-*-f90-*-
program point_example

  use nulib
  use inputparser
#if NUCLEI_HEMPEL
  use nuclei_hempel, only : set_up_Hempel
#endif
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

  character*200 :: parameters = "/projects/ceclub/gr1dnulib/GitHub/NuLib/parameters"

  !local variables
  real*8, allocatable,dimension(:,:) :: local_emissivity
  real*8, allocatable,dimension(:,:) :: local_absopacity
  real*8, allocatable,dimension(:,:) :: local_scatopacity
  real*8, allocatable,dimension(:,:) :: local_delta
  real*8, allocatable,dimension(:,:) :: local_Phi0, local_Phi1
  real*8, allocatable,dimension(:,:) :: blackbody_spectra
  real*8, allocatable,dimension(:) :: eos_variables
  real*8 :: matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe
  integer :: keytemp,keyerr
  real*8 :: precision = 1.0d-10
  real*8 :: xrho, xtemp, xye
  integer i,j
  real*8 dxfac,mindx

  !allocate the arrays for the point values
  allocate(local_emissivity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_absopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_scatopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_delta(mypoint_number_output_species,mypoint_number_groups))
  allocate(blackbody_spectra(mypoint_number_output_species,mypoint_number_groups))

  call input_parser(parameters)

  !this sets up many cooefficients and creates the energy grid (one
  !zone + log spacing) see nulib.F90:initialize_nulib
  call initialize_nulib(mypoint_neutrino_scheme,mypoint_number_species,mypoint_number_groups)
#if NUCLEI_HEMPEL
  call set_up_Hempel ! set's up EOS for nuclear abundances
#endif
  !read in EOS table & set reference mass
  call read_eos_table(eos_filename)
  m_ref = m_amu !for SFHo_EOS (Hempel)
  ! m_ref = m_n !for LS220

  !example point
  xrho = 1.0d12 !g/cm^3
  xtemp = 1.5d0 !MeV
  xye = 0.35d0 !dimensionless

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
  eos_variables(:) = 0.0d0
  eos_variables(rhoindex) = xrho
  eos_variables(tempindex) = xtemp
  eos_variables(yeindex) = xye

  !! EOS stuff
  call set_eos_variables(eos_variables)
  write(*,*) "Rho: ",eos_variables(rhoindex)," g/ccm"
  write(*,*) "T: ",eos_variables(tempindex)," MeV"
  write(*,*) "Ye: ",eos_variables(yeindex)
  write(*,*) "X_n: ",eos_variables(xnindex)
  write(*,*) "X_p: ",eos_variables(xpindex)
  write(*,*) "X_alpha: ",eos_variables(xaindex)
  
  !calculate the full emissivities and opacities, this averages the
  !values based on the mypoint_neutrino_scheme this assumes detailed
  !balance and uses the opacities to fullying calculate the
  !emissivities and vice versa, see single_point_return_all.
  call single_point_return_all(eos_variables, &
       local_emissivity,local_absopacity,local_scatopacity,local_delta, &
       mypoint_neutrino_scheme)

  write(*,*) "Example of a single point call with returning all emissivity, absorptive opacity, and scattering opacity"
  write(*,*) 
  do i=1,mypoint_number_output_species
     do j=1,mypoint_number_groups
        write(*,"(i4,i4,1P10E18.9)") i,j,energies(j),local_emissivity(i,j),local_absopacity(i,j), &
             local_scatopacity(i,j),local_delta(i,j)
     enddo
  enddo

  write(*,*) 
  write(*,*) 


  !these calls return all six neutrinos so arrays must match, also
  !this does not assume detailed balance values are pure emissivities
  !(i.e. no \nu_e emission from electron capture as that is calculated
  !based on opacities)

  !deallocate the arrays for the point values
  deallocate(local_emissivity)
  deallocate(local_absopacity)
  deallocate(local_scatopacity)
  deallocate(local_delta)
  !allocate the arrays for the point values
  allocate(local_emissivity(mypoint_number_species,mypoint_number_groups))
  allocate(local_absopacity(mypoint_number_species,mypoint_number_groups))
  allocate(local_scatopacity(mypoint_number_species,mypoint_number_groups))
  allocate(local_delta(mypoint_number_species,mypoint_number_groups))

  write(*,*) "Example of single point but only emissivity"
  write(*,*) 
  call return_emissivity_spectra_given_neutrino_scheme(local_emissivity,eos_variables)

  write(*,*) "Example of single point but only absorptive opacity"
  write(*,*) 
  call return_absorption_opacity_spectra_given_neutrino_scheme(local_absopacity,eos_variables)

  write(*,*) "Example of single point but only scattering opacity"
  write(*,*) 
  call return_scattering_opacity_spectra_given_neutrino_scheme(local_scatopacity,local_delta,eos_variables)

  do i=1,mypoint_number_output_species
     do j=1,mypoint_number_groups
        write(*,"(i4,i4,1P10E18.9)") i,j,energies(j),local_emissivity(i,j),local_absopacity(i,j), &
             local_scatopacity(i,j),local_delta(i,j)
     enddo
  enddo

  write(*,*)
  write(*,*)
  write(*,*) "black body function"
  
  call return_blackbody_spectra(blackbody_spectra,eos_variables)
  do i=1,mypoint_number_output_species
     do j=1,mypoint_number_groups
        write(*,"(i4,i4,1P10E18.9)") i,j,energies(j),blackbody_spectra(i,j),blackbody_spectra(i,j),blackbody_spectra(i,j)
     enddo
  enddo

  if (add_nue_Iscattering_electrons.or.add_anue_Iscattering_electrons.or. &
       add_numu_Iscattering_electrons.or.add_anumu_Iscattering_electrons.or. &
       add_nutau_Iscattering_electrons.or.add_anutau_Iscattering_electrons) then

     write(*,*)
     write(*,*)
     write(*,*) "inelastic call"
     
     allocate(local_Phi0(mypoint_number_output_species,mypoint_number_groups))
     allocate(local_Phi1(mypoint_number_output_species,mypoint_number_groups))
     
     call single_Ipoint_return_all(mypoint_number_groups,eos_variables(mueindex)/eos_variables(tempindex), &
          eos_variables(tempindex),local_Phi0,local_Phi1,mypoint_neutrino_scheme)              
     
     write(*,*) "sample scattering kernel for energy",energies(mypoint_number_groups)
     do i=1,mypoint_number_output_species
        do j=1,mypoint_number_groups
           write(*,*) energies(mypoint_number_groups),energies(j),local_Phi0(i,j),local_Phi1(i,j)
        enddo
     enddo
  endif


end program point_example
  
