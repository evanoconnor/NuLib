!-*-f90-*-
program point_example_with_HelmholtzEOS

  use nulib
#if HELMHOLTZ_EOS
  include 'other_eos/helmholtz/implno.dek'
  include 'other_eos/helmholtz/vector_eos.dek'
#else
  implicit none
#endif

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
  character*200 :: eos_filename = "./LS220.h5"

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

  !for the helmholtz EOS, only three species
  integer, parameter :: HELM_ionmax = 3
  real*8 HELM_xmass(HELM_ionmax),HELM_aion(HELM_ionmax),HELM_zion(HELM_ionmax),HELM_abar,HELM_zbar

  !allocate the arrays for the point values
  allocate(local_emissivity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_absopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_scatopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_delta(mypoint_number_output_species,mypoint_number_groups))
  allocate(blackbody_spectra(mypoint_number_output_species,mypoint_number_groups))

  !this sets up many cooefficients and creates the energy grid (one
  !zone + log spacing) see nulib.F90:initialize_nulib
  call initialize_nulib(mypoint_neutrino_scheme,mypoint_number_species,mypoint_number_groups)

  !read in EOS table & set reference mass
#if HELMHOLTZ_EOS
  call read_helm_table
  m_ref = m_amu
#else
  call readtable(eos_filename)
  m_ref = m_n !for LS220
#endif

  !example point
  xrho = 1.0d11 !g/cm^3
  xtemp = 10.0d0 !MeV
  xye = 0.5d0 !dimensionless

  !set up energies bins
  do_integrated_BB_and_emissivity = .false.
  mindx = 1.0d0
  bin_bottom(1) = 0.0d0 !MeV
  bin_bottom(2) = 4.0d0 !MeV
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
  eos_variables = 0.0d0
  eos_variables(rhoindex) = xrho
  eos_variables(tempindex) = xtemp
  eos_variables(yeindex) = xye

  !! EOS stuff
  
#if HELMHOLTZ_EOS
  ! set the mass fractions, z's and a's of the composition, there is a
  ! lot of manual stuff here depending on your application, the
  ! Helmholtz EOS is very general when it comes to compositions

  ! neutron, protons, alphas - for Helmholtz EOS
  HELM_xmass(1) = 0.50d0 ; HELM_aion(1)  = 1.0d0  ; HELM_zion(1)  = 0.0d0
  HELM_xmass(2) = 0.50d0 ; HELM_aion(2)  = 1.0d0  ; HELM_zion(2)  = 1.0d0
  HELM_xmass(3) = 0.0d0 ; HELM_aion(3)  = 4.0d0 ; HELM_zion(3)  = 2.0d0

  ! for NuLib - adjust!!
  eos_variables(xaindex) = HELM_xmass(3)
  eos_variables(xnindex) = HELM_xmass(1)
  eos_variables(xpindex) = HELM_xmass(2)
  eos_variables(xhindex) = 0.0d0 !note this

  ! average atomic weight and charge, used in Helmholtz EOS, not the same abar and zbar of NuLib 
  HELM_abar   = 1.0d0/sum(HELM_xmass(1:HELM_ionmax)/HELM_aion(1:HELM_ionmax))
  HELM_zbar   = HELM_abar*sum(HELM_xmass(1:HELM_ionmax)*HELM_zion(1:HELM_ionmax)/HELM_aion(1:HELM_ionmax))

  !NuLib does not include neutrons, proton, and alpha in abar and
  !zbar, unlike the Helmholtz EOS
  !setting these to 1 since eos_variables(xhindex) = 0.0d0, if you
  !have heavy nuclei other than neutrons,protons, and alphas, you mus
  !set these appropiately
  eos_variables(abarindex) = 1.0d0
  eos_variables(zbarindex) = 1.0d0

  !this overrides the ye specified above - here zbar and abar are an
  !average for ALL species, so defines ye, not true for
  !stellarcollapse.org EOS
  eos_variables(yeindex) = HELM_zbar/HELM_abar

  ! set the input vector. pipeline is only 1 element long
  temp_row(1) = eos_variables(tempindex)/kelvin_to_mev
  den_row(1)  = eos_variables(rhoindex)
  abar_row(1) = HELM_abar
  zbar_row(1) = HELM_zbar
  jlo_eos = 1 ; jhi_eos = 1

  ! call the eos
  call helmeos

  !set eos_variables
  eos_variables(energyindex) = etot_row(1)
  eos_variables(mueindex) = etaele_row(1)*eos_variables(tempindex) + 0.511d0 !add in electron rest mass

  !analytic mu's from EOSmaker on stellarcollapse.org, has correction
  !to mu_p for coulomb, following our convention, we add in the rest
  !mass difference of the neutron and proton into the chemical
  !potentials and into muhat.
  call mu_np(eos_variables(rhoindex),eos_variables(tempindex)/kelvin_to_mev, &
       eos_variables(xnindex),eos_variables(xpindex),eos_variables(yeindex), &
       eos_variables(munindex),eos_variables(mupindex))

  eos_variables(muhatindex) = eos_variables(munindex) - eos_variables(mupindex)

#else
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
#endif

  write(*,*) "Testing: ",eos_variables(xnindex),eos_variables(xpindex), &
       eos_variables(xaindex),eos_variables(munindex), &
       eos_variables(mupindex),eos_variables(mueindex),eos_variables(mueindex)-eos_variables(muhatindex)
  stop
  !! Done EOS stuff
  
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


end program point_example_with_HelmholtzEOS
  
