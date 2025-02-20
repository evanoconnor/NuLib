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
  integer :: mypoint_number_groups = 18 !18

  character*200 :: parameters = "./parameters"

  !local variables
  real*8, allocatable,dimension(:,:) :: local_emissivity
  real*8, allocatable,dimension(:,:) :: local_absopacity
  real*8, allocatable,dimension(:,:) :: local_scatopacity
  real*8, allocatable,dimension(:,:) :: local_delta
  real*8, allocatable,dimension(:,:) :: local_Phi0, local_Phi1
  real*8, allocatable,dimension(:,:,:) :: local_Phi0_epannihil, local_Phi1_epannihil
  real*8, allocatable,dimension(:,:,:) :: local_Phi0_bremsstrahlung
  real*8, allocatable,dimension(:,:,:) :: local_Phi0_bremsstrahlung2
  real*8, allocatable,dimension(:,:) :: blackbody_spectra
  real*8, allocatable,dimension(:) :: eos_variables
  real*8,allocatable, dimension(:) :: Phi0_brem,Phi0_brem2
  real*8,allocatable, dimension(:) :: Phi0_brem_ann,Phi0_brem2_ann
  real*8 :: matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe
  real*8 :: single_neutrino_emissivity_from_NNBrem_given_energypoint
  real*8 :: single_neutrino_emissivity_from_epannihil_given_energypoint
  integer :: keytemp,keyerr
  real*8 :: precision = 1.0d-10
  real*8 :: xrho, xtemp, xye
  real*8 :: n_N,eta_star
  integer :: i,j,inde
  real*8 :: dxfac,mindx
  real*8 :: Q,Q_ann,Q2,Q2_ann,M,dE_dEdt,dE_dEdt2,dE_dEdt3,k
  real*8 :: Kuro_inv_lamb,Hann_inv_lamb,epan_inv_lamb
  real*8 :: dn_dEdt,dn_dEdt2,dn_dEdt3,ratio,ratio_2,dn_dEdt_burrows
  real*8 :: dE_dEdt_2,dE_dEdt2_2,dE_dEdt3_2,dE_dEdt_burrows
  real*8 :: Bremsstrahlung_Phi0_Hannestad
  real*8 :: Bremsstrahlung_Phi0_Kuroda
  real*8 :: epannihil_Phi_Bruenn
!~   real*8 :: find_s
  real*8,dimension(100) :: integral1,integral2,integral3,temp_array,dens_array,ener_array
  real*8:: integral
  real*8 :: J_1,J_1_bar,phi_a,phi_p,energy,energy_2
  real*8,dimension(100) :: M1_mom,M1_mom_inv
  
  
  !fermi function
  real*8 :: fermidirac,fermidirac_dimensionless
  real*8 :: x
  real*8 :: find_s,find_s2
  real*8 :: find_g
  real*8 :: nue_absorption_on_n
  real*8 :: anue_absorption_on_p
  

  real*8,parameter :: nulib_emissivity_gf = 5.59424238d-55/(6.77140812d-06**3*2.03001708d+05)
  real*8,parameter :: nulib_opacity_gf = 1.0d0/6.77140812d-06
  real*8,parameter :: nulib_energy_gf = 1.60217733d-6*5.59424238d-55
  real*8,parameter :: nulib_kernel_gf = 6.77140812d-06**3/2.03001708d+05
  
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
  
  xrho =  4.0d14 !g/cm^3
  xtemp = 12.0d0 !MeV
  xye = 0.1d0! !dimensionless
  
  !set up energies bins
  do_integrated_BB_and_emissivity = .true.
  mindx = 2.0d0
  bin_bottom(1) = 0.0d0 !MeV
  bin_bottom(2) = 2.0d0 !MeV
  bin_bottom(3) = bin_bottom(2)+mindx
  bin_bottom(number_groups) = 200.0d0 ! MeV
  
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
  write(*,*) energies
  stop
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
  
  write(*,*) "eos", "xp",eos_variables(xpindex),"xn",eos_variables(xnindex), &
				"mue",eos_variables(mueindex),"xa",eos_variables(xaindex), "xh",&
				eos_variables(xhindex)

  write(*,*) eos_variables(xpindex)+eos_variables(xnindex)+eos_variables(xaindex)&
			+eos_variables(xhindex)

  write(*,*) "emi pair 10 bin nue",single_neutrino_emissivity_from_epannihil_given_energypoint( &
     1,energies(10),eos_variables)

  !calculate the full emissivities and opacities, this averages the
  !values based on the mypoint_neutrino_scheme this assumes detailed
  !balance and uses the opacities to fullying calculate the
  !emissivities and vice versa, see single_point_return_all.
  call single_point_return_all(eos_variables, &
       local_emissivity,local_absopacity,local_scatopacity,local_delta, &
       mypoint_neutrino_scheme)
       
write(*,*) local_emissivity(1,6), energies(6)
write(*,*) local_absopacity(1,6), energies(6)
write(*,*) nue_absorption_on_n(energies(6),eos_variables), energies(6)
write(*,*) anue_absorption_on_p(energies(6),eos_variables), energies(6)
stop
!~   write(*,*) "Example of a single point call with returning all emissivity, absorptive opacity, and scattering opacity"
!~   write(*,*) 
!~   do i=1,mypoint_number_output_species
!~      do j=1,mypoint_number_groups
!~         write(*,"(i4,i4,1P10E18.9)") i,j,energies(j),eos_variables(tempindex),local_emissivity(i,j),local_absopacity(i,j), &
!~              local_scatopacity(i,j),local_delta(i,j)
!~      enddo
!~   enddo

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
     
!~      write(*,*) "sample scattering kernel for energy",energies(mypoint_number_groups)
     do i=1,mypoint_number_output_species
        do j=1,mypoint_number_groups
!~            write(*,*) energies(mypoint_number_groups),energies(j),local_Phi0(i,j),local_Phi1(i,j)
        enddo
     enddo
  endif
  
  
  if (add_nue_kernel_epannihil.or.add_anue_kernel_epannihil.or. &
       add_numu_kernel_epannihil.or.add_anumu_kernel_epannihil.or. &
       add_nutau_kernel_epannihil.or.add_anutau_kernel_epannihil) then

     write(*,*)
     write(*,*)
     write(*,*) "epan call"
     
     allocate(local_Phi0_epannihil(mypoint_number_output_species,mypoint_number_groups,2))
     allocate(local_Phi1_epannihil(mypoint_number_output_species,mypoint_number_groups,2))
     
     call single_epannihil_kernel_point_return_all(mypoint_number_groups, &
          eos_variables(mueindex)/eos_variables(tempindex), &
          eos_variables(tempindex),local_Phi0_epannihil,local_Phi1_epannihil,mypoint_neutrino_scheme)              
   
	 write(*,*) "epan absorption for nu_x",local_Phi0_epannihil(3,:,1)


     deallocate(local_Phi0_epannihil)
     deallocate(local_Phi1_epannihil)
     
  endif
  
  if (add_nue_kernel_bremsstrahlung.or.add_anue_kernel_bremsstrahlung.or. &
       add_numu_kernel_bremsstrahlung.or.add_anumu_kernel_bremsstrahlung.or. &
       add_nutau_kernel_bremsstrahlung.or.add_anutau_kernel_bremsstrahlung) then

     write(*,*)
     write(*,*)
     write(*,*) "epan call"
     
     allocate(local_Phi0_bremsstrahlung(mypoint_number_output_species,mypoint_number_groups,2))     
          
	 call single_bremsstrahlung_kernel_point_return_all_Hannestad(i, &
	          n_N,eos_variables(tempindex),local_Phi0_bremsstrahlung,mypoint_neutrino_scheme) 
	          
	          
	write(*,*) "bremsstrahlung kernel absorption for nu_x", local_Phi0_bremsstrahlung 

     deallocate(local_Phi0_bremsstrahlung)
  endif
  
	


end program point_example
  
