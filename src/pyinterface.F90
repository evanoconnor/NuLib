!-*-f90-*-
module pynulib
  use class_ratelibrary
  use weakrates_interface, only : weakratelib

  ! the weak rate library object
!  type(RateLibrary) :: weakrate_library
  
  double precision, dimension(8140,18) :: global_emissivity
  double precision, dimension(18) :: global_emissivity_freep
  double precision, dimension(18) :: global_blocking_factor
  double precision :: global_dyedt_freep
  double precision, dimension(8140) :: global_dyedt
  
  integer, allocatable,dimension(:) :: nuclei_A
  integer, allocatable,dimension(:) :: nuclei_Z
  double precision, allocatable,dimension(:) :: number_densities
  double precision, allocatable,dimension(:) :: mass_fractions
  integer, allocatable,dimension(:,:) :: nucleus_index
  
  integer :: nspecies
  
  contains

    subroutine standard_nulib_init
      use nuclei_hempel
      use inputparser
      use nulib
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
      integer :: mytable_neutrino_scheme = 1

      !number of species for nulib to calculation interactions for, must
      !be six currently, average is done via above parameter
      integer :: mytable_number_species = 6

      !number of energy groups
      integer :: mytable_number_groups = 18

      !NuLib parameters file (weak rates and EOS)
      character(200) :: parameters_filename = "./parameters"
      double precision dxfac,mindx
      integer :: i

      call input_parser(parameters_filename)
      !this sets up many coefficients and creates the energy grid (one
      !zone + log spacing) see nulib.F90:initialize_nulib
      call initialize_nulib(mytable_neutrino_scheme,mytable_number_species,mytable_number_groups)
      call read_eos_table(eos_filename) !read in EOS table & set reference mass
      call set_up_Hempel !set's up EOS for nuclear abundances

      !!!!! The initialization of the weak rate library should be done through a seperate function call !!!!!
      ! initialize weak-rate library if it is turned on in requested_interactions.inc
      !if (add_nue_emission_weakinteraction_ecap.or.add_anue_emission_weakinteraction_poscap) then
      !   call initialize_weakratelib(parameters_filename)
      !endif

      adhoc_nux_factor = 0.0d0 !increase for adhoc nux heating (also set
      !add_nux_absorption_on_n_and_p to true)
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

    end subroutine standard_nulib_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine weakrate_init()
      use class_ratelibrary
      use nuclei_hempel
      character(200) :: parameters_filename = "./parameters"

      integer a, z
      
      call get_Hempel_number_of_species(nspecies)
      
      allocate(nuclei_A(nspecies))
      allocate(nuclei_Z(nspecies))
      allocate(number_densities(nspecies))
      allocate(mass_fractions(nspecies))
      allocate(nucleus_index(nspecies,nspecies))
      
      call get_Hempel_As_and_Zs(nuclei_A,nuclei_Z)
      
      weakratelib = new_RateLibrary(parameters_filename)
      
    end subroutine weakrate_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    subroutine get_dyedt(xrho,xtemp,xye)
      use nuclei_hempel 
      use inputparser
      use nulib, only : tempindex,  mueindex, rhoindex, yeindex, kelvin_to_mev,& 
           add_nue_emission_weakinteraction_ecap, single_point_return_all,neutrino_scheme,&
           mev_to_erg,energies,bin_widths, pi
      use sfho_frdm_composition_module, only : sfho_mass
      use weakrates_interface, only : emissivity_from_weak_interaction_rates
      use class_ratelibrary, only: in_table
      
      implicit none
      
      integer i
      double precision, intent(in) :: xrho,xtemp,xye
      double precision, dimension(15) :: eos_variables
      !          double precision, dimension(8140,18) :: emissivity,emissivity_temp
      double precision, dimension(3,18) :: local_emissivity
      double precision, dimension(3,18) :: local_absopacity
      double precision, dimension(3,18) :: local_scatopacity
      double precision, dimension(3,18) :: blackbody_spectra
      
      integer, dimension(5),save :: file_priority
      integer idxtable, A, Z 
      
      double precision :: q
      double precision :: logrhoYe,t9
      logical :: parameterized_rate  

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! integer j, k, myIndex
      ! !number of nuclei in high sensitivity region: 30, 74, 315
      ! integer, parameter :: nNuclei = 74
      ! double precision, dimension(nNuclei) :: hsA, hsZ
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      eos_variables = 0.0d0
      eos_variables(rhoindex) = xrho
      eos_variables(tempindex) = xtemp
      eos_variables(yeindex) = xye

      call set_eos_variables(eos_variables)
      
      !Hempel EOS and number of species are set up in readrates
      call nuclei_distribution_Hempel(&
           nspecies,nuclei_A,nuclei_Z,mass_fractions,number_densities,eos_variables)          
      global_emissivity = 0.0d0
      logrhoYe = log10(eos_variables(rhoindex)*eos_variables(yeindex))
      t9 = (eos_variables(tempindex)/kelvin_to_mev)*1.0d-9

      do i=1,nspecies
         
         parameterized_rate = .false.

         if(number_densities(i).eq.0.0d0)cycle

         A = nuclei_A(i) 
         Z = nuclei_Z(i)

         !if rate data from a table is not present and A>4 with iapprox nonzero,
         !use the parameterized rate function, else skip this nucleus
         ! check if rate exists in a table
         if (weakratelib%ntables.ne.0)then

            idxtable = in_table(weakratelib,A,Z,logrhoYe,t9)

         endif

         ! if no table contains the requested rate
         if(idxtable.eq.0) then
            ! and if the approximation is turned on
            if(weakratelib%priority(size(weakratelib%priority)).gt.0) then
               ! and the nucleus is above the A=4 isobars and below
               if (A.gt.4) then
               else
                  cycle
               end if
               ! check to make sure masses exist for both parent and daughter nucleus
               if(hempel_lookup_table(A,Z).eq.0.or.hempel_lookup_table(A,Z-1).eq.0) then
                  cycle
               end if

            else
               cycle
            end if
         end if

         !emissivity calculation
         global_emissivity(i,:) = emissivity_from_weak_interaction_rates(A,Z,&
              number_densities(i),eos_variables,1,idxtable)

         !calculate dyedt
         global_dyedt(i) = (4.0d0*pi/6.02214129d23/mev_to_erg/eos_variables(rhoindex))*Sum(bin_widths(:)*global_emissivity(i,:)*global_blocking_factor(:)/energies(:))

      end do
      
      !calculate emissivity for electron capture on free protons
      !note that in order to only calculate the above for free protons,electron
      !capture must be turned off in requested_interactions.inc or manually here
      add_nue_emission_weakinteraction_ecap = .false.    
      call single_point_return_all(eos_variables, &
           local_emissivity,local_absopacity,local_scatopacity,neutrino_scheme)
      global_emissivity_freep = local_emissivity(1,:)
      global_dyedt_freep = (4.0d0*pi/6.02214129d23/mev_to_erg/eos_variables(rhoindex))*Sum(bin_widths(:)*global_emissivity_freep(:)*global_blocking_factor(:)/energies(:))          

      ! resetting state
      add_nue_emission_weakinteraction_ecap = .true.    
      number_densities = 0.0d0
      mass_fractions = 0.0d0

      return

!      print *, " "
      
    end subroutine get_dyedt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine calculate_neutrino_blocking(xrho,xtemp,xye,xalp,f_nu)
      use nulib, only : mueindex, muhatindex, energies, rhoindex, tempindex, yeindex
      use nuclei_hempel
      
      double precision, intent(in) :: xrho,xtemp,xye,xalp
      double precision, intent(in),dimension(18) :: f_nu
      double precision, dimension(18) :: f_nu_eq
      double precision, dimension(15) :: eos_variables
      double precision :: xeta

      eos_variables = 0.0d0
      global_blocking_factor = 0.0d0

      eos_variables(rhoindex) = xrho
      eos_variables(tempindex) = xtemp
      eos_variables(yeindex) = xye
      
      call set_eos_variables(eos_variables)
      xeta = (eos_variables(mueindex) - eos_variables(muhatindex))/xtemp
      f_nu_eq(:) = 1.0d0/(exp(energies(:)/xtemp - xeta) + 1.0d0)
      global_blocking_factor(:) = xalp*(1-f_nu(:)/f_nu_eq(:))

    end subroutine calculate_neutrino_blocking
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine example
      use nuclei_hempel
      use class_ratetable
      use class_rateapproximation
      use nulib

      implicit none
      ! parameters file containing tables and loading priority
      character(200) :: parameters_filename = "./parameters"


      ! function parameters
      integer :: A,Z,table_index
      double precision :: T9,logrhoye
      double precision :: temp_mev, qvalue_mev, echempot_mev
      double precision :: density_gcm3, ye
      double precision :: rate

      double precision :: eos_variables(15)

      ! Initialization
      m_ref = m_amu !sets reference mass for NSE
      call set_up_Hempel !set's up EOS for nuclear abundances
      weakratelib = new_RateLibrary(parameters_filename)
      call readtable(weakratelib%eos_path) !read in EOS table

      ! ------------------------------------------------------------------------ !
      ! There are three ways to access the weak rates.                           !
      ! But in each case, the return_weakrate function interface is used         !
      ! and the difference lies in the function parameters that are passed in.   !
      ! These three methods are detailed below.                                  !
      ! ------------------------------------------------------------------------ !

      A = 56
      Z = 28  ! Ni56
      T9 = 10.0d0 ! 10 GK
      logrhoye = 12.0d0 ! log10(density*ye [g/cm3])
      table_index = in_table(weakratelib,A,Z,logrhoye,T9) ! retrieve table containing rate

      rate = return_weakrate(weakratelib,A,Z,T9,logrhoye,table_index,2)
      print *, "return_weakrate_from_table: ",log10(rate)

    end subroutine example



end module pynulib
