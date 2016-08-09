!-*-f90-*-
module pynulib
  use class_ratelibrary

  ! the weak rate library object
  type(RateLibrary) :: weakrate_library

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
      character*200 :: parameters_filename = "./parameters"
      real*8 dxfac,mindx
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


    subroutine example
      use nuclei_hempel
      use class_ratetable
      use class_rateapproximation
      use nulib

      implicit none
      ! parameters file containing tables and loading priority
      character*200 :: parameters_filename = "./parameters"


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
      weakrate_library = new_RateLibrary(parameters_filename)
      call readtable(weakrate_library%eos_path) !read in EOS table

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
      table_index = in_table(weakrate_library,A,Z,logrhoye,T9) ! retrieve table containing rate

      rate = return_weakrate(weakrate_library,A,Z,T9,logrhoye,table_index,2)
      print *, "return_weakrate_from_table: ",log10(rate)

    end subroutine example



end module pynulib
