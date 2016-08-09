!-*-f90-*-
module pynulib
  use class_ratelibrary

  ! the weak rate library object
  type(RateLibrary) :: weakrate_library

  contains

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
