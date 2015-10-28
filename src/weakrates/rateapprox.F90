!-*-f90-*- 
module class_rateapproximation

  implicit none
  private
  public :: RateApprox, new_RateApprox, return_hempel_qec, weakrates_approx

  ! constructor declaration
  interface new_RateApprox
     module procedure new_RateApproxDefault
  end interface new_RateApprox

  ! members
  type RateApprox
     integer nspecies
     integer, allocatable,dimension(:) :: nuclei_A ! Hempel EOS nuclei
     integer, allocatable,dimension(:) :: nuclei_Z ! Hempel EOS nuclei
     real*8, allocatable,dimension(:) :: number_densities
     real*8, allocatable,dimension(:) :: mass_fractions
  end type RateApprox

  
  ! methods
contains

!------------------------------------------------------------------------------------!
  
  function new_RateApproxDefault() result(this)
    !""" Default RateApprox constructor """

    use nuclei_hempel
    implicit none
    type(RateApprox) this
    
    ! Initialize Hempel EOS dependencies 
    call get_Hempel_number_of_species(this%nspecies) ! returns the total number of nuclei

    allocate(this%nuclei_A(this%nspecies))
    allocate(this%nuclei_Z(this%nspecies))
    allocate(this%number_densities(this%nspecies))
    allocate(this%mass_fractions(this%nspecies))

    this%nuclei_A = 0
    this%nuclei_Z = 0
    this%number_densities = 0
    this%mass_fractions = 0
    
    call get_Hempel_As_and_Zs(this%nuclei_A,this%nuclei_Z)    

    write(*,"(A25,I4,A9)") " Done loading masses for ", this%nspecies," species."
  end function new_RateApproxDefault
  
!------------------------------------------------------------------------------------!
  
  function return_hempel_qec(A,Z_p,Z_d) result(q)

    use sfho_frdm_composition_module, only : sfho_mass
    use nuclei_hempel, only : hempel_lookup_table

    implicit none
    real*8 :: q
    integer, intent(in) :: A
    integer, intent(in) :: Z_p
    integer, intent(in) :: Z_d

    q = sfho_mass(hempel_lookup_table(A,Z_p)) - sfho_mass(hempel_lookup_table(A,Z_d))

  end function return_hempel_qec

!------------------------------------------------------------------------------------!

  function weakrates_approx(n,temperature,q_gs,mue) result(rate)

    integer, intent(in) :: n ! 0 for electron capture rate, 1 for nuetrino energy loss rate
    real*8 :: rate
    real*8 :: temperature
    real*8 :: chi
    real*8 :: eta
    real*8 :: q_gs
    real*8 :: mue
    real*8 :: complete_fermi_integral
    include 'constants.inc'

    chi = (q_gs-2.5d0)/temperature
    eta = mue/temperature + chi

    rate = log(2.0d0)*4.6d0/6146.0d0*(temperature**(5.0d0+n)/m_e**(5.0d0))*&
         (complete_fermi_integral(4+n,eta) - &
         2.0d0*chi*complete_fermi_integral(3+n,eta) + &
         chi**2.0d0*complete_fermi_integral(2+n,eta))

    return

  end function weakrates_approx
!------------------------------------------------------------------------------------!
end module class_rateapproximation
