!-*-f90-*- 
module class_rateapproximation

  implicit none
  private
  public :: RateApprox, new_RateApprox, return_hempel_qec, weakrates_approx, print_approx_reference

  ! constructor declaration
  interface new_RateApprox
     module procedure new_RateApproxDefault
  end interface new_RateApprox

  type RateApprox
     integer nspecies
     integer, dimension(:),pointer :: nuclei_A ! Hempel EOS nuclei
     integer, dimension(:),pointer :: nuclei_Z ! Hempel EOS nuclei
     real*8, dimension(:),pointer :: number_densities
     real*8, dimension(:),pointer :: mass_fractions
  end type RateApprox

  ! class instance
  type(RateApprox), save, target :: this

  !$OMP THREADPRIVATE(this)
  
  ! methods
contains

!------------------------------------------------------------------------------------!
  
  function new_RateApproxDefault() result(approx)
    !""" Default RateApprox constructor """

    use nuclei_hempel
    implicit none
    type(RateApprox), pointer :: approx
       
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

    !$OMP SINGLE
    write(*,"(A29,I4,A9)") "    Done loading masses for ", this%nspecies," species."
    !$OMP END SINGLE
    
    approx => this
    !$OMP PARALLEL COPYIN(this)
    !$OMP END PARALLEL    
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
  subroutine print_approx_reference
    !$OMP SINGLE
    print *, "    Loading approximation. Make reference to: "
    print *, "    ------------------------------------------------------------------------------------------"
    print *, "    | iapprox | Langanke, K., & Mart\'{i}nez-Pinedo, G. (2003).                              |"
    print *, "    |         | Electron capture rates on nuclei and implications for stellar core collapse. |"
    print *, "    |         | Physical Review Letters 90, 241102.                                          |"
    print *, "    |         | http://prl.aps.org/abstract/PRL/v90/i24/e241102                              |"
    print *, "    ------------------------------------------------------------------------------------------"
    !$OMP END SINGLE
  end subroutine print_approx_reference
end module class_rateapproximation
