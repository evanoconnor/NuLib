!-*-f90-*-
module class_rateapproximation

  implicit none
  private
  public :: RateApprox, new_RateApprox, return_hempel_qec, weakrates_approx, weakrates_approx_raduta, print_approx_reference

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

  !------------------------------------------------------------------
  function weakrates_approx_raduta(n,temperature,q_gs,mue,rho,A,Z,model) result(rate)

    integer, intent(in) :: n ! 0 for electron capture rate, 1 for nuetrino energy loss rate
    real*8 :: rate
    real*8 :: temperature
    real*8 :: chi
    real*8 :: eta
    real*8 :: q_gs
    real*8 :: mue
    real*8 :: rho
    real*8 :: delta_E
    integer, intent(in) :: A
    integer, intent(in) :: Z
    integer, intent(in) :: model
    real*8 :: complete_fermi_integral
    real*8 :: get_delta_e_raduta
    include 'constants.inc'

    delta_E = get_delta_e_raduta(temperature, rho, A, Z, model)
    chi = (q_gs-delta_E)/temperature
    eta = mue/temperature + chi

    rate = log(2.0d0)*4.6d0/6146.0d0*(temperature**(5.0d0+n)/m_e**(5.0d0))*&
         (complete_fermi_integral(4+n,eta) - &
         2.0d0*chi*complete_fermi_integral(3+n,eta) + &
         chi**2.0d0*complete_fermi_integral(2+n,eta))

    return

  end function weakrates_approx_raduta
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
    print *, "    | Ad. R. Raduta (2017)                                                                   |"
    print *, "    | Stellar electron capture rates on neutron-rich nuclei and their impact on              |"
    print *, "    | stellar core collapse. Physical Review C, 95, 025805.                                  |"
    print *, "    | https://journals.aps.org/prc/abstract/10.1103/PhysRevC.95.025805                       |"
    print *, "    ------------------------------------------------------------------------------------------"
    !$OMP END SINGLE
  end subroutine print_approx_reference
  !------------------------------------------------------------------------------------!
end module class_rateapproximation
! !------------------------------------------------------------------------------------!
! function complete_fermi_integral(ifermi,eta)
!   implicit none
!   integer ifermi
!   real*8 complete_fermi_integral
!   real*8 eta
!   real*8 fermi_integral_analytical

!   fermi_integral_analytical = 0.0d0

!   ! Expressions for Fermi integrals given in Takahashi et al. 1978
!   if (eta.gt.1.D-3) then
!      select case (ifermi)
!      case (0)
!         fermi_integral_analytical = &
!              log10(1.0d0+exp(eta))
!      case (1)
!         fermi_integral_analytical = &
!              (eta**2/2.0D0 + 1.6449d0)/(1.0D0+EXP(-1.6855d0*eta))
!      case (2)
!         fermi_integral_analytical = &
!              (eta**3/3.0D0 + 3.2899d0*eta)/(1.0D0-EXP(-1.8246d0*eta))
!      case (3)
!         fermi_integral_analytical = &
!              (eta**4/4.0D0 + 4.9348d0*eta**2+11.3644d0) / &
!              (1.0D0+EXP(-1.9039d0*eta))
!      case (4)
!         fermi_integral_analytical = &
!              (eta**5/5.0D0 + 6.5797d0*eta**3+45.4576d0*eta) / &
!              (1.0D0-EXP(-1.9484d0*eta))
!      case (5)
!         fermi_integral_analytical = &
!              (eta**6/6.0D0 + 8.2247d0*eta**4 + 113.6439d0*eta**2 + &
!              236.5323d0)/(1.0D0+EXP(-1.9727d0*eta))
!      end select

!   else
!      select case (ifermi)
!      case (0)
!         fermi_integral_analytical = &
!              log10(1.0d0+exp(eta))
!      case (1)
!         fermi_integral_analytical = &
!              EXP(eta)/(1.0D0+0.2159d0*EXP(0.8857d0*eta))
!      case (2)
!         fermi_integral_analytical = &
!              2.0D0*EXP(eta)/(1.0D0+0.1092d0*EXP(0.8908d0*eta))
!      case (3)
!         fermi_integral_analytical = &
!              6.0D0*EXP(eta)/(1.0D0+0.0559d0*EXP(0.9069d0*eta))
!      case (4)
!         fermi_integral_analytical = &
!              24.0D0*EXP(eta)/(1.0D0+0.0287d0*EXP(0.9257d0*eta))
!      case (5)
!         fermi_integral_analytical = &
!              120.0D0*EXP(eta) / (1.0D0 + 0.0147d0*EXP(0.9431d0*eta))
!      end select

!   endif
!   complete_fermi_integral = fermi_integral_analytical

!   return
! end function complete_fermi_integral
! !------------------------------------------------------------------------------------!
