!-*-f90-*-
function Sion_beta(i,Gamma) result(beta_of_Gamma)

  use nulib
  implicit none

  !inputs
  real*8, intent(in) :: Gamma
  integer :: i !index of Gamma

  !output
  real*8 beta_of_Gamma

  !local
  integer :: j !counter
  real*8 :: sqrtGamma

  select case (i)
     case (0)
        beta_of_Gamma = log(0.300d0/(0.300d0+3.0d0*Gamma)) !ln
     case (1)
        beta_of_Gamma = 0.0d0
     case (2)
        beta_of_Gamma = 6.6666667d0
     case (3)
        sqrtGamma = sqrt(Gamma)
        beta_of_Gamma = Sion_coeffs(3,1)+ &
             Sion_coeffs(3,2)*sqrtGamma + &
             Sion_coeffs(3,3)*Gamma + &
             Sion_coeffs(3,4)*Gamma*sqrtGamma
     case (4)
        sqrtGamma = sqrt(Gamma)
        beta_of_Gamma = Sion_coeffs(4,1)+ &
             Sion_coeffs(4,2)*sqrtGamma + &
             Sion_coeffs(4,3)*Gamma + &
             Sion_coeffs(4,4)*Gamma*sqrtGamma
     case (5)
        sqrtGamma = sqrt(Gamma)
        beta_of_Gamma = Sion_coeffs(5,1)+ &
             Sion_coeffs(5,2)*sqrtGamma + &
             Sion_coeffs(5,3)*Gamma + &
             Sion_coeffs(5,4)*Gamma*sqrtGamma
     case (6)
        sqrtGamma = sqrt(Gamma)
        beta_of_Gamma = Sion_coeffs(6,1)+ &
             Sion_coeffs(6,2)*sqrtGamma + &
             Sion_coeffs(6,3)*Gamma + &
             Sion_coeffs(6,4)*Gamma*sqrtGamma
     case default
        stop "Sion_beta: unknown case"
     end select

end function Sion_beta

subroutine nulib_series2(nzones,xmin,xmax,mindx,dxfac)
! This routine is motivated by the "series2" subroutine of
! Cala resp. Prometheus. and copied from GR1D
!
! It solves for a factor dxfac by which each dx is slightly
! larger than the preceding dx.

  implicit none
  
  real*8 dxfac
  integer nzones
  real*8 xmin,xmax,mindx
  
  ! internal vars
  real*8 tol
  real*8 al,aold,ferror,sum,F,dsum,dFda,anew
  integer k,i,itermax
  
  tol = 1.0d-6
  itermax = 100
  
  ! solve for dxfac
  dxfac=0.0d0

  ! estimate
  al = log( (xmax-xmin)/mindx )/ ( real(nzones-2) )
  aold = exp(al)
  k = 1
  ferror = 1.0d0

  !-------------------------------------------------
  ! Solve: F = (xmax-xmin)/mindx - (Sum[ a^j],j=0,N)
  ! let x = a, y(x) = F
  !-------------------------------------------------

  ! evaluate F
  do while( (ferror.gt.tol).and.(k.lt.itermax))
     sum = 0.0d0
     do i=1,nzones-1
        sum = sum + aold**(i-1)
     enddo
     F = ( (xmax-xmin)/mindx) - sum
     
     ! evaluate dFDa
     dsum = 1.0d0
     do i=4,nzones
        dsum = dsum + (i-2)*(aold**(i-3))
     enddo
     dFda = -1.0d0*dsum
     ! next root
     anew = aold - F/dFda
     ferror = abs(anew-aold)/aold
     k = k + 1
     aold = anew
  enddo
  dxfac = anew
  
end subroutine  nulib_series2

subroutine mu_np(xrho,xtemp,xn,xp,xye,mu_n,mu_p)
  !taken from EOSmaker on stellarcollapse.org
  ! calculate chemical potentials of neutrons and protons at low densities
  use nulib, only: m_n,m_p,mev_to_erg,kelvin_to_mev,planck,pi,clight,avo
  implicit none
  
  real*8 :: xrho,xtemp,xn,xp,xye
  real*8 :: mu_n,mu_p
  real*8 :: n_n,n_p,a_e,Gamma_e,Gamma_p,mu_p_coul
  real*8 :: fac_1_erg,fac_1_mev,fac_2
  real*8,  parameter :: A_1 = -0.9052d0, &
       &                  A_2 = 0.6322d0,  &
       &                  A_3 = - sqrt(3.0d0)/2.0d0 - A_1/sqrt(A_2)
  real*8, parameter  :: m_n_cgs = m_n*mev_to_erg/clight**2,&
       &                  m_p_cgs = m_p*mev_to_erg/clight**2
  real*8, parameter  :: electron_charge = 4.803204184d-10
  real*8, parameter  :: tiny = 1.0d-20
  
  fac_1_erg = xtemp*kelvin_to_mev*mev_to_erg
  fac_1_mev = xtemp*kelvin_to_mev
  fac_2 = (planck**2/(2.0*pi*fac_1_erg))**1.5
  
  if(xn.le.0.0d0) then
     xn = tiny
  endif
  if(xp.le.0.0d0) then
     xp = tiny
  endif
  n_n = xn*xrho/m_n_cgs
  n_p = xp*xrho/m_p_cgs
  
  a_e = (4.0d0/3.0d0*pi*xrho*xye*avo)**(-1.0d0/3.0d0)
  Gamma_e = electron_charge**2/(fac_1_erg*a_e)
  Gamma_p = Gamma_e
  
  ! chemical potentials for Boltzmann gas (in unit of MeV)
  ! include rest mass but normalized by free neutron mass
  ! when n_n = 0 or n_p = 0, set mu_n = 0.0 and mu_p = 0, correspondingly
  ! 
  if(n_n.gt.0.0d0) then
     mu_n = fac_1_mev*log(n_n/2.0d0*fac_2/m_n_cgs**1.5) + m_n - m_n
  else 
     mu_n = 0.0d0
  endif
  if(n_p.gt.0.0d0) then
     mu_p = fac_1_mev*log(n_p/2.0d0*fac_2/m_p_cgs**1.5) + m_p - m_n
  else
     mu_p = 0.0d0
  endif
  ! coulomb correction for charged particles from Chabrier & Potkhin (1998)
  mu_p_coul = fac_1_mev * (A_1*(sqrt(Gamma_p*(A_2 + Gamma_p)) &
       - A_2 * log(sqrt(Gamma_p/A_2) + sqrt(1.0d0+Gamma_p/A_2))) &
       + 2.0d0*A_3*(sqrt(Gamma_p) - atan(sqrt(Gamma_p))))
  
  mu_p = mu_p + mu_p_coul
end subroutine mu_np

