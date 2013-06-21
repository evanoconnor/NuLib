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
