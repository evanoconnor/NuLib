!-*-f90-*-



function  Bremsstrahlung_Phi0_Hannestad(nu_energy_x,nubar_energy_x&
         ,matter_temperature,n_N,neutrino_species,pro_ann) result(Phi0_nn)
         
! Rates from Hannestad & Raffelt 1998
use nulib
implicit none 
real*8,intent(in) :: nu_energy_x,nubar_energy_x  ! dimensionless
real*8, intent(in) :: matter_temperature ! MeV
real*8, intent(in) :: n_N ! n/cm3
integer, intent(in) :: neutrino_species ! 1 to 6
integer,intent(in) :: pro_ann ! 0 or 1

real*8 :: eta_star
real*8 :: Gamma_sigma
real*8 :: S_sig
real*8 :: p_f
real*8 :: s_func
real*8 :: g 
real*8 :: x
real*8 :: y 
real*8 :: G_f  
real*8 :: G_rampp
real*8 :: GC2_raffelt
real*8,parameter :: C_a = gA/2.0d0
!~ real*8,parameter :: C_a = 0.5d0


!functions
real*8 :: find_g
real*8 :: find_s
real*8 :: find_s2

!output  
real*8 :: Phi0_nn
real*8 :: Phi0_nn_old
real*8 :: Phi0_nn_Raff



!Rampp
p_f = hbarc_mevcm/clight  * (3.0d0 * pi **2 * n_N)**(1.0d0/3.0d0) ! Mev s cm-1

eta_star = (p_f)**2 /(2.0d0 *m_amu/clight**2 * matter_temperature) !adimensional



Gamma_sigma = 8.0d0 * SQRT(2.0d0*pi) * 15.0d0**2 / (3.0d0*pi**2) * &
			  eta_star**(1.5d0) * &
			  matter_temperature**2/m_amu   ! MeV      !alpha_pi = 15
			 

G_f = Gfermi*(hbarc_mevcm)  ! Mev-1 cm
G_rampp = (G_f*hbarc_mevcm**2)**2/(pi*(hbarc_mevcm)**4) !cm2 MeV^-2
GC2_raffelt = 2.1d0 * 1d-44 !cm2  Mev-2

!Raffelt			  

!~ y = 1.94d0* 10.0d0/matter_temperature ! adim
y = m_pi**2/(matter_temperature*m_amu) !adim
x = nu_energy_x + nubar_energy_x  ! adimensionnal
s_func = find_s(eta_star,y,x)  ! adimenssional
g = find_g(eta_star,y)   ! adimenssional



!Phi 0

S_sig = Gamma_sigma / &
			( ( x *matter_temperature)**2  + &
			0.25d0 * ( Gamma_sigma * g) **2 ) &
			* s_func   ! MeV^-1


Phi0_nn = 3.0d0*G_rampp  * C_a**2 * n_N * S_sig &
			*hbarc_mevcm**3 *clight/(2.0d0*pi)**2 ! cm3   s-1
			
Phi0_nn_old = 6.0d0*G_f ** 2 * C_a**2 * n_N * S_sig &
			*(hbarc_mevcm)**3 *clight   ! cm3   s-1
	

Phi0_nn_Raff = GC2_raffelt * n_N * S_sig &
			*hbarc_mevcm**3 *clight/(2.0d0*pi) ! cm3   s-1


if (pro_ann .EQ. 1) then
       Phi0_nn = Phi0_nn_old ! absorption 
else if (pro_ann .EQ. 0) then 
       Phi0_nn = Phi0_nn_old * exp(-x) ! production
else
      write(*,*) "brem : annihilation or production wrongly inserted"
end if


if ( Phi0_nn .NE. Phi0_nn) then 
	write(*,*) "NaNing", s_func,x,y,matter_temperature,eta_star
	stop
endif
if ( Phi0_nn .LT. 0.0d0) then 
	write(*,*) "negative Han brem", s_func,x,y,matter_temperature,eta_star
	stop
endif
         
!~  write(*,*) S_sig
end function Bremsstrahlung_Phi0_Hannestad



function  Bremsstrahlung_Phi0_Kuroda(nu_energy_x,nubar_energy_x&
         ,matter_temperature,n_N,neutrino_species,pro_ann) result(Phi0_nn)

! Kuroda rates from Kuroda et al. 2016 ( Raffelt 1996) modified with factor from Raffelt & Seckel 1995
use nulib
implicit none 
real*8,intent(in) :: nu_energy_x,nubar_energy_x  ! dimensionless
real*8, intent(in) :: matter_temperature ! MeV
real*8, intent(in) :: n_N ! n/cm3
integer, intent(in) :: neutrino_species 
integer,intent(in) :: pro_ann ! 0 or 1

real*8 :: p_f
real*8 :: x
real*8 :: Phi_ND
real*8 :: Phi_ND_2
real*8 :: Phi_D
real*8 :: Phi_D_2
real*8 :: constante
real*8,parameter :: G = 1.55d-33 ! G**2 in Kuroda   Mev-2 cm3 s-1
real*8,parameter :: C_A = gA/2.0d0


!output 
real*8 :: Phi0_nn

!function
real*8 :: find_s2

x = nubar_energy_x + nu_energy_x  ! dimensionless
p_f = hbarc_mevcm/clight * ( 3.0d0 * pi**2 * n_N)**(1.0d0/3.0d0)  ! cm2 s-1  


Phi_D_2 =  G  / 4.0d0 * 15.0d0**2*(4.0d0*pi)**2 * pi**4 & 
		/ (12.0d0 * pi**9) * matter_temperature/x * (4.0d0 * pi**2 + x**2) &
		/ ( 1.0d0 - exp(-x)) * 3.0d0* C_A**2.0d0 * 2.0d0 * hbarc_mevcm  &
		* ( 3.0d0 *pi**2 * n_N) ** (1.0d0/3.0d0)


Phi_ND_2 = hbarc_mevcm**6 * G/4.0d0 * 15.0d0**2 *(4.0d0*pi)**2/m_amu**4&
		*pi**4 * 2.0d0 * m_amu**(1.5d0) &
		/ pi**(5.5d0) /( matter_temperature**(1.5d0) *x**2) &
		*3.0d0 * C_A**2.0d0 * n_N ** 2 * find_s2(0.0d0,0.0d0,x)! using the s_kl part as well 
		

if ( pro_ann .EQ. 0 ) then ! pro
	Phi0_nn = exp(-x)* MIN(Phi_D_2,Phi_ND_2) !cm3 s-1 
else if (pro_ann .EQ. 1 ) then !abs
	Phi0_nn =   MIN(Phi_D_2,Phi_ND_2)  ! cm3 s-1
else 
	stop " pro_ann not allowed"
endif


    
end function Bremsstrahlung_Phi0_Kuroda



function find_g(eta_star,y) result (g)
   implicit none
  real*8, intent(in) :: eta_star
  real*8, intent(in) :: y
  
  !output
  real*8 :: g
  
  real*8 :: alpha1
  real*8 :: alpha2
  real*8 :: alpha3
  real*8 :: p1
  real*8 :: p2
  
  alpha1 = (0.5d0 + 1.0d0/eta_star)/ &
           ((1.0d0 + 1.0d0/eta_star) * (25.0d0 * y**2 + 1.0d0)) + &
           (0.5d0 + eta_star/15.6d0)* (25.0d0 * y**2 )/(25.0d0*y**2 + 1.0d0)
			
			
  alpha2 = (0.63d0 + 0.04d0 * eta_star**(1.45d0) )/&
           (1.0d0 + 0.02d0 * eta_star**(2.5d0)) 

  alpha3 = 1.2d0 * exp(0.6d0*eta_star - 0.4d0 * eta_star**(1.5d0)) 

  p1= (1.8d0 + 0.45d0*eta_star) / ( 1.0d0 + 0.15d0 *eta_star**(1.5d0)) 

  p2 = 2.3d0 - 0.05d0*eta_star/(1.0d0 + 0.025*eta_star) 

  g = (alpha1 + alpha2 *y **(p1))/ &
      (1.0d0 + alpha3 * y **(p2) + alpha2 * y**(p1+2.0d0)/13.75d0)
      
  if (g .LT. 0.0d0 ) stop "brem : g < 0"
  
end function find_g



function find_s(eta, y,x) result(s)
use nulib
implicit none

real*8, intent(in) :: x
real*8,intent(in) :: y 
real*8,intent(in) :: eta

real*8 :: s_ND
real*8 :: s_D
real*8 :: u 
real*8 :: f 
real*8 :: h 
real*8 :: p
real*8 :: C
real*8 :: G
real*8 :: big_F

!result
real*8 :: s


s_ND = 2.0d0 * sqrt(pi) * ( x + 2.0d0 - exp(-y/12.0d0))**1.5d0 * &
			(x**2 + 2.0d0*x*y + 5.0d0/3.0d0 * y**2 +1.0d0) / ( sqrt(pi) &
			 + ( pi **(1.0d0/8.0d0) + x + y)**4) 

		
u = sqrt(y/(2.0d0*eta))

f = 1.0d0 - 5.0d0*u/6.0d0*atan(2.0d0/u) + u**2 /(3.0d0*(u**2 +4.0d0)) &
    + u**2 /(6.0d0*sqrt(2.*u**2 +4.0d0)) * atan( 2.0d0*sqrt(2.0d0*u**2+4.0d0)/u**2)

s_D =3.0d0*(pi/2.0d0)**2.5d0 * eta**(-2.5d0)* (x**2 + 4.0d0*pi**2)*x &
	/(4.0d0*pi**2 *(1.0d0 - exp(-x))) * f



h = 0.1d0 * eta /(2.39d0 + 0.1d0 * eta**1.1)

p = 0.67d0 + 0.18d0*y**0.4d0

C = 1.1d0* x**1.1d0 * h /(2.3 +h * x**0.93d0+0.0001d0*x**1.2)*30.0d0 &
	/(30.0d0+0.005d0*x**2.8d0)

G = 1.0d0 - 0.0044*x**1.1 * y/(0.8d0+0.06d0*y*1.05d0)* sqrt(eta)/(eta+0.2d0)

big_F = 1.0d0 + 1.0d0/(( 3.0d0 +( x-1.2d0)**2 + x**(-4.0d0)) * &
         (1.0d0+ eta**2.0d0)*(1+y**4.0d0))

s = ( s_ND**(-p) + ABS(s_D)**(-p)) ** (-1.0d0/p) * big_F *(1.00d0+ C*G)  ! adim

if ( s.NE.s) then 
	write(*,*) "find_s NaNing",s_ND,s_D, big_F,C,G,eta,y,x
	stop
endif
 
end function find_s

function find_s2(eta, y,x) result(s)
use nulib
implicit none

real*8, intent(in) :: x
real*8,intent(in) :: y 
real*8,intent(in) :: eta

real*8 :: s_ND
real*8 :: s_D
real*8 :: u 
real*8 :: f 
real*8 :: h 
real*8 :: p
real*8 :: C
real*8 :: G
real*8 :: big_F


!result
real*8 :: s


s_ND = sqrt(1.0d0+x*pi/4.0d0) - (1.0d0+(pi*x/4.0d0)**(5.0d0/4.0d0))**0.4d0&
		+ (1+ (64.0d0/(169.0d0*pi)*x)**(5.0d0/4.0d0))**(-0.4d0)   ! (s_0 - s_k.l) from B.12 & B.16 in Raffelt & Seckel 1995
		



s_D = (x**2 + 4.0d0 * pi**2 ) * x/(1.0d0-exp(-abs(x)))/(4.0d0*pi**2) &
		*3.0d0*(pi/2.0d0)**2.5d0/(eta**2.5d0)

s= s_ND   ! used in Kuroda rates


if ( s.NE.s) then 
	write(*,*) "find_s2 NaNing",s_ND,s_D, big_F,C,G,eta,y,x
	stop
endif
 
end function find_s2

function Bremsstrahlung_Phi0_gang(nu_energy_x,nubar_energy_x&
         ,matter_temperature,n_N,Ye,neutrino_species,pro_ann) result(Phi0_nn)
use nulib

implicit none 
real*8,intent(in) :: nu_energy_x,nubar_energy_x  ! dimensionless
real*8, intent(in) :: matter_temperature ! MeV
real*8, intent(in) :: n_N ! n/cm3
real*8, intent(in) :: Ye ! n/cm3
integer, intent(in) :: neutrino_species ! 1 to 6
integer,intent(in) :: pro_ann ! 0 or 1


real*8 :: Phi0_nn

! internal variables 
integer :: in_N,i
real*8 :: Itable_n_N(25)
real*8 :: temp_array(25) = (/(i, i=2,50, 2)/)
real*8 :: Ye_array(26) = (/(i, i=0,50,2)/)/100.0d0
real*8 :: om_array(40) = (/(i,i = 1,40,1)/)/10.0d0 -1.4d0
real*8 :: omega
real*8 :: S,n_N_fm
real*8 :: G_f  
real*8,parameter :: C_a = gA/2.0d0
real*8 ::  eas_brem(25)
real*8 :: dx,delt
integer :: indx
real*8 ::  t_try

 do in_N=1,25
	Itable_n_N(in_N) = &
		 (-4.0d0+dble(in_N-1)/dble(24)*(4.0d0))
 enddo

n_N_fm= log10(n_N*(1.0d-13)**3)

omega= log10((nu_energy_x+nubar_energy_x)*matter_temperature)
S=0.0d0
indx=0

if ( matter_temperature .GT. maxval(temp_array) .or. matter_temperature .LT. minval(temp_array) &
	.or.  Ye .GT. maxval(Ye_array) .or. Ye .LT. minval(Ye_array) &
	.or.  abs(n_N_fm) .GE. maxval(abs(Itable_n_N)) .or. abs(n_N_fm) .LE. minval(abs(Itable_n_N))  & 
	.or.  omega .GT. maxval(om_array) .or. omega .LT. minval(om_array) ) then
	
else 

	call intp3d_many (omega, Ye,n_N_fm, eas_brem, 1, gang_table,40,26,25,25&
						,om_array, Ye_array,Itable_n_N )
						
	
	
	dx = 0.5d0
	indx= 1+ INT((matter_temperature-temp_array(1))*dx)


	delt = (matter_temperature-temp_array(indx)) / 2.0d0 
	
	S = eas_brem(indx) + (eas_brem(indx+1)-eas_brem(indx))*delt
 

	
endif

            


!Rampp
G_f = Gfermi*(hbarc_mevcm)  ! Mev-1 cm

Phi0_nn = 6.0d0*G_f ** 2 * C_a**2 * n_N * S &
			*(hbarc_mevcm)**3 *clight ! cm3   s-1

if (pro_ann .EQ. 1) then
       Phi0_nn = Phi0_nn ! absorption 
else if (pro_ann .EQ. 0) then 
       Phi0_nn = Phi0_nn * exp(-(nu_energy_x+nubar_energy_x)) ! production
else
      write(*,*) "brem : annihilation or production wrongly inserted"
end if



if (Phi0_nn .LT. 0.0d0) then 
	write(*,*) " brem negative : problem"
	write(*,*) "S",S,"Phi0_nn :", Phi0_nn,"Ye :", Ye, "omega :",omega,"n :"&
				,n_N_fm,"n_cm : ",n_N,"T :", matter_temperature,"pro_ann :",pro_ann
	write(*,*) eas_brem
	write(*,*) "arg", shape(eas_brem), shape(gang_table)
	write(*,*)
	write(*,*) "nu,nubar,temp",nu_energy_x,nubar_energy_x&
         ,matter_temperature
    write(*,*) neutrino_species,pro_ann
    write(*,*) 
	write(*,*) om_array
	write(*,*) temp_array
	write(*,*) Ye_array
	write(*,*) Itable_n_N
	stop
endif


end function Bremsstrahlung_Phi0_gang
