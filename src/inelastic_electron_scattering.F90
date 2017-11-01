!-*-f90-*-
function NES_Phi0_ThompsonBruenn(nu_energy_in,nu_energy_out,matter_eta,matter_temperature,neutrino_species) result(Phi0)

  use nulib, only : Gfermi,clight,hbarc_mevcm,pi
  implicit none

  !inputs
  real*8, intent(in) :: nu_energy_in !incoming neutrino energy, MeV
  real*8, intent(in) :: nu_energy_out !outgoing neutirno energy, MeV
  real*8, intent(in) :: matter_eta !electron chemical potential / Temp, dimensionless
  real*8, intent(in) :: matter_temperature !temperature, MeV

  integer, intent(in) :: neutrino_species !neutrino species, will affect some coeffs

  !output
  real*8 :: Phi0 !dimensionless

  !local, GLQ integration
  integer :: i
  real*8 :: electron_x,fem,feother !the underscore x denotes /Temp
  real*8 :: nu_energy_in_x,nu_energy_out_x


  !function declarations
  real*8 :: NES_H0_ThompsonBruenn !function declaration
  real*8 :: fermidirac_dimensionless !function declaration

  real*8 :: bound,width
  integer :: npoints=128

  Phi0 = 0.0d0
  bound = max(0.0d0,matter_eta)+10.0d0 ! go 10 past the location of the FD exponential cutoff
  width = bound/float(npoints)
  nu_energy_in_x = nu_energy_in/matter_temperature
  nu_energy_out_x = nu_energy_out/matter_temperature 

  do i=1,npoints
     !must take interval change into account
     electron_x = width * (float(i)-0.5d0) !MeV/Temp
     fem = fermidirac_dimensionless(electron_x,matter_eta)
     feother = 1.0d0-fermidirac_dimensionless(electron_x+nu_energy_in_x - &
          nu_energy_out_x,matter_eta)
     Phi0 = Phi0 + &
          fem*feother*NES_H0_ThompsonBruenn(nu_energy_in_x,nu_energy_out_x, &
          electron_x,neutrino_species)
  enddo
  Phi0 = Phi0*width

  !constants that do not have to do with the integral over electron energy
  !need to restore temperatures, two factors, back to cgs units then
  !needs multiplying these appropriate constants
  !Gfermi**2*hbarc_mevcm**2*clight/pi, then units are cm^3/s (Bruenn C50)
  Phi0 = matter_temperature**2*Gfermi**2*hbarc_mevcm**2*clight*Phi0/ &
       (nu_energy_in_x**2*nu_energy_out_x**2*pi) 

end function NES_Phi0_ThompsonBruenn

function NES_Phi1_ThompsonBruenn(nu_energy_in,nu_energy_out,matter_eta, &
     matter_temperature,neutrino_species) result(Phi1)

  use nulib, only : Gfermi,clight,hbarc_mevcm,pi
  implicit none

  !inputs
  real*8, intent(in) :: nu_energy_in !incoming neutrino energy, MeV
  real*8, intent(in) :: nu_energy_out !outgoing neutirno energy, MeV
  real*8, intent(in) :: matter_eta !electron chemical potential / T, dimensionless
  real*8, intent(in) :: matter_temperature !temperature, MeV
  
  integer, intent(in) :: neutrino_species !neutrino species, will affect some coeffs

  !output
  real*8 :: Phi1 !dimensionless

  !local, GLQ integration
  integer :: i
  real*8 :: electron_x,fem,feother !the underscore x denotes /Temp
  real*8 :: nu_energy_in_x,nu_energy_out_x
  
  !function declarations
  real*8 :: NES_H1_ThompsonBruenn !function declaration
  real*8 :: fermidirac_dimensionless !function declaration

  real*8 :: bound,width
  integer :: npoints=128

  Phi1 = 0.0d0
  bound = max(0.0d0,matter_eta)+10.0d0 ! go 10 past the location of the FD exponential cutoff
  width = bound/float(npoints)
  nu_energy_in_x = nu_energy_in/matter_temperature
  nu_energy_out_x = nu_energy_out/matter_temperature 

  do i=1,npoints
     !must take interval change into account
     electron_x = bound/float(npoints) * (float(i)-0.5) !MeV/Temp
     fem = fermidirac_dimensionless(electron_x,matter_eta)
     feother = 1.0d0-fermidirac_dimensionless(electron_x+nu_energy_in_x-nu_energy_out_x,matter_eta)
     Phi1 = Phi1 + &
          fem*feother*NES_H1_ThompsonBruenn(nu_energy_in_x,nu_energy_out_x,electron_x,neutrino_species)
  enddo
  Phi1 = Phi1*width
  
  !constants that do not have to do with the integral over electron energy
  !need to restore temperatures, two factors, back to cgs units then
  !needs multiplying these appropriate constants
  !Gfermi**2*hbarc_mevcm**2*clight/pi, then units are cm^3/s
  Phi1 = matter_temperature**2*Gfermi**2*hbarc_mevcm**2*clight*Phi1/(nu_energy_in_x**2*nu_energy_out_x**2*pi) 

end function NES_Phi1_ThompsonBruenn

function NES_H0_ThompsonBruenn(nu_energy_in,nu_energy_out,electron_energy,neutrino_species) result(H0)
  
  use nulib, only : H0_constants
  implicit none
  
  !inputs
  real*8, intent(in) :: nu_energy_in !neutrino energy, MeV, or reduced energy
  real*8, intent(in) :: nu_energy_out !antineutrino energy, MeV, or reduced energy
  real*8, intent(in) :: electron_energy !electron energy, MeV, or reduced energy

  integer, intent(in) :: neutrino_species 

  !output
  real*8 :: H0 !final answer

  !local
  real*8 :: H0I,H0II !J functions
  real*8 :: Ea,Eb,Ec !rename energy variables for ease

  Ea = nu_energy_in
  Eb = nu_energy_out
  Ec = electron_energy

  if (Ec+Ea.lt.Eb) then
     H0 = 0.0d0
     return
  endif
  
  !Bruenn 1985, Yueh & Buchler 1977 (with corrections from Bruenn 85)
  if (Eb.ge.Ec) then
     H0I = L0I(Ea,Eb,Ec) 

     H0II = L0II(Ea,Eb,Ec)

  else
     H0I = a0(Ea,Eb) + b0(Ea,Eb)*Ec + c0(Ea,Eb)*Ec**2

     H0II = a0(Eb,Ea) - b0(Eb,Ea)*Ec + c0(Eb,Ea)*Ec**2 !invert Ea and Eb, flip sign on b0 term to use a0's

  endif

  H0 = H0_constants(1,neutrino_species)*H0I + H0_constants(2,neutrino_species)*H0II !same coefficients as e-e+ annhilation

contains
  function HeavySide(arguement) result(answer)

    implicit none
    
    !inputs
    real*8, intent(in) :: arguement

    !output
    real*8 :: answer

    if (arguement.gt.0.0d0) then
       answer = 1.0d0
    else if (arguement.eq.0.0d0) then
       answer = 0.5d0
    else
       answer = 0.0d0
    endif

  end function HeavySide
  
  function a0(E1,E2) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2 !energies
    
    !output
    real*8 :: answer !answer

    answer = (8.0d0/3.0d0*E1**2*E2**3-4.0d0*E1*E2**4+8.0d0/5.0d0*E2**5)*HeavySide(E1-E2) + &
         4.0d0/15.0d0*E1**5*HeavySide(E2-E1)

  end function a0

  function b0(E1,E2) result(answer)

    implicit none

    !inputs 
    real*8, intent(in) :: E1,E2 !energies
    
    !output
    real*8 :: answer !answer

    answer = (16.0d0/3.0d0*E1*E2**3-4.0d0*E2**4)*HeavySide(E1-E2) + &
         4.0d0/3.0d0*E1**4*HeavySide(E2-E1)

  end function b0

  function c0(E1,E2) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2 !energies
    
    !output
    real*8 :: answer !answer

    answer = 8.0d0/3.0d0*E2**3*HeavySide(E1-E2)+8.0d0/3.0d0*E1**3*HeavySide(E2-E1)

  end function c0

  function Gamma0(E1,E2,E3) result(answer)
    
    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies

    !output
    real*8 :: answer !answer

    answer = 8.0d0/3.0d0*E3**2*(E1**3-E2**3) + &
          4.0d0*E3*(E1-E2)**2*(E1**2/3.0d0+2.0d0*E1*E2/3.0d0+E2**2) + &
          4.0d0*(E1-E2)**3*(E1**2/15.0d0+E1*E2/5.0d0+2.0d0*E2**2/5.0d0)

  end function Gamma0

  function L0I(E1,E2,E3) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies

    !output
    real*8 :: answer !answer

    answer = 4.0d0/15.0d0*E3**5 + 4.0d0/3.0d0*E3**4*E1 + 8.0d0/3.0d0*E3**3*E1**2 + &
          HeavySide(E2-E1)*Gamma0(E1,E2,E3)

  end function L0I

  function L0II(E1,E2,E3) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies

    !output
    real*8 :: answer !answer

    answer = 4.0d0/15.0d0*E3**5 - 4.0d0/3.0d0*E3**4*E2 + 8.0d0/3.0d0*E3**3*E2**2 + &
          HeavySide(E2-E1)*Gamma0(-E2,-E1,E3) !note fix in Bruenn 85 to Ec**4*Eb term and Ec**3*Eb**2 term

  end function L0II

end function NES_H0_ThompsonBruenn
  
function NES_H1_ThompsonBruenn(nu_energy_in,nu_energy_out,electron_energy,neutrino_species) result(H1)
  
  use nulib, only : H0_constants !really H1 constants....
  implicit none
  
  !inputs
  real*8, intent(in) :: nu_energy_in !neutrino energy, MeV, or reduced energy
  real*8, intent(in) :: nu_energy_out !antineutrino energy, MeV, or reduced energy
  real*8, intent(in) :: electron_energy !electron energy, MeV, or reduced energy

  integer, intent(in) :: neutrino_species 

  !output
  real*8 :: H1 !final answer

  !local
  real*8 :: H1I,H1II !J functions
  real*8 :: Ea,Eb,Ec !rename energy variables for ease

  Ea = nu_energy_in
  Eb = nu_energy_out
  Ec = electron_energy

  if (Ec+Ea.lt.Eb) then
     H1 = 0.0d0
     return
  endif

  !Bruenn 1985, Yueh & Buchler 1977 (with corrections from Bruenn 85)
  if (Eb.ge.Ec) then
     H1I = 1.0d0/(Ea*Eb)*((Ea**2+Eb**2)/2.0d0*L0I(Ea,Eb,Ec) - 16.0d0/35.0d0*Ec**7 - &
          4.0d0/5.0d0*(3.0d0*Ea-Eb)*Ec**6 - 2.0d0/15.0d0*Ec**5*(37.0d0*Ea**2-26.0d0*Ea*Eb+Eb**2) - &
          2.0d0/3.0d0*Ec**4*Ea*(Ea-Eb)*(7.0d0*Ea-Eb) - 4.0d0/3.0d0*Ec**3*Ea**2*(Ea-Eb)**2 + &
          Gamma1(Ea,Eb,Ec)*HeavySide(Eb-Ea))

     H1II = 1.0d0/(Ea*Eb)*((Ea**2+Eb**2)/2.0d0*L0II(Ea,Eb,Ec) - 16.0d0/35.0d0*Ec**7 + &
          4.0d0/5.0d0*(3.0d0*Eb-Ea)*Ec**6 - 2.0d0/15.0d0*Ec**5*(37.0d0*Eb**2-26.0d0*Ea*Eb+Ea**2) - &
          2.0d0/3.0d0*Ec**4*Eb*(Ea-Eb)*(7.0d0*Eb-Ea) - 4.0d0/3.0d0*Ec**3*Eb**2*(Eb-Ea)**2 + &
          Gamma1(-Eb,-Ea,Ec)*HeavySide(Eb-Ea)) !error in YU77? did they forget to switch a w-w' to w'-w??

  else
     H1I = a1(Ea,Eb) + b1(Ea,Eb)*Ec + c1(Ea,Eb)*Ec**2

     H1II = a1(Eb,Ea) - b1(Eb,Ea)*Ec + c1(Eb,Ea)*Ec**2 !invert Ea and Eb, flip sign on b0 term

  endif

  H1 = H0_constants(1,neutrino_species)*H1I + H0_constants(2,neutrino_species)*H1II !same coefficients as e-e+ annhilation

contains
  function HeavySide(arguement) result(answer)

    implicit none
    
    !inputs
    real*8, intent(in) :: arguement

    !output
    real*8 :: answer

    if (arguement.gt.0.0d0) then
       answer = 1.0d0
    else if (arguement.eq.0.0d0) then
       answer = 0.5d0
    else
       answer = 0.0d0
    endif

  end function HeavySide
  
  function a1(E1,E2) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2 !energies
    
    !output
    real*8 :: answer !answer

    answer = (-4.0d0/3.0d0*E1**2*E2**3+16.0d0/5.0d0*E1*E2**4-16.0d0/5.0d0*E2**5+ &
         8.0d0/7.0d0*E2**6/E1)*HeavySide(E1-E2) + &
         (12.0d0/35.0d0*E1**6/E2-8.0d0/15.0d0*E1**5)*HeavySide(E2-E1)

  end function a1

  function b1(E1,E2) result(answer)

    implicit none

    !inputs 
    real*8, intent(in) :: E1,E2 !energies
    
    !output
    real*8 :: answer !answer

    answer = (-8.0d0/3.0d0*E1*E2**3+28.0d0/5.0d0*E2**4-16.0d0/5.0d0*E2**5/E1)* &
         HeavySide(E1-E2) + &
         (8.0d0/5.0d0*E1**5/E2-28.0d0/15.0d0*E1**4)*HeavySide(E2-E1) !corrected from YU77 to Bruenn 85

  end function b1

  function c1(E1,E2) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2 !energies
    
    !output
    real*8 :: answer !answer

    answer = (12.0d0/5.0d0*E2**4/E1-4.0d0/3.0d0*E2**3)*HeavySide(E1-E2) + &
         (12.0d0/5.0d0*E1**4/E2-4.0d0/3.0d0*E1**3)*HeavySide(E2-E1)

  end function c1

  function Gamma0(E1,E2,E3) result(answer)
    
    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies

    !output
    real*8 :: answer !answer

    answer = 8.0d0/3.0d0*E3**2*(E1**3-E2**3) + &
          4.0d0*E3*(E1-E2)**2*(E1**2/3.0d0+2.0d0*E1*E2/3.0d0+E2**2) + &
          4.0d0*(E1-E2)**3*(E1**2/15.0d0+E1*E2/5.0d0+2.0d0*E2**2/5.0d0)

  end function Gamma0

  function Gamma1(E1,E2,E3) result(answer)
    
    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies

    !output
    real*8 :: answer !answer

    answer = 4.0d0/15.0d0*E3**2*(E1-E2)**3*(4.0d0*E1**2+7.0d0*E1*E2+4.0d0*E2**2) + &
         2.0d0/15.0d0*E3*(E1-E2)**4*(7.0d0*E1**2+14.0d0*E1*E2+9.0d0*E2**2) + &
         2.0d0/105.0d0*(E1-E2)**5*(11.0d0*E1**2+27.0d0*E1*E2+18.0d0*E2**2)

  end function Gamma1

  function L0I(E1,E2,E3) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies

    !output
    real*8 :: answer !answer

    answer = 4.0d0/15.0d0*E3**5 + 4.0d0/3.0d0*E3**4*E1 + 8.0d0/3.0d0*E3**3*E1**2 + &
          HeavySide(E2-E1)*Gamma0(E1,E2,E3)

  end function L0I

  function L0II(E1,E2,E3) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies

    !output
    real*8 :: answer !answer

    answer = 4.0d0/15.0d0*E3**5 - 4.0d0/3.0d0*E3**4*E2 + 8.0d0/3.0d0*E3**3*E2**2 + &
          HeavySide(E2-E1)*Gamma0(-E2,-E1,E3) !note fix in Bruenn 85 to Ec**4*Eb term and Ec**3*Eb**2 term

  end function L0II

end function NES_H1_ThompsonBruenn


  
