!-*-f90-*-
function epannhil_dQdenu_BRT06(nu_energy_x,matter_eta,neutrino_species) result(dQdenu)

  use nulib, only : GLQ_n16_roots, GLQ_n16_weights
  implicit none

  !inputs
  real*8, intent(in) :: nu_energy_x !neutrino energy/Temp
  real*8, intent(in) :: matter_eta !electron chemical potential/Temp

  integer, intent(in) :: neutrino_species !neutrino species, will affect some coeffs

  !output
  real*8 :: dQdenu

  !local, GLQ integration
  integer :: i  
  real*8 :: nubar_energy_x

  !function declarations
  real*8 :: epannhil_Phi0_BRT06 !function declaration
  
  dQdenu = 0.0d0
  do i=1,16
     nubar_energy_x = GLQ_n16_roots(i)
     dQdenu = dQdenu + &
          nubar_energy_x**2*epannhil_Phi0_BRT06(nu_energy_x,nubar_energy_x, &
          matter_eta,neutrino_species)*GLQ_n16_weights(i)
  enddo

end function epannhil_dQdenu_BRT06

function epannhil_Phi0_BRT06(nu_energy_x,nubar_energy_x,matter_eta,neutrino_species) result(Phi)

  use nulib, only : GPQ_n16_roots, GPQ_n16_weights
  implicit none

  !inputs
  real*8, intent(in) :: nu_energy_x !neutrino energy/Temp
  real*8, intent(in) :: nubar_energy_x !antineutirno energy/Temp
  real*8, intent(in) :: matter_eta !electron chemical potential/Temp

  integer, intent(in) :: neutrino_species !neutrino species, will affect some coeffs

  !output
  real*8 :: Phi !dimensionless

  !local, GPQ integration
  integer :: i
  real*8 :: electron_x,positron_x,fep,fem !the underscore x denotes /Temp

  !function declarations
  real*8 :: epannhil_H0_BRT06 !function declaration
  real*8 :: fermidirac_dimensionless !function declaration

  Phi = 0.0d0

  do i=1,16
     !must take interval change into account
     electron_x = (nu_energy_x+nubar_energy_x)/2.0d0*(GPQ_n16_roots(i)+1.0d0) !MeV
     positron_x = max(nu_energy_x+nubar_energy_x-electron_x,1.0d-20) !MeV
     fem = fermidirac_dimensionless(electron_x,matter_eta)
     fep = fermidirac_dimensionless(positron_x,-matter_eta)
     Phi = Phi + &
          fep*fem*epannhil_H0_BRT06(nu_energy_x,nubar_energy_x,electron_x,neutrino_species)* &
          GPQ_n16_weights(i)
  enddo
  
  !interval change
  Phi = Phi*(nu_energy_x+nubar_energy_x)/2.0d0

end function epannhil_Phi0_BRT06

function epannhil_H0_BRT06(nu_energy,nubar_energy,electron_energy,neutrino_species) result(H0)
  
  use nulib, only : H0_constants
  implicit none
  
  !inputs
  real*8, intent(in) :: nu_energy !neutrino energy, MeV, or reduced energy
  real*8, intent(in) :: nubar_energy !antineutrino energy, MeV, or reduced energy
  real*8, intent(in) :: electron_energy !electron energy, MeV, or reduced energy

  integer, intent(in) :: neutrino_species 

  !output
  real*8 :: H0 !final answer

  !local
  real*8 :: J0I,J0II !J functions
  real*8 :: Ea,Eb,Ec !rename energy variables for ease

  Ea = nu_energy
  Eb = nubar_energy
  Ec = electron_energy

  !Bruenn 1985
  J0I = HeavySide(Ea+Eb-Ec)*( &
       (HeavySide(Ea-Ec)*HeavySide(Eb-Ea)+HeavySide(Eb-Ec)*HeavySide(Ea-Eb))*a0(Ea,Eb,Ec) + &
       (HeavySide(Ec-Ea)*HeavySide(Ea-Eb)+HeavySide(Ec-Eb)*HeavySide(Eb-Ea))*b0(Ea,Eb,Ec) + &
       HeavySide(Ec-Ea)*HeavySide(Eb-Ec)*c0(Ea,Eb,Ec) + &
       HeavySide(Ec-Eb)*HeavySide(Ea-Ec)*d0(Ea,Eb,Ec) &
       )/(Ea*Eb)

  Ea = nubar_energy
  Eb = nu_energy
  Ec = electron_energy

  !Bruenn 1985
  J0II = HeavySide(Ea+Eb-Ec)*( &
       (HeavySide(Ea-Ec)*HeavySide(Eb-Ea)+HeavySide(Eb-Ec)*HeavySide(Ea-Eb))*a0(Ea,Eb,Ec) + &
       (HeavySide(Ec-Ea)*HeavySide(Ea-Eb)+HeavySide(Ec-Eb)*HeavySide(Eb-Ea))*b0(Ea,Eb,Ec) + &
       HeavySide(Ec-Ea)*HeavySide(Eb-Ec)*c0(Ea,Eb,Ec) + &
       HeavySide(Ec-Eb)*HeavySide(Ea-Ec)*d0(Ea,Eb,Ec) &
       )/(Ea*Eb)
  

  H0 = H0_constants(1,neutrino_species)*J0I + H0_constants(2,neutrino_species)*J0II

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
  
  function a0(E1,E2,E3) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies
    
    !output
    real*8 :: answer !answer

    answer = (4.0d0*E3**5-20.0d0*E3**4*E2+40.0d0*E3**3*E2**2)/(15.0d0*E1*E2)

  end function a0

  !modified to correct bug in Bruenn 1985, units don't work for b0, multiply the a0 in this term by E1*E2, check via mathematica
  function b0(E1,E2,E3) result(answer)

    implicit none

    !inputs 
    real*8, intent(in) :: E1,E2,E3 !energies
    
    !output
    real*8 :: answer !answer

    answer = ( &
         -15.0d0*a0(E1,E2,E3)*E1*E2+ &
         40.0d0*E3**2*(E1**3+E2**3)- &
         20.0d0*E3*(E1+E2)**2*(E2**2-2.0d0*E1*E2+3.0d0*E1**2)+ &
         4.0d0*(E1+E2)**3*(E2**2-3.0d0*E1*E2+6.0d0*E1**2) &
         )/(15.0d0*E1*E2) 

  end function b0

  function c0(E1,E2,E3) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies
    
    !output
    real*8 :: answer !answer

    answer = E1**2*(40.0d0*E2**2+60.0d0*E1*E2+24.0d0*E1**2)/(15.0d0*E2) - &
         (E3*(16.0d0*E2*E1**2+12.0d0*E1**3) - 8.0d0*E3**2*E1**2)/(3.0d0*E2)

  end function c0

  function d0(E1,E2,E3) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies
    
    !output
    real*8 :: answer !answer

    answer = (4.0d0*E2**4-20.0d0*E2**3*E3+40.0d0*E2**2*E3**2)/(15.0d0*E1)

  end function d0

end function epannhil_H0_BRT06


  
