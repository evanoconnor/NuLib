!-*-f90-*-
function epannihil_dQdenu_BRT06(nu_energy_x,matter_eta,neutrino_species) result(dQdenu)

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
  real*8 :: epannihil_Phi_Bruenn !function declaration
  
  dQdenu = 0.0d0
  do i=1,16
     nubar_energy_x = GLQ_n16_roots(i)
     dQdenu = dQdenu + &
          nubar_energy_x**2*epannihil_Phi_Bruenn(nu_energy_x,nubar_energy_x, &
          matter_eta,neutrino_species,1,0)*GLQ_n16_weights(i)
  enddo

end function epannihil_dQdenu_BRT06

function epannihil_Phi_Bruenn(nu_energy_x,nubar_energy_x,matter_eta,neutrino_species,pro_or_ann,which_l) result(Phi)

  use nulib, only : GPQ_n16_roots, GPQ_n16_weights
  implicit none

  !inputs
  real*8, intent(in) :: nu_energy_x !neutrino energy/Temp
  real*8, intent(in) :: nubar_energy_x !antineutirno energy/Temp
  real*8, intent(in) :: matter_eta !electron chemical potential/Temp

  integer, intent(in) :: neutrino_species !neutrino species, will affect some coeffs
  integer, intent(in) :: pro_or_ann !to treat electron and positrons
                                    !as initial or final state. 1 is
                                    !production, i.e. we care about
                                    !how many electrons/postirons
                                    !there are. 2 is annihilation,
                                    !i.e. we care about how much space
                                    !is availabe
  integer, intent(in) :: which_l !i.e. 0 or 1

  !output
  real*8 :: Phi !dimensionless

  !local, GPQ integration
  integer :: i
  real*8 :: electron_x,positron_x,femtimesfep !the underscore x denotes /Temp

  !function declarations
  real*8 :: epannihil_H0_Bruenn !function declaration
  real*8 :: epannihil_H1_Bruenn !function declaration
  real*8 :: fermidirac_dimensionless !function declaration

  Phi = 0.0d0

  do i=1,16
     !must take interval change into account
     electron_x = (nu_energy_x+nubar_energy_x)/2.0d0*(GPQ_n16_roots(i)+1.0d0) !MeV
     positron_x = max(nu_energy_x+nubar_energy_x-electron_x,1.0d-20) !MeV
     if (pro_or_ann.eq.1) then
        if (electron_x-matter_eta.gt.35.0d0) then
           !1/(1+e) reduces to e^{-1} since e is big 
           if (positron_x+matter_eta.gt.35.0d0) then
              !1/(1+p) reduces to p^{-1} since p is big 
              femtimesfep = exp(-electron_x-positron_x)
           else
              !1/(1+p) cannot be reduced
              femtimesfep = exp(-electron_x+matter_eta)*fermidirac_dimensionless(positron_x,-matter_eta)
           endif
        else if (electron_x-matter_eta.lt.-35.0d0) then
           !1/(1+e) reduces to 1 since e is small
           if (positron_x+matter_eta.gt.35.0d0) then
              !1/(1+p) reduces to p^{-1} since p is big 
              femtimesfep = exp(-positron_x-matter_eta)
           else
              !1/(1+p) cannot be reduced
              femtimesfep = fermidirac_dimensionless(positron_x,-matter_eta)
           endif           
        else
           !1/(1+e) cannot be reduced
           if (positron_x+matter_eta.gt.35.0d0) then
              !1/(1+p) reduces to p^{-1} since p is big 
              femtimesfep = fermidirac_dimensionless(electron_x,matter_eta)*exp(-positron_x-matter_eta)
           else
              !1/(1+p) cannot be reduced
              femtimesfep = fermidirac_dimensionless(electron_x,matter_eta)* &
                   fermidirac_dimensionless(positron_x,-matter_eta)
           endif                      
        endif
     else if (pro_or_ann.eq.2) then
        if (electron_x-matter_eta.gt.35.0d0) then
           !1-1/(1+e) is e/(1+e) which reduces to 1 since e is big 
           if (positron_x+matter_eta.gt.35.0d0) then
              !1-1/(1+p) is p/(1+p) which reduces to 1 since p is big 
              femtimesfep = 1.0d0
           else
              !1-1/(1+p) is p/(1+p) which cannot be reduced
              femtimesfep = exp(positron_x+matter_eta)*fermidirac_dimensionless(positron_x,-matter_eta)
           endif
        else if (electron_x-matter_eta.lt.-35.0d0) then
           !1-1/(1+e) is e/(1+e) which reduces to e since e is small
           if (positron_x+matter_eta.gt.35.0d0) then
              !1-1/(1+p) is p/(1+p) which reduces to 1 since p is big 
              femtimesfep = exp(electron_x-matter_eta)
           else
              !1-1/(1+p) is p/(1+p) which cannot be reduced
              femtimesfep = exp(electron_x+positron_x)* &
                   fermidirac_dimensionless(positron_x,-matter_eta)
           endif           
        else
           !1-1/(1+e) is e/(1+e) which cannot be reduced 
           if (positron_x+matter_eta.gt.35.0d0) then
              !1-1/(1+p) is p/(1+p) which reduces to 1 since p is big 
              femtimesfep = fermidirac_dimensionless(electron_x,matter_eta)*exp(electron_x-matter_eta)
           else
              !1-1/(1+p) is p/(1+p) which cannot be reduced
              femtimesfep = exp(electron_x+positron_x)* &
                   fermidirac_dimensionless(electron_x,matter_eta)* &
                   fermidirac_dimensionless(positron_x,-matter_eta)
           endif                      
        endif

        if (femtimesfep.gt.1.0d0) then
           write(*,*) electron_x,positron_x,matter_eta,pro_or_ann,femtimesfep
           stop "femtimesfep terms are greater than 1, this is bad"
           
        endif
        if (femtimesfep.ne.femtimesfep) then
           stop "femtimesfep is NaN, this is really bad"
        endif

     else
        stop "correctly request whether you want production or annihilation"
     endif
     if (which_l.eq.0) then
        Phi = Phi + &
             femtimesfep*epannihil_H0_Bruenn(nu_energy_x,nubar_energy_x,electron_x,neutrino_species)* &
             GPQ_n16_weights(i)
     else if (which_l.eq.1) then
        Phi = Phi + &
             femtimesfep*epannihil_H1_Bruenn(nu_energy_x,nubar_energy_x,electron_x,neutrino_species)* &
             GPQ_n16_weights(i)
     else
        stop "requesting a legendre mode that is not coded"
     endif
  enddo
  
  !interval change
  Phi = Phi*(nu_energy_x+nubar_energy_x)/2.0d0

end function epannihil_Phi_Bruenn

function epannihil_H0_Bruenn(nu_energy,nubar_energy,electron_energy,neutrino_species) result(H0)
  
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

end function epannihil_H0_Bruenn

function epannihil_H1_Bruenn(nu_energy,nubar_energy,electron_energy,neutrino_species) result(H1)
  
  use nulib, only : H0_constants
  implicit none
  
  !inputs
  real*8, intent(in) :: nu_energy !neutrino energy, MeV, or reduced energy
  real*8, intent(in) :: nubar_energy !antineutrino energy, MeV, or reduced energy
  real*8, intent(in) :: electron_energy !electron energy, MeV, or reduced energy

  integer, intent(in) :: neutrino_species 

  !output
  real*8 :: H1 !final answer

  !local
  real*8 :: J1I,J1II !J functions
  real*8 :: Ea,Eb,Ec !rename energy variables for ease

  Ea = nu_energy
  Eb = nubar_energy
  Ec = electron_energy

  !Bruenn 1985
  J1I = HeavySide(Ea+Eb-Ec)*( &
       (HeavySide(Ea-Ec)*HeavySide(Eb-Ea)+HeavySide(Eb-Ec)*HeavySide(Ea-Eb))*a1(Ea,Eb,Ec) + &
       (HeavySide(Ec-Ea)*HeavySide(Ea-Eb)+HeavySide(Ec-Eb)*HeavySide(Eb-Ea))*b1(Ea,Eb,Ec) + &
       HeavySide(Ec-Ea)*HeavySide(Eb-Ec)*c1(Ea,Eb,Ec) + &
       HeavySide(Ec-Eb)*HeavySide(Ea-Ec)*d1(Ea,Eb,Ec) &
       )/(Ea*Eb)

  Ea = nubar_energy
  Eb = nu_energy
  Ec = electron_energy

  !Bruenn 1985
  J1II = HeavySide(Ea+Eb-Ec)*( &
       (HeavySide(Ea-Ec)*HeavySide(Eb-Ea)+HeavySide(Eb-Ec)*HeavySide(Ea-Eb))*a1(Ea,Eb,Ec) + &
       (HeavySide(Ec-Ea)*HeavySide(Ea-Eb)+HeavySide(Ec-Eb)*HeavySide(Eb-Ea))*b1(Ea,Eb,Ec) + &
       HeavySide(Ec-Ea)*HeavySide(Eb-Ec)*c1(Ea,Eb,Ec) + &
       HeavySide(Ec-Eb)*HeavySide(Ea-Ec)*d1(Ea,Eb,Ec) &
       )/(Ea*Eb)
  

  H1 = H0_constants(1,neutrino_species)*J1I + H0_constants(2,neutrino_species)*J1II

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
  
  function a1(E1,E2,E3) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies
    
    !output
    real*8 :: answer !answer

    answer = (16.0d0/35.0d0*E3**7 - 4.0d0/5.0d0*E3**6*(E1+3.0d0*E2) + 4.0d0/15.0d0*E3**5*E2*(13.0d0*E1+18.0d0*E2) - &
         4.0d0/3.0d0*E3**4*E2**2*(4.0d0*E1+3.0d0*E2) + 8.0d0/3.0d0*E3**3*E1*E2**3 )/(E1**2*E2**2)

  end function a1

  !modified to correct bug in Bruenn 1985, units don't work for b0, multiply the a0 in this term by E1**2*E2**2, check via mathematica
  function b1(E1,E2,E3) result(answer)

    implicit none

    !inputs 
    real*8, intent(in) :: E1,E2,E3 !energies
    
    !output
    real*8 :: answer !answer

    answer = ( -a1(E1,E2,E3)*E1**2*E2**2 - E3**2*(12.0d0/5.0d0*E1**5 + 4.0d0/3.0d0*E1**4*E2+4.0d0/3.0d0*E1*E2**4 + &
         12.0d0/5.0d0*E2**5) + E3*(16.0d0/5.0d0*E1**6 + 28.0d0/5.0d0*E1**5*E2 + 8.0d0/3.0d0*E1**4*E2**2 + &
         28.0d0/15.0d0*E1*E2**5 + 8.0d0/5.0d0*E2**6) - (8.0d0/7.0d0*E1**7 + 16.0d0/5.0d0*E1**6*E2 + &
         16.0d0/5.0d0*E1**5*E2**2 + 4.0d0/3.0d0*E1**4*E2**3 + 8.0d0/15.0d0*E1*E2**6 + 12.0d0/35.0d0*E2**7))/(E1**2*E2**2) 

  end function b1

  function c1(E1,E2,E3) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies
    
    !output
    real*8 :: answer !answer

    answer = -E1**2*(8.0d0/7.0d0*E1**3 + 16.0d0/5.0d0*E1**2*E2 + 16.0d0/5.0d0*E1*E2**2 + 4.0d0/3.0d0*E2**3)/E2**2 + &
         E3*E1**2*(16.0d0/5.0d0*E1**2 + 28.0d0/5.0d0*E1*E2 + 8.0d0/3.0d0*E2**2)/E2**2 - &
         E3**2*E1**2*(12.0d0/5.0d0*E1+4.0d0/3.0d0*E2)/E2**2

  end function c1

  function d1(E1,E2,E3) result(answer)

    implicit none

    !inputs
    real*8, intent(in) :: E1,E2,E3 !energies
    
    !output
    real*8 :: answer !answer

    answer = -E2**4*(8.0d0/15.0d0*E1 + 12.0d0/35.0d0*E2)/E1**2 + &
         E3*E2**3*(28.0d0/15.0d0*E1 + 8.0d0/5.0d0*E2)/E1**2 - &
         E3**2*E2**2*(4.0d0/3.0d0*E1+12.0d0/5.0d0*E2)/E1**2

  end function d1

end function epannihil_H1_Bruenn

