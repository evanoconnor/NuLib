!-*-f90-*-
function weak_mag_correction_absorption(nu_energy,type) result(correction)

  use nulib
  implicit none
  real*8 :: correction !answer, dimensionless
  real*8, intent(in)  :: nu_energy !in MeV
  integer, intent(in) :: type !0 for neutrino absorption on neutrons, 1 for antineutrino on proton

  real*8 :: c_v,c_a,F2,e

  !from Horowitz Phys. Rev. D. 65 043001, 2002, critical for proper behaviour at large antineutrino energies
  c_v = 1.0d0
  c_a = gA
  F2 = 3.706d0

  if (type.eq.0) then !neutrino capture crosssection on neutron
     e = nu_energy/m_n
     correction = (c_v**2*(1.0d0+4.0d0*e+16.0d0/3.0d0*e**2)+ &
          3.0d0*c_a**2*(1.0d0+4.0d0/3.0d0*e)**2+4.0d0*(c_v+F2)*c_a*e* &
          (1.0d0+4.0d0/3.0d0*e)+8.0d0/3.0d0*c_v*F2*e**2+&
          5.0d0/3.0d0*e**2*(1.0d0+2.0d0/5.0d0*e)*F2**2)/ &
          ((c_v**2+3.0d0*c_a**2)*(1.0d0+2.0d0*e)**3)
  else if (type.eq.1) then !antineutrino capture crosssection on proton
     e = nu_energy/m_p
     correction = (c_v**2*(1.0d0+4.0d0*e+16.0d0/3.0d0*e**2)+ &
          3.0d0*c_a**2*(1.0d0+4.0d0/3.0d0*e)**2-4.0d0*(c_v+F2)*c_a*e* &
          (1.0d0+4.0d0/3.0d0*e)+8.0d0/3.0d0*c_v*F2*e**2+&
          5.0d0/3.0d0*e**2*(1.0d0+2.0d0/5.0d0*e)*F2**2)/ &
          ((c_v**2+3.0d0*c_a**2)*(1.0d0+2.0d0*e)**3)
  else
     stop "weak_mag_correction_absorption: Wrong type"
  endif

end function weak_mag_correction_absorption

function weak_mag_correction_scattering_differential(nu_energy,x,type,lepton) result(correction)
  
  use nulib
  implicit none

  real*8 :: correction !answer, dimensionless
  real*8, intent(in)  :: nu_energy !in MeV
  real*8, intent(in)  :: x !angle, cos(theta)
  integer, intent(in) :: type !0 for (anti)neutrino scattering on neutrons, 1 for (anti)neutrino scattering on proton
  integer, intent(in) :: lepton !-1 for antineutrino, 1 for neutrino

  real*8 :: c_v,c_a,F2,sign !constants for correction term
  real*8 :: e,eta,eta2 !reduced energies

  if (.not.(lepton.eq.-1.or.lepton.eq.1)) then
     stop "weak_mag_correction_scattering: Wrong lepton number, 1 for neutrinos, -1 for antineutrinos"
  endif

  if (type.eq.0) then 
     !scattering on neutrons
     c_v = -0.5d0
     c_a = -gA/2.0d0
     F2 = -0.963d0
     e = nu_energy/m_n
     eta = 1.0d0/(1.0d0+e*(1.0d0-x))
     eta2 = eta*eta
     sign = real(lepton)
     
  else if (type.eq.1) then
     !scattering on protons
     c_v = 0.5d0-2.0d0*sin2thetaW
     c_a = gA/2.0d0
     F2 = 1.019d0
     e = nu_energy/m_p
     eta = 1.0d0/(1.0d0+e*(1.0d0-x))
     eta2 = eta*eta
     sign = real(lepton)

  else 
     stop "wrong type"
  endif

  !lepton causes sign change for antineutrino (lepton=-1) and neutrino (lepton=1)
  correction = (eta2*(c_v**2*(1.0d0+eta2-eta*(1.0d0-x)) + &
       c_a**2*(1.0d0+eta2+eta*(1.0d0-x)) + &
       sign*2.0d0*c_a*(c_v+F2)*(1.0d0-eta**2) + &
       eta2*e*e/2.0d0*(1.0d0-x)*(F2**2*(3.0d0-x)+4.0d0*c_v*F2*(1.0d0-x)))) / &
       (c_v**2*(1.0d0+x)+c_a**2*(3.0d0-x)) 

end function weak_mag_correction_scattering_differential

function weak_mag_correction_scattering_transport(nu_energy,type,lepton) result(correction)

  use nulib
  implicit none

  real*8 :: correction !answer, dimensionless
  real*8, intent(in)  :: nu_energy !in MeV
  integer, intent(in) :: type !0 for (anti)neutrino scattering on neutrons, 1 for (anti)neutrino scattering on proton
  integer, intent(in) :: lepton !-1 for antineutrino, 1 for neutrino

  real*8 :: c_v,c_a,F2,sign !constants for correction term
  real*8 :: e,e2,e3 !reduced energies

  if (.not.(lepton.eq.-1.or.lepton.eq.1)) then
     stop "weak_mag_correction_scattering: Wrong lepton number, 1 for neutrinos, -1 for antineutrinos"
  endif

  if (type.eq.0) then
     !scattering on neutrons
     c_v = -0.5d0
     c_a = -gA/2.0d0
     F2 = -0.963d0
     e = nu_energy/m_n
     e2 = e*e
     e3 = e2*e
     sign = real(lepton)
     !lepton causes sign change for antineutrino (lepton=-1) and neutrino (lepton=1)
     correction = ( &
          c_v**2*((e-1.0d0)/(2.0d0*e3)*log(1.0d0+2.0d0*e)+ & 
          (3.0d0+12.0d0*e+9.0*e2-10.0d0*e3)/(3.0d0*e2*(1.0d0+2.0d0*e)**3)) + &
          c_a**2*((1.0d0+e)/(2.0d0*e3)*log(1.0d0+2.0d0*e)- &
          (10.d0*e3+27.0d0*e2+18.0d0*e+3.0d0)/(3.0d0*e2*(1.0d0+2.0d0*e)**3)) + &
          sign*(c_v+F2)*c_a*(log(1.0d0+2.0d0*e)/e2- &
          (2.0d0+10.d0*e+28.0d0*e2/3.0d0)/(e*(1.0d0+2.0d0*e)**3)) + &
          c_v*F2*(log(1.0d0+2.0d0*e)/e2-2.0d0* &
          (3.0d0+15.0d0*e+22.0d0*e2)/(3.0d0*e*(1.0d0+2.0d0*e)**3)) + &
          F2**2*(log(1.0d0+2.0d0*e)/(4.0d0*e2)+ &
          (8.0d0*e3-22.0d0*e2-15.0d0*e-3.0d0)/(6.0d0*e*(1.0d0+2.0d0*e)**3)))/ &
          (2.0d0*(c_v**2+5.0d0*c_a**2)/3.0d0)
  else if (type.eq.1) then
     !scattering on protons
     c_v = 0.5d0-2.0d0*sin2thetaW
     c_a = gA/2.0d0
     F2 = 1.019d0
     e = nu_energy/m_p
     e2 = e*e
     e3 = e2*e
     sign = real(lepton)
     !lepton causes sign change for antineutrino (lepton=-1) and neutrino (lepton=1)
     correction = ( &
          c_v**2*((e-1.0d0)/(2.0d0*e3)*log(1.0d0+2.0d0*e)+ & 
          (3.0d0+12.0d0*e+9.0*e2-10.0d0*e3)/(3.0d0*e2*(1.0d0+2.0d0*e)**3)) + &
          c_a**2*((1.0d0+e)/(2.0d0*e3)*log(1.0d0+2.0d0*e)- &
          (10.d0*e3+27.0d0*e2+18.0d0*e+3.0d0)/(3.0d0*e2*(1.0d0+2.0d0*e)**3)) + &
          sign*(c_v+F2)*c_a*(log(1.0d0+2.0d0*e)/e2- &
          (2.0d0+10.d0*e+28.0d0*e2/3.0d0)/(e*(1.0d0+2.0d0*e)**3)) + &
          c_v*F2*(log(1.0d0+2.0d0*e)/e2-2.0d0* &
          (3.0d0+15.0d0*e+22.0d0*e2)/(3.0d0*e*(1.0d0+2.0d0*e)**3)) + &
          F2**2*(log(1.0d0+2.0d0*e)/(4.0d0*e2)+ &
          (8.0d0*e3-22.0d0*e2-15.0d0*e-3.0d0)/(6.0d0*e*(1.0d0+2.0d0*e)**3)))/ &
          (2.0d0*(c_v**2+5.0d0*c_a**2)/3.0d0)

  else
     stop "weak_mag_correction_scattering: Wrong type"
  endif

end function weak_mag_correction_scattering_transport

function weak_mag_correction_scattering_total(nu_energy,type,lepton) result(correction)

  use nulib
  implicit none

  real*8 :: correction !answer, dimensionless
  real*8, intent(in)  :: nu_energy !in MeV
  integer, intent(in) :: type !0 for (anti)neutrino scattering on neutrons, 1 for (anti)neutrino scattering on proton
  integer, intent(in) :: lepton !-1 for antineutrino, 1 for neutrino

  real*8 :: c_v,c_a,F2,sign !constants for correction term
  real*8 :: e,e2,e3 !reduced energies

  if (.not.(lepton.eq.-1.or.lepton.eq.1)) then
     stop "weak_mag_correction_scattering: Wrong lepton number, 1 for neutrinos, -1 for antineutrinos"
  endif

  if (type.eq.0) then
     !scattering on neutrons
     c_v = -0.5d0
     c_a = -gA/2.0d0
     F2 = -0.963d0
     e = nu_energy/m_n
     e2 = e*e
     sign = real(lepton)

  else if (type.eq.1) then
     !scattering on protons
     c_v = 0.5d0-2.0d0*sin2thetaW
     c_a = gA/2.0d0
     F2 = 1.019d0
     e = nu_energy/m_p
     e2 = e*e
     sign = real(lepton)
  else
     stop "weak_mag_correction_scattering: Wrong type"
  endif

  !lepton causes sign change for antineutrino (lepton=-1) and neutrino (lepton=1)
  correction = (c_v**2*(1.0d0+4.0d0*e+16.0d0/3.0d0*e2)+ &
       3.0d0*c_a**2*(1.0d0+4.0d0/3.0d0*e)**2+sign*4.0d0*(c_v+F2)*c_a*e* &
       (1.0d0+4.0d0/3.0d0*e)+8.0d0/3.0d0*c_v*F2*e2+&
       5.0d0/3.0d0*e2*(1.0d0+2.0d0/5.0d0*e)*F2**2)/ &
       ((c_v**2+3.0d0*c_a**2)*(1.0d0+2.0d0*e)**3)

end function weak_mag_correction_scattering_total

