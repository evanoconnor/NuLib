!-*-f90-*-
!always keep term dimensionless unless it is the units of the final answer...
!total and transport scattering crosssections
subroutine nu_scatter_elastic_e_total(neutrino_energy,type,lepton,eos_variables,crosssection)

  use nulib
  implicit none

  !Thompson Thesis + Bower's & Wilson 1982 (ApJS) here we use the
  !interpolated cross section between the degenerate and non
  !degenerate limits. There is actually an energy loss with this
  !crosssection, that will come eventually

  !inputs
  real*8, intent(in)  :: eos_variables(total_eos_variables)
  real*8, intent(in) :: neutrino_energy !in MeV
  integer, intent(in) :: type !transport (type.eq1) or absorption (type.eq.0) crosssection
  integer, intent(in) :: lepton !lepton number, 1 for neutrino,  -1 for antineutrino

  !output
  real*8, intent(out) :: crosssection !in cm^2

  !locals
  real*8 :: c_v_prime
  real*8 :: c_a_prime
  
  if (lepton.eq.1) then
     c_v_prime = 0.5d0+2.0d0*sin2thetaW
     c_a_prime = 0.5d0
  else if(lepton.eq.-1) then
     c_v_prime = 0.5d0+2.0d0*sin2thetaW
     c_a_prime = -0.5d0
  else if(lepton.ge.2) then
     c_v_prime = -0.5d0+2.0d0*sin2thetaW
     c_a_prime = -0.5d0
  else if(lepton.le.-2) then
     c_v_prime = -0.5d0+2.0d0*sin2thetaW
     c_a_prime = 0.5d0
  else
     stop "nu_scatter_elastic_e_total: lepton number weird...."
  endif

  if (type.eq.1) then
     crosssection = 3.0d0/8.0d0*sigma0/m_e**2*neutrino_energy* &
          (eos_variables(tempindex)+eos_variables(mueindex)/4.0d0)* &
          ((c_v_prime+c_a_prime)**2+(c_v_prime-c_a_prime)**2/3.0d0)
  else
     stop "no transport crosssection?"
  endif


end subroutine nu_scatter_elastic_e_total

subroutine nu_scatter_elastic_alpha_total(neutrino_energy,type,lepton,crosssection)

  use nulib
  implicit none
  
  !inputs
  real*8, intent(in) :: neutrino_energy !in MeV
  integer, intent(in) :: type !transport (type.eq1) or absorption (type.eq.0) crosssection
  integer, intent(in) :: lepton !lepton number, >0 for neutrino,  <0 for antineutrino

  !output
  real*8,intent(out) :: crosssection !in cm^2

  crosssection = 4.0d0*sigma0*(neutrino_energy/m_e)**2*sin2thetaW**2

  !transport crosssection
  if (type.eq.1) then
     crosssection = 2.0d0/3.0d0*crosssection
  end if

end subroutine nu_scatter_elastic_alpha_total

subroutine nu_scatter_elastic_p_total(neutrino_energy,type,lepton,eos_variables,crosssection)

  use nulib
  implicit none
  
  !inputs
  real*8, intent(in) :: neutrino_energy !in MeV
  real*8, intent(in)  :: eos_variables(total_eos_variables)
  integer, intent(in) :: type !transport (type.eq1) or absorption (type.eq.0) crosssection
  integer, intent(in) :: lepton !lepton number, >0 for neutrino,  <0 for antineutrino

  !output
  real*8 :: crosssection !in cm^2
  
  !local
  real*8 :: c_v_prime
  real*8 :: c_a_prime
  real*8 :: delta_p
  real*8 :: weak_mag
  integer :: lepton_local
  
  !strange virial correction via S_A
  real*8 :: gAeff,S_A,S_V
  real*8 :: virialA,virialB,virialC
  real*8 :: virial_nb,virial_ye,virial_temp

  !function declaration
  real*8 :: weak_mag_correction_scattering_transport
  real*8 :: weak_mag_correction_scattering_total

  !for scattering on p, only need sign of lepton number
  if (lepton.gt.0) then
     lepton_local = 1
  else
     lepton_local = -1
  endif

  !a la Horowitz et al. (2016)
  !note, this includes a BS (1998) estimate at high densities
  if (do_nc_virial_correction) then
     virial_nb = eos_variables(rhoindex)*1.0d-39/(m_ref*mev_to_gram)
     virial_ye = eos_variables(yeindex)
     virial_temp = eos_variables(tempindex)
     
     virialA = 920.0d0*virial_nb*(1.0d0-virial_ye+virial_ye**2)/virial_temp**1.22d0
     virialB = 3.05d0/virial_temp**0.75d0
     virialC = 6140.0d0*virial_nb*virial_ye*(1.0d0-virial_ye)/virial_temp**0.5d0 + &
          1.5d13*virial_nb**4/virial_temp**6
     S_A = 1.0d0/(1.0d0+virialA*(1.d0+virialB*exp(-virialC)))
     S_V = 1.0d0

  else
     S_A = 1.0d0
     S_V = 1.0d0
  endif

  !Horowitz (2002)
  if (do_strange_coupling) then
     gAeff = gA-gAs ! - because protons
  else
     gAeff = gA
  endif

  c_v_prime = 0.5d0+2.0d0*sin2thetaW
  c_a_prime = 0.5d0

  !BRT06 Eq.15
  crosssection = sigma0/4.0d0 * & !cm^2
       (neutrino_energy/m_e)**2 * & !dimensionless
       (S_V*(c_v_prime-1.0d0)**2+3.0d0*gAeff**2*S_A*(c_a_prime-1.0d0)**2)

  !actually want transport crosssection, in cm^2
  if (type.eq.1) then
     delta_p = (S_V*(c_v_prime-1.0d0)**2-gAeff**2*S_A*(c_a_prime-1.0d0)**2)/ &
          (S_V*(c_v_prime-1.0d0)**2+3.0d0*gAeff**2*S_A*(c_a_prime-1.0d0)**2)
     if (do_weak_mag_corrections) then
        weak_mag = weak_mag_correction_scattering_transport(neutrino_energy,1,lepton_local) !1 for proton
     else
        weak_mag = 1.0d0
     endif
     crosssection = crosssection*(1.0d0-delta_p/3.0d0)*weak_mag !implicit integration over \Omega
  else if (type.eq.0) then
     if (do_weak_mag_corrections) then
        weak_mag = weak_mag_correction_scattering_total(neutrino_energy,1,lepton_local) !1 for proton
     else
        weak_mag = 1.0d0
     endif
     crosssection = crosssection*weak_mag !implicit integration over \Omega
  endif
  
end subroutine nu_scatter_elastic_p_total

subroutine nu_scatter_elastic_n_total(neutrino_energy,type,lepton,eos_variables,crosssection)

  use nulib
  implicit none

  !inputs
  real*8, intent(in)  :: neutrino_energy !in MeV
  real*8, intent(in)  :: eos_variables(total_eos_variables)
  integer, intent(in) :: type !transport (type.eq1) or absorption (type.eq.0) crosssection
  integer, intent(in) :: lepton !lepton number, >0 for neutrino,  <0 for antineutrino

  !output
  real*8,intent(out) :: crosssection !in cm^2

  !local
  real*8 :: delta_n,weak_mag
  integer :: lepton_local

  !strange virial correction via S_A
  real*8 :: gAeff,S_A,S_V
  real*8 :: virialA,virialB,virialC
  real*8 :: virial_nb,virial_ye,virial_temp

  !function declaration
  real*8 :: weak_mag_correction_scattering_transport
  real*8 :: weak_mag_correction_scattering_total

  !for scattering on n, only need sign of lepton number
  if (lepton.gt.0) then
     lepton_local = 1
  else
     lepton_local = -1
  endif

  !a la Horowitz et al. (2016)
  !note, this includes a BS (1998) estimate at high densities
  if (do_nc_virial_correction) then
     virial_nb = eos_variables(rhoindex)*1.0d-39/(m_ref*mev_to_gram)
     virial_ye = eos_variables(yeindex)
     virial_temp = eos_variables(tempindex)
     
     virialA = 920.0d0*virial_nb*(1.0d0-virial_ye+virial_ye**2)/virial_temp**1.22d0
     virialB = 3.05d0/virial_temp**0.75d0
     virialC = 6140.0d0*virial_nb*virial_ye*(1.0d0-virial_ye)/virial_temp**0.5d0 + &
          1.5d13*virial_nb**4/virial_temp**6
     S_A = 1.0d0/(1.0d0+virialA*(1.d0+virialB*exp(-virialC)))
     S_V = 1.0d0

  else
     S_A = 1.0d0
     S_V = 1.0d0
  endif

  !Horowitz (2002)
  if (do_strange_coupling) then
     gAeff = gA+gAs ! + because neutrons
  else
     gAeff = gA
  endif


  !BRT06 Eq.20
  crosssection = sigma0/4.0d0 * & !cm^2
       (neutrino_energy/m_e)**2 * & !dimensionless
       (S_V+3.0d0*gAeff**2*S_A)/4.0d0 ! dimensionless

  !actually want transport crosssection, in cm^2
  if (type.eq.1) then
     delta_n = (S_V-gAeff**2*S_A)/(S_V+3.0d0*gAeff**2*S_A)
     if (do_weak_mag_corrections) then
        weak_mag = weak_mag_correction_scattering_transport(neutrino_energy,0,lepton_local) !0 for neutron
     else
        weak_mag = 1.0d0
     endif
     crosssection = crosssection*(1.0d0-delta_n/3.0d0)*weak_mag !implicit integration over \Omega
  else if (type.eq.0) then
     if (do_weak_mag_corrections) then
        weak_mag = weak_mag_correction_scattering_total(neutrino_energy,0,lepton_local) !0 for neutron
     else
        weak_mag = 1.0d0
     endif
     crosssection = crosssection*weak_mag !implicit integration over \Omega
  endif

end subroutine nu_scatter_elastic_n_total

!differential scattering crosssections
function nu_scatter_elastic_p_differential(neutrino_energy,cosine_angle,lepton,eos_variables) result(dcrosssection)

  use nulib
  implicit none

  !inputs
  real*8, intent(in)  :: neutrino_energy !in MeV
  real*8, intent(in)  :: eos_variables(total_eos_variables)
  real*8, intent(in)  :: cosine_angle !i.e. \mu = cos(scattering_angle)
  integer, intent(in) :: lepton !-1 for antineutrino, 1 for neutrino

  !output
  real*8 :: dcrosssection !cm^2/srad
  
  !local
  real*8  :: delta_p
  real*8  :: c_v_prime = 0.5d0+2.0d0*sin2thetaW
  real*8  :: c_a_prime = 0.5d0
  integer :: type = 0 !if using BRT, we take total cross section and
                      !the \delta defintion to get the differential
                      !cross section, for this we set type = 0
  integer :: lepton_local
  real*8 :: c_v = 0.5d0-2.0d0*sin2thetaW
  real*8 :: c_a = gA/2.0d0
  real*8 :: weak_mag !correction factor

  !function declaration
  real*8 :: weak_mag_correction_scattering_differential

  !for scattering on p, only need sign of lepton number
  if (lepton.gt.0) then
     lepton_local = 1
  else
     lepton_local = -1
  endif

!  stop "This differential routine does not integrate over \phi, is that what you want?"

  if (do_nc_virial_correction) stop "virial correction not in differential cross section"
  if (do_strange_coupling) stop "strange coupling not in differential cross section"

  if (do_weak_mag_corrections) then
     
     !note that the integral of this dcrosssection over {\phi,0,2pi} and {\theta,0,pi} gives nu_scatter_elastic_n_total
     weak_mag = weak_mag_correction_scattering_differential(neutrino_energy,cosine_angle,0,lepton)
     dcrosssection = sigma0/(16.0d0*pi) * & !cm^2/sterad
          (neutrino_energy/m_e)**2 * & !dimensionless
          (c_v**2*(1.0d0+cosine_angle)+c_a**2*(3.0d0-cosine_angle)) * &
          weak_mag !correction

  else
     !BRT way for ease
     delta_p = ((c_v_prime-1.0d0)**2-gA**2*(c_a_prime-1.0d0)**2)/ &
          ((c_v_prime-1.0d0)**2+3.0d0*gA**2*(c_a_prime-1.0d0)**2)
     !BRT06 Eq.21
     !note that the integral of this dcrosssection over {\phi,0,2pi} and {\theta,0,pi} gives nu_scatter_elastic_p_total
     call nu_scatter_elastic_p_total(neutrino_energy,type,lepton_local,eos_variables,dcrosssection)
     dcrosssection = dcrosssection/(4.0d0*pi) * & !total crosssection cm^2/sterad
          (1.0d0+cosine_angle*delta_p) !angular part, dimensionless
  endif
       
end function nu_scatter_elastic_p_differential

function nu_scatter_elastic_n_differential(neutrino_energy,cosine_angle,lepton,eos_variables) result(dcrosssection)

  use nulib
  implicit none

  !inputs
  real*8, intent(in)  :: neutrino_energy !in MeV
  real*8, intent(in)  :: eos_variables(total_eos_variables)
  real*8, intent(in)  :: cosine_angle !i.e. \mu = cos(scattering_angle)
  integer, intent(in) :: lepton !-1 for antineutrino, 1 for neutrino

  !output
  real*8 :: dcrosssection !cm^2/srad
  
  !local
  real*8 :: delta_n
  integer :: type = 0 !we do not want transport differential crosssection, not sure what this is....
  integer :: lepton_local
  real*8 :: c_v = -0.5d0 
  real*8 :: c_a = -gA/2.0d0
  real*8 :: weak_mag !correction factor

  !function declaration
  real*8 :: weak_mag_correction_scattering_differential

  !for scattering on n, only need sign of lepton number
  if (lepton.gt.0) then
     lepton_local = 1
  else
     lepton_local = -1
  endif

  stop "This differential routine does not integrate over \phi, is that what you want?"

  if (do_nc_virial_correction) stop "virial correction not in differential cross section"
  if (do_strange_coupling) stop "strange coupling not in differential cross section"

  if (do_weak_mag_corrections) then
     
     !note that the integral of this dcrosssection over {\phi,0,2pi} and {\theta,0,pi} gives nu_scatter_elastic_n_total
     weak_mag = weak_mag_correction_scattering_differential(neutrino_energy,cosine_angle,1,lepton)
     dcrosssection = sigma0/(16.0d0*pi) * & !cm^2/sterad
          (neutrino_energy/m_e)**2 * & !dimensionless
          (c_v**2*(1.0d0+cosine_angle)+c_a**2*(3.0d0-cosine_angle)) * &
          weak_mag !correction

  else
     !BRT way for ease
     delta_n = (1.0d0-gA**2)/(1.0d0+3.0d0*gA**2)
     !BRT06 Eq.21
     !note that the integral of this dcrosssection over {\phi,0,2pi} and {\theta,0,pi} gives nu_scatter_elastic_n_total
     call nu_scatter_elastic_n_total(neutrino_energy,type,lepton_local,eos_variables,dcrosssection)
     dcrosssection = dcrosssection/(4.0d0*pi) * & !total crosssection cm^2/sterad
          (1.0d0+cosine_angle*delta_n) !angular part, dimensionless
  endif

       
end function nu_scatter_elastic_n_differential

!total crossection for neutrino nucleo scattering (integral of differential)
subroutine nu_scatter_elastic_heavy_total(neutrino_energy,transport,lepton,eos_variables,crosssection)

  use nulib
  implicit none

  !inputs
  real*8, intent(in)  :: eos_variables(total_eos_variables)
  real*8, intent(in)  :: neutrino_energy !in MeV
  integer, intent(in) :: transport !1 for transport, 0 for total
  integer, intent(in) :: lepton !-1 for antineutrino, 1 for neutrino

  !outputs
  real*8, intent(out) :: crosssection !total crosssection for neutrino nuclei scattering

  !local
  integer :: i !index for integral
  real*8 :: cosine !angle
  integer :: lepton_local

  !function declaration
  real*8 :: nu_scatter_elastic_heavy_differential
  
  if (transport.eq.0) then
     stop "crossection assumes transport (i.e. (1.0d0-cosine) term)"
  endif

  !for scattering on heavy, only need sign of lepton number
  if (lepton.gt.0) then
     lepton_local = 1
  else
     lepton_local = -1
  endif

  !integrate over theta
  crosssection = 0.0d0
  do i=1,16
     cosine = GPQ_n16_roots(i)
     crosssection = crosssection + (1.0d0-cosine)* &
          nu_scatter_elastic_heavy_differential(neutrino_energy, &
          cosine,transport,lepton_local,eos_variables)*GPQ_n16_weights(i)
  enddo
  
  !integrate over phi
  crosssection = 2.0d0*pi*crosssection

end subroutine nu_scatter_elastic_heavy_total
  

!differential crosssection for neutrino nuclei scatters
function nu_scatter_elastic_heavy_differential(neutrino_energy, &
     cosine_angle,transport,lepton,eos_variables) result(dcrosssection)
  
  use nulib
  implicit none

  !inputs
  real*8, intent(in)  :: eos_variables(total_eos_variables)
  real*8, intent(in)  :: neutrino_energy !in MeV
  real*8, intent(in)  :: cosine_angle !i.e. \mu = cos(scattering_angle)
  integer, intent(in) :: transport !1 for transport, 0 for total
  integer, intent(in) :: lepton !-1 for antineutrino, 1 for neutrino

  !outputs
  real*8 :: dcrosssection !differential crosssection for neutrino nuclei scattering

  !local
  integer :: i !counter
  real*8 :: CFF,y !form factor term BRT06 Eq. 31
  real*8 :: W !cal W from BRT06 Eq. 26
  real*8 :: CLOS !electron polarization correction BRT06 Eq. 29-30
  integer :: lepton_local

  !ion-ion correlation BRT06 Eq. 27-28 Horowitz 1997 PRD 55 4577, 
  !note that Sion from Chuck is angle averaged, any couple between the terms over theta is ignored, 
  !however, the correction term of 3/4 in Eq. 5 of Chuck's paper corrects for the 
  !extra factors of (1-cos) and (1+cos).  Therefore Sion must only be used for transport crosssections
  real*8 :: Sion

  real*8 :: electron_number_density ! electron number density, number/cm^3
  real*8 :: fermi_energy !fermi energy of electrons, MeV
  real*8 :: fermi_momentum !fermi momentum of electrons, MeV
  real*8 :: k !momentum transfer, 1/cm
  real*8 :: r_debye ! debye radius, cm, BRT06 Eq. 30
  real*8 :: exparguement !for an exponetial
  real*8 :: Ebar, Estar, Gamma, electric_term, ion_numden,interparticle_spacing !Sion variables

  !function declarations
  real*8 :: Sion_beta

  real*8 ktimesrdsquared

  !for scattering on heavies, only need sign of lepton number
  if (lepton.gt.0) then
     lepton_local = 1
  else
     lepton_local = -1
  endif

  if (transport.eq.0) then
     stop "Sion correction is determined for transport differential cross section"
  endif

  W = 1.0d0 - 2.0d0*eos_variables(zbarindex)/eos_variables(abarindex)*(1.0d0-2.0d0*sin2thetaW)
  
  !Sion correction calculation
  ion_numden = eos_variables(xhindex)*eos_variables(rhoindex)/(eos_variables(abarindex)*m_ref*mev_to_gram) !ion density, divide heavy nuclei mass fraction by average mass of heavy nuclei
  interparticle_spacing = (3.0d0/(4.0d0*pi*ion_numden))**(1.0d0/3.0d0) !in cm
  Ebar = neutrino_energy*interparticle_spacing/hbarc_mevcm !dimensionless
  electric_term = 1.43996439d-13 !e^2/(4*pi*\epsilon_0) in MeV cm
  Gamma = min(max(1.0d0,eos_variables(zbarindex)**2*electric_term/interparticle_spacing/eos_variables(tempindex)),150.0d0) !dimensionless
  Estar = 3.0d0 + 4.0d0/sqrt(Gamma)
  
  if (Ebar.gt.Estar) then
     Sion = 1.0d0
  else
     exparguement = 0.0d0
     do i = 0,6
        exparguement = exparguement - Sion_beta(i,Gamma)*Ebar**i
     enddo
     Sion =  1.0d0/(1.0d0+exp(exparguement))
  endif
  
  electron_number_density = eos_variables(rhoindex)*eos_variables(yeindex)/(m_ref*mev_to_gram) !number /cm^3
  fermi_energy = eos_variables(mueindex)!I think this is the chemical potential,  MeV
  !when far out of equilibrium, matter_mue < m_e.  This is not good, unphysical, and causes errors.
  fermi_energy = max(eos_variables(mueindex),m_e+1.0d-10)

  fermi_momentum = sqrt(fermi_energy**2-m_e**2) !MeV
  k = sqrt(2.0d0*(neutrino_energy/hbarc_mevcm)**2*(1.0d0-cosine_angle)) !cm^-1
  r_debye = sqrt(pi*hbarc_mevcm**2/(4.0d0*FSC*fermi_momentum*fermi_energy)) !cm
  !does this matter for mu ans tau neutrino scattering on heavies???, 
  CLOS = eos_variables(zbarindex)/eos_variables(abarindex)*((1.0d0+4.0d0*sin2thetaW)/(1.0d0+(k*r_debye)**2))
  y = (neutrino_energy/56.0d0)**2*(eos_variables(abarindex)/100.0d0)**(2.0d0/3.0d0)
  CFF = exp(-y*(1.0d0-cosine_angle)/2.0d0)

  if (.not.do_ionion_correlation) then
     Sion = 1.0d0
  endif

  if (.not.do_heavyscat_formfactor) then
     CFF = 1.0d0
  endif
  
  if (.not.do_electronpolarization_correction) then
     CLOS = 0.0d0
  endif

  dcrosssection = sigma0/(64.0d0*pi)*(neutrino_energy/m_e)**2* &
       eos_variables(abarindex)**2*(W*CFF+CLOS)**2*Sion*(1.0d0+cosine_angle)
  
end function nu_scatter_elastic_heavy_differential

subroutine total_scattering_opacity(neutrino_species,neutrino_energy,scattering_opacity,eos_variables)

  use nulib
  implicit none
  
  !inputs
  real*8, intent(in) :: eos_variables(total_eos_variables)  
  integer, intent(in) :: neutrino_species  !integer 1 through 9
  real*8, intent(in) :: neutrino_energy !MeV
  
  !outputs
  real*8 :: scattering_opacity !cm^-1
  
  !temporary variables
  real*8 :: crosssection
  
  scattering_opacity = 0.0d0

  !electron neutrino
  if (neutrino_species.eq.1) then
     !scattering (transport cross section) on neutrons
     if (add_nue_scattering_n) then
        !function call takes neutrino energy, transport=1, lepton number = 1
        call nu_scatter_elastic_n_total(neutrino_energy,1,1,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xnindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section) on protons
     if (add_nue_scattering_p) then
        !function call takes neutrino energy, transport=1, lepton number = 1
        call nu_scatter_elastic_p_total(neutrino_energy,1,1,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xpindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section) on heavies
     if (add_nue_scattering_heavies.and.eos_variables(xhindex).ne.0.0d0) then
        !function call takes neutrino energy, transport=1, lepton = 1
        call nu_scatter_elastic_heavy_total(neutrino_energy,1,1,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xhindex)/eos_variables(abarindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

     !scattering (transport cross section??) on electrons
     if (add_nue_scattering_electrons) then
        ! matter temperature, density, transport=1, lepton = 1
        call nu_scatter_elastic_e_total(neutrino_energy,1,1,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(yeindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

     !scattering on alpha
     if (add_nue_scattering_alphas) then
        ! matter temperature, density, transport=1, lepton = 1
        call nu_scatter_elastic_alpha_total(neutrino_energy,1,1,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * (eos_variables(xaindex)/4.0d0)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

  !electron antineutrino
  else if (neutrino_species.eq.2) then
     !scattering (transport cross section) on neutrons
     if (add_anue_scattering_n) then
        !function call takes neutrino energy, transport=1, lepton number = -1
        call nu_scatter_elastic_n_total(neutrino_energy,1,-1,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xnindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section) on protons
     if (add_anue_scattering_p) then
        !function call takes neutrino energy, transport=1, lepton number = -1
        call nu_scatter_elastic_p_total(neutrino_energy,1,-1,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xpindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif
     
     !scattering (transport cross section) on heavies
     if (add_anue_scattering_heavies.and.eos_variables(xhindex).ne.0.0d0) then
        !function call takes neutrino energy, transport=1, lepton = -1
        call nu_scatter_elastic_heavy_total(neutrino_energy,1,-1,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xhindex)/eos_variables(abarindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section??) on electrons
     if (add_anue_scattering_electrons) then
        ! matter temperature, density, transport=1, lepton = -1
        call nu_scatter_elastic_e_total(neutrino_energy,1,-1,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(yeindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

     !scattering on alpha
     if (add_anue_scattering_alphas) then
        ! matter temperature, density, transport=1, lepton = -1
        call nu_scatter_elastic_alpha_total(neutrino_energy,1,-1,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * (eos_variables(xaindex)/4.0d0)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 
     
  !mu neutrino
  else if (neutrino_species.eq.3) then
     !scattering (transport cross section) on neutrons
     if (add_numu_scattering_n) then
        !function call takes neutrino energy, transport=1, lepton number = 2
        call nu_scatter_elastic_n_total(neutrino_energy,1,2,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xnindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section) on protons
     if (add_numu_scattering_p) then
        !function call takes neutrino energy, transport=1, lepton number = 2
        call nu_scatter_elastic_p_total(neutrino_energy,1,2,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xpindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section) on heavies
     if (add_numu_scattering_heavies.and.eos_variables(xhindex).ne.0.0d0) then
        !function call takes neutrino energy, transport=1, lepton = 2
        call nu_scatter_elastic_heavy_total(neutrino_energy,1,2,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xhindex)/eos_variables(abarindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section??) on electrons
     if (add_numu_scattering_electrons) then
        ! matter temperature, density, transport=1, lepton = 2
        call nu_scatter_elastic_e_total(neutrino_energy,1,2,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(yeindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

     !scattering on alpha
     if (add_numu_scattering_alphas) then
        ! matter temperature, density, transport=1, lepton = 2
        call nu_scatter_elastic_alpha_total(neutrino_energy,1,2,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * (eos_variables(xaindex)/4.0d0)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

  !mu antineutrino
  else if (neutrino_species.eq.4) then
     !scattering (transport cross section) on neutrons
     if (add_anumu_scattering_n) then
        !function call takes neutrino energy, transport=1, lepton number = -2
        call nu_scatter_elastic_n_total(neutrino_energy,1,-2,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xnindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section) on protons
     if (add_anumu_scattering_p) then
        !function call takes neutrino energy, transport=1, lepton number = -2
        call nu_scatter_elastic_p_total(neutrino_energy,1,-2,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xpindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif
     
     !scattering (transport cross section) on heavies
     if (add_anumu_scattering_heavies.and.eos_variables(xhindex).ne.0.0d0) then
        !function call takes neutrino energy, transport=1, lepton = -2
        call nu_scatter_elastic_heavy_total(neutrino_energy,1,-2,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xhindex)/eos_variables(abarindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section??) on electrons
     if (add_anumu_scattering_electrons) then
        ! matter temperature, density, transport=1, lepton = -2
        call nu_scatter_elastic_e_total(neutrino_energy,1,-2,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(yeindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

     !scattering on alpha
     if (add_anumu_scattering_alphas) then
        ! matter temperature, density, transport=1, lepton = -2
        call nu_scatter_elastic_alpha_total(neutrino_energy,1,-2,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * (eos_variables(xaindex)/4.0d0)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

  !tau neutrino
  else if (neutrino_species.eq.5) then
     !scattering (transport cross section) on neutrons
     if (add_nutau_scattering_n) then
        !function call takes neutrino energy, transport=1, lepton number = 3
        call nu_scatter_elastic_n_total(neutrino_energy,1,3,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xnindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section) on protons
     if (add_nutau_scattering_p) then
        !function call takes neutrino energy, transport=1, lepton number = 3
        call nu_scatter_elastic_p_total(neutrino_energy,1,3,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xpindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif
     
     !scattering (transport cross section) on heavies
     if (add_nutau_scattering_heavies.and.eos_variables(xhindex).ne.0.0d0) then
        !function call takes neutrino energy, transport=1, lepton = 3
        call nu_scatter_elastic_heavy_total(neutrino_energy,1,3,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xhindex)/eos_variables(abarindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section??) on electrons
     if (add_nutau_scattering_electrons) then
        ! matter temperature, density, transport=1, lepton = 3
        call nu_scatter_elastic_e_total(neutrino_energy,1,3,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(yeindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

     !scattering on alpha
     if (add_nutau_scattering_alphas) then
        ! matter temperature, density, transport=1, lepton = 3
        call nu_scatter_elastic_alpha_total(neutrino_energy,1,3,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * (eos_variables(xaindex)/4.0d0)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

  !tau antineutrino
  else if (neutrino_species.eq.6) then
     !scattering (transport cross section) on neutrons
     if (add_anutau_scattering_n) then
        !function call takes neutrino energy, transport=1, lepton number = -3
        call nu_scatter_elastic_n_total(neutrino_energy,1,-3,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xnindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section) on protons
     if (add_anutau_scattering_p) then
        !function call takes neutrino energy, transport=1, lepton number = -3
        call nu_scatter_elastic_p_total(neutrino_energy,1,-3,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xpindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif
     
     !scattering (transport cross section) on heavies
     if (add_anutau_scattering_heavies.and.eos_variables(xhindex).ne.0.0d0) then
        !function call takes neutrino energy, transport=1, lepton = -3
        call nu_scatter_elastic_heavy_total(neutrino_energy,1,-3,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(xhindex)/eos_variables(abarindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif

     !scattering (transport cross section??) on electrons
     if (add_anutau_scattering_electrons) then
        ! matter temperature, density, transport=1, lepton = -3
        call nu_scatter_elastic_e_total(neutrino_energy,1,-3,eos_variables,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * eos_variables(yeindex)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

     !scattering on alpha
     if (add_anutau_scattering_alphas) then
        ! matter temperature, density, transport=1, lepton = -3
        call nu_scatter_elastic_alpha_total(neutrino_energy,1,-3,crosssection)
        scattering_opacity = scattering_opacity + &
             crosssection * (eos_variables(xaindex)/4.0d0)*eos_variables(rhoindex)/(m_ref*mev_to_gram)
     endif 

  else
     stop "total_scattering_opacity: How did you get in here??"
  endif
     

end subroutine total_scattering_opacity


subroutine return_scattering_opacity_spectra_given_neutrino_scheme(scattering_opacity_spectra,eos_variables)
  
  use nulib
  implicit none
  
  !inputs & outputs
  real*8, intent(in) :: eos_variables(total_eos_variables)
  real*8, intent(out) :: scattering_opacity_spectra(number_species,number_groups)
  
  !locals
  integer :: ns,ng
  real*8 :: opacity
  
  if (size(scattering_opacity_spectra,1).ne.number_species) then
     stop "return_scatter_opacity_spectra_given_neutrino_scheme:provided array has wrong number of species"
  endif
  if (size(scattering_opacity_spectra,2).ne.number_groups) then
     stop "return_scatter_opacity_spectra_given_neutrino_scheme:provided array has wrong number of groups"
  endif
  
  do ns=1,number_species
     do ng=1,number_groups
        call total_scattering_opacity(ns,energies(ng),opacity,eos_variables)
        scattering_opacity_spectra(ns,ng) = opacity
     enddo
  enddo
  
end subroutine return_scattering_opacity_spectra_given_neutrino_scheme
