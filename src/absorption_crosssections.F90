!-*-f90-*-
!always keep term dimensionless unless it is the units of the final answer...
function nue_absorption_on_n(neutrino_energy,eos_variables) result(crosssection)
  
  use nulib
  implicit none

  real*8, intent(in) :: eos_variables(total_eos_variables)  
  real*8, intent(in) :: neutrino_energy !MeV
  real*8 :: crosssection !final answer in cm^2

  !local variables
  real*8 :: final_electron_energy !MeV
  real*8 :: feminus_exp_log10ed,SA_exp_log10ed,feminus_over_SA_exp_log10ed !dimensionless
  real*8 :: one_plus_feminus_exp_log10ed,one_plus_SA_exp_log10ed !dimensionless
  real*8 :: weak_mag !dimensionless
  real*8 :: mu_nu_eq !MeV
  real*8 :: logterm,expterm !dimensionless

  !function declarations
  real*8 :: weak_mag_correction_absorption
  real*8 :: fermidirac_exptermonly_log10ed

  !based on Todd Thompson PhD Appendix B, Eq. B8, final state proton
  !blocking is done outside of this subroutine (bottom of
  !absorption_crosssections.F90)

  final_electron_energy = neutrino_energy + delta_np
  mu_nu_eq = eos_variables(mueindex)-eos_variables(muhatindex)
  feminus_exp_log10ed = fermidirac_exptermonly_log10ed(final_electron_energy, &
       eos_variables(mueindex),eos_variables(tempindex)) !final state electron blocking
  SA_exp_log10ed = fermidirac_exptermonly_log10ed(neutrino_energy,mu_nu_eq,eos_variables(tempindex))
  if (do_weak_mag_corrections) then
     weak_mag = weak_mag_correction_absorption(neutrino_energy,0) !to first order 1.0d0+1.01d0*neutrino_energy/m_ref, Horowitz 2002
  else
     weak_mag = 1.0d0
  endif

  !Note on the exp's.  The Fermi functions and simulated absorption
  !terms can be huge/small.  The best way to deal with this is to
  !play tricks.  We write them all in terms of the log of the exp, not
  !the fermi functions, and then combine appropiately (taking into
  !account the size of the exp to deal with the +1's appropiately

  !Note, we can also combine the feminus_exp/SA_exp term
  !the arguement of feminus_exp is neutrino_energy + delta_np - matter_mue
  !the arguement of SA_exp is neutrino_energy - matter_nue + matter_muhat
  !feplus_exp/SA_exp = exponential with arguement of
  !delta_np-matter_muhat
  feminus_over_SA_exp_log10ed = (delta_np-eos_variables(muhatindex))/eos_variables(tempindex)*log10exp

  !The full term we are lumping together is:
  !(1-f_{e-})/(1-f_{\nu_e}^{eq}) =
  !fexp_{e-}*(1+fexp_{SA})/(1+fexp_{e-})/fexp_{SA} =
  !10.0d0**(feminus_over_SA_exp_log10ed - one_plus_feminus_exp_log10ed
  !+ one_plus_SA_exp_log10ed)

  !deal with fermi functions +1's
  if (feminus_exp_log10ed.gt.30.0d0) then !exp >> 1
     one_plus_feminus_exp_log10ed = feminus_exp_log10ed 
  else if (feminus_exp_log10ed.lt.-30.0d0) then !exp << 1
     one_plus_feminus_exp_log10ed = 0.0d0
  else 
     one_plus_feminus_exp_log10ed = log10(1.0d0+10.0d0**(feminus_exp_log10ed))
  endif

  if (SA_exp_log10ed.gt.30.0d0) then  !exp >> 1
     one_plus_SA_exp_log10ed = SA_exp_log10ed
  else if (SA_exp_log10ed.lt.-30.0d0) then !exp << 1
     one_plus_SA_exp_log10ed = 0.0d0
  else
     one_plus_SA_exp_log10ed = log10(1.0d0+10.0d0**(SA_exp_log10ed))
  endif

  logterm = min(200.0d0,max(-200.0d0,feminus_over_SA_exp_log10ed - one_plus_feminus_exp_log10ed + one_plus_SA_exp_log10ed))
  expterm = 10.0d0**(logterm)

  crosssection = sigma0 * & !cm^2
          (1.0d0+3.0d0*gA**2) / 4.0d0 * & !dimensionless
          (final_electron_energy/m_e)**2 * & !dimensionless
          (1.0d0-(m_e/final_electron_energy)**2)**0.5d0 * & !dimensionless
          expterm*weak_mag !dimensionless

end function nue_absorption_on_n

function anue_absorption_on_p(neutrino_energy,eos_variables) result(crosssection)
  
  use nulib
  implicit none

  real*8, intent(in) :: eos_variables(total_eos_variables)  
  real*8, intent(in) :: neutrino_energy !MeV
  real*8  :: crosssection !final answer in cm^2

  !local variables
  real*8 :: final_positron_energy !MeV
  real*8 :: feplus_exp_log10ed,SA_exp_log10ed,feplus_over_SA_exp_log10ed !dimensionless
  real*8 :: one_plus_SA_exp_log10ed,one_plus_feplus_exp_log10ed !dimensionless
  real*8 :: weak_mag !dimensionless
  real*8 :: mu_nu_eq !MeV
  real*8 :: expterm, logterm !dimensionless

  !function declarations
  real*8 :: fermidirac_exptermonly_log10ed
  real*8 :: weak_mag_correction_absorption

  !based on Todd Thompson PhD Appendix B, Eq. B11, final state neutron
  !blocking is done outside of this subroutine (bottom of
  !absorption_crosssections.F90)
  final_positron_energy = neutrino_energy - delta_np
  if (final_positron_energy.lt.2.0d0*m_e) then
     !only happens if neutrino energy < ~2.3 MeV
     !we choose 2m_e because GR1D has trouble with zero opacity
     final_positron_energy = 2.0d0*m_e
  endif

  mu_nu_eq = -(eos_variables(mueindex)-eos_variables(muhatindex))
  feplus_exp_log10ed = fermidirac_exptermonly_log10ed(final_positron_energy,-eos_variables(mueindex),eos_variables(tempindex))
  SA_exp_log10ed = fermidirac_exptermonly_log10ed(neutrino_energy,mu_nu_eq,eos_variables(tempindex))
  if (do_weak_mag_corrections) then
     weak_mag = weak_mag_correction_absorption(neutrino_energy,1) !to first order 1.0d0-7.1d0*neutrino_energy/m_p, Horowitz 2002
  else
     weak_mag = 1.0d0
  endif

  !Note on the exp's.  The Fermi functions and simulated absorbtion
  !terms can be huges/small.  The best way to deal with this is to
  !play tricks.  We write them all in terms of the log of the exp, not
  !the fermi functions, and then combine appropiately (taking into
  !account the size of the exp to deal with the +1 appropiately

  !Note, we can also combine the feplus_exp/SA_exp term
  !the arguement of feplus_exp is neutrino_energy - delta_np + matter_mue
  !the arguement of SA_exp is neutrino_energy + matter_nue - matter_muhat
  !feplus_exp/SA_exp = exponential with arguement of
  !-delta_np+matter_muhat
  feplus_over_SA_exp_log10ed = (-delta_np+eos_variables(muhatindex))/eos_variables(tempindex)*log10exp

  !The full term we are lumping together is:
  !(1-f_{e+})/(1-f_{\bar{nu}_e}^{eq}) =
  !fexp_{e+}*(1+fexp_{SA})/(1+fexp_{e+})/fexp_{SA} =
  !10.0d0**(feplus_over_SA_exp_log10ed - one_plus_feplus_exp_log10ed
  !+ one_plus_SA_exp_log10ed)

  !deal with fermi functions +1's
  if (feplus_exp_log10ed.gt.30.0d0) then !exp >> 1
     one_plus_feplus_exp_log10ed = feplus_exp_log10ed 
  else if (feplus_exp_log10ed.lt.-30.0d0) then !exp << 1
     one_plus_feplus_exp_log10ed = 0.0d0
  else 
     one_plus_feplus_exp_log10ed = log10(1.0d0+10.0d0**(feplus_exp_log10ed))
  endif

  if (SA_exp_log10ed.gt.30.0d0) then  !exp >> 1
     one_plus_SA_exp_log10ed = SA_exp_log10ed
  else if (SA_exp_log10ed.lt.-30.0d0) then !exp << 1
     one_plus_SA_exp_log10ed = 0.0d0
  else
     one_plus_SA_exp_log10ed = log10(1.0d0+10.0d0**(SA_exp_log10ed))
  endif

  !add all log terms up together, and also apply limits for extreme cases!!!
  logterm = min(200.0d0,max(-200.0d0,feplus_over_SA_exp_log10ed - one_plus_feplus_exp_log10ed + one_plus_SA_exp_log10ed))

  expterm = 10.0d0**(logterm)

  crosssection = sigma0 * & !cm^2
       (1.0d0+3.0d0*gA**2) / 4.0d0 * & !dimensionless
       (final_positron_energy/m_e)**2 * & !dimensionless
       (1.0d0-(m_e/final_positron_energy)**2)**0.5d0 * & !dimensionless
       expterm*weak_mag !dimensionless

end function anue_absorption_on_p

function nux_absorption_on_n_and_p(neutrino_energy,eos_variables) result(crosssection)

  use nulib
  implicit none

  real*8, intent(in) :: eos_variables(total_eos_variables)  
  real*8, intent(in) :: neutrino_energy !MeV
  real*8  :: crosssection !final answer in cm^2

  real*8 :: entfactor,rhofactor

  !only want heating in certain regions, smoothly cut it off outside
  !that region.  We'll follow Fischer's choice; s>6, rho<1e10, we'll
  !smooth it.

  if (eos_variables(entropyindex).gt.6.0d0) then
     entfactor = 1.0d0
  else if (eos_variables(entropyindex).gt.4.0d0) then
     entfactor = (eos_variables(entropyindex)-4.0d0)/2.0d0
  else
     entfactor = 0.0d0
  endif

  if (eos_variables(rhoindex).lt.1.0d10) then
     rhofactor = 1.0d0
  else if (eos_variables(rhoindex).lt.2.0d10) then
     rhofactor = 2.0d0-eos_variables(rhoindex)/1.0d10
  else
     rhofactor = 0.0d0
  endif
  
  crosssection = sigma0 * & !cm^2
       (1.0d0+3.0d0*gA**2) / 4.0d0 * & !dimensionless
       (neutrino_energy/m_e)**2 !dimensionless

  crosssection = crosssection*entfactor*rhofactor*adhoc_nux_factor

end function nux_absorption_on_n_and_p

function nue_absorption_on_A(neutrino_energy,eos_variables) result(crosssection)

  use nulib
  implicit none

  real*8, intent(in) :: eos_variables(total_eos_variables)  
  real*8, intent(in) :: neutrino_energy !MeV
  real*8 :: crosssection !final answer in cm^2

  !local variables
  real*8 :: N_p_Z, N_n_N !dimensionless
  real*8 :: final_electron_energy !MeV, energy of electron in products of interaction
  real*8 :: feminus_exp_log10ed,SA_exp_log10ed !dimensionless
  real*8 :: one_plus_feminus_exp_log10ed,one_plus_SA_exp_log10ed !dimensionless
  real*8 :: exponential_term_log10ed !dimensionless, place holder
  real*8 :: mu_nu_eq
  real*8 :: logterm, expterm

  !function declarations
  real*8 :: fermidirac_exptermonly_log10ed

  !based on Todd Thompson PhD Appendix B, Eq. B13, and BRT06 Eq. 12,
  !BRT06 say better expressions out there...  nu_energy + Qprime, as
  !in Thompson, mus have rest mass difference included
  final_electron_energy = neutrino_energy  + eos_variables(muhatindex) + 3.0d0 
  !At low energies, muhat can be quite negative, therefore
  !final_energy_energy < m_e... Rampp & Janka set crosssection to zero
  if (final_electron_energy.lt.m_e+1.d-10) then
     crosssection = 0.0d0
     return
  endif

  mu_nu_eq = eos_variables(mueindex)-eos_variables(muhatindex)
  feminus_exp_log10ed = fermidirac_exptermonly_log10ed(final_electron_energy,eos_variables(mueindex),eos_variables(tempindex))
  SA_exp_log10ed =fermidirac_exptermonly_log10ed(neutrino_energy,mu_nu_eq,eos_variables(tempindex))
  !Boltzmann term for neutron to be in 3.0MeV exicted state (see Thompson and Bruenn 1985) 
  exponential_term_log10ed = -3.0d0/eos_variables(tempindex)*log10exp

  !Zbar_Nterm 
  if (eos_variables(zbarindex).lt.20.0d0) then
     N_p_Z = 0.0d0
  else if (eos_variables(zbarindex).lt.28.0d0) then
     N_p_Z = eos_variables(zbarindex) - 20.0d0
  else 
     N_p_Z = 8.0d0
  endif
  
  !Nbar_Nterm 
  if (eos_variables(abarindex)-eos_variables(zbarindex).lt.34.0d0) then
     N_n_N = 6.0d0
  else if (eos_variables(abarindex)-eos_variables(zbarindex).lt.40.0d0) then
     N_n_N = 40.0d0 - (eos_variables(abarindex)-eos_variables(zbarindex))
  else 
     N_n_N = 0.0d0
  endif

  if (SA_exp_log10ed.gt.30.0d0) then  !exp >> 1
     one_plus_SA_exp_log10ed = SA_exp_log10ed
  else if (SA_exp_log10ed.lt.-30.0d0) then !exp << 1
     one_plus_SA_exp_log10ed = 0.0d0
  else
     one_plus_SA_exp_log10ed = log10(1.0d0+10.0d0**(SA_exp_log10ed))
  endif

  if (feminus_exp_log10ed.gt.30.0d0) then  !exp >> 1
     one_plus_feminus_exp_log10ed = feminus_exp_log10ed
  else if (feminus_exp_log10ed.lt.-30.0d0) then !exp << 1
     one_plus_feminus_exp_log10ed = 0.0d0
  else
     one_plus_feminus_exp_log10ed = log10(1.0d0+10.0d0**(feminus_exp_log10ed))
  endif

  !note, in principle we do this
  !       feminus_exp/(1.0d0+feminus_exp)*(1.0d0+SA_exp)/SA_exp * & !dimensionless
  !       exponential_term !dimensionless
  !but in practice, the expontential_term*feminus_exp/SA_exp = 1, therefore
  logterm = one_plus_SA_exp_log10ed - one_plus_feminus_exp_log10ed
  
  expterm = 10.0d0**(logterm)

  crosssection = sigma0 / 14.0d0 * gA**2 * & !cm^2
       N_p_Z * N_n_N * & ! dimensionless
       (final_electron_energy/m_e)**2 * & !dimensionless
       (1.0d0 - (m_e/final_electron_energy)**2)**0.5d0 * & !dimensionless
       expterm

end function nue_absorption_on_A


!neutrino species
!1 = electron neutrino
!2 = electron anti-neutrino
!3 = muon neutrino
!4 = muon anti-neutrino
!5 = tau neutrino
!6 = tau anti-neutrino
!7 = x neutrino = (3+4+5+6)
!8 = y neutrino = (3+5)
!9 = z anti-neutrino = (4+6)

subroutine total_absorption_opacities(neutrino_species,neutrino_energy,absorption_opacity,eos_variables)

  use nulib
  implicit none

  !inputs
  real*8, intent(in) :: eos_variables(total_eos_variables)  
  integer, intent(in) :: neutrino_species  !integer 1 through 9
  real*8, intent(in) :: neutrino_energy !MeV
  
  !outputs
  real*8, intent(out) :: absorption_opacity !cm^-1
  
  !local
  real*8 :: proton_number_density !proton number density, # neutron/cm^3
  real*8 :: neutron_number_density !neutron number density, # neutron/cm^3
  real*8 :: heavy_number_density !heaviy number density, # neutron/cm^3
  real*8 :: matter_muhat0 !matter_muhat0 is difference without rest
                          !masses included muhat =
                          !(mun+mn-mn)-(mup+mn-mp) = mun-mup-Delta_np
  real*8 :: mu_nu_eq !neutrino equilibrium chemical potential
  real*8 :: expterm,inverse_bottom !to calculate stimulated absorption

  !function declarations
  real*8 :: fermidirac !function declaration
  real*8 :: fermidirac_exptermonly !function declaration
  real*8 :: nue_absorption_on_n !function declaration
  real*8 :: anue_absorption_on_p !function declaration
  real*8 :: nue_absorption_on_A !function declaration
  real*8 :: nux_absorption_on_n_and_p !function declaration

  neutron_number_density = max(1.0d-100,eos_variables(xnindex))* &
       eos_variables(rhoindex)/(m_ref*mev_to_gram) !# neutrons/cm^3
  proton_number_density = max(1.0d-100,eos_variables(xpindex))* &
       eos_variables(rhoindex)/(m_ref*mev_to_gram) !# protons/cm^3
  if (eos_variables(abarindex).eq.0.0d0) then
     if (eos_variables(xhindex).gt.1.0d-10) write(*,*) "Warning!! matter_xh>0 but matter_abar=0.0"
     heavy_number_density = 0.0d0
  else
     heavy_number_density = eos_variables(xhindex)*eos_variables(rhoindex)/ &
          (eos_variables(abarindex)*m_ref*mev_to_gram) !# heavies/cm^3
  endif

  !matter_muhat0 is difference without rest masses included muhat =
  !(mun+mn-mn)-(mup+mn-mp) = mun-mup-Delta_np
  matter_muhat0 = eos_variables(muhatindex) - delta_np 
  absorption_opacity = 0.0d0

  !add in the electron neutrino absorption on free neutrons
  if (add_nue_absorption_on_n.and.neutrino_species.eq.1) then

     !some low density, low temperature regions give issues with nucleon blocking, so if rho<1e11, assume no nulceon blocking
     if (eos_variables(rhoindex).ge.1.0d11) then
        absorption_opacity = absorption_opacity + & !total opacity, dimensions cm^-1
             nue_absorption_on_n(neutrino_energy,eos_variables)* & !crosssection, cm^2
             neutron_number_density* &! # neutons/cm^3
             max(0.0d0,(proton_number_density/neutron_number_density-1.0d0)/ & !final state proton blocking, =1 if blocking is irrelevant
             (exp(-matter_muhat0/eos_variables(tempindex))-1.0d0)) !Bruenn 1985, Eq. C14
     else
        absorption_opacity = absorption_opacity + & !total opacity, dimensions cm^-1
             nue_absorption_on_n(neutrino_energy,eos_variables)* & !crosssection, cm^2
             neutron_number_density! # neutons/cm^3
     endif

     !just a check
     if (absorption_opacity.ne.absorption_opacity) then
        write(*,*) eos_variables(rhoindex),eos_variables(tempindex),eos_variables(yeindex)
        stop "NaN's in absorption_crosssections.F90: add_nue_absorption_on_n:"
     endif

  endif

  !add in the electron antineutrino absorption on free protons
  if (add_anue_absorption_on_p.and.neutrino_species.eq.2) then

     !use minus electron chemical potential, equals positron chemical potential
     !some low density, low temperature regions give issues with nucleon blocking, so if rho<1e11, assume no nulceon blocking
     if (eos_variables(rhoindex).ge.1.0d11) then
        absorption_opacity = absorption_opacity + & ! total opacity, dimensions cm^-1
             anue_absorption_on_p(neutrino_energy,eos_variables)* & !crosssection, cm^2
             proton_number_density* & ! # protons/cm^3
             max(0.0d0,(neutron_number_density/proton_number_density-1.0d0)/ & !final state neutron blocking, =1 if blocking is irrelevant
             (exp(matter_muhat0/eos_variables(tempindex))-1.0d0)) !Bruenn 1985, Eq. C14, with n and p switched
     else
        absorption_opacity = absorption_opacity + & ! total opacity, dimensions cm^-1
             anue_absorption_on_p(neutrino_energy,eos_variables)* & !crosssection, cm^2
             proton_number_density ! # protons/cm^3
     endif

     !just a check
     if (absorption_opacity.ne.absorption_opacity) then
        write(*,*) eos_variables(rhoindex),eos_variables(tempindex),eos_variables(yeindex)
        stop "NaN's in absorption_crosssections.F90: add_anue_absorption_on_p:"
     endif
  endif
  
  !add in the electron antineutrino absorption on heavy nuclei
  if (add_nue_absorption_on_A.and.neutrino_species.eq.1.and.heavy_number_density.ne.0.0d0) then

     absorption_opacity = absorption_opacity + &  !total opacity, dimensions cm^-1
          nue_absorption_on_A(neutrino_energy,eos_variables)* &  
          heavy_number_density ! # heavies/cm^3
  endif

  !add in the adhoc nux absorption on neutrons and protons
  if (add_nux_absorption_on_n_and_p.and.(neutrino_species.eq.3.or.neutrino_species.eq.4 &
       .or.neutrino_species.eq.5.or.neutrino_species.eq.6)) then

     absorption_opacity = absorption_opacity + &  !total opacity, dimensions cm^-1
          nux_absorption_on_n_and_p(neutrino_energy,eos_variables)* & !cm^2  
          (neutron_number_density+proton_number_density) ! # nucleons/cm^3

  endif

end subroutine total_absorption_opacities

subroutine return_absorption_opacity_spectra_given_neutrino_scheme(absorption_opacity_spectra,eos_variables)
  
  use nulib
  implicit none
  
  !inputs & outputs
  real*8, intent(in) :: eos_variables(total_eos_variables)  
  real*8, intent(out) :: absorption_opacity_spectra(number_species,number_groups)
  
  !locals
  integer :: ns,ng
  real*8 opacity
  
  if (size(absorption_opacity_spectra,1).ne.number_species) then
     write(*,*) size(absorption_opacity_spectra,1), number_species
     stop "return_absorption_opacity_spectra_given_neutrino_scheme:provided array has wrong number of species"
  endif
  if (size(absorption_opacity_spectra,2).ne.number_groups) then
     stop "return_absorption_opacity_spectra_given_neutrino_scheme:provided array has wrong number of groups"
  endif
  
  do ns=1,number_species
     do ng=1,number_groups
        call total_absorption_opacities(ns,energies(ng),opacity,eos_variables)
        absorption_opacity_spectra(ns,ng) = opacity
     enddo
  enddo
  
end subroutine return_absorption_opacity_spectra_given_neutrino_scheme
