!-*-f90-*-
function fermidirac(energy,mu,T) result(f)

  implicit none
  
  real*8, intent(in) :: energy !MeV
  real*8, intent(in) :: mu !includes rest mass, MeV
  real*8, intent(in) :: T !MeV
  real*8 :: f !dimensionless

  f = 1.0d0/(1.0d0+exp((energy-mu)/T))

end function fermidirac

function fermidirac_exptermonly(energy,mu,T) result(expterm)

  implicit none

  real*8, intent(in) :: energy !MeV
  real*8, intent(in) :: mu !includes rest mass, MeV
  real*8, intent(in) :: T !MeV
  real*8 :: expterm !dimensionless

  expterm = exp((energy-mu)/T)

end function fermidirac_exptermonly

function fermidirac_exptermonly_log10ed(energy,mu,T) result(expterm_log10ed)

  use nulib, only : log10exp
  implicit none

  real*8, intent(in) :: energy !MeV
  real*8, intent(in) :: mu !includes rest mass, MeV
  real*8, intent(in) :: T !MeV
  real*8 :: expterm_log10ed !dimensionless

  expterm_log10ed = (energy-mu)/T*log10exp

end function fermidirac_exptermonly_log10ed

function fermidirac_dimensionless(x,eta) result(f)

  implicit none
  
  real*8, intent(in) :: x !like energy/T
  real*8, intent(in) :: eta !like chemical potential/T
  real*8 :: f !dimensionless

  f = 1.0d0/(1.0d0+exp(x-eta))

end function fermidirac_dimensionless

function incomplete_fermidirac_integral(n,eta,a,b) result(incomplete_integral)

  use nulib, only : GPQ_n32_roots,GPQ_n32_weights
  implicit none
  
  !inputs
  real*8, intent(in) :: eta
  real*8, intent(in) :: a,b
  integer, intent(in) :: n

  !output
  real*8 :: incomplete_integral

  !function_declaration
  real*8 :: fermidirac_dimensionless

  !local
  integer :: i
  real*8 :: x

  incomplete_integral = 0.0d0
  do i=1,32
     x = (b-a)/2.0d0*GPQ_n32_roots(i) + (b+a)/2.0d0
     incomplete_integral = incomplete_integral + &
          GPQ_n32_weights(i)*x**n*fermidirac_dimensionless(x,eta)
  enddo
  incomplete_integral = (b-a)*incomplete_integral/2.0d0

end function incomplete_fermidirac_integral

subroutine return_blackbody_spectra(blackbody_spectra,eos_variables) 
  
  use nulib
  implicit none

  !inputs & outputs
  real*8, intent(in) :: eos_variables(total_eos_variables)
  real*8, intent(out) :: blackbody_spectra(3,number_groups)

  !local
  integer :: ns,ng,i
  real*8 :: energy_top_x,energy_bottom_x, bin_width_x! dimensionless
  real*8 :: energy_x, eta

  !function declaration
  real*8 :: fermidirac_dimensionless
  real*8 :: complete_fermi_integral

  eta = (eos_variables(mueindex)-eos_variables(muhatindex))/eos_variables(tempindex)
  blackbody_spectra = 0.0d0

  if (do_integrated_BB_and_emissivity) then
     do ng=1,number_groups
    
        energy_bottom_x = bin_bottom(ng)/eos_variables(tempindex)
        energy_top_x = bin_top(ng)/eos_variables(tempindex)
        bin_width_x = bin_widths(ng)/eos_variables(tempindex)
        
        do i=1,32
           energy_x = (energy_top_x-energy_bottom_x)/2.0d0*GPQ_n32_roots(i)+ &
                (energy_top_x+energy_bottom_x)/2.0d0
           blackbody_spectra(1,ng) = blackbody_spectra(1,ng) + &
                energy_x**3*fermidirac_dimensionless(energy_x,eta)*GPQ_n32_weights(i)
           blackbody_spectra(2,ng) = blackbody_spectra(2,ng) +  &
                energy_x**3*fermidirac_dimensionless(energy_x,-eta)*GPQ_n32_weights(i)
           blackbody_spectra(3,ng) = blackbody_spectra(3,ng) + &
                energy_x**3*fermidirac_dimensionless(energy_x,0.0d0)*GPQ_n32_weights(i)
           
        enddo
        
        blackbody_spectra(1,ng) = blackbody_spectra(1,ng)* &
             (energy_top_x-energy_bottom_x)/2.0d0/bin_width_x !dimensionless
        blackbody_spectra(2,ng) = blackbody_spectra(2,ng)* &
             (energy_top_x-energy_bottom_x)/2.0d0/bin_width_x !dimensionless
        blackbody_spectra(3,ng) = blackbody_spectra(3,ng)* &
             (energy_top_x-energy_bottom_x)/2.0d0/bin_width_x !dimensionless
        
     enddo

     blackbody_spectra = blackbody_spectra*eos_variables(tempindex)**3*clight/ &
          (2.0d0*pi*hbarc_mevcm)**3*mev_to_erg ! ergs/cm^2/s/MeV/srad

  else
     !this is the pointwise definition of the BB function, better
     !matching for the source term in MGFLD There can be a big
     !difference (i.e. factor of 4) between the integrate BB and the
     !pointwise BB, especially if the exponential regime
     do ng=1,number_groups
        energy_x = energies(ng)/eos_variables(tempindex)
        blackbody_spectra(1,ng) = clight*energies(ng)**3* &
             (mev_to_erg/(2.0d0*pi*hbarc_mevcm)**3)*fermidirac_dimensionless(energy_x,eta)
        blackbody_spectra(2,ng) = clight*energies(ng)**3* &
             (mev_to_erg/(2.0d0*pi*hbarc_mevcm)**3)*fermidirac_dimensionless(energy_x,-eta)
        blackbody_spectra(3,ng) = clight*energies(ng)**3* &
             (mev_to_erg/(2.0d0*pi*hbarc_mevcm)**3)*fermidirac_dimensionless(energy_x,0.0d0)
     enddo
     
  endif

end subroutine return_blackbody_spectra

  
!######################################################################
function complete_fermi_integral(ifermi,eta)
  implicit none
  integer ifermi
  real*8 complete_fermi_integral
  real*8 eta
  real*8 fermi_integral_analytical
  
  fermi_integral_analytical = 0.0d0
  
  ! Expressions for Fermi integrals given in Takahashi et al. 1978 
  if (eta.gt.1.D-3) then  
     select case (ifermi)
     case (0)
        fermi_integral_analytical = &
             log10(1.0d0+exp(eta))
     case (1)
        fermi_integral_analytical = &
             (eta**2/2.0D0 + 1.6449d0)/(1.0D0+EXP(-1.6855d0*eta))
     case (2)
        fermi_integral_analytical = &
             (eta**3/3.0D0 + 3.2899d0*eta)/(1.0D0-EXP(-1.8246d0*eta))
     case (3)
        fermi_integral_analytical = & 
             (eta**4/4.0D0 + 4.9348d0*eta**2+11.3644d0) / &
             (1.0D0+EXP(-1.9039d0*eta))        
     case (4)
        fermi_integral_analytical = &
             (eta**5/5.0D0 + 6.5797d0*eta**3+45.4576d0*eta) / &
             (1.0D0-EXP(-1.9484d0*eta))        
     case (5)
        fermi_integral_analytical = &
             (eta**6/6.0D0 + 8.2247d0*eta**4 + 113.6439d0*eta**2 + &
             236.5323d0)/(1.0D0+EXP(-1.9727d0*eta))
     end select
     
  else
     select case (ifermi)
     case (0)
        fermi_integral_analytical = &
             log10(1.0d0+exp(eta))
     case (1)
        fermi_integral_analytical = &
             EXP(eta)/(1.0D0+0.2159d0*EXP(0.8857d0*eta))
     case (2)
        fermi_integral_analytical = & 
             2.0D0*EXP(eta)/(1.0D0+0.1092d0*EXP(0.8908d0*eta))
     case (3)
        fermi_integral_analytical = & 
             6.0D0*EXP(eta)/(1.0D0+0.0559d0*EXP(0.9069d0*eta))
     case (4)
        fermi_integral_analytical = & 
             24.0D0*EXP(eta)/(1.0D0+0.0287d0*EXP(0.9257d0*eta))
     case (5)
        fermi_integral_analytical = &
             120.0D0*EXP(eta) / (1.0D0 + 0.0147d0*EXP(0.9431d0*eta))
     end select
     
  endif
  complete_fermi_integral = fermi_integral_analytical

  return
end function complete_fermi_integral
  
