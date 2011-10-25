!-*-f90-*-
module nulib

  implicit none

  !variables initialization
  integer :: neutrino_scheme
  !many people use different number of species, this is to denote how they are devided up.
  ! neutrino_scheme = 1 (three species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon+tau neutrino+antineutrino
  ! neutrino_scheme = 2 (four species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon+tau neutrino             #4 mu and tau antineutrino
  ! neutrino_scheme = 3 (six species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon neutrino                 #4 mu antineutrino
  !         #5: tau neutrino                  #6 tau antineutrino

  integer :: number_species
  integer :: number_groups

  real*8, allocatable,dimension(:) :: energies! MeV, energy of bin
  real*8, allocatable,dimension(:) :: bin_widths ! MeV, energy width of bin
  real*8, allocatable,dimension(:) :: bin_bottom ! MeV, energy at bottom of bin
  real*8, allocatable,dimension(:) :: bin_top ! MeV, energy at top of bin
  
  !EOS reference mass (Shen, m_amu, LS, m_n)
  real*8 :: m_ref !MeV

  !constants for routines
  real*8 :: H0_constants(2,6)
  real*8 :: Sion_coeffs(6,4)

  !Gauss Laguerre (L) and Legendre (P) weights and roots
  real*8 :: GLQ_n4_weights(4), GLQ_n4_roots(4) !weights and roots for n=4
  real*8 :: GLQ_n16_weights(16), GLQ_n16_roots(16) !weights and roots for n=16
  real*8 :: GLQ_n32_weights(32), GLQ_n32_roots(32) !weights and roots for n=32
  real*8 :: GLQ_n64_weights(64), GLQ_n64_roots(64) !weights and roots for n=64

  real*8 :: GPQ_n4_weights(4), GPQ_n4_roots(4) !weights and roots for n=4
  real*8 :: GPQ_n16_weights(16), GPQ_n16_roots(16) !weights and roots for n=16
  real*8 :: GPQ_n32_weights(32), GPQ_n32_roots(32) !weights and roots for n=32
  real*8 :: GPQ_n64_weights(64), GPQ_n64_roots(64) !weights and roots for n=64

  !EOS variables index holders, we carry this around with us to
  !instead of globally setting it for easy parallization
  integer :: total_eos_variables = 14
  integer :: rhoindex =1
  integer :: tempindex = 2
  integer :: yeindex = 3
  integer :: energyindex = 4
  integer :: xaindex = 5
  integer :: xhindex = 6
  integer :: xnindex = 7
  integer :: xpindex = 8
  integer :: abarindex = 9
  integer :: zbarindex = 10
  integer :: mueindex = 11
  integer :: munindex = 12
  integer :: mupindex = 13
  integer :: muhatindex = 14

  logical :: debug = .false.

  include 'constants.inc'
  include 'requested_interactions.inc'

  contains
!#################################################################
    subroutine initialize_nulib(neutrino_scheme_in,number_species_in,number_groups_in)
      
      integer :: neutrino_scheme_in
      integer :: number_species_in
      integer :: number_groups_in

      integer :: i
      real*8 dxfac,mindx

      neutrino_scheme = neutrino_scheme_in
      number_species = number_species_in
      number_groups = number_groups_in

      !double check we know what we are doing
      if (neutrino_scheme.eq.1) then
         if (number_species.eq.6) then
            write(*,*) "You choose to use neutrino scheme #",neutrino_scheme, &
                 ".  This will average nu_mu, anu_mu, nu_tau, anu_tau together, for 3 independant species"
         else
            write(*,*) "nuLib can only handle 6 species, the scheme may average after, depending on your neutrino scheme"
            stop "initialize: please adjust neutrino species/neutrio scheme"
         endif
      else if (neutrino_scheme.eq.2) then
         if (number_species.eq.6) then
            write(*,*) "You choose to use neutrino scheme #",neutrino_scheme, &
                 ".  This will average nu_mu and nu_tau, and also average anu_mu, anu_tau, for 4 independant species"
         else
            write(*,*) "nuLib can only handle 6 species, the scheme may average after, depending on your neutrino scheme"
            stop "initialize: please adjust neutrino species/neutrio scheme"
         endif
      else if (neutrino_scheme.eq.3) then
         if (number_species.eq.6) then
            write(*,*) "You choose to use neutrino scheme #",neutrino_scheme, &
                 ".  This will maintain 6 independant species"
         else
            write(*,*) "nuLib can only handle 6 species, the scheme may average after, depending on your neutrino scheme"
            stop "initialize: please adjust neutrino species/neutrio scheme"
         endif
      else
         stop "initialize_nulib: neutrino scheme unknown"
      endif

      !allocate some arrays
      allocate(energies(number_groups))
      allocate(bin_widths(number_groups))
      allocate(bin_bottom(number_groups))
      allocate(bin_top(number_groups))

      !set up energies bins
      mindx = 1.0d0
      bin_bottom(1) = 0.0d0 !MeV
      bin_bottom(2) = 4.0d0 !MeV
      bin_bottom(3) = bin_bottom(2)+mindx
      bin_bottom(number_groups) = 250.0d0

      call nulib_series2(number_groups-1,bin_bottom(2),bin_bottom(number_groups),mindx,dxfac)
      do i=4,number_groups
         bin_bottom(i) = bin_bottom(i-1)+(bin_bottom(i-1)-bin_bottom(i-2))*dxfac
      enddo

      !calculate bin widths & energies from the bottom of the bin & energy at top on bin
      do i=1,number_groups-1
         energies(i) = (bin_bottom(i)+bin_bottom(i+1))/2.0d0
         bin_widths(i) = bin_bottom(i+1)-bin_bottom(i)
         bin_top(i) = bin_bottom(i+1)
      enddo
      energies(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups-1)*dxfac/2.0d0
      bin_widths(number_groups) = 2.0*(energies(number_groups)-bin_bottom(number_groups))
      bin_top(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups)

      !setup H0_constants for electron positron annihilation
      !standard
      H0_constants(1,1) =  (0.5d0+2.0d0*sin2thetaW + 0.5d0)**2 !(Cv+Ca)^2, electron neutrino
      H0_constants(2,1) =  (0.5d0+2.0d0*sin2thetaW - 0.5d0)**2 !(Cv-Ca)^2, electron neutrino

      !Ca flips sign from standard
      H0_constants(1,2) =  (0.5d0+2.0d0*sin2thetaW + (-0.5d0))**2 !(Cv+Ca)^2, electron antineutrino
      H0_constants(2,2) =  (0.5d0+2.0d0*sin2thetaW - (-0.5d0))**2 !(Cv-Ca)^2, electron antineutrino

      !Cv -> Cv-1, Ca -> Ca-1  from standard
      H0_constants(1,3) =  ((0.5d0+2.0d0*sin2thetaW - 1.0d0) + (0.5d0 - 1.0d0))**2 !(Cv+Ca)^2, mu neutrino
      H0_constants(2,3) =  ((0.5d0+2.0d0*sin2thetaW - 1.0d0) - (0.5d0 - 1.0d0))**2 !(Cv-Ca)^2, mu neutrino

      !Cv -> Cv-1, Ca -> Ca-1  from standard
      H0_constants(1,5) =  ((0.5d0+2.0d0*sin2thetaW - 1.0d0) + (0.5d0 - 1.0d0))**2 !(Cv+Ca)^2, tau neutrino
      H0_constants(2,5) =  ((0.5d0+2.0d0*sin2thetaW - 1.0d0) - (0.5d0 - 1.0d0))**2 !(Cv-Ca)^2, tau neutrino

      !Cv -> Cv-1, Ca -> Ca-1 , then Ca-1 flips sign 
      H0_constants(1,4) =  ((0.5d0+2.0d0*sin2thetaW - 1.0d0) + (-(0.5d0 - 1.0d0)))**2 !(Cv+Ca)^2, mu antineutrino
      H0_constants(2,4) =  ((0.5d0+2.0d0*sin2thetaW - 1.0d0) - (-(0.5d0 - 1.0d0)))**2 !(Cv-Ca)^2, mu antineutrino

      !Cv -> Cv-1, Ca -> Ca-1 , then Ca-1 flips sign 
      H0_constants(1,6) =  ((0.5d0+2.0d0*sin2thetaW - 1.0d0) + (-(0.5d0 - 1.0d0)))**2 !(Cv+Ca)^2, tau antineutrino
      H0_constants(2,6) =  ((0.5d0+2.0d0*sin2thetaW - 1.0d0) - (-(0.5d0 - 1.0d0)))**2 !(Cv-Ca)^2, tau antineutrino

      !setup ion-ion correlation cooefficients for Chuck's fits
      Sion_coeffs = 0.0d0

      Sion_coeffs(3,1) = -7.362056d0
      Sion_coeffs(3,2) = 0.5371365d0
      Sion_coeffs(3,3) = -0.1078845d0
      Sion_coeffs(3,4) = 4.189612d-3
      
      Sion_coeffs(4,1) = 3.4489581d0
      Sion_coeffs(4,2) = -0.40251656d0
      Sion_coeffs(4,3) = 9.0877878d-2
      Sion_coeffs(4,4) = -3.4353581d-3

      Sion_coeffs(5,1) =  -0.74128645d0
      Sion_coeffs(5,2) =  0.11019855d0
      Sion_coeffs(5,3) =  -2.5359361d-2
      Sion_coeffs(5,4) =  9.0487744d-4

      Sion_coeffs(6,1) =  5.9573285d-2
      Sion_coeffs(6,2) =  -1.0186552d-2
      Sion_coeffs(6,3) =  2.2791369d-3
      Sion_coeffs(6,4) =  -7.4614597d-5  

      !setup Gauss-Laguerre weights and roots
      call GaussLaguerreQuadrature_roots_and_weights(4,GLQ_n4_roots,GLQ_n4_weights)
      call GaussLaguerreQuadrature_roots_and_weights(16,GLQ_n16_roots,GLQ_n16_weights)
      call GaussLaguerreQuadrature_roots_and_weights(32,GLQ_n32_roots,GLQ_n32_weights)
      call GaussLaguerreQuadrature_roots_and_weights(64,GLQ_n64_roots,GLQ_n64_weights)

      call GaussLegendreQuadrature_weights_and_roots(4,GPQ_n4_roots,GPQ_n4_weights)
      call GaussLegendreQuadrature_weights_and_roots(16,GPQ_n16_roots,GPQ_n16_weights)
      call GaussLegendreQuadrature_weights_and_roots(32,GPQ_n32_roots,GPQ_n32_weights)
      call GaussLegendreQuadrature_weights_and_roots(64,GPQ_n64_roots,GPQ_n64_weights)

    end subroutine initialize_nulib

    subroutine single_point_return_all(eos_variables, &
         emissivities,absorption_opacity,scattering_opacity,neutrino_local_scheme)
      
      implicit none

      !inputs & outputs

      real*8, intent(in) :: eos_variables(total_eos_variables)
      integer, intent(in) :: neutrino_local_scheme

      !neutrino species
      !1 = electron neutrino !local scheme 1,2,3
      !2 = electron anti-neutrino !scheme 1,2,3
      !3 = muon neutrino !scheme 1,2,3
      !4 = muon anti-neutrino !scheme 2,3
      !5 = tau neutrino !scheme 3
      !6 = tau anti-neutrino !scheme 3

      !x neutrino = (3+4+5+6) !scheme 1
      !y neutrino = (3+5) !scheme 2
      !z anti-neutrino = (4+6) !scheme 2

      real*8, dimension(:,:), intent(out) :: emissivities
      real*8, dimension(:,:), intent(out) :: absorption_opacity
      real*8, dimension(:,:), intent(out) :: scattering_opacity

      !local
      real*8 :: temporary_spectra(number_species,number_groups)
      real*8 :: blackbody_spectra(3,number_groups)
      integer :: number_local_species

      integer ns,ng,i

      if (neutrino_local_scheme.ne.neutrino_scheme) then
         write(*,*) neutrino_local_scheme, neutrino_scheme
         stop "you are requesting different schemes"
      endif

      if (neutrino_scheme.eq.1) then
         number_local_species = 3
      else if (neutrino_scheme.eq.2) then
         number_local_species = 4
      else if (neutrino_scheme.eq.3) then
         number_local_species = 6
      else
         stop "single_point_return_all:incorrect neutrino scheme"
      endif

      if (size(emissivities,1).ne.number_local_species) then
         stop "single_point_return_all:provided array has wrong number of species"
      endif
      if (size(emissivities,2).ne.number_groups) then
         stop "single_point_return_all:provided array has wrong number of groups"
      endif

      if (size(absorption_opacity,1).ne.number_local_species) then
         stop "single_point_return_all:provided array has wrong number of species"
      endif
      if (size(absorption_opacity,2).ne.number_groups) then
         stop "single_point_return_all:provided array has wrong number of groups"
      endif

      if (size(scattering_opacity,1).ne.number_local_species) then
         stop "single_point_return_all:provided array has wrong number of species"
      endif
      if (size(scattering_opacity,2).ne.number_groups) then
         stop "single_point_return_all:provided array has wrong number of groups"
      endif

      emissivities = 0.0d0
      absorption_opacity = 0.0d0
      scattering_opacity = 0.0d0

      !first get black body spectra, ergs/cm^2/s/MeV/srad
      call return_blackbody_spectra(blackbody_spectra,eos_variables)

      !then get absorption opacities
      temporary_spectra = 0.0d0
      call return_absorption_opacity_spectra_given_neutrino_scheme(temporary_spectra,eos_variables)

      !to fix minimum here to ensure emissivity agrees
      temporary_spectra = max(temporary_spectra,1.0d-30)
      absorption_opacity(1,1:number_groups) = absorption_opacity(1,1:number_groups) + &
           temporary_spectra(1,1:number_groups)
      absorption_opacity(2,1:number_groups) = absorption_opacity(2,1:number_groups) + &
           temporary_spectra(2,1:number_groups)
      emissivities(1,1:number_groups) = emissivities(1,1:number_groups) + &
           temporary_spectra(1,1:number_groups)*blackbody_spectra(1,1:number_groups)
      emissivities(2,1:number_groups) = emissivities(2,1:number_groups) + &
           temporary_spectra(2,1:number_groups)*blackbody_spectra(2,1:number_groups)

      if (number_local_species.eq.3) then
         absorption_opacity(3,1:number_groups) = absorption_opacity(3,1:number_groups) + &
              (temporary_spectra(3,1:number_groups) + temporary_spectra(4,1:number_groups) + &
              temporary_spectra(5,1:number_groups) + temporary_spectra(6,1:number_groups))/4.0d0
         emissivities(3,1:number_groups) = emissivities(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups)*blackbody_spectra(3,1:number_groups) + & 
              temporary_spectra(4,1:number_groups)*blackbody_spectra(3,1:number_groups) + & 
              temporary_spectra(5,1:number_groups)*blackbody_spectra(3,1:number_groups) + &
              temporary_spectra(6,1:number_groups)*blackbody_spectra(3,1:number_groups)

      else if (number_local_species.eq.4) then
         absorption_opacity(3,1:number_groups) = absorption_opacity(3,1:number_groups) + &
              (temporary_spectra(3,1:number_groups)+temporary_spectra(5,1:number_groups))/2.0d0
         absorption_opacity(4,1:number_groups) = absorption_opacity(4,1:number_groups) + &
              (temporary_spectra(4,1:number_groups)+temporary_spectra(6,1:number_groups))/2.0d0
         emissivities(3,1:number_groups) = emissivities(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups)*blackbody_spectra(3,1:number_groups) + &
              temporary_spectra(5,1:number_groups)*blackbody_spectra(3,1:number_groups)
         emissivities(4,1:number_groups) = emissivities(4,1:number_groups) + &
              temporary_spectra(4,1:number_groups)*blackbody_spectra(3,1:number_groups) + &
              temporary_spectra(6,1:number_groups)*blackbody_spectra(3,1:number_groups)

      else if (number_local_species.eq.6) then
         absorption_opacity(3,1:number_groups) = absorption_opacity(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups)
         absorption_opacity(4,1:number_groups) = absorption_opacity(4,1:number_groups) + &
              temporary_spectra(4,1:number_groups)
         absorption_opacity(5,1:number_groups) = absorption_opacity(5,1:number_groups) + &
              temporary_spectra(5,1:number_groups)
         absorption_opacity(6,1:number_groups) = absorption_opacity(6,1:number_groups) + &
              temporary_spectra(6,1:number_groups)
         emissivities(3,1:number_groups) = emissivities(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups)*blackbody_spectra(3,1:number_groups)
         emissivities(4,1:number_groups) = emissivities(4,1:number_groups) + &
              temporary_spectra(4,1:number_groups)*blackbody_spectra(3,1:number_groups)
         emissivities(5,1:number_groups) = emissivities(5,1:number_groups) + &
              temporary_spectra(5,1:number_groups)*blackbody_spectra(3,1:number_groups)
         emissivities(6,1:number_groups) = emissivities(6,1:number_groups) + &
              temporary_spectra(6,1:number_groups)*blackbody_spectra(3,1:number_groups)

      endif
      
      if(debug) then
         write(*,*) "debug #2: emissivities and absorption opacities from \kappa for species 1:", &
              absorption_opacity(1,1:number_groups),"em",emissivities(1,1:number_groups)
         write(*,*) "debug #2: emissivities and absorption opacities from \kappa for species 2:", &
              absorption_opacity(2,1:number_groups),"em",emissivities(2,1:number_groups)
         write(*,*) "debug #2: emissivities and absorption opacities from \kappa for species 3:", &
              absorption_opacity(3,1:number_groups),"em",emissivities(3,1:number_groups)

      endif

      !now get scattering opacites
      temporary_spectra = 0.0d0
      call return_scattering_opacity_spectra_given_neutrino_scheme(temporary_spectra,eos_variables)

      !to fix minimum, 1.0e-30 in units of cm^-1 is effectivly zero
      temporary_spectra = max(temporary_spectra,1.0d-30)

      scattering_opacity(1,1:number_groups) = scattering_opacity(1,1:number_groups) + &
           temporary_spectra(1,1:number_groups)
      scattering_opacity(2,1:number_groups) = scattering_opacity(2,1:number_groups) + &
           temporary_spectra(2,1:number_groups)
      
      if (number_local_species.eq.3) then
         scattering_opacity(3,1:number_groups) = scattering_opacity(3,1:number_groups) + &
              (temporary_spectra(3,1:number_groups) + temporary_spectra(4,1:number_groups) + &
              temporary_spectra(5,1:number_groups) + temporary_spectra(6,1:number_groups))/4.0d0

      else if (number_local_species.eq.4) then
         scattering_opacity(3,1:number_groups) = scattering_opacity(3,1:number_groups) + &
              (temporary_spectra(3,1:number_groups)+temporary_spectra(5,1:number_groups))/2.0d0
         scattering_opacity(4,1:number_groups) = scattering_opacity(4,1:number_groups) + &
              (temporary_spectra(4,1:number_groups)+temporary_spectra(6,1:number_groups))/2.0d0

      else if (number_local_species.eq.6) then
         scattering_opacity(3,1:number_groups) = scattering_opacity(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups)
         scattering_opacity(4,1:number_groups) = scattering_opacity(4,1:number_groups) + &
              temporary_spectra(4,1:number_groups)
         scattering_opacity(5,1:number_groups) = scattering_opacity(5,1:number_groups) + &
              temporary_spectra(5,1:number_groups)
         scattering_opacity(6,1:number_groups) = scattering_opacity(6,1:number_groups) + &
              temporary_spectra(6,1:number_groups)
      endif

      temporary_spectra = 0.0d0
      call return_emissivity_spectra_given_neutrino_scheme(temporary_spectra,eos_variables)

      if(debug) then
         write(*,*) "debug #3: emissivities and black bodies for species 1:", &
              temporary_spectra(1,1:number_groups),"bb",blackbody_spectra(1,1:number_groups)
         write(*,*) "debug #3: emissivities and black bodies for species 2:", &
              temporary_spectra(2,1:number_groups),"bb",blackbody_spectra(2,1:number_groups)
         write(*,*) "debug #3: emissivities and black bodies for species 3:", &
              temporary_spectra(3,1:number_groups),"bb",blackbody_spectra(3,1:number_groups)

      endif
     
      !in units of the BBS, 1.0d-30 is ergs/cm^2/s/MeV/srad, i.e. 10km
      !object would radiate in a 1MeV bin 1e-17 ergs/s (1e-26 less than
      !a light bulb...).  This limits the BB to a finite value, with
      !the expontials in the definition it can get out of control and
      !be 0.0d0 which is not good for the division in the next line
      blackbody_spectra(:,:) = max(1.0d-30,blackbody_spectra(:,:))
     
      absorption_opacity(1,1:number_groups) = absorption_opacity(1,1:number_groups) + &
           temporary_spectra(1,1:number_groups)/blackbody_spectra(1,1:number_groups)
      absorption_opacity(2,1:number_groups) = absorption_opacity(2,1:number_groups) + &
           temporary_spectra(2,1:number_groups)/blackbody_spectra(2,1:number_groups)
      emissivities(1,1:number_groups) = emissivities(1,1:number_groups) + &
           temporary_spectra(1,1:number_groups)
      emissivities(2,1:number_groups) = emissivities(2,1:number_groups) + &
           temporary_spectra(2,1:number_groups)

      if(debug) then
         write(*,*) "debug #4: absorption opacities after \eta/B_\nu term for species 1:", &
              absorption_opacity(1,1:number_groups)
         write(*,*) "debug #4: absorption opacities after \eta/B_\nu term for species 2:", &
              absorption_opacity(2,1:number_groups)
         write(*,*) "debug #4: absorption opacities after \eta/B_\nu term for species 3:", &
              absorption_opacity(3,1:number_groups)
      endif
      if (number_local_species.eq.3) then
         absorption_opacity(3,1:number_groups) = absorption_opacity(3,1:number_groups) + &
              (temporary_spectra(3,1:number_groups)/blackbody_spectra(3,1:number_groups) + &
              temporary_spectra(4,1:number_groups)/blackbody_spectra(3,1:number_groups) + & 
              temporary_spectra(5,1:number_groups)/blackbody_spectra(3,1:number_groups) + &
              temporary_spectra(6,1:number_groups)/blackbody_spectra(3,1:number_groups))/4.0d0
         emissivities(3,1:number_groups) = emissivities(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups) + temporary_spectra(4,1:number_groups) + &
              temporary_spectra(5,1:number_groups) + temporary_spectra(6,1:number_groups)

      !average neutrinos and antineutrinos individually
      else if (number_local_species.eq.4) then
         absorption_opacity(3,1:number_groups) = absorption_opacity(3,1:number_groups) + &
              (temporary_spectra(3,1:number_groups)/blackbody_spectra(3,1:number_groups) + &
              temporary_spectra(5,1:number_groups)/blackbody_spectra(3,1:number_groups))/2.0d0
         absorption_opacity(4,1:number_groups) = absorption_opacity(4,1:number_groups) + &
              (temporary_spectra(4,1:number_groups)/blackbody_spectra(3,1:number_groups) + &
              temporary_spectra(6,1:number_groups)/blackbody_spectra(3,1:number_groups))/2.0d0
         emissivities(3,1:number_groups) = emissivities(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups) + temporary_spectra(5,1:number_groups)
         emissivities(4,1:number_groups) = emissivities(4,1:number_groups) + &
              temporary_spectra(4,1:number_groups) + temporary_spectra(6,1:number_groups)

      !no averaging at all, what six different species
      else if (number_local_species.eq.6) then
         absorption_opacity(3,1:number_groups) = absorption_opacity(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups)/blackbody_spectra(3,1:number_groups)
         absorption_opacity(4,1:number_groups) = absorption_opacity(4,1:number_groups) + &
              temporary_spectra(4,1:number_groups)/blackbody_spectra(3,1:number_groups)
         absorption_opacity(5,1:number_groups) = absorption_opacity(5,1:number_groups) + &
              temporary_spectra(5,1:number_groups)/blackbody_spectra(3,1:number_groups)
         absorption_opacity(6,1:number_groups) = absorption_opacity(6,1:number_groups) + &
              temporary_spectra(6,1:number_groups)/blackbody_spectra(3,1:number_groups)
         emissivities(3,1:number_groups) = emissivities(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups)
         emissivities(4,1:number_groups) = emissivities(4,1:number_groups) + &
              temporary_spectra(4,1:number_groups)
         emissivities(5,1:number_groups) = emissivities(5,1:number_groups) + &
              temporary_spectra(5,1:number_groups)
         emissivities(6,1:number_groups) = emissivities(6,1:number_groups) + &
              temporary_spectra(6,1:number_groups)

      endif

      !check for NaNs
      do ng=1,number_groups
         do i=1,number_local_species
            if(emissivities(i,ng).ne.emissivities(i,ng)) stop "NaN in emissivities"
            if(absorption_opacity(i,ng).ne.absorption_opacity(i,ng)) stop "NaN in absorption_opacity"
            if(scattering_opacity(i,ng).ne.scattering_opacity(i,ng)) stop "NaN in scattering_opacity"
         enddo      
      enddo

      emissivities(:,:) = max(1.0d-300,emissivities(:,:)) !ergs/cm^3/s/MeV/srad
      absorption_opacity(:,:) = max(1.0d-300,absorption_opacity(:,:)) !cm^-1
      scattering_opacity(:,:) = max(1.0d-300,scattering_opacity(:,:)) !cm^-1


    end subroutine single_point_return_all

 end module nulib


