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

  !EOS table
  character*200 :: eos_filename

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
  real*8 :: GLQ_n128_weights(128), GLQ_n128_roots(128) !weights and roots for n=128

  real*8 :: GPQ_n4_weights(4), GPQ_n4_roots(4) !weights and roots for n=4
  real*8 :: GPQ_n16_weights(16), GPQ_n16_roots(16) !weights and roots for n=16
  real*8 :: GPQ_n32_weights(32), GPQ_n32_roots(32) !weights and roots for n=32
  real*8 :: GPQ_n64_weights(64), GPQ_n64_roots(64) !weights and roots for n=64
  real*8 :: GPQ_n128_weights(128), GPQ_n128_roots(128) !weights and roots for n=128

  !EOS variables index holders, we carry this around with us to
  !instead of globally setting it for easy parallization
  integer :: total_eos_variables = 15
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
  integer :: entropyindex = 15

#ifdef DEBUG
  logical :: debug = .true.
#else
  logical :: debug = .false.
#endif
  logical :: do_integrated_BB_and_emissivity

  !tabulated weak rate table bounds
  real*8, dimension(4) :: table_bounds

  !special terms
  real*8 :: adhoc_nux_factor = 0.0d0

  include 'constants.inc'
  include 'requested_interactions.inc'

  contains
!#################################################################
    subroutine initialize_nulib(neutrino_scheme_in,number_species_in,number_groups_in)

      integer :: neutrino_scheme_in
      integer :: number_species_in
      integer :: number_groups_in

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

      if (add_nue_Iscattering_electrons.and.add_nue_scattering_electrons) then
         stop "initialize_nulib: you can't both inelastically and elastically scatter off of electrons: nue"
      endif
      if (add_anue_Iscattering_electrons.and.add_anue_scattering_electrons) then
         stop "initialize_nulib: you can't both inelastically and elastically scatter off of electrons: anue"
      endif
      if (add_numu_Iscattering_electrons.and.add_numu_scattering_electrons) then
         stop "initialize_nulib: you can't both inelastically and elastically scatter off of electrons: numu"
      endif
      if (add_anumu_Iscattering_electrons.and.add_anumu_scattering_electrons) then
         stop "initialize_nulib: you can't both inelastically and elastically scatter off of electrons: anumu"
      endif
      if (add_nutau_Iscattering_electrons.and.add_nutau_scattering_electrons) then
         stop "initialize_nulib: you can't both inelastically and elastically scatter off of electrons: nutau"
      endif
      if (add_anutau_Iscattering_electrons.and.add_anutau_scattering_electrons) then
         stop "initialize_nulib: you can't both inelastically and elastically scatter off of electrons: anutau"
      endif
      if (add_nue_kernel_epannihil.and.add_nue_emission_epannihil) then
         stop "initialize_nulib: you can't both thermally emit nu with/without a kernel: nue"
      endif
      if (add_anue_kernel_epannihil.and.add_anue_emission_epannihil) then
         stop "initialize_nulib: you can't both thermally emit nu with/without a kernel: anue"
      endif
      if (add_numu_kernel_epannihil.and.add_numu_emission_epannihil) then
         stop "initialize_nulib: you can't both thermally emit nu with/without a kernel: numu"
      endif
      if (add_anumu_kernel_epannihil.and.add_anumu_emission_epannihil) then
         stop "initialize_nulib: you can't both thermally emit nu with/without a kernel: anumu"
      endif
      if (add_nutau_kernel_epannihil.and.add_nutau_emission_epannihil) then
         stop "initialize_nulib: you can't both thermally emit nu with/without a kernel: nutau"
      endif
      if (add_anutau_kernel_epannihil.and.add_anutau_emission_epannihil) then
         stop "initialize_nulib: you can't both thermally emit nu with/without a kernel: anutau"
      endif

      if (neutrino_scheme.eq.1) then
         if ((add_numu_Iscattering_electrons.neqv.add_anumu_Iscattering_electrons).or. &
              (add_numu_Iscattering_electrons.neqv.add_nutau_Iscattering_electrons).or. &
              (add_numu_Iscattering_electrons.neqv.add_anutau_Iscattering_electrons)) then
            stop "With neutrino scheme 1, all 4 heavy lepton inelastic scattering must be the same"
         endif
         if ((add_numu_scattering_n.neqv.add_anumu_scattering_n).or. &
              (add_numu_scattering_n.neqv.add_nutau_scattering_n).or. &
              (add_numu_scattering_n.neqv.add_anutau_scattering_n)) then
            stop "With neutrino scheme 1, all 4 heavy lepton n-scattering must be the same"
         endif
         if ((add_numu_scattering_p.neqv.add_anumu_scattering_p).or. &
              (add_numu_scattering_p.neqv.add_nutau_scattering_p).or. &
              (add_numu_scattering_p.neqv.add_anutau_scattering_p)) then
            stop "With neutrino scheme 1, all 4 heavy lepton p-scattering must be the same"
         endif
         if ((add_numu_scattering_heavies.neqv.add_anumu_scattering_heavies).or. &
              (add_numu_scattering_heavies.neqv.add_nutau_scattering_heavies).or. &
              (add_numu_scattering_heavies.neqv.add_anutau_scattering_heavies)) then
            stop "With neutrino scheme 1, all 4 heavy lepton havy-scattering must be the same"
         endif
         if ((add_numu_scattering_electrons.neqv.add_anumu_scattering_electrons).or. &
              (add_numu_scattering_electrons.neqv.add_nutau_scattering_electrons).or. &
              (add_numu_scattering_electrons.neqv.add_anutau_scattering_electrons)) then
            stop "With neutrino scheme 1, all 4 heavy lepton electron-scattering must be the same"
         endif
         if ((add_numu_scattering_alphas.neqv.add_anumu_scattering_alphas).or. &
              (add_numu_scattering_alphas.neqv.add_nutau_scattering_alphas).or. &
              (add_numu_scattering_alphas.neqv.add_anutau_scattering_alphas)) then
            stop "With neutrino scheme 1, all 4 heavy lepton alpha-scattering must be the same"
         endif
         if ((add_numu_emission_epannihil.neqv.add_anumu_emission_epannihil).or. &
              (add_numu_emission_epannihil.neqv.add_nutau_emission_epannihil).or. &
              (add_numu_emission_epannihil.neqv.add_anutau_emission_epannihil)) then
            stop "With neutrino scheme 1, all 4 heavy lepton epannihils must be the same"
         endif
         if ((add_numu_emission_NNBrems.neqv.add_anumu_emission_NNBrems).or. &
              (add_numu_emission_NNBrems.neqv.add_nutau_emission_NNBrems).or. &
              (add_numu_emission_NNBrems.neqv.add_anutau_emission_NNBrems)) then
            stop "With neutrino scheme 1, all 4 heavy lepton NBrems must be the same"
         endif
         if ((add_numu_kernel_epannihil.neqv.add_anumu_kernel_epannihil).or. &
              (add_numu_kernel_epannihil.neqv.add_nutau_kernel_epannihil).or. &
              (add_numu_kernel_epannihil.neqv.add_anutau_kernel_epannihil)) then
            stop "With neutrino scheme 1, all 4 heavy lepton kernel epannihils must be the same"
         endif
      else if (neutrino_scheme.eq.2) then
         if ((add_numu_Iscattering_electrons.neqv.add_nutau_Iscattering_electrons).or. &
              (add_anumu_Iscattering_electrons.neqv.add_anutau_Iscattering_electrons)) then
            stop "With neutrino scheme 2, each heavy lepton nu/anu inelastic scattering must be the same"
         endif
         if ((add_numu_scattering_n.neqv.add_nutau_scattering_n).or. &
              (add_anumu_scattering_n.neqv.add_anutau_scattering_n)) then
            stop "With neutrino scheme 2, each heavy lepton nu/anu n-scattering must be the same"
         endif
         if ((add_numu_scattering_p.neqv.add_nutau_scattering_p).or. &
              (add_anumu_scattering_p.neqv.add_anutau_scattering_p)) then
            stop "With neutrino scheme 2, each heavy lepton nu/anu p-scattering must be the same"
         endif
         if ((add_numu_scattering_heavies.neqv.add_nutau_scattering_heavies).or. &
              (add_anumu_scattering_heavies.neqv.add_anutau_scattering_heavies)) then
            stop "With neutrino scheme 2, each heavy lepton nu/anu havy-scattering must be the same"
         endif
         if ((add_numu_scattering_electrons.neqv.add_nutau_scattering_electrons).or. &
              (add_anumu_scattering_electrons.neqv.add_anutau_scattering_electrons)) then
            stop "With neutrino scheme 2, each heavy lepton nu/anu electron-scattering must be the same"
         endif
         if ((add_numu_scattering_alphas.neqv.add_nutau_scattering_alphas).or. &
              (add_anumu_scattering_alphas.neqv.add_anutau_scattering_alphas)) then
            stop "With neutrino scheme 2, each heavy lepton nu/anu alpha-scattering must be the same"
         endif
         if ((add_numu_emission_epannihil.neqv.add_nutau_emission_epannihil).or. &
              (add_anumu_emission_epannihil.neqv.add_anutau_emission_epannihil)) then
            stop "With neutrino scheme 2, each heavy lepton nu/anu epannihils must be the same"
         endif
         if ((add_numu_emission_NNBrems.neqv.add_nutau_emission_NNBrems).or. &
              (add_anumu_emission_NNBrems.neqv.add_anutau_emission_NNBrems)) then
            stop "With neutrino scheme 2, each heavy lepton nu/anu NBrems must be the same"
         endif
         if ((add_numu_kernel_epannihil.neqv.add_nutau_kernel_epannihil).or. &
              (add_anumu_kernel_epannihil.neqv.add_anutau_kernel_epannihil)) then
            stop "With neutrino scheme 2, each heavy lepton nu/anu kernel epannihils must be the same"
         endif
      endif

      !allocate some arrays
      allocate(energies(number_groups))
      allocate(bin_widths(number_groups))
      allocate(bin_bottom(number_groups))
      allocate(bin_top(number_groups))

      !setup H0_constants for electron positron annihilation and inelastic electron scattering, also constants for the H1 term...
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
      call GaussLaguerreQuadrature_roots_and_weights(128,GLQ_n128_roots,GLQ_n128_weights)

      call GaussLegendreQuadrature_weights_and_roots(4,GPQ_n4_roots,GPQ_n4_weights)
      call GaussLegendreQuadrature_weights_and_roots(16,GPQ_n16_roots,GPQ_n16_weights)
      call GaussLegendreQuadrature_weights_and_roots(32,GPQ_n32_roots,GPQ_n32_weights)
      call GaussLegendreQuadrature_weights_and_roots(64,GPQ_n64_roots,GPQ_n64_weights)
      call GaussLegendreQuadrature_weights_and_roots(128,GPQ_n128_roots,GPQ_n128_weights)

    end subroutine initialize_nulib

    subroutine single_point_return_all(eos_variables, &
         emissivities,absorption_opacity,scattering_opacity,delta,neutrino_local_scheme)
      
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
      real*8, dimension(:,:), intent(out) :: delta

      !local
      real*8 :: temporary_spectra(number_species,number_groups)
      real*8 :: temporary_delta(number_species,number_groups)
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

      if (size(delta,1).ne.number_local_species) then
         stop "single_point_return_all:provided array has wrong number of species"
      endif
      if (size(delta,2).ne.number_groups) then
         stop "single_point_return_all:provided array has wrong number of groups"
      endif

      emissivities = 0.0d0
      absorption_opacity = 0.0d0
      scattering_opacity = 0.0d0
      delta = 0.0d0

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
      temporary_delta = 0.0
      call return_scattering_opacity_spectra_given_neutrino_scheme(temporary_spectra,temporary_delta,eos_variables)

      !to fix minimum, 1.0e-30 in units of cm^-1 is effectivly zero
      temporary_spectra = max(temporary_spectra,1.0d-30)

      scattering_opacity(1,1:number_groups) = scattering_opacity(1,1:number_groups) + &
           temporary_spectra(1,1:number_groups)
      scattering_opacity(2,1:number_groups) = scattering_opacity(2,1:number_groups) + &
           temporary_spectra(2,1:number_groups)
      delta(1,1:number_groups) = delta(1,1:number_groups) + &
           temporary_spectra(1,1:number_groups) * temporary_delta(1,1:number_groups)
      delta(2,1:number_groups) = delta(2,1:number_groups) + &
           temporary_spectra(2,1:number_groups) * temporary_delta(2,1:number_groups)

      if (number_local_species.eq.3) then
         scattering_opacity(3,1:number_groups) = scattering_opacity(3,1:number_groups) + &
              (temporary_spectra(3,1:number_groups) + temporary_spectra(4,1:number_groups) + &
              temporary_spectra(5,1:number_groups) + temporary_spectra(6,1:number_groups))/4.0d0
         delta(3,1:number_groups) = delta(3,1:number_groups) + ( &
              temporary_spectra(3,1:number_groups)*temporary_delta(3,1:number_groups) + &
              temporary_spectra(4,1:number_groups)*temporary_delta(4,1:number_groups) + &
              temporary_spectra(5,1:number_groups)*temporary_delta(5,1:number_groups) + &
              temporary_spectra(6,1:number_groups)*temporary_delta(6,1:number_groups))/4.0d0

      else if (number_local_species.eq.4) then
         scattering_opacity(3,1:number_groups) = scattering_opacity(3,1:number_groups) + &
              (temporary_spectra(3,1:number_groups)+temporary_spectra(5,1:number_groups))/2.0d0
         scattering_opacity(4,1:number_groups) = scattering_opacity(4,1:number_groups) + &
              (temporary_spectra(4,1:number_groups)+temporary_spectra(6,1:number_groups))/2.0d0
         delta(3,1:number_groups) = delta(3,1:number_groups) + ( &
              temporary_spectra(3,1:number_groups)*temporary_delta(3,1:number_groups) + &
              temporary_spectra(5,1:number_groups)*temporary_delta(5,1:number_groups))/2.0d0
         delta(4,1:number_groups) = delta(4,1:number_groups) + ( &
              temporary_spectra(4,1:number_groups)*temporary_delta(4,1:number_groups) + &
              temporary_spectra(6,1:number_groups)*temporary_delta(6,1:number_groups))/2.0d0

      else if (number_local_species.eq.6) then
         scattering_opacity(3,1:number_groups) = scattering_opacity(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups)
         scattering_opacity(4,1:number_groups) = scattering_opacity(4,1:number_groups) + &
              temporary_spectra(4,1:number_groups)
         scattering_opacity(5,1:number_groups) = scattering_opacity(5,1:number_groups) + &
              temporary_spectra(5,1:number_groups)
         scattering_opacity(6,1:number_groups) = scattering_opacity(6,1:number_groups) + &
              temporary_spectra(6,1:number_groups)
         delta(3,1:number_groups) = delta(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups)*temporary_delta(3,1:number_groups)
         delta(4,1:number_groups) = delta(4,1:number_groups) + &
              temporary_spectra(4,1:number_groups)*temporary_delta(4,1:number_groups)
         delta(5,1:number_groups) = delta(5,1:number_groups) + &
              temporary_spectra(5,1:number_groups)*temporary_delta(5,1:number_groups)
         delta(6,1:number_groups) = delta(6,1:number_groups) + &
              temporary_spectra(6,1:number_groups)*temporary_delta(6,1:number_groups)
      endif
      delta = delta / scattering_opacity !pull out the weights to get the opacity-weighted average

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
     
      if (apply_kirchoff_to_pair_creation) then
         absorption_opacity(1,1:number_groups) = absorption_opacity(1,1:number_groups) + &
              temporary_spectra(1,1:number_groups)/blackbody_spectra(1,1:number_groups)
         absorption_opacity(2,1:number_groups) = absorption_opacity(2,1:number_groups) + &
              temporary_spectra(2,1:number_groups)/blackbody_spectra(2,1:number_groups)
      end if
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

         if (apply_kirchoff_to_pair_creation) then
            absorption_opacity(3,1:number_groups) = absorption_opacity(3,1:number_groups) + &
                 (temporary_spectra(3,1:number_groups)/blackbody_spectra(3,1:number_groups) + &
                 temporary_spectra(4,1:number_groups)/blackbody_spectra(3,1:number_groups) + &
                 temporary_spectra(5,1:number_groups)/blackbody_spectra(3,1:number_groups) + &
                 temporary_spectra(6,1:number_groups)/blackbody_spectra(3,1:number_groups))/4.0d0
         end if
         emissivities(3,1:number_groups) = emissivities(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups) + temporary_spectra(4,1:number_groups) + &
              temporary_spectra(5,1:number_groups) + temporary_spectra(6,1:number_groups)

      !average neutrinos and antineutrinos individually
      else if (number_local_species.eq.4) then
         if (apply_kirchoff_to_pair_creation) then
            absorption_opacity(3,1:number_groups) = absorption_opacity(3,1:number_groups) + &
                 (temporary_spectra(3,1:number_groups)/blackbody_spectra(3,1:number_groups) + &
                 temporary_spectra(5,1:number_groups)/blackbody_spectra(3,1:number_groups))/2.0d0
            absorption_opacity(4,1:number_groups) = absorption_opacity(4,1:number_groups) + &
                 (temporary_spectra(4,1:number_groups)/blackbody_spectra(3,1:number_groups) + &
                 temporary_spectra(6,1:number_groups)/blackbody_spectra(3,1:number_groups))/2.0d0
         end if
         emissivities(3,1:number_groups) = emissivities(3,1:number_groups) + &
              temporary_spectra(3,1:number_groups) + temporary_spectra(5,1:number_groups)
         emissivities(4,1:number_groups) = emissivities(4,1:number_groups) + &
              temporary_spectra(4,1:number_groups) + temporary_spectra(6,1:number_groups)

      !no averaging at all, what six different species
      else if (number_local_species.eq.6) then
         if (apply_kirchoff_to_pair_creation) then
            absorption_opacity(3,1:number_groups) = absorption_opacity(3,1:number_groups) + &
                 temporary_spectra(3,1:number_groups)/blackbody_spectra(3,1:number_groups)
            absorption_opacity(4,1:number_groups) = absorption_opacity(4,1:number_groups) + &
                 temporary_spectra(4,1:number_groups)/blackbody_spectra(3,1:number_groups)
            absorption_opacity(5,1:number_groups) = absorption_opacity(5,1:number_groups) + &
                 temporary_spectra(5,1:number_groups)/blackbody_spectra(3,1:number_groups)
            absorption_opacity(6,1:number_groups) = absorption_opacity(6,1:number_groups) + &
                 temporary_spectra(6,1:number_groups)/blackbody_spectra(3,1:number_groups)
         end if
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
            if(.not. do_transport_opacities) then
               if(delta(i,ng).ne.delta(i,ng)) stop "NaN in scattering delta"
            endif
         enddo      
      enddo

      emissivities(:,:) = max(1.0d-300,emissivities(:,:)) !ergs/cm^3/s/MeV/srad
      absorption_opacity(:,:) = max(1.0d-300,absorption_opacity(:,:)) !cm^-1
      scattering_opacity(:,:) = max(1.0d-300,scattering_opacity(:,:)) !cm^-1


    end subroutine single_point_return_all

    !calcualates the expansion of the scattering kernal, This is
    !\Phi_{out}_{0/1}, it does not contain the
    !exp(-\beta(\omega-\omega^\prime)), see make_table_example for symmetries
    subroutine single_Ipoint_return_all(iin,eta,temperature, &
         Phi0s,Phi1s,neutrino_local_scheme)

      implicit none

      !input
      real*8, dimension(:,:), intent(out) :: Phi0s
      real*8, dimension(:,:), intent(out) :: Phi1s
      integer, intent(in) :: iin
      real*8, intent(in) :: eta,temperature
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

      !local
      integer :: ns,ng
      integer :: number_local_species
      real*8  :: nuenergyin 

      !functions
      real*8 :: NES_Phi0_ThompsonBruenn
      real*8 :: NES_Phi1_ThompsonBruenn

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
         stop "single_Ipoint_return_all:incorrect neutrino scheme"
      endif

      if (size(Phi0s,1).ne.number_local_species) then
         stop "single_Ipoint_return_all:provided array has wrong number of species"
      endif
      if (size(Phi0s,2).ne.number_groups) then
         stop "single_Ipoint_return_all:provided array has wrong number of groups"
      endif

      if (size(Phi1s,1).ne.number_local_species) then
         stop "single_Ipoint_return_all:provided array has wrong number of species"
      endif
      if (size(Phi1s,2).ne.number_groups) then
         stop "single_Ipoint_return_all:provided array has wrong number of groups"
      endif
      
      nuenergyin = energies(iin)

      Phi0s = 0.0d0
      Phi1s = 0.0d0

      !best way to calculate a kernel with E_iout > E_in is to
      !calculate it with E_in^\prime = E_out and E_out^\prime = E_in,
      !and use that fact that Rout(e2,e1) = Rin(e1,e2) =
      !exp(-(e2-e1)/T * Rout(e1,e2). Numerically this is desired
      !because E_out > E_in doesn't work so well in these routines (in
      !fact it may be assumed in their derivation that Ein>=Eout [or
      !at least Ein >= Eout-Ee] although I cannot find the place)
      do ng=1,iin
         !electron neutrinos
         if (add_nue_Iscattering_electrons) then
            Phi0s(1,ng) = NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng), &
                 eta,temperature,1)
            Phi1s(1,ng) = NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng), &
                 eta,temperature,1)
         endif

         !electron antineutrinos
         if (add_anue_Iscattering_electrons) then
            Phi0s(2,ng) = NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng), &
                 eta,temperature,2)
            Phi1s(2,ng) = NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng), &
                 eta,temperature,2)
         endif

         if (number_local_species.eq.3) then
            !already ensure that all 4 are the all included or not
            if (add_numu_Iscattering_electrons) then
               Phi0s(3,ng) = (NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,3) + NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,4) + NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,5) + NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,6))/4.0d0
               Phi1s(3,ng) = (NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,3) + NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,4) + NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,5) + NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,6))/4.0d0
            endif

         else if (number_local_species.eq.4) then

            !already ensure that mu and tau the same
            if (add_numu_Iscattering_electrons) then
               Phi0s(3,ng) = (NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,3) + NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,5))/2.0d0
               Phi1s(3,ng) = (NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,3) + NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,5))/2.0d0
            endif

            !already ensure that amu and atau the same
            if (add_anumu_Iscattering_electrons) then
               Phi0s(4,ng) = (NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,4) + NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,6))/2.0d0
               Phi1s(4,ng) = (NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,4) + NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng), &
                    eta,temperature,6))/2.0d0
            endif

         else if (number_local_species.eq.6) then

            if (add_numu_Iscattering_electrons) then
               Phi0s(3,ng) = NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng),eta,temperature,3)
               Phi1s(3,ng) = NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng),eta,temperature,3)
            endif

            if (add_anumu_Iscattering_electrons) then
               Phi0s(4,ng) = NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng),eta,temperature,4)
               Phi1s(4,ng) = NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng),eta,temperature,4)
         endif

            if (add_nutau_Iscattering_electrons) then
               Phi0s(5,ng) = NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng),eta,temperature,5)
               Phi1s(5,ng) = NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng),eta,temperature,5)
            endif

            if (add_anutau_Iscattering_electrons) then
               Phi0s(6,ng) = NES_Phi0_ThompsonBruenn(nuenergyin,energies(ng),eta,temperature,6)
               Phi1s(6,ng) = NES_Phi1_ThompsonBruenn(nuenergyin,energies(ng),eta,temperature,6)
            endif

         else
            stop "shouldn't be here"
         endif

      enddo

      if(debug) then
         ! check that phis are consistent
         do ns=1,number_local_species
            do ng=1,iin
               if(abs(Phi1s(ns,ng)/Phi0s(ns,ng)) > 1.0d0) then
                  write(*,*) ns, ng, Phi1s(ns,ng), Phi0s(ns,ng), Phi1s(ns,ng)/Phi0s(ns,ng)
               end if
            enddo
         enddo
      endif

    end subroutine single_Ipoint_return_all

    !calcualates the expansion of the epannihilation kernal, These are
    !|phi^{prod/annihl}_{0/1}
    subroutine single_epannihil_kernel_point_return_all(iin,eta,temperature, &
         Phi0s,Phi1s,neutrino_local_scheme)

      implicit none

      !input
      real*8, dimension(:,:,:), intent(out) :: Phi0s !species,other neutrino's energy, prod/annihilation
      real*8, dimension(:,:,:), intent(out) :: Phi1s !species,other neutrino's energy, prod/annihilation
      integer, intent(in) :: iin !this neutrinos energy
      real*8, intent(in) :: eta,temperature
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

      !local
      integer :: ns,ng
      integer :: number_local_species
      real*8  :: nu_energy_x
      real*8  :: nuother_energy_x

      !functions
      real*8 :: epannihil_Phi_Bruenn

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
         stop "single_epannihil_kernel_return_all:incorrect neutrino scheme"
      endif

      if (size(Phi0s,1).ne.number_local_species) then
         stop "single_epannihil_kernel_return_all:provided array has wrong number of species"
      endif
      if (size(Phi0s,2).ne.number_groups) then
         stop "single_epannihil_kernel_return_all:provided array has wrong number of groups"
      endif
      if (size(Phi0s,3).ne.2) then
         stop "single_epannihil_kernel_return_all:provided array has wrong number of kernels"
      endif

      if (size(Phi1s,1).ne.number_local_species) then
         stop "single_epannihil_kernel_return_all:provided array has wrong number of species"
      endif
      if (size(Phi1s,2).ne.number_groups) then
         stop "single_epannihil_kernel_return_all:provided array has wrong number of groups"
      endif
      if (size(Phi1s,3).ne.2) then
         stop "single_epannihil_kernel_return_all:provided array has wrong number of kernels"
      endif
      
      nu_energy_x = energies(iin)/temperature

      Phi0s = 0.0d0
      Phi1s = 0.0d0

      do ng=1,number_groups
         nuother_energy_x = energies(ng)/temperature
         !electron neutrinos
         if (add_nue_kernel_epannihil) then
            Phi0s(1,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,1,1,0)!production of Phi_0
            Phi0s(1,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,1,2,0)!annihilation of Phi_0
            Phi1s(1,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,1,1,1)!production of Phi_1
            Phi1s(1,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,1,2,1)!annihilation of Phi_1
         endif

         !electron antineutrinos
         if (add_anue_kernel_epannihil) then
            Phi0s(2,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,2,1,0)!production of Phi_0
            Phi0s(2,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,2,2,0)!annihilation of Phi_0
            Phi1s(2,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,2,1,1)!production of Phi_1
            Phi1s(2,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,2,2,1)!annihilation of Phi_1
         endif

         if (number_local_species.eq.3) then
            !already ensure that all 4 are the all included or not
            if (add_numu_kernel_epannihil) then
               Phi0s(3,ng,1) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,1,0) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,1,0) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,1,0) + & 
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,1,0))/4.0d0 !production of Phi_0
               Phi0s(3,ng,2) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,2,0) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,2,0) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,2,0) + & 
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,2,0))/4.0d0 !annihilation of Phi_0
               Phi1s(3,ng,1) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,1,1) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,1,1) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,1,1) + & 
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,1,1))/4.0d0 !production of Phi_1
               Phi1s(3,ng,2) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,2,1) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,2,1) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,2,1) + & 
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,2,1))/4.0d0 !annihilation of Phi_1

            endif

         else if (number_local_species.eq.4) then

            !already ensure that mu and tau the same
            if (add_numu_kernel_epannihil) then
               Phi0s(3,ng,1) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,1,0) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,1,0))/2.0d0 !production of Phi_0
               Phi0s(3,ng,2) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,2,0) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,2,0))/2.0d0 !annihilation of Phi_0
               Phi1s(3,ng,1) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,1,1) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,1,1))/2.0d0 !production of Phi_1
               Phi1s(3,ng,2) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,2,1) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,2,1))/2.0d0 !annihilation of Phi_1

            endif

            !already ensure that amu and atau the same
            if (add_anumu_kernel_epannihil) then
               Phi0s(3,ng,1) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,1,0) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,1,0))/2.0d0 !production of Phi_0
               Phi0s(3,ng,2) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,2,0) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,2,0))/2.0d0 !annihilation of Phi_0
               Phi1s(3,ng,1) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,1,1) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,1,1))/2.0d0 !production of Phi_1
               Phi1s(3,ng,2) = (epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,2,1) + &
                    epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,2,1))/2.0d0 !annihilation of Phi_1

            endif

         else if (number_local_species.eq.6) then

            if (add_numu_Iscattering_electrons) then
               Phi0s(3,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,1,0) !production of Phi_0
               Phi0s(3,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,2,0) !annihilation of Phi_0
               Phi1s(3,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,1,1) !production of Phi_1
               Phi1s(3,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,3,2,1) !annihilation of Phi_1
            endif

            if (add_anumu_Iscattering_electrons) then
               Phi0s(4,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,1,0) !production of Phi_0
               Phi0s(4,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,2,0) !annihilation of Phi_0
               Phi1s(4,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,1,1) !production of Phi_1
               Phi1s(4,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,4,2,1) !annihilation of Phi_1
            endif

            if (add_nutau_Iscattering_electrons) then
               Phi0s(5,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,1,0) !production of Phi_0
               Phi0s(5,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,2,0) !annihilation of Phi_0
               Phi1s(5,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,1,1) !production of Phi_1
               Phi1s(5,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,5,2,1) !annihilation of Phi_1
            endif

            if (add_anutau_Iscattering_electrons) then
               Phi0s(6,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,1,0) !production of Phi_0
               Phi0s(6,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,2,0) !annihilation of Phi_0
               Phi1s(6,ng,1) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,1,1) !production of Phi_1
               Phi1s(6,ng,2) = epannihil_Phi_Bruenn(nu_energy_x,nuother_energy_x,eta,6,2,1) !annihilation of Phi_1
            endif

         else
            stop "shouldn't be here"
         endif
         
      enddo

      !comming from epannihil_Phi_Bruenn, the kernels have no
      !units.  we need cm^3/s and need to restore 2 factors of
      !temperature, Bruenn C62
      Phi0s = Phi0s*2.0d0*Gfermi**2/(2.0d0*pi)*hbarc_mevcm**2*temperature**2*clight !cm^3/s
      Phi1s = Phi1s*2.0d0*Gfermi**2/(2.0d0*pi)*hbarc_mevcm**2*temperature**2*clight !cm^3/s



    end subroutine single_epannihil_kernel_point_return_all


 end module nulib



