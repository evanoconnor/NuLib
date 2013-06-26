  !requested interactions

  !corrections
  logical :: do_weak_mag_corrections = .false.
  logical :: do_ionion_correlation = .false.
  logical :: do_heavyscat_formfactor = .false.
  logical :: do_electronpolarization_correction = .false.

  !absorptions
  logical :: add_nue_absorption_on_n = .false.
  logical :: add_anue_absorption_on_p = .false.
  logical :: add_nue_absorption_on_A = .false.

  !scatterings
  logical :: add_nue_scattering_n = .false.
  logical :: add_nue_scattering_p = .false.
  logical :: add_nue_scattering_heavies = .false.
  logical :: add_nue_scattering_electrons = .false.
  logical :: add_nue_scattering_alphas = .false.
  logical :: add_anue_scattering_n = .false.
  logical :: add_anue_scattering_p = .false.
  logical :: add_anue_scattering_heavies = .false.
  logical :: add_anue_scattering_electrons = .false.
  logical :: add_anue_scattering_alphas = .false.
  logical :: add_numu_scattering_n = .false.
  logical :: add_numu_scattering_p = .false.
  logical :: add_numu_scattering_heavies = .false.
  logical :: add_numu_scattering_electrons = .false.
  logical :: add_numu_scattering_alphas = .false.
  logical :: add_anumu_scattering_n = .false.
  logical :: add_anumu_scattering_p = .false.
  logical :: add_anumu_scattering_heavies = .false.
  logical :: add_anumu_scattering_electrons = .false.
  logical :: add_anumu_scattering_alphas = .false.
  logical :: add_nutau_scattering_n = .false.
  logical :: add_nutau_scattering_p = .false.
  logical :: add_nutau_scattering_heavies = .false.
  logical :: add_nutau_scattering_electrons = .false.
  logical :: add_nutau_scattering_alphas = .false.
  logical :: add_anutau_scattering_n = .false.
  logical :: add_anutau_scattering_p = .false.
  logical :: add_anutau_scattering_heavies = .false.
  logical :: add_anutau_scattering_electrons = .false.
  logical :: add_anutau_scattering_alphas = .false.

  !inelastic scatterings
  logical :: add_nue_Iscattering_electrons = .false.
  logical :: add_anue_Iscattering_electrons = .false.
  logical :: add_numu_Iscattering_electrons = .false.
  logical :: add_anumu_Iscattering_electrons = .false.
  logical :: add_nutau_Iscattering_electrons = .false.
  logical :: add_anutau_Iscattering_electrons = .false.

  !emissivities
  logical :: add_nue_emission_epannihil = .false.
  logical :: add_anue_emission_epannihil = .false.
  logical :: add_numu_emission_epannihil = .false.
  logical :: add_anumu_emission_epannihil = .false.
  logical :: add_nutau_emission_epannihil = .false.
  logical :: add_anutau_emission_epannihil = .false.

  logical :: add_nue_emission_NNBrems = .false.
  logical :: add_anue_emission_NNBrems = .false.
  logical :: add_numu_emission_NNBrems = .false.
  logical :: add_anumu_emission_NNBrems = .false.
  logical :: add_nutau_emission_NNBrems = .false.
  logical :: add_anutau_emission_NNBrems = .false.	

  logical :: add_nue_emission_weakinteraction_ecap = .true.
  logical :: add_anue_emission_weakinteraction_poscap = .false.
