  !requested interactions

  !corrections
  logical :: do_weak_mag_corrections = .true.
  logical :: do_ionion_correlation = .true.
  logical :: do_heavyscat_formfactor = .true.
  logical :: do_electronpolarization_correction = .true.
  logical :: do_nc_virial_correction = .false.
  logical :: do_strange_coupling = .false.
  logical :: do_transport_opacities = .true.
  real*8, parameter :: gAs = -0.2d0 !dimensionless, for the strange coupling

  !final state blocking options
  logical :: do_breunn_final_state_nucleon_blocking_abs = .true.
  logical :: do_breunn_final_state_nucleon_blocking_scat = .true.

  !absorptions
  logical :: add_nue_absorption_on_n = .true.
  logical :: add_anue_absorption_on_p = .true.
  logical :: add_nue_absorption_on_A = .false.
  logical :: add_nux_absorption_on_n_and_p = .false.

  !scatterings
  logical :: add_nue_scattering_n = .true.
  logical :: add_nue_scattering_p = .true.
  logical :: add_nue_scattering_heavies = .true.
  logical :: add_nue_scattering_electrons = .false.
  logical :: add_nue_scattering_alphas = .true.
  logical :: add_anue_scattering_n = .true.
  logical :: add_anue_scattering_p = .true.
  logical :: add_anue_scattering_heavies = .true.
  logical :: add_anue_scattering_electrons = .false.
  logical :: add_anue_scattering_alphas = .true.
  logical :: add_numu_scattering_n = .true.
  logical :: add_numu_scattering_p = .true.
  logical :: add_numu_scattering_heavies = .true.
  logical :: add_numu_scattering_electrons = .false.
  logical :: add_numu_scattering_alphas = .true.
  logical :: add_anumu_scattering_n = .true.
  logical :: add_anumu_scattering_p = .true.
  logical :: add_anumu_scattering_heavies = .true.
  logical :: add_anumu_scattering_electrons = .false.
  logical :: add_anumu_scattering_alphas = .true.
  logical :: add_nutau_scattering_n = .true.
  logical :: add_nutau_scattering_p = .true.
  logical :: add_nutau_scattering_heavies = .true.
  logical :: add_nutau_scattering_electrons = .false.
  logical :: add_nutau_scattering_alphas = .true.
  logical :: add_anutau_scattering_n = .true.
  logical :: add_anutau_scattering_p = .true.
  logical :: add_anutau_scattering_heavies = .true.
  logical :: add_anutau_scattering_electrons = .false.
  logical :: add_anutau_scattering_alphas = .true.

  !inelastic scatterings
  logical :: add_nue_Iscattering_electrons = .true.
  logical :: add_anue_Iscattering_electrons = .true.
  logical :: add_numu_Iscattering_electrons = .true.
  logical :: add_anumu_Iscattering_electrons = .true.
  logical :: add_nutau_Iscattering_electrons = .true.
  logical :: add_anutau_Iscattering_electrons = .true.

  !emissivities
  logical :: add_nue_emission_epannihil = .false.
  logical :: add_anue_emission_epannihil = .false.
  logical :: add_numu_emission_epannihil = .true.
  logical :: add_anumu_emission_epannihil = .true.
  logical :: add_nutau_emission_epannihil = .true.
  logical :: add_anutau_emission_epannihil = .true.

  logical :: add_nue_emission_NNBrems = .false.
  logical :: add_anue_emission_NNBrems = .false.
  logical :: add_numu_emission_NNBrems = .true.
  logical :: add_anumu_emission_NNBrems = .true.
  logical :: add_nutau_emission_NNBrems = .true.
  logical :: add_anutau_emission_NNBrems = .true.	

  logical :: add_nue_emission_weakinteraction_ecap = .true.
  logical :: add_anue_emission_weakinteraction_poscap = .true.

  logical :: apply_kirchoff_to_pair_creation = .true.

  !kernels for full calculation of thermal processes
  logical :: add_nue_kernel_epannihil = .false.
  logical :: add_anue_kernel_epannihil = .false.
  logical :: add_numu_kernel_epannihil = .false.
  logical :: add_anumu_kernel_epannihil = .false.
  logical :: add_nutau_kernel_epannihil = .false.
  logical :: add_anutau_kernel_epannihil = .false.