module CosmoRec

# CosmoRec.jl: a source-backed, native-Julia reimplementation (in progress) of
# the default recombination calculation performed by CosmoRec
# (https://github.com/cmbant/CosmoRec), targeting full compatibility with the
# SciML automatic-differentiation stack (ForwardDiff, Mooncake).
#
# Design and source-backed physics/derivative analysis live outside this
# repository for now (see the project's `local_analysis` tree:
# VERIFIED_ANALYSIS.md, sciml_chunk0/REPORT.md, sciml_chunk0/CHUNK0_REVIEW.md).
#
# Chunk 6a adds src/HIPDECoefficients.jl: population splines, get_rates_all, pd / Rp-Rm and the pd / Dnem coefficient splines of the HI diffusion PDE.
# Chunk 5a adds src/RecombinationODE.jl: the 12-/7-state ODE right-hand side dy/dz (native fcn_effective(int*, ...)) and the explicit production-table loader.
# Chunk 4d adds src/HeliumSwitch.jl: sampled-HeI state switch (discrete decision/reset/packing, frozen branch).
# Chunk 4c adds src/RecfastPP.jl: preliminary Recfast++ history (production terms, node grid, Saha segments; solver supplied by the caller).
# Chunk 4b adds src/CosmosAccessors.jl: GSL-compatible natural cubic splines and the explicit Cosmos accessors (history nodes, Hubble table supplied).
# Chunk 4a adds src/SahaInit.jl: explicit-input Saha per-level initialization and copy_LI_to_ysol packing.
# Chunk 3e adds src/DPescCoh.jl: explicit coherent-scattering DPesc integral (native off-table fallback of the H-I absorber).
# Chunk 3d adds src/FcnEffective.jl: assembled pointwise first-pass default H/He RHS (native fcn_effective layout).
# Chunk 3c adds src/HIAbsorption.jl: default approximate H-I absorption of HeI photons (explicit DP/fcorr/Bitot tables).
# Chunk 2 adds src/RateTable.jl: explicit-table port of the native HI effective-rate lookup get_rates().
#
# Chunk 1a: package skeleton only. This module intentionally
# exports no public API and contains no recombination physics yet. Chunk 1a's
# deliverable is a toolchain compatibility probe -- proving that a small stiff
# SciML ODE solve can be differentiated end-to-end with both ForwardDiff and
# Mooncake (via DifferentiationInterface + SciMLSensitivity) -- which lives
# under `test/` as a test-only helper, not as part of this module's public
# surface. See test/chunk1a_stiff_ad_probe.jl and docs/fixture_contract.md.

include("RateTable.jl")
include("RHS.jl")
include("HeRHS.jl")
include("HIAbsorption.jl")
include("DPescCoh.jl")
include("FcnEffective.jl")
include("SahaInit.jl")
include("CosmosAccessors.jl")
include("RecfastPP.jl")
include("HeliumSwitch.jl")
include("RecombinationODE.jl")
include("HIPDECoefficients.jl")
include("HIPDEProfiles.jl")
include("HIPDEDefine.jl")
include("HIPDESolver.jl")
include("HIPDEIntegrals.jl")
include("HIDiffusionFeedback.jl")
include("RecombinationDiffusion.jl")

export AtomicRateTable, ResolvedStates, RateConstants, NATIVE_CONSTANTS, RateResult, RateTableDomainError,
    get_rates, get_rates!, log_qnl_qe, stencil_start, read_native_rate_file, load_rate_table,
    RHSConstants, NATIVE_RHS_CONSTANTS, HydrogenAtom, NATIVE_HYDROGEN_ATOM, rho_g_fac, matter_temperature_rate, two_photon!, lyman_channel!, lyman!,
    continuum!, interlevel!, hydrogen_rhs!, proton_fraction, electron_fraction,
    HeliumConstants, NATIVE_HELIUM_CONSTANTS, HeliumAtom, NATIVE_HELIUM_ATOM, HeliumRateTable, read_native_helium_rate_file, load_helium_rate_table,
    helium_A, helium_stencil_start, get_helium_rates, helium_continuum!, helium_interlevel!, helium_two_photon!, helium_lyman!, helium_base_rhs!,
    DPTableDomainError, HIAbsConstants, NATIVE_HIABS_CONSTANTS, HIAbsLine, NATIVE_HIABS_SINGLET, NATIVE_HIABS_TRIPLET, NATIVE_HIABS_PROFILE_A21, DPSheet, DPTable, read_native_dp_table,
    FcorrSpline, read_native_fcorr, fcorr, BitotSeries, read_native_helium_bitot, helium_Bitot, pd_singlet, sobolev_p, tau_S_gw, exp_nu, dp_sheet_start, dp_lookup, dp_correction,
    hi_abs_eta, hi_abs_singlet, hi_abs_triplet, hi_absorption_rhs!,
    VoigtLine, HLycCross, NATIVE_HLYC, DPescModel, NATIVE_DPESC, PATTERSON_NPOINTS, PATTERSON_X, PATTERSON_W, dawson_jc, erf_jc, voigt_phi, voigt_xi_int, voigt_a, voigt_DnuT,
    lyc_ratio, patterson_integrate, dpesc_appr_I_sym, dpesc_coh, dpesc_fallback, dp_T_in_table,
    SahaConstants, NATIVE_SAHA_CONSTANTS, SahaLevels, NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM, saha_lte_hydrogen, saha_lte_helium, SahaInputs, saha_initial_state,
    NATIVE_RESHI, NATIVE_RESHE, pack_ysol,
    NaturalCubicSpline, natural_cubic_spline, SplineDomainError, spline_eval, spline_eval_native, CosmosConstants, RecfastSplines, recfast_splines, HubbleTable, hubble_table,
    CosmosAccessors, cosmos_H, cosmos_TCMB, cosmos_Ttoz, cosmos_Nb, cosmos_NH, cosmos_fHe, cosmos_sigT, cosmos_rho_g_1_cm3, cosmos_Xe_Seager, cosmos_Xe_b, cosmos_dXe_dz, cosmos_Ne,
    cosmos_Ntot, cosmos_kappa_cool, cosmos_Te_Tg, cosmos_Te, cosmos_SahaBoltz_HeIII, cosmos_SahaBoltz_HeII, cosmos_SahaBoltz_HeII1s, cosmos_SahaBoltz_HeI1s, cosmos_SahaBoltz_HI1s,
    RecfastConstants, NATIVE_RECFAST_CONSTANTS, RECFAST_THETA_NAMES, rf_fHe, rf_H0, rf_NH, rf_TCMB, recfast_rhs, recfast_rhs3!, rf_Te_QS, rf_SahaBoltz_HeIII, rf_SahaBoltz_HeII,
    load_hi_bitot_table, get_rates_all, HIPopulationSplines, hi_Xe, hi_rho, hi_Xi, HIPDELevels, NATIVE_HI_PDE_LEVELS, hi_rp_rm_pd, HIPDECoefficientSplines, hi_pde_coefficients, hi_pd, hi_Dnem,
    NATIVE_HYDROGEN_RESOLVED, load_native_effective_model, RecombinationModel, ode_nstate, ode_background, ode_unpack, recombination_rhs!, recombination_rhs, init_xarr_linear, RecombinationODEParams, recombination_ode!, ode_initial_state, ode_observables, recombination_pass, init_xarr_log, RecfastTailParams, recfast_tail_rhs!, recfast_tail, recfast_tail_inputs, recombination_output_rows, return_solution_to_grid, recombination_history,
    NATIVE_XI_HEI_SWITCH, NATIVE_HE_SWITCH_ZOFF, NATIVE_HE_FLOOR, helium_switch_condition, helium_switch_reset, helium_switch, helium_switch_ysize,
    rf_alphaH, rf_alphaHe, rf_boltzmann, rf_inv_boltzmann, rf_sahaboltz, rf_inv_sahaboltz, RecfastGrid, recfast_grid, recfast_history,
    cosmos_X1s, cosmos_Xp, cosmos_XHeII1s, cosmos_XHeI1s, cosmos_NHeI, cosmos_NHeII, cosmos_NHeIII, saha_inputs_at,
    EffectiveBackground, EffectiveModel, effective_fractions, fcn_effective!, fcn_effective,
    pde_log10factorial, hydrogen_A_SH, HITransitionData, hi_transition_data, gamma_np, A_npks, A_npkd, hi_R1snp, hi_Rksnp, hi_Rkdnp, HIProfileData, load_hi_profile_data,
    sigma_2s1s_2gamma, sigma_ns1s_2gamma, sigma_nd1s_2gamma, sigma_ns1s_2gamma_ratio, sigma_nd1s_2gamma_ratio, sigma_2s1s_Raman_ratio, hi_pde_grid, HIPDESetup, hi_pde_setup,
    HIPDEConstants, NATIVE_HI_PDE_CONSTANTS, HIPDEAtom, NATIVE_HI_PDE_ATOM, voigt_dphi_dx, HIPDEModel, HIPDEState, hi_pde_rhs_coefficients!, hi_pde_rhs_coefficients,
    LagrangeO2, lagrange_o2, HIPDEStepper, hi_pde_step!, polint_jc, hi_pde_march,
    sobolev_p_ij, hi_tau_S, hi_pde_integrals, hi_pde_corrections, PattersonLevels, replay,
    HIDiffusionFeedback, hi_diffusion_feedback, hi_DI1_2s, hi_DF_2g, hi_DF_R, hi_diffusion_rhs!, with_diffusion,
    pass_solution_rows, ScaledBackground, HIDiffusionInputs, hi_diffusion_stage, recombination_history_diffusion

end # module CosmoRec
