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
    EffectiveBackground, EffectiveModel, effective_fractions, fcn_effective!, fcn_effective

end # module CosmoRec
