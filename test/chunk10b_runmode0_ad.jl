# Chunk 10b: derivatives of the COMPOSED production (runmode 0) iteration: 3 ODE passes and 2 HI diffusion PDE stages with feedback, Recfast tail and
# output assembly, w.r.t. p = [F, A2s1s (Recfast), hscale, nbscale] (H(z) and NH(z) scales, applied consistently to the ODE passes AND the PDE stages).
# ForwardDiff through the whole composition vs central finite differences with step convergence (elasticity error, Chunk 5e metric).
# FROZEN (no derivative claimed): the helium-switch node (k_switch, all passes), the Saha initial state and the preliminary Recfast history inside the
# Cosmos accessors (Xe_Seager, rho splines, the loaded Hubble table shape), the PDE frequency grid and coefficient-grid size, the polint_JC stencils,
# the Patterson stopping decisions (primal; FD is evaluated on the same frozen decisions, recorded at p0), the feedback range ends (500, 2000) and the Recfast tail grid. A prepared Mooncake gradient of the composed
# objective is NOT part of this file (see docs/ACCEPTANCE_CHUNK10A.md). REQUIRES COSMOREC_NATIVE_DATA_DIR.
using Test
using Random
using ForwardDiff
using CosmoRec
@isdefined(SOLVEP5E) || include("chunk5e_helpers.jl")
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const OBS10B = Dict{String,Float64}()
rec10b(k, v) = (OBS10B[k] = max(get(OBS10B, k, 0.0), v))
const D10B = HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6)
function full_out10(p; levels = nothing)
    h = recombination_history_diffusion(RM5, THETA5E_vec(p), D10B, SOLVEP5E, solve_tail5, ZG5E; hscale = p[3], nbscale = p[4], k_switch = K5E, patterson_levels = levels)
    return vcat(h.Xe, h.Te)
end
# the Patterson stopping decisions of both PDE stages recorded at P0 (native decisions), then replayed (frozen branch)
const LEV10B = [PattersonLevels(), PattersonLevels()]
const F0REC10B = full_out10(P0_5E; levels = LEV10B)
full_out10_frozen(p) = full_out10(p; levels = replay.(LEV10B))

@testset "Chunk 10b: gradients of the composed runmode-0 iteration (ForwardDiff vs FD)" begin
    f0 = full_out10(P0_5E)
    @test all(isfinite, f0) && length(f0) == 2 * length(ZG5E)
    @test full_out10_frozen(P0_5E) == f0 == F0REC10B          # replaying the recorded decisions reproduces the native-decision evaluation bitwise
    J = ForwardDiff.jacobian(full_out10, P0_5E)
    @test all(isfinite, J)
    @test ForwardDiff.jacobian(full_out10_frozen, P0_5E) == J  # ForwardDiff differentiates the frozen-decision map (decisions are primal)
    errs = Float64[]; errs_free = Float64[]
    for hrel in (1e-2, 1e-3, 1e-4)
        push!(errs, elast_err(J, fd_jac(full_out10_frozen, P0_5E, hrel), f0, P0_5E))
        push!(errs_free, elast_err(J, fd_jac(full_out10, P0_5E, hrel), f0, P0_5E))
    end
    rec10b("FD_frozen_elasticity_best", minimum(errs)); rec10b("FD_frozen_elasticity_worst", maximum(errs))
    # DIAGNOSTIC: with free (re-decided) quadrature stopping rules the composed map has small jumps; FD then does not converge to any derivative
    rec10b("diag:FD_free_decisions_elasticity_best", minimum(errs_free)); rec10b("diag:FD_free_decisions_elasticity_worst", maximum(errs_free))
    @test minimum(errs) < 1.0e-4                     # same gate as the single pass (5e): FD is limited by the solver tolerance
    # the diffusion iterations change the derivative: compare with the single-pass Jacobian (5e route) at the same parameters
    J1 = ForwardDiff.jacobian(full_out, P0_5E)
    rec10b("jacobian_change_vs_single_pass_rel", maximum(abs.(J .- J1)) / maximum(abs.(J1)))
    @test maximum(abs.(J .- J1)) > 0
    rng = Random.Xoshiro(10)
    v = randn(rng, 4) .* abs.(P0_5E)
    dd = ForwardDiff.derivative(t -> full_out10(P0_5E .+ t .* v), 0.0)
    es = Float64[]
    for h in (1e-2, 1e-3, 1e-4)
        ref = (full_out10_frozen(P0_5E .+ h .* v) .- full_out10_frozen(P0_5E .- h .* v)) ./ (2h)
        push!(es, maximum(abs.(dd .- ref) ./ abs.(f0)))
    end
    rec10b("JVP_dir_best", minimum(es))
    @test minimum(es) < 1e-5
    foreach(kv -> println("Chunk10b ", kv[1], " = ", kv[2]), sort!(collect(OBS10B)))
end
