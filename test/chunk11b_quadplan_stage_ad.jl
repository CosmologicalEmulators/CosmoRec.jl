# Chunk 11b: full PDE-stage prepared Mooncake gradient through the opt-in quadrature plan (reverse rule of ext/CosmoRecChainRulesCoreExt.jl via the
# reverse-only Mooncake bridge) vs the default path (Mooncake-derived), under the REFERENCE-BASED gate approved by M. Bonici (2026-10-06):
#   primary   E_new <= E_old, E_x = |g_x - g_ref|_inf / |g_ref|_inf, g_ref = 256-bit frozen-Patterson-level same-stored-factor central-difference
#             reference (worker mode `ref`), required PER SEED (2705, 4242) and PER POINT ([1,1], [0.999,1.001]); all E values are printed;
#   diagnostic old-vs-new relative difference, REPORTED, with a noise-scale sanity bound 1e-5 (the Float64 derivative of this objective is
#             ill-conditioned at ~1e-6..1e-5: old/new agreement at 1e-10 is unattainable between Float64 realizations, see
#             cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk29/SWEEP_LOCALIZATION_REDUCED.md);
#   plus Patterson record/replay bookkeeping (worker mode `record`).
# Kernel / contract / primal-ForwardDiff parity / stale-plan hard gates live in chunk11_quadplan_ad.jl.
# Resource protocol: each full-stage preparation retains ~10-13 GB, so each configuration runs in its OWN fresh Julia process, sequentially, with
# the same julia flags, active project and environment as this process; results come back as text (exact Float64 bit patterns).
# The strict MC-vs-ForwardDiff 1e-6 comparison of this objective is a pre-existing failure of BOTH paths and is NOT asserted here.
# REQUIRES COSMOREC_NATIVE_DATA_DIR.
using Test

function run_worker11b(mode, dir, seed = 2705)
    out = joinpath(dir, "chunk11b_$(mode)_$(seed).tsv")
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $(joinpath(@__DIR__, "chunk11b_stage_worker.jl")) $mode $out $seed`
    ok = success(pipeline(cmd; stdout = stdout, stderr = stderr))
    d = Dict{String,String}()
    if ok && isfile(out)
        for l in eachline(out)
            k, v = split(l, '\t'; limit = 2); d[k] = v
        end
        foreach(l -> println("Chunk11b worker[$mode,$seed] ", l), filter(l -> !occursin("_decimal", l) && !startswith(l, "grad_"), readlines(out)))
    end
    return ok, d
end
frombits11b(s) = [reinterpret(Float64, parse(UInt64, x; base = 16)) for x in split(s, ',')]

@testset "Chunk 11b: opt-in quadrature plan, full-stage gradient, reference-based gate (fresh process per configuration)" begin
    haskey(ENV, "COSMOREC_NATIVE_DATA_DIR") || error("COSMOREC_NATIVE_DATA_DIR is not set (never skipped)")
    dir = mktempdir()
    for seed in (2705, 4242)
        ok_old, old = run_worker11b("old", dir, seed)       # sequential: each worker has exited before the next starts
        @test ok_old
        ok_new, new = run_worker11b("plan", dir, seed)
        @test ok_new
        ok_ref, ref = run_worker11b("ref", dir, seed)
        @test ok_ref
        (ok_old && ok_new && ok_ref) || continue
        @test old["mode"] == "old" && new["mode"] == "plan" && ref["mode"] == "ref"
        for k in ("seed", "julia", "mooncake", "cosmorec_src", "check_bounds", "threads", "F0_sha256", "W_sha256", "W_length", "points")
            @test old[k] == new[k] == ref[k]                 # identical inputs, objective, weights and toolchain
        end
        for k in 1:2
            @test ref["replay_bitwise_$k"] == "true"
            go = frombits11b(old["grad_$k"]); gn = frombits11b(new["grad_$k"]); gr = frombits11b(ref["grad_$k"])
            @test length(go) == length(gn) == length(gr) == 2 && all(isfinite, go) && all(isfinite, gn) && all(isfinite, gr)
            Eold = maximum(abs, go .- gr) / maximum(abs, gr); Enew = maximum(abs, gn .- gr) / maximum(abs, gr)
            dnew = maximum(abs, gn .- go) / maximum(abs, go)
            println("Chunk11b GATE seed=", seed, " point=", k, " ref=", ref["grad_$(k)_decimal"], " old=", old["grad_$(k)_decimal"],
                    " new=", new["grad_$(k)_decimal"], " E_old=", Eold, " E_new=", Enew, " old_vs_new(diagnostic)=", dnew)
            @test Enew <= Eold                               # primary reference-based criterion, per seed and point
            @test dnew <= 1.0e-5                             # diagnostic sanity bound (gross breakage only)
        end
    end
    ok_rec, rec = run_worker11b("record", dir)
    @test ok_rec
    if ok_rec
        @test !isempty(rec["levels_old"])
        @test rec["levels_old"] == rec["levels_plan"]        # primal record sequences identical on both paths
        @test rec["record_grad_equals_replay_grad_old"] == "true"
        @test rec["record_grad_equals_replay_grad_plan"] == "true"
    end
end
