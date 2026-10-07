# Chunk 11b worker: ONE configuration per fresh Julia process (spawned by chunk11b_quadplan_stage_ad.jl; full-stage Mooncake preparations retain
# ~10-13 GB each, so the default-path and the plan-path preparations must never share a process). Usage:
#   julia --project=<active test env> chunk11b_stage_worker.jl <mode> <outfile> <seed>     mode = old | plan | ref | record
# Writes `key<TAB>value` text: provenance (versions, flags, input/weight checksums) and the results; Float64 values as exact bit patterns.
#   old / plan : prepared Mooncake gradient of the objective dot(W, stage outputs) at both points (default path / opt-in quadrature plan).
#   ref        : 256-bit reference of the same objective: central difference h = 1e-30 in BigFloat (precision 256) with the Patterson levels
#                REPLAYED from the Float64 run at each point, on the default (generic) path, so the natural-spline stored factors stay Float64
#                and the Patterson nodes/weights are the same Float64 constants (same-stored-factor, frozen-branch reference). No Mooncake.
#   record     : Patterson record/replay bookkeeping on the short production march, default vs plan path.
# Objective: W = randn(Xoshiro(seed)) ./ (per-block max of the default-path outputs at [1,1]) ./ length, mask 500 < z < 2000.
# REQUIRES COSMOREC_NATIVE_DATA_DIR (inherited from the parent test process).
using Random
using LinearAlgebra
using SHA
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
using CosmoRec: _hi_quad_plans
include(joinpath(@__DIR__, "chunk7_helpers.jl"))

const MODE11W = ARGS[1]
const OUT11W = ARGS[2]
const SEED11W = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 2705
const MC11W = AutoMooncake(; config = nothing)
const D11W = HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6)
const R11W = NODES5[:, 1:9]
const PTS11W = ([1.0, 1.0], [0.999, 1.001])
bits11(v) = join(string.(reinterpret(UInt64, collect(Float64, v)); base = 16, pad = 16), ",")
sha11(v) = bytes2hex(sha256(reinterpret(UInt8, collect(Float64, v))))
function out11w(p; plan = nothing, levels = nothing)
    o, _ = hi_diffusion_stage(RM5, D11W, R11W; hscale = p[1], nbscale = p[2], quadplan = plan, levels = levels)
    mask = (o.z .< 2000.0) .& (o.z .> 500.0)
    vcat(o.DI1_2s[mask], o.DF_2g[1][mask], o.DF_2g[2][mask], o.DF_R[1][mask])
end
function weights11w(seed)
    F0 = out11w([1.0, 1.0]); NB = length(F0) ÷ 4      # default path, so W cannot differ between configurations
    SC = vcat([fill(maximum(abs, F0[((k - 1) * NB + 1):(k * NB)]), NB) for k in 1:4]...)
    return F0, randn(Random.Xoshiro(seed), length(F0)) ./ SC ./ length(F0)
end

res = Pair{String,String}[]
push!(res, "mode" => MODE11W, "seed" => string(SEED11W), "julia" => string(VERSION), "mooncake" => string(pkgversion(Mooncake)),
      "cosmorec_src" => pathof(CosmoRec), "check_bounds" => string(Base.JLOptions().check_bounds), "threads" => string(Threads.nthreads()))
if MODE11W in ("old", "plan", "ref")
    F0, W = weights11w(SEED11W)
    push!(res, "F0_sha256" => sha11(F0), "W_sha256" => sha11(W), "W_length" => string(length(W)), "points" => repr(PTS11W))
    if MODE11W == "ref"
        H = 1.0e-30
        for (k, p) in enumerate(PTS11W)
            rec = PattersonLevels(); f64 = dot(W, out11w(p; levels = rec))
            ok = dot(W, out11w(p; levels = replay(rec))) === f64       # quality: Float64 replay reproduces the adaptive objective bitwise
            g = setprecision(BigFloat, 256) do
                Wb = big.(W); pb = big.(p)
                [begin
                     e = zeros(BigFloat, 2); e[i] = BigFloat(H)
                     (dot(Wb, out11w(pb .+ e; levels = replay(rec))) - dot(Wb, out11w(pb .- e; levels = replay(rec)))) / (2 * BigFloat(H))
                 end for i in 1:2]
            end
            push!(res, "levels_$(k)_sha256" => bytes2hex(sha256(reinterpret(UInt8, rec.levels))), "replay_bitwise_$k" => string(ok),
                  "grad_$k" => bits11(Float64.(g)), "grad_$(k)_decimal" => repr(Float64.(g)), "grad_$(k)_big" => join(string.(g), ","))
        end
    else
        plan = MODE11W == "plan" ? _hi_quad_plans(SETUP7) : nothing
        obj(p, w) = dot(w, out11w(p; plan = plan))
        prep = prepare_gradient(obj, MC11W, [1.0, 1.0], Constant(W))
        for (k, p) in enumerate(PTS11W)
            g = gradient(obj, prep, MC11W, p, Constant(W))
            push!(res, "grad_$k" => bits11(g), "grad_$(k)_decimal" => repr(g))
        end
    end
elseif MODE11W == "record"
    # Patterson bookkeeping, short production march, default vs plan path. Mooncake's reverse pass UNDOES the forward push! into rec.levels
    # (rule contract: pullbacks restore mutated primal state; chunk29/jobs/record_repro), so the levels are compared from PRIMAL records, and the
    # record-mode gradient must equal the gradient that replays exactly those levels (the forward pass realised them).
    plan = _hi_quad_plans(SETUP7)
    x0 = vec(R11W[:, 2:9])
    function lev(x, rec, qp)
        rows = hcat(R11W[:, 1], reshape(x, size(R11W, 1), 8)); pops = HIPopulationSplines(rows; zs = ZS6, ze = ZE6)
        m = HIPDEModel(SETUP7, pops, hi_pde_coefficients(rows, pops, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6), ACC5)
        return sum(hi_pde_corrections(m; zs = 2500.0, ze = 2470.0, T = eltype(x), levels = rec, quadplan = qp).DI1_2s)
    end
    for (lab, qp) in (("old", nothing), ("plan", plan))
        r0 = PattersonLevels(); lev(x0, r0, qp); L = copy(r0.levels)
        ra = PattersonLevels(); fa(x) = lev(x, ra, qp); pa = prepare_gradient(fa, MC11W, x0)
        empty!(ra.levels); ra.pos = 0; ga = gradient(fa, pa, MC11W, x0)
        rr = PattersonLevels(:replay, copy(L), 0); fr(x) = lev(x, rr, qp); pr = prepare_gradient(fr, MC11W, x0)
        rr.pos = 0; gr = gradient(fr, pr, MC11W, x0)
        push!(res, "levels_$lab" => join(L, ","), "record_grad_equals_replay_grad_$lab" => string(ga == gr))
    end
else
    error("chunk11b_stage_worker: unknown mode $MODE11W")
end
push!(res, "maxrss_MB" => string(Sys.maxrss() / 2^20))
open(OUT11W, "w") do io
    foreach(kv -> println(io, kv.first, '\t', kv.second), res)
end
