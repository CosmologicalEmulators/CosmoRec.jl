# Phase 5 helpers: explicit native data location, production tables, accessors. The native tables are NEVER bundled: the tests REQUIRE the environment variable
# COSMOREC_NATIVE_DATA_DIR (the native `Rec_database` directory; read-only) and fail loudly if it is unset or the files are missing (no skip).
using CosmoRec
using SHA
@isdefined(HIST4B) || include("chunk4b_helpers.jl")
@isdefined(FX4A) || include("chunk4a_helpers.jl")
@isdefined(FX4D) || include("chunk4d_helpers.jl")

function native_data_dir5()
    d = get(ENV, "COSMOREC_NATIVE_DATA_DIR", "")
    isempty(d) && error("COSMOREC_NATIVE_DATA_DIR is not set: the Phase 5 tests require the native Rec_database directory (e.g. /home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/tools/CosmoRec/Rec_database); they never skip")
    isdir(d) || error("COSMOREC_NATIVE_DATA_DIR = $d is not a directory")
    for sub in ("Effective_Rates.HI", "Effective_Rates.HeI", "Pesc_Data")
        isdir(joinpath(d, sub)) || error("COSMOREC_NATIVE_DATA_DIR = $d lacks $sub")
    end
    return d
end
fcorr_path5(d = native_data_dir5()) = (p = joinpath(dirname(d), "Development", "Recombination", "Data.fcorr", "f.corr.dat"); isfile(p) || error("missing native f.corr.dat at $p"); p)

const DATADIR5 = native_data_dir5()
const EFF5 = load_native_effective_model(DATADIR5, fcorr_path5(DATADIR5))
const ACC5 = accessors4b()
const RM5 = RecombinationModel(EFF5, ACC5)

function read_fixture5a(path)
    rows = NamedTuple[]
    for l in eachline(path)
        (isempty(l) || startswith(l, "#") || !startswith(l, "F5A\t")) && continue
        f = split(l, '\t')
        iF = findfirst(==("F"), f)
        push!(rows, (id = parse(Int, f[2]), flag = f[3] == "1", z = parse(Float64, f[4]), pat = parse(Int, f[5]), Tg = parse(Float64, f[7]), NH = parse(Float64, f[9]), Hz = parse(Float64, f[11]),
                     Xp = parse(Float64, f[13]), y = parse.(Float64, f[15:(iF - 1)]), f = parse.(Float64, f[(iF + 1):end])))
    end
    return rows
end

# ---- SciML stiff solve callback for the recombination pass (test-only: the package has no ODE dependency) ----
@isdefined(primal_norm) || include("chunk4c_helpers.jl")
using ForwardDiff
using SciMLBase: ODEFunction, ODEProblem, solve
import SciMLBase
using OrdinaryDiffEqRosenbrock: Rodas5P, Rodas4P

function jac5!(Jm, u, p, z)
    Jm .= ForwardDiff.jacobian(uu -> (du = similar(uu); recombination_ode!(du, uu, p, z); du), u); return nothing
end
function tgrad5!(dT, u, p, z)
    dT .= ForwardDiff.derivative(zz -> (du = zeros(promote_type(eltype(u), typeof(zz)), length(u)); recombination_ode!(du, u, p, zz); du), z); return nothing
end
# explicit FullSpecialize (caller callback only; the library has no solver default): avoids SciMLBase's AutoSpecialize FunctionWrappersWrapper around the RHS
const ODEFUN5 = ODEFunction{true, SciMLBase.FullSpecialize}(recombination_ode!; jac = jac5!, tgrad = tgrad5!)

"""abstol per component for the packed 12-/7-state: rho, X1s: `a1`; H excited: `aex`; He 1s: `a1`; He excited: `aex`."""
abstol5(n; a1 = 1.0e-12, aex = 1.0e-30) = n == 12 ? [a1, a1, fill(aex, 5)..., a1, fill(aex, 4)...] : [a1, a1, fill(aex, 5)...]

# Per-solve prepared state Jacobian (caller callback only). Built from the ACTUAL u0/p/z0 of one solve (Float64, Dual, nested Dual, BigFloat element types as
# supplied); every call writes the LIVE p and z of that call into the callable before differentiating (nothing frozen). Calls whose types differ from the
# prepared ones (e.g. a Dual z into a Float64-z cache) use the unprepared oracle jac5! (counted in `nfallback`), never a conversion. One instance per
# solve: no global workspace, no sharing across solves/chains. Same chunk size ForwardDiff picks by default (n <= 12).
mutable struct RHSAt5{P,Z}
    p::P
    z::Z
end
(r::RHSAt5)(du, u) = (recombination_ode!(du, u, r.p, r.z); nothing)
struct PreparedJac5{R<:RHSAt5,C,Y}
    f::R
    cfg::C
    y::Y
    nfallback::Base.RefValue{Int}
end
function PreparedJac5(u0, p, z0)
    f = RHSAt5(p, z0); y = similar(u0)
    return PreparedJac5(f, ForwardDiff.JacobianConfig(f, y, u0, ForwardDiff.Chunk{length(u0)}()), y, Ref(0))
end
function (j::PreparedJac5{RHSAt5{P,Z}})(Jm, u, p, z) where {P,Z}
    if p isa P && z isa Z && typeof(u) === typeof(j.y)
        j.f.p = p; j.f.z = z
        ForwardDiff.jacobian!(Jm, j.f, j.y, u, j.cfg)
    else
        j.nfallback[] += 1
        jac5!(Jm, u, p, z)
    end
    return nothing
end
# Per-solve PRIMAL RHS workspace (Step 3 increment 1; private CosmoRec._rhs_workspace / _recombination_rhs_ws!, same result as recombination_ode!).
# Built from the solve's actual u0/p/z0 types; a call with other argument types runs the public allocating path (counted in `nfallback`). The state
# Jacobian (PreparedJac5) keeps using the allocating recombination_ode!.
struct RHSWS5{W}
    ws::W
    nfallback::Base.RefValue{Int}
end
RHSWS5(u0, p, z0) = RHSWS5(CosmoRec._rhs_workspace(u0, z0, p.rm; hscale = p.hscale, nbscale = p.nbscale), Ref(0))
function (r::RHSWS5)(du, u, p, z)
    CosmoRec._recombination_rhs_ws!(du, z, u, p.rm, r.ws; flag_He = p.flag, hscale = p.hscale, nbscale = p.nbscale) || (r.nfallback[] += 1)
    return nothing
end
# Per-solve BUFFERED state Jacobian (Step 3 increment 2, docs/ODE_CORE_WORKSPACE_DESIGN.md §8–§9): the differentiated functor writes through a private
# Dual-typed RHS workspace (CosmoRec._rhs_workspace / _recombination_rhs_ws!) with the LIVE p and z of each call. Its tag is the FAMILY tag
# Tag{BufRHS5,V}: it passes ForwardDiff's DEFAULT tag check (f isa BufRHS5) and avoids the recursive type a concrete Tag{BufRHS5{W},V} would need, since
# W holds Duals of that tag. Ordering: ForwardDiff's internal `tagcount` is triggered once at construction, exactly as `ForwardDiff.Tag(f, V)` does
# for concrete callables (an explicit dependency on unexported ForwardDiff internals, isolated in `bufrhs5_tag`). Built per solve from the actual
# u0/p/z0 types; a call with other types uses the unprepared oracle jac5! (counted), never a conversion. PreparedJac5 and jac5! remain the oracles.
mutable struct BufRHS5{W,P,Z}
    ws::W
    p::P
    z::Z
    nfallback::Base.RefValue{Int}
end
function (b::BufRHS5)(du, u)
    CosmoRec._recombination_rhs_ws!(du, b.z, u, b.p.rm, b.ws; flag_He = b.p.flag, hscale = b.p.hscale, nbscale = b.p.nbscale) || (b.nfallback[] += 1)
    return nothing
end
function bufrhs5_tag(::Type{V}) where {V}
    TT = ForwardDiff.Tag{BufRHS5,V}
    ForwardDiff.tagcount(TT)                         # internal ForwardDiff API: what Tag(f, V) does for a concrete f
    return TT()
end
struct WSJac5{B<:BufRHS5,C,Y}
    f::B
    cfg::C
    y::Y
    nfallback::Base.RefValue{Int}
end
function prepare_wsjac5(u0, p, z0)                   # factory (distinct name: no clash with the default constructor)
    V = eltype(u0); n = length(u0); tag = bufrhs5_tag(V)
    ws = CosmoRec._rhs_workspace(ForwardDiff.Dual{typeof(tag),V,n}.(u0), z0, p.rm; hscale = p.hscale, nbscale = p.nbscale)
    f = BufRHS5(ws, p, z0, Ref(0)); y = similar(u0)
    return WSJac5(f, ForwardDiff.JacobianConfig(f, y, u0, ForwardDiff.Chunk{n}(), tag), y, Ref(0))
end
function (j::WSJac5{BufRHS5{W,P,Z}})(Jm, u, p, z) where {W,P,Z}
    if p isa P && z isa Z && typeof(u) === typeof(j.y)
        j.f.p = p; j.f.z = z
        ForwardDiff.jacobian!(Jm, j.f, j.y, u, j.cfg)      # default tag check
    else
        j.nfallback[] += 1
        jac5!(Jm, u, p, z)
    end
    return nothing
end
wsjac5_fallbacks(j::WSJac5) = j.nfallback[] + j.f.nfallback[]
# Per-solve PREPARED time gradient (caller callback only; docs/ODE_CORE_WORKSPACE_DESIGN.md pattern of §8 applied to d f/d z): ONE Dual-redshift RHS workspace and a scalar
# ForwardDiff.DerivativeConfig built from the solve's ACTUAL u0/p/z0 types; the mutating `ForwardDiff.derivative!` differentiates the functor `BufTG5`, which writes the LIVE u and p
# of each call (z arrives as the Dual argument). The Dual output element type is read from `promote_type(eltype(u0), Dual{tag,typeof(z0),1})`, the type the unprepared tgrad5! gives
# its output vector (a promotion that is not of the form Dual{tag,Y,1} makes every call take the oracle). Tag: FAMILY tag Tag{BufTG5,V} (as 5h's BufRHS5), which keeps the
# workspace type free of itself; V = promote of ALL input element types (state, z, hscale, nbscale, feedback), so any outer Dual's tag type is part of V and the tagcount ForwardDiff assigns
# at first construction of Tag{BufTG5,V} orders above every Dual already present in u0/p/z0 (a V that omitted e.g. the state's tag would reuse a stale, lower count). For a scalar x ForwardDiff's default `checktag` falls to its permissive method (no throw); checktag(cfg, ...) is kept as the pinned default.
# A call whose u/p/z types differ from the prepared ones runs the unprepared oracle tgrad5! (counted in `nfallback`), never a conversion. No state is shared across solves.
mutable struct BufTG5{W,U,P}
    ws::W
    u::U
    p::P
    nfallback::Base.RefValue{Int}
end
function (b::BufTG5)(du, zd)
    CosmoRec._recombination_rhs_ws!(du, zd, b.u, b.p.rm, b.ws; flag_He = b.p.flag, hscale = b.p.hscale, nbscale = b.p.nbscale) || (b.nfallback[] += 1)
    return nothing
end
function buftg5_tag(::Type{V}) where {V}
    TT = ForwardDiff.Tag{BufTG5,V}
    ForwardDiff.tagcount(TT)                         # internal ForwardDiff API: what Tag(f, V) does for a concrete f
    return TT()
end
struct WSTGrad5{B,C,Y,Z}
    f::B
    cfg::C
    y::Y
    nfallback::Base.RefValue{Int}
    ncalls::Base.RefValue{Int}
end
function prepare_wstgrad5(u0, p, z0)                 # factory (distinct name: no clash with the default constructor)
    Z = typeof(z0)
    V = promote_type(eltype(u0), Z, typeof(p.hscale), typeof(p.nbscale), CosmoRec.feedback_eltype(p.rm.diffusion))      # the family tag's V carries EVERY input element type (see above)
    tag = buftg5_tag(V); T = typeof(tag)
    PT = promote_type(eltype(u0), ForwardDiff.Dual{T,Z,1})
    if !(PT <: ForwardDiff.Dual{T} && PT === ForwardDiff.Dual{T,PT.parameters[2],1})
        return WSTGrad5{Nothing,Nothing,Nothing,Z}(nothing, nothing, nothing, Ref(0), Ref(0))      # unsupported promotion: every call takes the oracle
    end
    Y = PT.parameters[2]
    ws = CosmoRec._rhs_workspace(u0, ForwardDiff.Dual{T}(z0, one(z0)), p.rm; hscale = p.hscale, nbscale = p.nbscale)
    f = BufTG5(ws, u0, p, Ref(0)); y = Vector{Y}(undef, length(u0))
    cfg = ForwardDiff.DerivativeConfig(f, y, z0, tag)
    return WSTGrad5{typeof(f),typeof(cfg),typeof(y),Z}(f, cfg, y, Ref(0), Ref(0))
end
function (g::WSTGrad5{BufTG5{W,U,P},C,Y,Z})(dT, u, p, z) where {W,U,P,C,Y,Z}
    g.ncalls[] += 1
    if typeof(u) === U && p isa P && z isa Z
        g.f.u = u; g.f.p = p
        ForwardDiff.derivative!(dT, g.f, g.y, z, g.cfg)      # default tag check
    else
        g.nfallback[] += 1
        tgrad5!(dT, u, p, z)
    end
    return nothing
end
function (g::WSTGrad5{Nothing})(dT, u, p, z)
    g.ncalls[] += 1; g.nfallback[] += 1
    tgrad5!(dT, u, p, z)
    return nothing
end
wstgrad5_fallbacks(g::WSTGrad5) = g.nfallback[] + (g.f === nothing ? 0 : g.f.nfallback[])
odefun5(u0, p, z0) = ODEFunction{true, SciMLBase.FullSpecialize}(RHSWS5(u0, p, z0); jac = prepare_wsjac5(u0, p, z0), tgrad = prepare_wstgrad5(u0, p, z0))
prob_phase5(u0, z0, znodes, p) = ODEProblem{true, SciMLBase.FullSpecialize}(odefun5(u0, p, z0), u0, (z0, znodes[end]), p)
# Defaults = the measured forward config (-50% stage-B wall vs the legacy tight config; CMB Tier A verified vs native-history CAMB, l = 2..10000): Rodas4P, reltol 1e-10,
# a1 1e-18, aex 1e-12 (the literal), solver steps forced onto the save nodes (tstops = znodes). The tight legacy config is alg = Rodas5P(), reltol 1e-12, a1 1e-18,
# aex 1e-14, tstops = false; every fixture-gated regression test pins it explicitly (the trajectories of those tests do not depend on these defaults).
function solve_phase5(f!, u0, z0, znodes, p; alg = Rodas4P(), reltol = 1.0e-10, a1 = 1.0e-18, aex = 1.0e-12, tstops = true, sensealg = nothing)
    prob = prob_phase5(u0, z0, znodes, p)
    kw = (; reltol = reltol, abstol = abstol5(length(u0); a1 = a1, aex = aex), saveat = znodes, internalnorm = primal_norm)
    tstops && (kw = merge(kw, (; tstops = znodes)))
    sol = sensealg === nothing ? solve(prob, alg; kw...) : solve(prob, alg; kw..., sensealg = sensealg)
    SciMLBase.successful_retcode(sol) || error("ODE solve failed: retcode = $(sol.retcode) after $(length(sol.u)) saved nodes (z0 = $z0)")
    return reduce(hcat, sol.u)
end

# ---- Chunk 5b fixtures: one native pass (runmode 1) ----
const FIX5B = Dict("nodes" => joinpath(@__DIR__, "fixtures", "native_ode_pass_nodes.txt"), "output" => joinpath(@__DIR__, "fixtures", "native_ode_pass_output.txt"),
                   "grid" => joinpath(@__DIR__, "fixtures", "native_ode_pass_grid.txt"), "summary" => joinpath(@__DIR__, "fixtures", "native_ode_pass_summary.txt"))
const FIX5B_SHA = Dict("nodes" => "e3f06ae0ca9f517c65adb9051cae8f576e4188d30cd9c4999e511a38fa39762a", "output" => "2a96e12d2ae215a08aaa6bf66d4579e4e205e62290cb44d8a45c0de25fe93695",
                       "grid" => "d9bf4d2af4deddaccc47acd86e2b9fe58683b8e2f2061c0c148b36da6eb9e70a", "summary" => "533a10ca0bbb31fe7cfc898118ac29242684d85346abf3cc159bc1e95d02a392")
read_table5(path) = reduce(hcat, [parse.(Float64, split(l, '\t')) for l in eachline(path) if !isempty(l) && !startswith(l, "#")])'
function read_summary5(path)
    d = Dict{String,Any}()
    for l in eachline(path)
        (isempty(l) || startswith(l, "#")) && continue
        f = split(l, '\t')
        if f[1] == "FINAL_X"; d[f[1]] = parse.(Float64, f[2:end])
        elseif f[1] == "RECFAST_INPUT"; d[f[1]] = Dict(String(f[k]) => parse(Float64, f[k + 1]) for k in 2:2:length(f))
        elseif f[1] == "HE_OFF"; d[f[1]] = parse(Float64, f[3])
        elseif length(f) == 2; d[f[1]] = parse(Float64, f[2])
        end
    end
    return d
end
const NODES5 = read_table5(FIX5B["nodes"])      # 3000 x 11: z, Xe, X1s, 2s, 2p, 3s, 3p, 3d, rho, 0, 0
const OUT5 = read_table5(FIX5B["output"])       # 3199 x 3: z, Xe, Te
const GRID5 = read_table5(FIX5B["grid"])        # 10000 x 3: z, Xe, Te
const SUM5 = read_summary5(FIX5B["summary"])

# ---- Chunk 5c: Recfast tail solve callback ----
using ForwardDiff
const RFCONST5 = NATIVE_RECFAST_CONSTANTS
function jac_tail5!(Jm, u, p, z)
    Jm .= ForwardDiff.jacobian(uu -> (du = similar(uu); recfast_tail_rhs!(du, uu, p, z); du), u); return nothing
end
function tgrad_tail5!(dT, u, p, z)
    dT .= ForwardDiff.derivative(zz -> (du = zeros(promote_type(eltype(u), typeof(zz)), length(u)); recfast_tail_rhs!(du, u, p, zz); du), z); return nothing
end
const ODEFUN_TAIL5 = ODEFunction(recfast_tail_rhs!; jac = jac_tail5!, tgrad = tgrad_tail5!)
function solve_tail5(f!, u0, z0, znodes, p; alg = Rodas5P(), reltol = 1.0e-12, abstol = [1.0e-16, 1.0e-16, 1.0e-12], sensealg = nothing)
    prob = ODEProblem(ODEFUN_TAIL5, u0, (z0, znodes[end]), p)
    kw = (; reltol = reltol, abstol = abstol, saveat = znodes, internalnorm = primal_norm)
    sol = sensealg === nothing ? solve(prob, alg; kw...) : solve(prob, alg; kw..., sensealg = sensealg)
    SciMLBase.successful_retcode(sol) || error("tail ODE solve failed: retcode = $(sol.retcode) after $(length(sol.u)) saved nodes")
    return reduce(hcat, sol.u)
end
# production Recfast parameters theta = [F, A2s1s, Yp, Omega_b, h100, T0] from the 4c fixture
theta5() = theta4c()
hfun5() = z -> cosmos_H(ACC5, z)
