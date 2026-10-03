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
using OrdinaryDiffEqRosenbrock: Rodas5P

function jac5!(Jm, u, p, z)
    Jm .= ForwardDiff.jacobian(uu -> (du = similar(uu); recombination_ode!(du, uu, p, z); du), u); return nothing
end
function tgrad5!(dT, u, p, z)
    dT .= ForwardDiff.derivative(zz -> (du = zeros(promote_type(eltype(u), typeof(zz)), length(u)); recombination_ode!(du, u, p, zz); du), z); return nothing
end
const ODEFUN5 = ODEFunction(recombination_ode!; jac = jac5!, tgrad = tgrad5!)

"""abstol per component for the packed 12-/7-state: rho, X1s: `a1`; H excited: `aex`; He 1s: `a1`; He excited: `aex`."""
abstol5(n; a1 = 1.0e-12, aex = 1.0e-30) = n == 12 ? [a1, a1, fill(aex, 5)..., a1, fill(aex, 4)...] : [a1, a1, fill(aex, 5)...]

function solve_phase5(f!, u0, z0, znodes, p; alg = Rodas5P(), reltol = 1.0e-9, a1 = 1.0e-12, aex = 1.0e-30, sensealg = nothing)
    prob = ODEProblem(ODEFUN5, u0, (z0, znodes[end]), p)
    kw = (; reltol = reltol, abstol = abstol5(length(u0); a1 = a1, aex = aex), saveat = znodes, internalnorm = primal_norm)
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
