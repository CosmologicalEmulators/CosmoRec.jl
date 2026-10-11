# CosmoRec.jl <-> CAMB bridge (chunk36 scratch). Production process model: Julia -> PyCall -> Python camb -> ctypes camblib.so (Fortran)
# -> C shim cosmorec_calc_h_cpp_ (bridge/shim) -> this Julia @cfunction, in ONE process, same thread (OMP_NUM_THREADS = 1).
# The Julia entry receives exactly the CAMB arguments (cosmorec.f90:211-213: Omega_c, Omega_b, Omega_k, N_eff(), H0, tcmb, Yhe, zrec, Hz, Nz) and
# writes Xe/Tm back into CAMB's arrays. Science: CosmoRec.jl cosmorec_calc_h (own closure, own preliminary history, runmode-0 history) with the
# package's existing production callbacks (production_callbacks.jl: preliminary Rodas5P dense 1e-12; passes solve_phase5 = Rodas4P 1e-10/1e-18/1e-12
# tstops; tail solve_tail5). Tables: user-provided original CosmoRec data tree (COSMOREC_TABLES_ROOT; read-only, never copied).
# Any exception -> status 1 -> the shim aborts the process (no fallback, no zero output).
using CosmoRec, PyCall, Libdl
include(joinpath(@__DIR__, "production_callbacks.jl"))

const BRIDGE_DIR = @__DIR__
const CAMB_BRIDGE_PY = joinpath(BRIDGE_DIR, "camb_bridge")                 # python package dir (contains camb/)
const CAMBLIB = joinpath(CAMB_BRIDGE_PY, "camb", "camblib.so")
const SRF_BRIDGE = (a...) -> solve_rodas_dense(a...; reltol = 1.0e-12, abstol = [1.0e-16, 1.0e-16, 1.0e-12])
const TABLES_BRIDGE = Ref{Any}(nothing)
const NCALLS_BRIDGE = Ref(0)
const BRIDGE_KEEP_CALLS = 16        # diagnostic capture bound (production safety: no unbounded retention of full input/output arrays)
const CALLS_BRIDGE = Any[]          # per call: inputs and returned arrays (copies); only the LAST BRIDGE_KEEP_CALLS are retained

function camb_cosmorec_entry(runmode::Ptr{Cint}, runpars::Ptr{Cdouble}, oc::Ptr{Cdouble}, ob::Ptr{Cdouble}, ok::Ptr{Cdouble}, nnu::Ptr{Cdouble},
                             h0::Ptr{Cdouble}, tcmb::Ptr{Cdouble}, yhe::Ptr{Cdouble}, zHz::Ptr{Cdouble}, Hz::Ptr{Cdouble}, nz::Ptr{Cint},
                             za::Ptr{Cdouble}, xe::Ptr{Cdouble}, tb::Ptr{Cdouble}, len::Ptr{Cint}, label::Ptr{Cint})::Cint
    # Python ctypes RELEASES the GIL around the CAMB call; Julia GC finalizers of PyCall objects (pydecref) may run inside this callback, so the GIL
    # must be re-acquired for its whole duration (observed: segfault in PyObject_Free from a finalizer without it, bridge_test_nogil.log).
    gil = ccall((:PyGILState_Ensure, PyCall.libpython), Cint, ())
    try
        n = Int(unsafe_load(nz)); m = Int(unsafe_load(len))
        zH = copy(unsafe_wrap(Array, zHz, n)); H = copy(unsafe_wrap(Array, Hz, n)); zo = copy(unsafe_wrap(Array, za, m))
        rp = copy(unsafe_wrap(Array, runpars, 21)); rm = Int(unsafe_load(runmode))
        inp = (; runmode = rm, runpars = rp, Omega_c = unsafe_load(oc), Omega_b = unsafe_load(ob), Omega_k = unsafe_load(ok), N_eff = unsafe_load(nnu),
               H0 = unsafe_load(h0), T_CMB = unsafe_load(tcmb), Y_He = unsafe_load(yhe), label = Int(unsafe_load(label)))
        t = time()
        r = cosmorec_calc_h(rm, rp, inp.Omega_c, inp.Omega_b, inp.Omega_k, inp.N_eff, inp.H0, inp.T_CMB, inp.Y_He, zH, H, zo, TABLES_BRIDGE[];
                            solve_recfast = SRF_BRIDGE, solve_pass = solve_phase5, solve_tail = solve_tail5)
        X = Vector{Float64}(r.Xe); Tm = Vector{Float64}(r.Tm)
        (length(X) == m && length(Tm) == m) || error("CosmoRec.jl returned $(length(X))/$(length(Tm)) values for $m output redshifts")
        all(x -> isfinite(x) && x > 0, X) && all(x -> isfinite(x) && x > 0, Tm) || error("CosmoRec.jl returned non-finite or non-positive Xe/Tm")
        unsafe_copyto!(xe, pointer(X), m); unsafe_copyto!(tb, pointer(Tm), m)
        NCALLS_BRIDGE[] += 1
        length(CALLS_BRIDGE) < BRIDGE_KEEP_CALLS || deleteat!(CALLS_BRIDGE, 1)
        push!(CALLS_BRIDGE, (; inp..., z_Hz = zH, Hz = H, z_out = zo, Xe = X, Tm = Tm, seconds = time() - t, thread = Threads.threadid()))
        return Cint(0)
    catch e
        showerror(stderr, e, catch_backtrace()); println(stderr)
        return Cint(1)
    finally
        ccall((:PyGILState_Release, PyCall.libpython), Cvoid, (Cint,), gil)
    end
end

"""Import the bridge CAMB through PyCall, check that the loaded camblib.so IS the shim build, register the Julia entry. Returns the camb module."""
function setup_bridge(; tables_root = ENV["COSMOREC_TABLES_ROOT"])
    get(ENV, "OMP_NUM_THREADS", "") == "1" || error("bridge runs require OMP_NUM_THREADS=1 (production setting)")
    TABLES_BRIDGE[] === nothing && (TABLES_BRIDGE[] = load_runmode0_tables(tables_root))
    sys = pyimport("sys"); pushfirst!(PyVector(sys."path"), CAMB_BRIDGE_PY)
    camb = pyimport("camb")
    realpath(camb.__file__) == realpath(joinpath(CAMB_BRIDGE_PY, "camb", "__init__.py")) || error("imported camb is not the bridge build: $(camb.__file__)")
    h = Libdl.dlopen(CAMBLIB, Libdl.RTLD_NOLOAD)          # must already be loaded by Python's ctypes (same handle, no second copy)
    h === nothing && error("camblib.so not loaded by camb import")
    ccall(Libdl.dlsym(h, :cosmorecjl_is_shim), Cint, ()) == 20261011 || error("camblib.so is not the shim build")
    fp = @cfunction(camb_cosmorec_entry, Cint, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}))
    ccall(Libdl.dlsym(h, :cosmorecjl_register), Cvoid, (Ptr{Cvoid},), fp)
    return camb
end
shim_ncalls() = ccall(Libdl.dlsym(Libdl.dlopen(CAMBLIB, Libdl.RTLD_NOLOAD), :cosmorecjl_ncalls), Clong, ())
shim_last_status() = ccall(Libdl.dlsym(Libdl.dlopen(CAMBLIB, Libdl.RTLD_NOLOAD), :cosmorecjl_last_status), Cint, ())
mapped_libs() = unique([split(l)[end] for l in eachline("/proc/self/maps") if occursin(".so", l)])

"""Run the EXISTING runner ../camb/camb_forward_c36.py in this process (runpy; sys.argv as on the command line). Returns the runner's globals."""
function run_forward_c36(pfile, out, preset)
    sys = pyimport("sys"); runpy = pyimport("runpy")
    sys."argv" = [joinpath(BRIDGE_DIR, "..", "camb", "camb_forward_c36.py"), pfile, out, "native", preset]   # "native" = runner mode without dump/inject
    return runpy.run_path(joinpath(BRIDGE_DIR, "..", "camb", "camb_forward_c36.py"), run_name = "__main__")
end
