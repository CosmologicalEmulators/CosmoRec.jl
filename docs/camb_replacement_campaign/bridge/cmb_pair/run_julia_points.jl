# Julia route for the paired CMB comparison: same runner (../camb/camb_forward_c36.py, preset worker) through the bridge CAMB in THIS process.
include(joinpath(@__DIR__, "..", "cosmorecjl_bridge.jl"))
setup_bridge()
for pt in ARGS
    n0 = NCALLS_BRIDGE[]
    run_forward_c36(joinpath(@__DIR__, "$pt.json"), joinpath(@__DIR__, "$(pt)_julia"), "worker")
    shim_ncalls() == NCALLS_BRIDGE[] || error("shim/Julia call counters differ")
    open(joinpath(@__DIR__, "$(pt)_julia_bridge.txt"), "w") do io
        println(io, "julia_calls_this_point $(NCALLS_BRIDGE[] - n0)\nshim_calls_total $(shim_ncalls())\njulia_seconds $(join([round(c.seconds, digits=2) for c in CALLS_BRIDGE[n0+1:end]], ' '))")
        c = CALLS_BRIDGE[end]; println(io, "inputs Omega_c=$(c.Omega_c) Omega_b=$(c.Omega_b) Omega_k=$(c.Omega_k) N_eff=$(c.N_eff) H0=$(c.H0) T_CMB=$(c.T_CMB) Y_He=$(c.Y_He) nz=$(length(c.z_Hz))")
    end
    println("JULIA POINT $pt done; calls $(NCALLS_BRIDGE[] - n0)")
end
