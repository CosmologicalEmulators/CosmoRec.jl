# Focused bridge feasibility test (cheap CAMB settings: lmax 2000, default accuracy; recombination_model CosmoRec, runmode-0 defaults).
# Proves: real CAMB -> shim -> Julia -> CAMB in one PyCall process at OMP 1; counter; returned arrays are what CAMB uses; native unlinked/not invoked;
# unregistered call aborts (no fallback); REPEATED calls in one process are safe (GIL re-acquisition + PyCall finalizers under GC pressure --
# the no-GIL variant segfaulted in PyObject_Free from a finalizer) and deterministic (same inputs -> bitwise-identical Xe/Tm).
# NOT a CMB accuracy gate (see run_cmb_pair.sh for that).
using Test
include(joinpath(@__DIR__, "cosmorecjl_bridge.jl"))
@testset "CosmoRec.jl CAMB bridge feasibility" begin
    @test Threads.nthreads() == 1
    camb = setup_bridge()
    @test shim_ncalls() == 0 && NCALLS_BRIDGE[] == 0
    pars = camb.set_params(H0 = 67.36, ombh2 = 0.02237, omch2 = 0.12, mnu = 0.06, num_massive_neutrinos = 1, nnu = 3.046, omk = 0.0, TCMB = 2.7255,
                           YHe = 0.2454, As = 2.1e-9, ns = 0.9649, tau = 0.0544, lmax = 2000, recombination_model = "CosmoRec")
    @test pybuiltin("type")(pars."Recomb")."__name__" == "CosmoRec" && pars."Recomb"."runmode" == 0
    t = time(); res = camb.get_results(pars); wall = time() - t
    @test NCALLS_BRIDGE[] >= 1 && shim_ncalls() == NCALLS_BRIDGE[] && shim_last_status() == 0
    c = CALLS_BRIDGE[end]
    # CAMB's actual inputs reached Julia (cosmorec.f90:165-213 conversions; nothing native-produced)
    @test c.runmode == 0 && c.H0 == 67.36 && c.T_CMB == 2.7255 && c.Y_He == 0.2454 && c.Omega_b == 0.02237 / (67.36 / 100)^2 && length(c.z_Hz) == 10000
    @test c.z_out == c.z_Hz && all(isfinite, c.Xe) && all(isfinite, c.Tm)
    # CAMB uses the returned arrays: CAMB's x_e/T_b at the returned nodes 200 < z < 2500 (no reionization contribution there) equal Julia's output
    sel = findall(z -> 200 < z < 2500, c.z_out)[1:7:end]
    ev = res.get_background_redshift_evolution(c.z_out[sel], ["x_e", "T_b"], format = "array")
    rx = maximum(abs.(ev[:, 1] ./ c.Xe[sel] .- 1)); rt = maximum(abs.(ev[:, 2] ./ c.Tm[sel] .- 1))
    println("BRIDGE calls=$(NCALLS_BRIDGE[]) shim_calls=$(shim_ncalls()) julia_s=$(round.(getfield.(CALLS_BRIDGE, :seconds), digits=1)) camb_wall_s=$(round(wall, digits=1)) maxrel CAMB x_e vs returned=$rx T_b=$rt")
    # first run (a-priori 1e-9) measured 1.38e-8: CAMB re-tabulates x_e on its own thermo grid (T_b matched 9.3e-12); 1e-7 still pins the returned array
    @test rx < 1e-7
    @test rt < 1e-6
    cl = res.get_lensed_scalar_cls(CMB_unit = "muK", raw_cl = true)
    @test all(isfinite, cl) && cl[1001, 1] > 0
    # repeated calls in ONE process (A-B-A): the GIL must be re-acquired around every callback (PyCall finalizers run inside it), and
    # the same cosmology must return bitwise-identical output (single thread, deterministic solver)
    x1, t1 = c.Xe, c.Tm
    GC.gc()
    parsB = camb.set_params(H0 = 50.0, ombh2 = 0.025, omch2 = 0.16, mnu = 0.0, num_massive_neutrinos = 0, nnu = 3.046, omk = 0.0, TCMB = 2.7255,
                            YHe = 0.2163021707042025, As = 2.1e-9, ns = 1.05, tau = 0.15, lmax = 2000, recombination_model = "CosmoRec")
    camb.get_results(parsB)
    GC.gc()
    camb.get_results(pars)
    @test NCALLS_BRIDGE[] == 3 && shim_ncalls() == 3 && shim_last_status() == 0
    cB, c3 = CALLS_BRIDGE[end-1], CALLS_BRIDGE[end]
    @test cB.H0 == 50.0 && cB.Xe != x1                  # the B call really was the other cosmology
    @test c3.Xe == x1 && c3.Tm == t1                    # A repeated: bitwise-identical output
    # native CosmoRec neither linked nor loaded; original camblib not mapped
    libs = mapped_libs()
    @test !any(l -> occursin("libgsl", l) || occursin("libCosmoRec", l) || occursin("tools/CAMB-cosmorec/camb/camblib.so", l), libs)
    @test any(l -> l == CAMBLIB, libs)
    proof = read(joinpath(@__DIR__, "link_proof.txt"), String)
    @test occursin("(expect none)\n0\n", proof) && !occursin("libgsl", split(proof, "(expect no libgsl)")[2])
    # no fallback: the same camblib without a registered Julia entry aborts in a plain Python process
    py = PyCall.python
    err = IOBuffer()
    p = run(pipeline(ignorestatus(addenv(`$py -c "import camb; p = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.12, lmax=100, recombination_model='CosmoRec'); camb.get_background(p)"`,
                                  "PYTHONPATH" => CAMB_BRIDGE_PY, "OMP_NUM_THREADS" => "1")), stderr = err, stdout = devnull))
    msg = String(take!(err))
    @test !success(p) && occursin("no Julia entry is registered", msg)
    println("BRIDGE unregistered subprocess: exit=$(p.exitcode) termsignal=$(p.termsignal) msg=", first(split(msg, '\n')))
end
