"""One CAMB run with CosmoRec recombination for paired native-vs-Julia comparisons (REWRITTEN 2026-10-11; replaces the lost chunk14 camb_forward.py,
whose full text was not retained -- output FORMAT kept compatible with cl_compare.py: <out>_cls.txt, <out>_xe.txt, <out>_meta.json).
usage: camb_forward_c36.py <params.json> <out_prefix> <native | dump:<file> | inject:<file>> <preset: chunk14 | worker>
  The CAMB package is whatever `import camb` resolves (PYTHONPATH selects the ORIGINAL tools/CAMB-cosmorec or the chunk36 scratch copy).
  dump/inject need the scratch build with camb_diag.diff (env CAMB_COSMOREC_DUMP / CAMB_COSMOREC_INJECT, read by fortran/cosmorec.f90 CosmoRec_diag_io).
  params.json keys: H0 ombh2 omch2 mnu omk TCMB nnu num_massive_neutrinos YHe ln10As ns tau w wa [neutrino_hierarchy standard_neutrino_neff Alens];
  YHe = null -> CAMB BBN consistency (worker classes degenerate/alens).
Presets (sources):
  chunk14: lmax out 3000, computed lmax 5000 (set_for_lmax(5000, lens_potential_accuracy=8)), kmax 10, k_per_logint 130, lAccuracyBoost 1.2,
           min_l_logl_sampling 6000, DoLateRadTruncation False, NonLinear_both + mead2020  (chunk14 settings recorded in BOLTZMANN_FORWARD_VALIDATION.md)
  worker : emulator-zoo/Capse.jl/*/camb_worker.py _build_params: lmax 9500, lens_output_margin 2050, lens_potential_accuracy 8, kmax 10,
           k_per_logint 130, AccuracyBoost 1, lSampleBoost 1, lAccuracyBoost 1.2, min_l_logl_sampling 6000, DoLateRadTruncation False, mead2020, PPF."""
import sys, os, json, time, hashlib, numpy as np
pfile, out, mode, preset = sys.argv[1:5]
if mode.startswith("dump:"): os.environ["CAMB_COSMOREC_DUMP"] = mode[5:]
elif mode.startswith("inject:"): os.environ["CAMB_COSMOREC_INJECT"] = mode[7:]
elif mode != "native": raise SystemExit("mode must be native | dump:<f> | inject:<f>")
import camb
P = json.load(open(pfile))
pars = camb.CAMBparams()
cos = dict(H0=P["H0"], ombh2=P["ombh2"], omch2=P["omch2"], omk=P["omk"], TCMB=P["TCMB"], mnu=P["mnu"],
           num_massive_neutrinos=P["num_massive_neutrinos"], nnu=P["nnu"])
if P.get("YHe") is not None: cos["YHe"] = P["YHe"]
for k in ("neutrino_hierarchy", "standard_neutrino_neff"):
    if k in P: cos[k] = P[k]
pars.set_cosmology(**cos)
pars.InitPower.set_params(As=np.exp(P["ln10As"]) * 1e-10, ns=P["ns"], pivot_scalar=0.05)
pars.Reion.set_tau(P["tau"])
pars.set_dark_energy(w=P["w"], wa=P["wa"], dark_energy_model="ppf")
if "Alens" in P: pars.Alens = float(P["Alens"])
pars.NonLinear = camb.model.NonLinear_both
pars.NonLinearModel.set_params(halofit_version="mead2020")
common = dict(cp=pars, recombination_model="CosmoRec", kmax=10, k_per_logint=130, lens_potential_accuracy=8, AccuracyBoost=1.0, lSampleBoost=1.0,
              lAccuracyBoost=1.2, min_l_logl_sampling=6000, DoLateRadTruncation=False, halofit_version="mead2020")
if preset == "worker":
    LMAX = 9500; pars = camb.set_params(lens_output_margin=2050, lmax=LMAX, **common)
elif preset == "chunk14":
    LMAX = 3000; pars = camb.set_params(lmax=LMAX + 2000, **common)   # computed lmax 5000, output 3000 (chunk14 doc); NOT byte-verified vs the lost script
else: raise SystemExit("preset must be chunk14 | worker")
assert type(pars.Recomb).__name__ == "CosmoRec", type(pars.Recomb).__name__
cfg = {k: getattr(pars.Recomb, k) for k in ("runmode", "accuracy", "fdm", "diff_iteration_max", "n_shells", "n_shells_hei", "flag_hi_absorption", "n_s_2gamma", "n_s_raman")}
assert cfg["runmode"] == 0 and cfg["accuracy"] == 0 and cfg["fdm"] == 0 and all(cfg[k] == -1 for k in list(cfg)[3:]), cfg
t = time.time(); res = camb.get_results(pars); elapsed = time.time() - t
ul = res.get_unlensed_scalar_cls(CMB_unit="muK", raw_cl=True, lmax=LMAX)
le = res.get_lensed_scalar_cls(CMB_unit="muK", raw_cl=True, lmax=LMAX)
pp = res.get_lens_potential_cls(lmax=LMAX, raw_cl=True)
ell = np.arange(LMAX + 1)
np.savetxt(out + "_cls.txt", np.column_stack([ell, ul[:, 0], ul[:, 1], ul[:, 3], le[:, 0], le[:, 1], le[:, 2], le[:, 3], pp[:, 0]]),
           header="ell  unlensed_TT unlensed_EE unlensed_TE  lensed_TT lensed_EE lensed_BB lensed_TE  phiphi   (raw C_ell; CMB in muK^2)", fmt="%.17g")
zz = np.concatenate([np.linspace(0, 50, 501), np.linspace(50.5, 3000, 5900)])
ev = res.get_background_redshift_evolution(zz, ["x_e", "T_b"], format="array")
np.savetxt(out + "_xe.txt", np.column_stack([zz, ev[:, 0], ev[:, 1]]), header="z  x_e  T_b  (CAMB's interpolated history incl. reionization)", fmt="%.17g")
meta = {"camb_file": camb.__file__, "camb_version": camb.__version__,
        "camblib_sha256": hashlib.sha256(open(os.path.join(os.path.dirname(camb.__file__), "camblib.so"), "rb").read()).hexdigest(),
        "mode": mode, "preset": preset, "params": P, "cosmorec_cfg": cfg, "derived": {"YHe": pars.YHe, "N_eff": pars.N_eff},
        "derived_params": res.get_derived_params(), "elapsed_s": elapsed, "lmax_out": LMAX, "omp_num_threads": os.environ.get("OMP_NUM_THREADS")}
json.dump(meta, open(out + "_meta.json", "w"), indent=1, sort_keys=True, default=str)
print(f"{out}: done in {elapsed:.1f} s; camb {camb.__file__}; mode {mode}; preset {preset}")
