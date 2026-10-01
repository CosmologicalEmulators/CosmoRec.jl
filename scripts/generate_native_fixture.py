#!/usr/bin/env python3
"""Capture the FINAL public CAMB thermodynamic history produced by the original CosmoRec solver.

Run inside the existing environment (never installs or builds anything):
    source /home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/tools/activate-camb-cosmorec.sh
    python scripts/generate_native_fixture.py generate OUT.txt
    python scripts/generate_native_fixture.py compare A.txt B.txt

Outputs are OUTPUTS ONLY (x_e, T_b, H after CAMB's own resampling); no CosmoRec internal arrays are exposed by the
public CAMB API. Parameter construction is the production `camb_worker._build_params` (read-only import).
"""
import hashlib, os, subprocess, sys
import numpy as np

TOOLS = "/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/tools"
CAMB_ROOT = f"{TOOLS}/CAMB-cosmorec"
COSMOREC_ROOT = f"{TOOLS}/CosmoRec"
WORKER_DIR = "/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/emulator-zoo/Capse.jl/camb_mnuw0wacdm"
PINNED = {CAMB_ROOT: "fa3f097343fbbe427cc04b4f5f0041c22c6ec764", COSMOREC_ROOT: "086769055f61ae0c244a53dd381ee65b624d0ac3"}
PARAMS = dict(ln10As=3.044, ns=0.9649, tau=0.0568, H0=67.36, omega_b=0.02237, omega_c=0.12, Mnu=0.06, w0=-1.0, wa=0.0)
C_KMS = 299792.458
COLUMNS = ["z", "role", "x_e", "T_b", "H", "x_e_noreion", "T_b_noreion"]


def query_grid():
    """Ascending z in [0, 1e4]; role 'grid' = interpolation nodes, 'holdout' = off-node native queries."""
    def rng(a, b, step):
        return list(np.round(np.arange(a, b + 0.5 * step, step), 10))
    grid = set(rng(0, 3, 0.25)) | set(rng(3.5, 10, 0.5)) | set(rng(12.5, 50, 2.5)) | set([15.0, 20.0, 30.0, 40.0])
    grid |= set(rng(60, 200, 10)) | set([250.0, 300.0, 400.0, 500.0])
    grid |= set(rng(550, 1500, 25))            # hydrogen recombination (dense)
    grid |= set(rng(1500, 2500, 50))           # He I recombination region
    grid |= set(rng(2500, 4000, 100))
    grid |= set(rng(4000, 7500, 100))          # He II -> He I transition region (z ~ 6000)
    grid |= set(rng(7500, 10000, 250))
    grid |= {3000.0, 200.0, 50.0, 9999.0, 10000.0}
    grid = sorted(grid)
    hold = [0.05, 0.17, 0.875, 1.625, 2.125, 4.25, 8.75, 11.0, 26.0, 36.0, 46.0, 55.0, 75.0, 105.0, 175.0, 225.0, 275.0, 350.0, 450.0, 520.0,
            562.5, 612.5, 662.5, 712.5, 762.5, 812.5, 862.5, 912.5, 962.5, 1012.5, 1062.5, 1112.5, 1162.5, 1212.5, 1262.5,
            1312.5, 1362.5, 1412.5, 1462.5, 1525.0, 1675.0, 1825.0, 1975.0, 2125.0, 2275.0, 2425.0, 2550.0, 2950.0, 3450.0,
            4050.0, 4550.0, 5050.0, 5550.0, 5950.0, 6050.0, 6550.0, 7050.0, 7450.0, 7625.0, 8125.0, 8875.0, 9625.0, 9990.0]
    assert not set(hold) & set(grid)
    z = np.array(sorted(set(grid) | set(hold)))
    role = ["grid" if v in set(grid) else "holdout" for v in z]
    return z, role


def sha_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def git(root, *args):
    return subprocess.run(["git", "-C", root, *args], check=True, capture_output=True, text=True).stdout.strip()


def provenance():
    import camb
    for root, sha in PINNED.items():
        got = git(root, "rev-parse", "HEAD")
        assert got == sha, f"{root}: HEAD {got} != pinned {sha}"
        assert git(root, "status", "--porcelain") == "", f"{root} has uncommitted changes"
    assert os.path.realpath(camb.__file__).startswith(os.path.realpath(CAMB_ROOT) + "/camb/"), camb.__file__
    assert os.environ.get("CAMB_COSMOREC_ROOT") == CAMB_ROOT and os.environ.get("COSMOREC_ROOT") == COSMOREC_ROOT
    lib = os.path.join(os.path.dirname(camb.__file__), "camblib.so")
    return {"camb_version": camb.__version__, "camb_file": os.path.realpath(camb.__file__), "camblib_sha256": sha_file(lib),
            "camb_cosmorec_sha": PINNED[CAMB_ROOT], "cosmorec_sha": PINNED[COSMOREC_ROOT],
            "python": sys.version.split()[0], "numpy": np.__version__,
            "camb_worker_sha256": sha_file(f"{WORKER_DIR}/camb_worker.py"), "omp_num_threads": os.environ.get("OMP_NUM_THREADS", "unset")}


def run_config(z, reionization):
    import camb
    sys.path.insert(0, WORKER_DIR)
    import camb_worker as cw
    pars = cw._build_params(PARAMS, cw.OUTPUT_LMAX)
    if not reionization:
        pars.Reion.Reionization = False
    rec = pars.Recomb
    assert type(rec).__name__ == "CosmoRec", type(rec).__name__
    cfg = {k: getattr(rec, k) for k, *_ in rec._fields_}
    assert cfg == {"runmode": 0, "fdm": 0.0, "accuracy": 0.0, "diff_iteration_max": -1, "n_shells": -1, "n_shells_hei": -1,
                   "flag_hi_absorption": -1, "n_s_2gamma": -1, "n_s_raman": -1}, cfg
    res = camb.get_background(pars)
    # ONE call with the whole grid: CAMB initialises its thermo table from min(times), so the query set matters.
    d = res.get_background_redshift_evolution(z, vars=["x_e", "T_b"])
    H = res.hubble_parameter(z)
    # independent H: H^2[1/Mpc^2] = (8 pi G a^4 rho_tot)/(3 a^4) from CAMB's densities; c in km/s
    a = 1.0 / (1.0 + z)
    tot = res.get_background_densities(a, vars=["tot"])["tot"]
    H_ind = C_KMS * np.sqrt(tot / 3.0) / a**2
    assert np.max(np.abs(H / H_ind - 1)) < 1e-10, np.max(np.abs(H / H_ind - 1))
    assert abs(res.hubble_parameter(0.0) / PARAMS["H0"] - 1) < 1e-12
    return dict(x_e=np.array(d["x_e"]), T_b=np.array(d["T_b"]), H=np.array(H), cfg=cfg, YHe=pars.YHe,
                bbn=type(pars.bbn_predictor).__name__ + ":" + str(getattr(pars.bbn_predictor, "file", "")),
                Hdiff=float(np.max(np.abs(H / H_ind - 1))), tau=pars.Reion.optical_depth if reionization else None)


def generate(out):
    prov = provenance()
    z, role = query_grid()
    A = run_config(z, True)
    B = run_config(z, False)
    for X in (A, B):
        for k in ("x_e", "T_b", "H"):
            assert np.all(np.isfinite(X[k])) and np.all(X[k] > 0), k
    rows = [f"{z[i]:.17g} {role[i]} {A['x_e'][i]:.17g} {A['T_b'][i]:.17g} {A['H'][i]:.17g} {B['x_e'][i]:.17g} {B['T_b'][i]:.17g}"
            for i in range(len(z))]
    body = "\n".join(rows) + "\n"
    head = [
        "# CosmoRec.jl native fixture: FINAL public CAMB thermodynamic history from the original CosmoRec solver (OUTPUTS ONLY)",
        "# format: 1; columns whitespace separated; one header line then rows; '#' lines are metadata",
        f"# source_camb_cosmorec_sha: {prov['camb_cosmorec_sha']}", f"# source_cosmorec_sha: {prov['cosmorec_sha']}",
        f"# camb_version: {prov['camb_version']}", f"# camb_file: {prov['camb_file']}", f"# camblib_sha256: {prov['camblib_sha256']}",
        f"# python: {prov['python']}", f"# numpy: {prov['numpy']}", f"# omp_num_threads_env: {prov['omp_num_threads']}",
        f"# params_builder: camb_worker._build_params (sha256 {prov['camb_worker_sha256']}) lmax={9500}",
        "# params: " + " ".join(f"{k}={v!r}" for k, v in PARAMS.items()) + " omk=0.0 TCMB=2.7255 nnu=3.046 num_massive_neutrinos=1 dark_energy=ppf",
        f"# YHe: {A['YHe']!r} (BBN consistency: {A['bbn']})",
        "# recombination_model: CosmoRec (type(pars.Recomb).__name__ asserted) fields: " + " ".join(f"{k}={v}" for k, v in A["cfg"].items()),
        "# cosmorec_settings: all CAMB-adapter defaults (runmode 0 with diffusion; accuracy 0.0; -1 overrides = CosmoRec batch/preset defaults, not overridden here)",
        "# cosmorec_native_zgrid_in_adapter: 10000 uniform z from 10000 down to 0; H array from CAMB dtauda; native Hubble endpoint behavior is whatever the compiled build does (not patched, not asserted)",
        "# api: camb.get_background(pars).get_background_redshift_evolution(z, vars=[x_e,T_b]) in ONE call; this is CAMB's resampled thermo history, not raw CosmoRec samples",
        "# order: z ascending; role=grid are interpolation nodes, role=holdout are independent off-node native queries (same call)",
        "# units: z dimensionless; x_e = n_e/n_H as returned by CAMB (He included, >1 when He fully ionised or reionised); T_b in K; H in km/s/Mpc",
        f"# H_check: hubble_parameter agrees with c*sqrt(8piG a^4 rho_tot/3)/a^2 from get_background_densities to {A['Hdiff']:.3e} relative; H(0)=H0 asserted to 1e-12",
        "# columns x_e,T_b,H: fiducial incl. reionization (tau=0.0568, tanh); x_e_noreion,T_b_noreion: same run with Reion.Reionization=False (pure recombination/adiabatic history at low z)",
        "# NOT captured: any CosmoRec internal arrays (level populations, rates, PDE fields, corrections); only these final outputs",
        f"# n_rows: {len(rows)}", f"# sha256_rows: {hashlib.sha256(body.encode()).hexdigest()}",
        " ".join(COLUMNS)]
    with open(out, "w") as f:
        f.write("\n".join(head) + "\n" + body)
    print(f"wrote {out}: {len(rows)} rows sha256_rows={hashlib.sha256(body.encode()).hexdigest()}")


def read(path):
    meta, rows, cols = {}, [], None
    for line in open(path):
        if line.startswith("#"):
            k, _, v = line[1:].strip().partition(":")
            meta[k.strip()] = v.strip()
        elif cols is None:
            cols = line.split()
        else:
            rows.append(line.split())
    return meta, cols, rows


def compare(pa, pb):
    ma, ca, ra = read(pa); mb, cb, rb = read(pb)
    assert ca == cb and len(ra) == len(rb)
    worst = (0.0, None); nbad = 0
    for i, (x, y) in enumerate(zip(ra, rb)):
        for j, name in enumerate(ca):
            if x[j] != y[j]:
                nbad += 1
                if name != "role":
                    rel = abs(float(x[j]) - float(y[j])) / max(abs(float(x[j])), 1e-300)
                    if rel > worst[0]: worst = (rel, (i, name, x[j], y[j]))
    print(f"rows {len(ra)}; differing fields {nbad}; worst relative difference {worst}")
    print("metadata identical:", ma == mb, "; files byte-identical:", open(pa, "rb").read() == open(pb, "rb").read(),
          "; sha256:", sha_file(pa), sha_file(pb))
    return nbad == 0 and ma == mb


if __name__ == "__main__":
    if sys.argv[1] == "generate":
        generate(sys.argv[2])
    elif sys.argv[1] == "compare":
        sys.exit(0 if compare(sys.argv[2], sys.argv[3]) else 1)
    else:
        raise SystemExit("usage: generate OUT | compare A B")
