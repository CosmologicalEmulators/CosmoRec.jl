"""Functional replacement gate (Marco 2026-10-11): pointwise |C_Julia/C_native - 1| < 5e-4 for LENSED TT and EE at EVERY 2 <= ell <= 9500.
Unlensed TT/EE are reported against the same threshold (informational). TE/BB/PP: diagnostics only, NO threshold (TE is sign-changing:
correlation-normalized |dTE|/sqrt(TT*EE), Marco's Knox-style preference; never dTE/TE).
EXIT CONTRACT (run_cmb_pair.sh relies on it): 0 = pass, 1 = gate FAIL, 2 = broken comparison (missing/mismatched/non-finite outputs).
usage: gate.py <dir> <tagJulia> <tagNative>"""
import sys, os, json, datetime
import numpy as np

def broken(msg):
    print(f"GATE BROKEN: {msg}")
    sys.exit(2)

d, A, B = sys.argv[1:4]
fa, fb = os.path.join(d, f"{A}_cls.txt"), os.path.join(d, f"{B}_cls.txt")
for f in (fa, fb):
    if not os.path.isfile(f):
        broken(f"missing {f}")
a = np.loadtxt(fa); b = np.loadtxt(fb)
if a.ndim != 2 or a.shape[1] != 9 or b.shape != a.shape:
    broken(f"bad shape {a.shape} vs {b.shape} (expect (9501, 9))")
if not np.array_equal(a[:, 0], b[:, 0]):
    broken("ell grids differ")
if a[0, 0] != 0 or a[-1, 0] != 9500 or not np.array_equal(a[:, 0], np.arange(9501)):
    broken("ell grid must be exactly 0..9500 step 1")
m = a[:, 0] >= 2
ell = a[m, 0]
out = {"A": A, "B": B, "when": datetime.datetime.now().isoformat(),
       "gate": "lensed TT and EE pointwise max|C_J/C_N-1| < 5e-4 at every ell 2..9500"}
fail = False
for k, c in (("unlensed_TT", 1), ("unlensed_EE", 2), ("lensed_TT", 4), ("lensed_EE", 5)):
    x, y = a[m, c], b[m, c]
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        broken(f"non-finite values in {k}")
    if np.any(y <= 0):
        broken(f"non-positive reference values in {k}")
    r = np.abs(x / y - 1.0)
    i = int(np.argmax(r))
    ok = bool(r.max() < 5e-4)
    if k.startswith("lensed") and not ok:
        fail = True
    out[k] = {"max_frac": float(r.max()), "argmax_ell": int(ell[i]), "pass_5e-4": ok}
    print(f"GATE {A} {k}: max|C_J/C_N-1| = {r.max():.3e} at ell {int(ell[i])}  -> {'PASS' if ok else 'FAIL'}")
# diagnostics only, no threshold
for k, c, tt, ee in (("unlensed_TE", 3, 1, 2), ("lensed_TE", 7, 4, 5)):
    dte = np.abs(a[m, c] - b[m, c])
    scale = np.sqrt(np.abs(b[m, tt] * b[m, ee]))
    out[k] = {"max_corrnorm": float((dte / scale).max()), "max_abs_muK2": float(dte.max())}
    print(f"GATE {A} {k} (diagnostic): max|dTE|/sqrt(TT*EE) = {out[k]['max_corrnorm']:.3e}  max|dTE| = {out[k]['max_abs_muK2']:.3e} muK^2")
for k, c in (("lensed_BB", 6), ("phiphi", 8)):
    x, y = a[m, c], b[m, c]
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        broken(f"non-finite values in {k}")
    rec = {"max_abs": float(np.abs(x - y).max())}
    if np.all(y > 0):
        rec["max_frac"] = float(np.abs(x / y - 1.0).max())
    out[k] = rec
    print(f"GATE {A} {k} (diagnostic): {rec}")
out["GATE_lensed_TT_EE"] = not fail
print(f"GATE {A} lensed TT&EE all ell 2..9500: {'PASS' if not fail else 'FAIL'}")
with open(os.path.join(d, "gate.jsonl"), "a") as fh:
    fh.write(json.dumps(out) + "\n")
sys.exit(1 if fail else 0)
