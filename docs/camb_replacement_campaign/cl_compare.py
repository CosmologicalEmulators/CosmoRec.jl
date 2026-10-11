"""Compare two camb_forward.py outputs. usage: cl_compare.py <dir> <tagA> <tagB> [label]   (A = test, B = reference)
Per spectrum over 2 <= ell <= 3000 (unlensed TT EE TE; lensed TT EE BB TE; phiphi):
  max_frac   = max |dC|/|C_B|            (TT, EE, BB, phiphi; not reported for TE, which has zeros)
  max_frac_TE_corr = max |dC_TE|/sqrt(C_TT C_EE)  (B spectra, same lensing state): correlation-normalized, NOT dTE/TE (TE crosses zero)
  max_abs    = max |dC|  [muK^2, raw C_ell]
  max_cv / chi2_cv = max |dC|/sigma_CV and sum (dC/sigma_CV)^2, full-sky cosmic variance only (no noise, beam = 1, fsky = 1, every ell):
                     sigma = sqrt(2/(2l+1)) |C| for TT EE BB phiphi;  sigma_TE = sqrt((C_TE^2 + C_TT C_EE)/(2l+1)) (finite at TE zeros).
  chi2_cv is a DIAGONAL per-spectrum diagnostic sum_l (dC/sigma)^2 -- NOT a joint TT/TE/EE (or lensing) likelihood chi^2: TT/EE/TE covariances
  and the non-Gaussian lensed-BB covariance are ignored.
Also CAMB's interpolated x_e, T_b (z grid of camb_forward.py): max relative difference per z range.
Prints one block; also appends a JSON line to <dir>/compare.jsonl."""
import sys, json, numpy as np
d, A, B = sys.argv[1], sys.argv[2], sys.argv[3]; label = sys.argv[4] if len(sys.argv) > 4 else f"{A} vs {B}"
a = np.loadtxt(f"{d}/{A}_cls.txt"); b = np.loadtxt(f"{d}/{B}_cls.txt"); assert np.array_equal(a[:, 0], b[:, 0])
m = a[:, 0] >= 2; ell = a[m, 0]
cols = {"uTT": 1, "uEE": 2, "uTE": 3, "lTT": 4, "lEE": 5, "lBB": 6, "lTE": 7, "PP": 8}
out = {"label": label, "A": A, "B": B, "bitwise_cls": bool(np.array_equal(a, b))}
print(f"=== {label}: bitwise identical C_ell arrays: {out['bitwise_cls']}   (* TE: correlation-normalized, not dTE/TE; diag sum = per-spectrum diagonal, not a joint chi^2)")
for k, c in cols.items():
    x, y = a[m, c], b[m, c]; dC = x - y
    if k.endswith("TE"):
        tt, ee = b[m, cols[k[0] + "TT"]], b[m, cols[k[0] + "EE"]]
        sig = np.sqrt((y**2 + tt * ee) / (2 * ell + 1)); frac = np.abs(dC) / np.sqrt(tt * ee); fl = "|dTE|/sqrt(TTEE)*"
    else:
        sig = np.sqrt(2 / (2 * ell + 1)) * np.abs(y); frac = np.abs(dC) / np.abs(y); fl = "max|dC|/|C|"
    cv = np.abs(dC) / sig
    r = {"max_frac": float(frac.max()), "ell_frac": int(ell[np.argmax(frac)]), "max_abs": float(np.abs(dC).max()),
         "max_cv": float(cv.max()), "ell_cv": int(ell[np.argmax(cv)]), "chi2_cv": float((cv**2).sum())}
    out[k] = r
    print(f"  {k:4s} {fl:20s} {r['max_frac']:.3e} (l={r['ell_frac']:4d})  max|dC| {r['max_abs']:.3e}  max|dC|/sigma_CV {r['max_cv']:.3e} (l={r['ell_cv']:4d})  diag sum(dC/sigma)^2 {r['chi2_cv']:.3e}")
xa = np.loadtxt(f"{d}/{A}_xe.txt"); xb = np.loadtxt(f"{d}/{B}_xe.txt"); z = xb[:, 0]
for nm, c in (("x_e", 1), ("T_b", 2)):
    rr = np.abs(xa[:, c] / xb[:, c] - 1); s = []
    for lo, hi in ((0, 10), (10, 50), (50, 500), (500, 1600), (1600, 3000.01)):
        mm = (z >= lo) & (z < hi); s.append(f"[{lo},{hi:.0f}):{rr[mm].max():.2e}")
        out[f"{nm}_maxrel_{lo}_{int(hi)}"] = float(rr[mm].max())
    print(f"  CAMB {nm}: max rel " + " ".join(s))
open(f"{d}/compare.jsonl", "a").write(json.dumps(out) + "\n")
