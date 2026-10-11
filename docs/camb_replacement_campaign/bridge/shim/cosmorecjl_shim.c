/* CosmoRec.jl bridge shim (chunk36 scratch). Replaces libCosmoRec.a in a SCRATCH CAMB build: exports the exact CAMB-facing symbol
 * cosmorec_calc_h_cpp_ (tools/CosmoRec/CosmoRec.h:153-158; called from CAMB fortran/cosmorec.f90:211-213) and forwards every argument
 * pointer, unchanged, to a Julia @cfunction registered at runtime, OR (input-contract mode, 2026-10-11) to the ORIGINAL native
 * CosmoRec loaded via dlopen (env COSMORECJL_NATIVE_LIB), so both providers can run behind the byte-identical camblib.
 * In BOTH modes env COSMORECJL_DUMP=<file> writes every CAMB->provider input (runmode, runpars[21], 7 scalars, z_Hz[nz],
 * Hz[nz], nz, za_in[len], len, label) as one raw dump, written by this same C code on both routes (bitwise-comparable).
 * No native CosmoRec code and no GSL are linked into camblib. There is NO fallback:
 * an unregistered call or a nonzero Julia status aborts the process (no native call, no zero/garbage output reaches CAMB). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <dlfcn.h>

typedef int (*cosmorecjl_fn)(const int *runmode, const double *runpars, const double *omega_c, const double *omega_b, const double *omega_k,
                             const double *num_nu, const double *h0, const double *t_cmb, const double *y_he, const double *z_Hz, const double *Hz,
                             const int *nz, double *za_in, double *xe_out, double *tb_out, const int *len, const int *label);

static cosmorecjl_fn g_fn = NULL;
static long g_ncalls = 0;
static int g_last_status = -1;

void cosmorecjl_register(cosmorecjl_fn fn) { g_fn = fn; }
long cosmorecjl_ncalls(void) { return g_ncalls; }
int cosmorecjl_last_status(void) { return g_last_status; }
int cosmorecjl_is_shim(void) { return 20261011; }

/* Input dump: every value CAMB hands the recombination provider. One file per call (CAMB calls once per process). */
static void cosmorecjl_dump_inputs(const int *runmode, const double *runpars, const double *omega_c, const double *omega_b,
                                   const double *omega_k, const double *num_nu, const double *h0, const double *t_cmb,
                                   const double *y_he, const double *z_Hz, const double *Hz, const int *nz,
                                   double *za_in, const int *len, const int *label)
{
    const char *dp = getenv("COSMORECJL_DUMP");
    if (!dp || !dp[0]) return;
    FILE *f = fopen(dp, "wb");
    if (!f) { fprintf(stderr, "cosmorecjl shim: cannot open COSMORECJL_DUMP=%s; aborting\n", dp); abort(); }
    int32_t ver = 1;
    fwrite("CRJLDUMP1", 1, 9, f);
    fwrite(&ver, sizeof(int32_t), 1, f);
    int32_t rm = (int32_t)*runmode;    fwrite(&rm, sizeof(int32_t), 1, f);
    fwrite(runpars, sizeof(double), 21, f);
    fwrite(omega_c, sizeof(double), 1, f); fwrite(omega_b, sizeof(double), 1, f); fwrite(omega_k, sizeof(double), 1, f);
    fwrite(num_nu, sizeof(double), 1, f);  fwrite(h0, sizeof(double), 1, f);      fwrite(t_cmb, sizeof(double), 1, f);
    fwrite(y_he, sizeof(double), 1, f);
    int32_t n = (int32_t)*nz;          fwrite(&n, sizeof(int32_t), 1, f);
    fwrite(z_Hz, sizeof(double), (size_t)n, f); fwrite(Hz, sizeof(double), (size_t)n, f);
    int32_t l = (int32_t)*len;         fwrite(&l, sizeof(int32_t), 1, f);
    fwrite(za_in, sizeof(double), (size_t)l, f);
    int32_t lb = (int32_t)*label;      fwrite(&lb, sizeof(int32_t), 1, f);
    if (fflush(f) != 0 || fclose(f) != 0) { fprintf(stderr, "cosmorecjl shim: dump write failed; aborting\n"); abort(); }
}

void cosmorec_calc_h_cpp_(const int *runmode, const double *runpars, const double *omega_c, const double *omega_b, const double *omega_k,
                          const double *num_nu, const double *h0, const double *t_cmb, const double *y_he, const double *z_Hz, const double *Hz,
                          const int *nz, double *za_in, double *xe_out, double *tb_out, const int *len, const int *label)
{
    cosmorecjl_dump_inputs(runmode, runpars, omega_c, omega_b, omega_k, num_nu, h0, t_cmb, y_he, z_Hz, Hz, nz, za_in, len, label);
    const char *nlib = getenv("COSMORECJL_NATIVE_LIB");
    if (nlib && nlib[0]) {              /* input-contract mode: forward to the ORIGINAL native CosmoRec via dlopen */
        static void *h = NULL;
        static cosmorecjl_fn realfn = NULL;
        if (!h) {
            h = dlopen(nlib, RTLD_NOW | RTLD_LOCAL);
            if (!h) { fprintf(stderr, "cosmorecjl shim: dlopen(%s) failed: %s; aborting\n", nlib, dlerror()); abort(); }
            realfn = (cosmorecjl_fn)dlsym(h, "cosmorec_calc_h_cpp_");
            if (!realfn) { fprintf(stderr, "cosmorecjl shim: dlsym(cosmorec_calc_h_cpp_) in %s failed; aborting\n", nlib); abort(); }
        }
        g_ncalls++;
        g_last_status = realfn(runmode, runpars, omega_c, omega_b, omega_k, num_nu, h0, t_cmb, y_he, z_Hz, Hz, nz, za_in, xe_out, tb_out, len, label);
        if (g_last_status != 0) { fprintf(stderr, "cosmorecjl shim: native CosmoRec returned status %d; aborting\n", g_last_status); abort(); }
        return;
    }
    if (g_fn == NULL) { fprintf(stderr, "cosmorecjl shim: cosmorec_calc_h_cpp_ called but no Julia entry is registered; aborting (no fallback)\n"); abort(); }
    g_ncalls++;
    g_last_status = g_fn(runmode, runpars, omega_c, omega_b, omega_k, num_nu, h0, t_cmb, y_he, z_Hz, Hz, nz, za_in, xe_out, tb_out, len, label);
    if (g_last_status != 0) { fprintf(stderr, "cosmorecjl shim: Julia CosmoRec.jl returned status %d; aborting (no fallback)\n", g_last_status); abort(); }
}

/* The other native entry points are not used by CAMB 2.0.4 (cosmorec.f90 calls only the _h_ variant); trap them. */
void cosmorec_calc_cpp_(void) { fprintf(stderr, "cosmorecjl shim: cosmorec_calc_cpp_ is not provided; aborting\n"); abort(); }
void cosmorec_calc_hptr_cpp_(void) { fprintf(stderr, "cosmorecjl shim: cosmorec_calc_hptr_cpp_ is not provided; aborting\n"); abort(); }
