/* CosmoRec.jl bridge shim (chunk36 scratch). Replaces libCosmoRec.a in a SCRATCH CAMB build: exports the exact CAMB-facing symbol
 * cosmorec_calc_h_cpp_ (tools/CosmoRec/CosmoRec.h:153-158; called from CAMB fortran/cosmorec.f90:211-213) and forwards every argument
 * pointer, unchanged, to a Julia @cfunction registered at runtime. No native CosmoRec code and no GSL are linked. There is NO fallback:
 * an unregistered call or a nonzero Julia status aborts the process (no native call, no zero/garbage output reaches CAMB). */
#include <stdio.h>
#include <stdlib.h>

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

void cosmorec_calc_h_cpp_(const int *runmode, const double *runpars, const double *omega_c, const double *omega_b, const double *omega_k,
                          const double *num_nu, const double *h0, const double *t_cmb, const double *y_he, const double *z_Hz, const double *Hz,
                          const int *nz, double *za_in, double *xe_out, double *tb_out, const int *len, const int *label)
{
    if (g_fn == NULL) { fprintf(stderr, "cosmorecjl shim: cosmorec_calc_h_cpp_ called but no Julia entry is registered; aborting (no fallback)\n"); abort(); }
    g_ncalls++;
    g_last_status = g_fn(runmode, runpars, omega_c, omega_b, omega_k, num_nu, h0, t_cmb, y_he, z_Hz, Hz, nz, za_in, xe_out, tb_out, len, label);
    if (g_last_status != 0) { fprintf(stderr, "cosmorecjl shim: Julia CosmoRec.jl returned status %d; aborting (no fallback)\n", g_last_status); abort(); }
}

/* The other native entry points are not used by CAMB 2.0.4 (cosmorec.f90 calls only the _h_ variant); trap them. */
void cosmorec_calc_cpp_(void) { fprintf(stderr, "cosmorecjl shim: cosmorec_calc_cpp_ is not provided; aborting\n"); abort(); }
void cosmorec_calc_hptr_cpp_(void) { fprintf(stderr, "cosmorecjl shim: cosmorec_calc_hptr_cpp_ is not provided; aborting\n"); abort(); }
