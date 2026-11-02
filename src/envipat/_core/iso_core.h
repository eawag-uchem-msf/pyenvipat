//
//  iso_core.h
//  enviPat
//
//  Copyright (c) 2014 cg. All rights reserved.
//  Extracted and modified from main.c in 2026 for the Python port
//  (pyEnvipat); distributed under the GNU GPL v2 (GPL-2.0-only), same
//  as the original R package. See LICENSE, AUTHORS and NOTICE.
//
//  Pure-C core for the 5 algorithm entry points (pattern algo 1/4, ppm
//  trace, profile with trace, centroid). No R dependency in this file or
//  in iso_core.c -- this is the boundary the CPython extension
//  (enviPat_py.c) calls into.
//
//  Error handling convention: every iso_core_* entry point returns an
//  ISO_CORE_OK / ISO_CORE_ERR_* status code and, on failure, fills the
//  caller-supplied errmsg buffer (ISO_CORE_ERRMSG_SIZE bytes) with a
//  human-readable message (no trailing newline); reporting it is up to
//  the caller.
//

#ifndef ISO_CORE_H
#define ISO_CORE_H

#include <stddef.h>

#include "preferences.h"
#include "isotope.h"
#include "element.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------- */
/* status codes                                                           */
/* ---------------------------------------------------------------------- */

#define ISO_CORE_OK                       0
#define ISO_CORE_ERR_ALLOC                1
#define ISO_CORE_ERR_INVALID_PEAK_LIMIT   2
#define ISO_CORE_ERR_INVALID_RTM          3
#define ISO_CORE_ERR_PARSE                4
#define ISO_CORE_ERR_INVALID_ISO_AMOUNT   5
#define ISO_CORE_ERR_COMBINATION          6
#define ISO_CORE_ERR_TOO_MANY_TRACE_PTS   7
#define ISO_CORE_ERR_INVALID_LENGTH       8

#define ISO_CORE_ERRMSG_SIZE 256

/* ---------------------------------------------------------------------- */
/* shared utility (moved verbatim from main.c, already pure C there)      */
/* ---------------------------------------------------------------------- */

int alloc_peaks(int p_l, size_t iso_amount, double **m_, double **a_, int **cc_);

/* ---------------------------------------------------------------------- */
/* iso_ppm_Call core: generates a ppm-spaced mass trace                   */
/* ---------------------------------------------------------------------- */

int iso_core_ppm_trace(double start,
                        double end,
                        double ppm,
                        double **trace_out,
                        size_t *trace_len,
                        char errmsg[ISO_CORE_ERRMSG_SIZE]);

/* ---------------------------------------------------------------------- */
/* iso_centroid_Call core: centroid / intensoid / valley detection        */
/* type: 0 = centroid, 1 = intensoid (local maxima), 2 = valley           */
/* ---------------------------------------------------------------------- */

int iso_core_centroid(double *profile_mass,
                       double *profile_abundance,
                       int n,
                       int type,
                       double **out_mass,
                       double **out_abundance,
                       int *out_len,
                       char errmsg[ISO_CORE_ERRMSG_SIZE]);

/* ---------------------------------------------------------------------- */
/* iso_profile_with_trace_Call core                                       */
/* ---------------------------------------------------------------------- */

int iso_core_profile_with_trace(int profile_type,
                                 double *mass,
                                 size_t mass_len,
                                 double *abundance,
                                 double *trace,
                                 size_t trace_len,
                                 int resolution,
                                 double threshold,
                                 int filter,
                                 double **out_mass,
                                 double **out_abundance,
                                 unsigned int *out_len,
                                 char errmsg[ISO_CORE_ERRMSG_SIZE]);

/* ---------------------------------------------------------------------- */
/* iso_pattern (algo=1) / iso_pattern_4 (algo=2) core                     */
/* ---------------------------------------------------------------------- */

typedef struct {
    Element *elements;       /* length element_amount */
    unsigned short element_amount;
    unsigned short iso_amount;
    double *mass;             /* length peak_amount */
    double *abundance;        /* length peak_amount */
    int    *compound;         /* length peak_amount * iso_amount */
    size_t  peak_amount;
    char   *names;             /* (iso_amount + 2) * MAX_NAME_SIZE packed
                                   fixed-width strings: "mass", "abundance",
                                   then one per isotope */
} IsoPatternResult;

void iso_pattern_result_init(IsoPatternResult *result);
void iso_pattern_result_free(IsoPatternResult *result);

int iso_core_pattern_algo1(char *sum,
                            int peak_limit,
                            double threshold,
                            int rel_to_mono,
                            char *iso_list_elem,
                            char *iso_list_iso,
                            double *iso_list_mass,
                            double *iso_list_abu,
                            unsigned short iso_amount_global,
                            IsoPatternResult *result,
                            char errmsg[ISO_CORE_ERRMSG_SIZE]);

int iso_core_pattern_algo4(char *sum,
                            int peak_limit,
                            double threshold,
                            int rel_to_mono,
                            char *iso_list_elem,
                            char *iso_list_iso,
                            double *iso_list_mass,
                            double *iso_list_abu,
                            unsigned short iso_amount_global,
                            IsoPatternResult *result,
                            char errmsg[ISO_CORE_ERRMSG_SIZE]);

#ifdef __cplusplus
}
#endif

#endif /* ISO_CORE_H */
