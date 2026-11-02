//
//  iso_core.c
//  enviPat
//
//  Copyright (c) 2014 cg. All rights reserved.
//  Extracted and modified from main.c in 2026 for the Python port
//  (pyEnvipat); distributed under the GNU GPL v2 (GPL-2.0-only), same
//  as the original R package. See LICENSE, AUTHORS and NOTICE.
//
//  Pure-C core algorithms for the 5 algorithm entry points. This file has
//  no R dependency. It is a structural extraction from the original
//  main.c: the floating point operations below run in exactly the same
//  order as the original code they were lifted from -- this is a pure
//  refactor, not a rewrite.
//
//  See iso_core.h for the error-code convention used uniformly across all
//  entry points here.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "preferences.h"
#include "isotope.h"
#include "parse.h"
#include "element.h"
#include "combination.h"
#include "peak.h"
#include "profile.h"
#include "iso_core.h"

/* ---------------------------------------------------------------------- */
/* shared utility (moved verbatim from main.c)                            */
/* ---------------------------------------------------------------------- */

int alloc_peaks(int p_l, size_t iso_amount, double **m_, double **a_, int **cc_){

#if USE_REALLOC == 1
    if (ALLOC_START_PL > p_l) {
        size_t alloc_ = sizeof(double[p_l]);
        *m_ = (double*)malloc(alloc_);
        if(*m_ == NULL){
            return 1;
        }
        *a_ = (double*)malloc(alloc_);
        if(*a_ == NULL){
            free(*m_);
            return 1;
        }
        alloc_ = sizeof(int[p_l * iso_amount]);
        *cc_ = (int*)malloc(alloc_);
        if(*cc_ == NULL){
            free(*m_);
            free(*a_);
            return 1;
        }
    }else{
        size_t alloc_ = ALLOC_START_PL * sizeof(double);
        *m_ = (double*)malloc(alloc_);
        if(*m_ == NULL){
            return 1;
        }
        alloc_ = ALLOC_START_PL * sizeof(double);
        *a_ = (double*)malloc(alloc_);
        if(*a_ == NULL){
            free(*m_);
            return 1;
        }
        alloc_ = ALLOC_START_PL * iso_amount * sizeof(int);
        *cc_ = (int*)malloc(alloc_);
        if(*cc_ == NULL){
            free(*m_);
            free(*a_);
            return 1;
        }
    }

#else
    size_t alloc_ = sizeof(double[p_l]);
    *m_ = (double*)malloc(alloc_);
    if(*m_ == NULL){
        return 1;
    }
    *a_ = (double*)malloc(alloc_);
    if(*a_ == NULL){
        free(*m_);
        return 1;
    }
    alloc_ = sizeof(int[p_l * iso_amount]);
    *cc_ = (int*)malloc(alloc_);
    if(*cc_ == NULL){
        free(*m_);
        free(*a_);
        return 1;
    }
#endif
    return 0;
}

/* ---------------------------------------------------------------------- */
/* iso_ppm_Call core                                                      */
/* ---------------------------------------------------------------------- */

int iso_core_ppm_trace(double start,
                        double end,
                        double ppm,
                        double **trace_out,
                        size_t *trace_len,
                        char errmsg[ISO_CORE_ERRMSG_SIZE])
{
    unsigned int storage_limit = 2e7;
    double *tr = (double*)malloc((size_t)storage_limit * sizeof(double));
    if (tr == NULL) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot allocate memory for ppm trace");
        return ISO_CORE_ERR_ALLOC;
    }

    *tr = start;

    int n = 0;
    while (*(tr + n) < end) {
        double tmp = *(tr + n) + *(tr + n) * ppm / 1e6;
        n++;
        *(tr + n) = tmp;
        if (n > (int)storage_limit - 1) {
            snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "too many mass points for ppm trace");
            free(tr);
            return ISO_CORE_ERR_TOO_MANY_TRACE_PTS;
        }
    }

    /* Trim the over-allocated buffer down to the values actually written
       (mirrors the original code's resize to n + 1 elements). */
    double *trimmed = (double*)realloc(tr, (size_t)(n + 1) * sizeof(double));
    if (trimmed != NULL) {
        tr = trimmed;
    }

    *trace_out = tr;
    *trace_len = (size_t)(n + 1);
    return ISO_CORE_OK;
}

/* ---------------------------------------------------------------------- */
/* iso_centroid_Call core                                                 */
/* type: 0 = centroid, 1 = intensoid (local maxima), 2 = valley           */
/* ---------------------------------------------------------------------- */

int iso_core_centroid(double *profile_mass,
                       double *profile_abundance,
                       int n,
                       int type,
                       double **out_mass,
                       double **out_abundance,
                       int *out_len,
                       char errmsg[ISO_CORE_ERRMSG_SIZE])
{
    if (n <= 0) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "invalid profile length");
        return ISO_CORE_ERR_INVALID_LENGTH;
    }

    double *p_m = profile_mass;
    double *p_a = profile_abundance;
    int t = type;

    size_t cent_size_alloc = (size_t)n * sizeof(double);

    double* c_m = (double*)malloc(cent_size_alloc);
    double* c_a = (double*)malloc(cent_size_alloc);
    if (c_m == NULL || c_a == NULL) {
        free(c_m);
        free(c_a);
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot allocate memory for centroid output");
        return ISO_CORE_ERR_ALLOC;
    }

    double centroid = 0.0;
    double sum_sticks = 0.0;

    double upper_sum_a = 0.0;
    double lower_sum_a = 0.0;
    double max_centroid = 0.0;
    double max_intensoid = 0.0;
    int centroid_count = 0;

    int j = 0;
    int  step = 1;

    for (ptrdiff_t  i = step; i < n - step; i++) {
        if (t == 1)
        {
            if (    (
                     *(p_a + i) > *(p_a + i + step) && *(p_a + i) > *(p_a + i - step)
                     )
                &&
                (
                 (*(p_m + i - step) < *(p_m + i)) && (*(p_m + i) < *(p_m + i + step))
                 )

                )
            {
                *(c_a + j) = *(p_a + i);
                *(c_m + j) = *(p_m + i);
                j++;

                if (*(c_a + j) > max_intensoid) {
                    max_intensoid = *(c_a + j);
                }
            }
        }
        else{

            centroid += *(p_a + i) * *(p_m + i);
            sum_sticks += *(p_a + i);

            double diff_mass = *(p_m + i) - *(p_m + i - 1);
            upper_sum_a += diff_mass * *(p_a + i);
            lower_sum_a += diff_mass * *(p_a + i - 1);

            if (( (*(p_a + i) <= *(p_a + i + step) && *(p_a + i) < *(p_a + i - step)) || i == n - step - 1)) {

                double centroid_temp = (upper_sum_a + lower_sum_a)/2.0;
                if (t == 2) {
                    *(c_a + j) = *(p_a + i);
                    *(c_m + j) = *(p_m + i);
                    j++;
                }

                if (max_centroid < centroid_temp) {
                    max_centroid = centroid_temp;
                }

                if (t == 0) {
                    if (sum_sticks > 0.0 && centroid_temp > 0.0 ) {
                        *(c_a + j) = centroid_temp;
                        *(c_m + j) = centroid/sum_sticks;
                        j++;
                        centroid_count++;
                    }
                }
                centroid = 0.0;
                sum_sticks = 0.0;

                upper_sum_a = 0.0;
                lower_sum_a = 0.0;
            }
        }
    }

    if(t == 0){
        for (ptrdiff_t l = 0; l < centroid_count; l++) {
            double temp = *(c_a +  l) * (100/max_centroid);
            *(c_a + l) = temp;
        }
    }

    *out_mass = c_m;
    *out_abundance = c_a;
    *out_len = j;
    return ISO_CORE_OK;
}

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
                                 char errmsg[ISO_CORE_ERRMSG_SIZE])
{
    if (trace_len >= (size_t)INT_LIMIT || trace_len < 1) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "invalid amount of trace values");
        return ISO_CORE_ERR_INVALID_LENGTH;
    }
    if (mass_len >= (size_t)INT_LIMIT || mass_len < 1) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "invalid peak limit");
        return ISO_CORE_ERR_INVALID_LENGTH;
    }

    double* p_m = (double*)malloc(sizeof(double[trace_len]));
    double* p_a = (double*)malloc(sizeof(double[trace_len]));
    if (p_m == NULL || p_a == NULL) {
        free(p_m);
        free(p_a);
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot allocate memory for profile output");
        return ISO_CORE_ERR_ALLOC;
    }

    unsigned int p_n = 0;

    int msg = calc_profile_with_trace(mass_len, mass, abundance, trace_len, trace, p_m, p_a, &p_n, resolution, profile_type, threshold, filter);
    if (msg) {
        free(p_m);
        free(p_a);
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "could not calculate profile");
        return ISO_CORE_ERR_COMBINATION;
    }

    *out_mass = p_m;
    *out_abundance = p_a;
    *out_len = p_n;
    return ISO_CORE_OK;
}

/* ---------------------------------------------------------------------- */
/* iso_pattern (algo=1) / iso_pattern_4 (algo=2) core                     */
/* ---------------------------------------------------------------------- */

void iso_pattern_result_init(IsoPatternResult *result)
{
    if (result == NULL) return;
    result->elements = NULL;
    result->element_amount = 0;
    result->iso_amount = 0;
    result->mass = NULL;
    result->abundance = NULL;
    result->compound = NULL;
    result->peak_amount = 0;
    result->names = NULL;
}

void iso_pattern_result_free(IsoPatternResult *result)
{
    if (result == NULL) return;
    free(result->elements);
    free(result->mass);
    free(result->abundance);
    free(result->compound);
    free(result->names);
    iso_pattern_result_init(result);
}

/* Builds the (iso_amount + 2) * MAX_NAME_SIZE packed name buffer
   ("mass", "abundance", then one entry per isotope in element order) that
   both algo1 and algo4 need identically. */
static char* build_pattern_names(Element *elements, unsigned short element_amount, unsigned short iso_amount)
{
    char *l_n = (char*)calloc((size_t)(iso_amount + 2) * MAX_NAME_SIZE, sizeof(char));
    if (l_n == NULL) {
        return NULL;
    }
    strncpy(l_n, "mass\0", 5 * sizeof(char));
    strncpy(l_n + MAX_NAME_SIZE, "abundance\0", 10 * sizeof(char));
    size_t index = 0;
    for (ptrdiff_t b = 0; b < element_amount; b++) {
        for (ptrdiff_t bb = 0; bb < (elements+b)->iso_amount; bb++) {
            strncpy(l_n + (index + 2) * MAX_NAME_SIZE, ((elements + b)->isotopes + bb)->isotope, MAX_NAME_SIZE*sizeof(char));
            index++;
        }
    }
    return l_n;
}

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
                            char errmsg[ISO_CORE_ERRMSG_SIZE])
{
    iso_pattern_result_init(result);

    unsigned short element_amount = 0;
    unsigned short iso_amount = 0;

    Element* elements = (Element*)calloc(MAX_ELEMENTS, sizeof(Element));
    if (elements == NULL) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot allocate memory for elements pointer");
        return ISO_CORE_ERR_ALLOC;
    }

    if ( parse_sum_formula_vector_R(
                                    elements,
                                    sum,
                                    &element_amount,
                                    &iso_amount,
                                    iso_amount_global,
                                    iso_list_elem,
                                    iso_list_iso,
                                    iso_list_mass,
                                    iso_list_abu
                                    )
        ) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot parse sum formula with the given isolist");
        free(elements);
        return ISO_CORE_ERR_PARSE;
    }
    qsort(elements, element_amount, sizeof(Element), elements_sort_by_isoamount_inc);

    if (peak_limit >= INT_LIMIT || peak_limit < 1) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "invalid peak limit");
        free(elements);
        return ISO_CORE_ERR_INVALID_PEAK_LIMIT;
    }

    if ( rel_to_mono > 4 || rel_to_mono < 0) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "wrong value for rtm");
        free(elements);
        return ISO_CORE_ERR_INVALID_RTM;
    }

    if (iso_amount >= MAX_ISO_SIZE
        || iso_amount <= 0
        || element_amount <= 0
        || element_amount >= MAX_ELEMENTS
        ) {
        free(elements);
        /* Matches the original code: this specific check failed silently,
           without a message. errmsg is left empty so callers can tell
           this case apart (enviPat_py.c falls back to a generic message). */
        errmsg[0] = '\0';
        return ISO_CORE_ERR_INVALID_ISO_AMOUNT;
    }

    double* m_ = NULL;
    double* a_ = NULL;
    int* cc_ = NULL;
    int msg_peaks = alloc_peaks(peak_limit, iso_amount, &m_, &a_, &cc_);
    if (msg_peaks) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "pointer allocation, error code: %d", msg_peaks);
        free(elements);
        return ISO_CORE_ERR_ALLOC;
    }

    char* l_n = build_pattern_names(elements, element_amount, iso_amount);
    if(l_n == NULL){
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot allocate memory for list name pointer");
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        return ISO_CORE_ERR_ALLOC;
    }

    CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
    if(monoisotopic == NULL){
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot allocate memory for list name pointer");
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        free(l_n);
        return ISO_CORE_ERR_ALLOC;
    }

    calc_monoisotopic(elements, element_amount, monoisotopic);
    long double mono_abundance = monoisotopic->abundance;
    free(monoisotopic);

    size_t peak_amount = 0;
    int msg = calc_pattern_algo_1(elements, &peak_amount, threshold, iso_amount,element_amount, mono_abundance,peak_limit,rel_to_mono, &m_, &a_, &cc_);
    if(msg || peak_amount < 1){
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "unable to calculate isotope pattern, exit code %d",msg);
        free(cc_);
        free(l_n);
        free(m_);
        free(a_);
        free(elements);
        return ISO_CORE_ERR_COMBINATION;
    }

    result->elements = elements;
    result->element_amount = element_amount;
    result->iso_amount = iso_amount;
    result->mass = m_;
    result->abundance = a_;
    result->compound = cc_;
    result->peak_amount = peak_amount;
    result->names = l_n;

    return ISO_CORE_OK;
}

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
                            char errmsg[ISO_CORE_ERRMSG_SIZE])
{
    iso_pattern_result_init(result);

    unsigned short element_amount = 0;
    unsigned short iso_amount = 0;

    Element* elements = (Element*)calloc(MAX_ELEMENTS, sizeof(Element));
    if (elements == NULL) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot allocate memory for elements pointer");
        return ISO_CORE_ERR_ALLOC;
    }

    if ( parse_sum_formula_vector_R(
                                    elements,
                                    sum,
                                    &element_amount,
                                    &iso_amount,
                                    iso_amount_global,
                                    iso_list_elem,
                                    iso_list_iso,
                                    iso_list_mass,
                                    iso_list_abu
                                    )
        ) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot parse sum formula with the given isolist");
        free(elements);
        return ISO_CORE_ERR_PARSE;
    }
    qsort(elements, element_amount, sizeof(Element), elements_sort_by_isoamount_inc);

    if (peak_limit >= INT_LIMIT || peak_limit < 1) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "invalid peak limit");
        free(elements);
        return ISO_CORE_ERR_INVALID_PEAK_LIMIT;
    }

    if ( rel_to_mono > 4 || rel_to_mono < 0) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "wrong value for rtm");
        free(elements);
        return ISO_CORE_ERR_INVALID_RTM;
    }

    if (iso_amount >= MAX_ISO_SIZE
        || iso_amount <= 0
        || element_amount <= 0
        || element_amount >= MAX_ELEMENTS
        ) {
        free(elements);
        /* Matches the original code: this specific check failed silently,
           without a message. errmsg is left empty so callers can tell
           this case apart (enviPat_py.c falls back to a generic message). */
        errmsg[0] = '\0';
        return ISO_CORE_ERR_INVALID_ISO_AMOUNT;
    }

    double* m_ = NULL;
    double* a_ = NULL;
    int* cc_ = NULL;
    int msg_peaks = alloc_peaks(peak_limit, iso_amount, &m_, &a_, &cc_);
    if (msg_peaks) {
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "pointer allocation, error code: %d", msg_peaks);
        free(elements);
        return ISO_CORE_ERR_ALLOC;
    }

    char* l_n = build_pattern_names(elements, element_amount, iso_amount);
    if(l_n == NULL){
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot allocate memory column names");
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        return ISO_CORE_ERR_ALLOC;
    }

    long double mono_abundance = -1.0;
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
        mono_abundance *= pow((elements + i)->isotopes[0].abundance, (elements + i)->amount);
    }

    size_t peak_amount = 0;
    long double max_a = 0;
    int msg = calc_pattern_algo_4(elements, &peak_amount, threshold, iso_amount, element_amount, mono_abundance, peak_limit, l_n, rel_to_mono, &m_, &a_, &cc_, &max_a);
    if(msg || peak_amount == 0){
        snprintf(errmsg, ISO_CORE_ERRMSG_SIZE, "cannot combine combinations, exit code: %d", msg);
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        free(l_n);
        return ISO_CORE_ERR_COMBINATION;
    }

    result->elements = elements;
    result->element_amount = element_amount;
    result->iso_amount = iso_amount;
    result->mass = m_;
    result->abundance = a_;
    result->compound = cc_;
    result->peak_amount = peak_amount;
    result->names = l_n;

    return ISO_CORE_OK;
}
