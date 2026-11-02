/*
 * enviPat_py.c
 *
 * CPython + NumPy C-API extension exposing the 5 algorithm entry points of
 * the pure-C core (iso_core.c/h, extracted from the R package's src/main.c
 * in issue #4) directly to Python, with NumPy arrays in/out.
 *
 * This is the Python-side analogue of the R package's SEXP wrappers in
 * src/main.c -- but it deliberately does NOT replicate their R-specific
 * result-shaping logic (e.g. the `rica`/return_iso_calc_amount branch that
 * returns a completely different SEXP type, or algo4's zero-count element
 * name concatenation within that branch). Every pattern_algo* call here
 * always returns the same dict shape, with the raw per-element
 * all_iso_calc_amount values included; a future, higher-level Python port
 * of the R wrapper functions themselves (isopattern(), etc.) can build
 * whatever R-equivalent views it needs on top of this.
 *
 * Error handling: iso_core_* functions return an ISO_CORE_OK/ISO_CORE_ERR_*
 * status code and fill a caller-supplied errmsg buffer on failure (see
 * iso_core.h). Every failure here is surfaced as a raised
 * envipat._core.enviPat_py.EnviPatCoreError (falling back to a generic
 * message when errmsg is left empty -- iso_core_pattern_algo1/4 do this
 * deliberately for one specific, silent-in-R error case; see
 * ISO_CORE_ERR_INVALID_ISO_AMOUNT in iso_core.c).
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <string.h>

#include "iso_core.h"
#include "preferences.h" /* MAX_NAME_SIZE */

static PyObject *EnviPatCoreError;

/* ------------------------------------------------------------------- */
/* small helpers                                                       */
/* ------------------------------------------------------------------- */

/* Copies `n` doubles out of a malloc'd buffer into a new, NumPy-owned
   1-D float64 array. Does not free `buf` -- the caller does that, mirroring
   the R SEXP wrappers' own memcpy-then-free pattern in main.c. */
static PyObject *
double_array_from_buffer(const double *buf, npy_intp n)
{
    PyObject *arr = PyArray_SimpleNew(1, &n, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    if (n > 0) {
        memcpy(PyArray_DATA((PyArrayObject *)arr), buf, (size_t)n * sizeof(double));
    }
    return arr;
}

/* Copies one MAX_NAME_SIZE-wide packed fixed-width name out of `names` at
   slot `index` into a new Python str (NUL-terminated defensively, in case
   a name fills the whole MAX_NAME_SIZE buffer with no trailing NUL -- see
   how build_pattern_names()/Element.name are filled in iso_core.c). */
static PyObject *
name_from_slot(const char *names, size_t index)
{
    char buf[MAX_NAME_SIZE + 1];
    memcpy(buf, names + index * MAX_NAME_SIZE, MAX_NAME_SIZE);
    buf[MAX_NAME_SIZE] = '\0';
    return PyUnicode_FromString(buf);
}

/* Raises EnviPatCoreError with `errmsg` when non-empty, or `fallback`
   otherwise (iso_core_pattern_algo1/4 leave errmsg empty on one specific
   error path -- ISO_CORE_ERR_INVALID_ISO_AMOUNT -- that the original SEXP
   wrapper treats as a silent R_NilValue return; every other error path
   always sets errmsg). */
static void
set_core_error(const char errmsg[ISO_CORE_ERRMSG_SIZE], const char *fallback)
{
    PyErr_SetString(EnviPatCoreError, errmsg[0] != '\0' ? errmsg : fallback);
}

/* ------------------------------------------------------------------- */
/* iso_ppm_Call core: envipat._core.enviPat_py.ppm_trace                */
/* ------------------------------------------------------------------- */

static PyObject *
envipat_ppm_trace(PyObject *self, PyObject *args, PyObject *kwargs)
{
    (void)self;
    static char *kwlist[] = {"start", "end", "ppm", NULL};
    double start, end, ppm;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ddd", kwlist, &start, &end, &ppm)) {
        return NULL;
    }

    double *trace = NULL;
    size_t trace_len = 0;
    char errmsg[ISO_CORE_ERRMSG_SIZE] = {0};

    int status = iso_core_ppm_trace(start, end, ppm, &trace, &trace_len, errmsg);
    if (status != ISO_CORE_OK) {
        set_core_error(errmsg, "iso_core_ppm_trace failed");
        return NULL;
    }

    PyObject *result = double_array_from_buffer(trace, (npy_intp)trace_len);
    free(trace);
    return result;
}

/* ------------------------------------------------------------------- */
/* iso_centroid_Call core: envipat._core.enviPat_py.centroid            */
/* type: 0 = centroid, 1 = intensoid (local maxima), 2 = valley         */
/* ------------------------------------------------------------------- */

static PyObject *
envipat_centroid(PyObject *self, PyObject *args, PyObject *kwargs)
{
    (void)self;
    static char *kwlist[] = {"profile_mass", "profile_abundance", "type", NULL};
    PyObject *mass_obj, *abundance_obj;
    int type;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOi", kwlist, &mass_obj, &abundance_obj, &type)) {
        return NULL;
    }

    PyArrayObject *mass_arr = (PyArrayObject *)PyArray_FROM_OTF(
        mass_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (mass_arr == NULL) {
        return NULL;
    }
    PyArrayObject *abundance_arr = (PyArrayObject *)PyArray_FROM_OTF(
        abundance_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (abundance_arr == NULL) {
        Py_DECREF(mass_arr);
        return NULL;
    }

    PyObject *ret = NULL;
    if (PyArray_NDIM(mass_arr) != 1 || PyArray_NDIM(abundance_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "profile_mass and profile_abundance must be 1-D arrays");
        goto cleanup;
    }
    if (PyArray_DIM(mass_arr, 0) != PyArray_DIM(abundance_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "profile_mass and profile_abundance must have the same length");
        goto cleanup;
    }

    {
        int n = (int)PyArray_DIM(mass_arr, 0);
        double *out_mass = NULL, *out_abundance = NULL;
        int out_len = 0;
        char errmsg[ISO_CORE_ERRMSG_SIZE] = {0};

        int status = iso_core_centroid(
            (double *)PyArray_DATA(mass_arr), (double *)PyArray_DATA(abundance_arr),
            n, type, &out_mass, &out_abundance, &out_len, errmsg);
        if (status != ISO_CORE_OK) {
            set_core_error(errmsg, "iso_core_centroid failed");
            goto cleanup;
        }

        PyObject *mass_out = double_array_from_buffer(out_mass, (npy_intp)out_len);
        PyObject *abundance_out = mass_out ? double_array_from_buffer(out_abundance, (npy_intp)out_len) : NULL;
        free(out_mass);
        free(out_abundance);
        if (mass_out == NULL || abundance_out == NULL) {
            Py_XDECREF(mass_out);
            Py_XDECREF(abundance_out);
            goto cleanup;
        }
        ret = Py_BuildValue("(NN)", mass_out, abundance_out);
    }

cleanup:
    Py_DECREF(mass_arr);
    Py_DECREF(abundance_arr);
    return ret;
}

/* ------------------------------------------------------------------- */
/* iso_profile_with_trace_Call core:                                    */
/* envipat._core.enviPat_py.profile_with_trace                          */
/* ------------------------------------------------------------------- */

static PyObject *
envipat_profile_with_trace(PyObject *self, PyObject *args, PyObject *kwargs)
{
    (void)self;
    static char *kwlist[] = {"profile_type", "mass", "abundance", "trace",
                              "resolution", "threshold", "filter", NULL};
    int profile_type;
    PyObject *mass_obj, *abundance_obj, *trace_obj;
    int resolution;
    double threshold;
    int filter;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iOOOidi", kwlist, &profile_type, &mass_obj,
                                      &abundance_obj, &trace_obj, &resolution, &threshold, &filter)) {
        return NULL;
    }

    PyArrayObject *mass_arr = (PyArrayObject *)PyArray_FROM_OTF(
        mass_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (mass_arr == NULL) {
        return NULL;
    }
    PyArrayObject *abundance_arr = (PyArrayObject *)PyArray_FROM_OTF(
        abundance_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (abundance_arr == NULL) {
        Py_DECREF(mass_arr);
        return NULL;
    }
    PyArrayObject *trace_arr = (PyArrayObject *)PyArray_FROM_OTF(
        trace_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (trace_arr == NULL) {
        Py_DECREF(mass_arr);
        Py_DECREF(abundance_arr);
        return NULL;
    }

    PyObject *ret = NULL;
    if (PyArray_NDIM(mass_arr) != 1 || PyArray_NDIM(abundance_arr) != 1 || PyArray_NDIM(trace_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "mass, abundance and trace must be 1-D arrays");
        goto cleanup;
    }
    if (PyArray_DIM(mass_arr, 0) != PyArray_DIM(abundance_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "mass and abundance must have the same length");
        goto cleanup;
    }

    {
        size_t mass_len = (size_t)PyArray_DIM(mass_arr, 0);
        size_t trace_len = (size_t)PyArray_DIM(trace_arr, 0);
        double *out_mass = NULL, *out_abundance = NULL;
        unsigned int out_len = 0;
        char errmsg[ISO_CORE_ERRMSG_SIZE] = {0};

        int status = iso_core_profile_with_trace(
            profile_type,
            (double *)PyArray_DATA(mass_arr), mass_len, (double *)PyArray_DATA(abundance_arr),
            (double *)PyArray_DATA(trace_arr), trace_len,
            resolution, threshold, filter,
            &out_mass, &out_abundance, &out_len, errmsg);
        if (status != ISO_CORE_OK) {
            set_core_error(errmsg, "iso_core_profile_with_trace failed");
            goto cleanup;
        }

        PyObject *mass_out = double_array_from_buffer(out_mass, (npy_intp)out_len);
        PyObject *abundance_out = mass_out ? double_array_from_buffer(out_abundance, (npy_intp)out_len) : NULL;
        free(out_mass);
        free(out_abundance);
        if (mass_out == NULL || abundance_out == NULL) {
            Py_XDECREF(mass_out);
            Py_XDECREF(abundance_out);
            goto cleanup;
        }
        ret = Py_BuildValue("(NN)", mass_out, abundance_out);
    }

cleanup:
    Py_DECREF(mass_arr);
    Py_DECREF(abundance_arr);
    Py_DECREF(trace_arr);
    return ret;
}

/* ------------------------------------------------------------------- */
/* iso_pattern (algo=1) / iso_pattern_4 (algo=2) core:                  */
/* envipat._core.enviPat_py.pattern_algo1 / pattern_algo2                */
/* ------------------------------------------------------------------- */

typedef int (*PatternAlgoFn)(char *sum, int peak_limit, double threshold, int rel_to_mono,
                              char *iso_list_elem, char *iso_list_iso,
                              double *iso_list_mass, double *iso_list_abu,
                              unsigned short iso_amount_global,
                              IsoPatternResult *result, char errmsg[ISO_CORE_ERRMSG_SIZE]);

/* Packs an IsoPatternResult into a Python dict:
 *   mass                        : float64[peak_amount]
 *   abundance                   : float64[peak_amount]
 *   compound                    : int32[peak_amount, iso_amount] -- compound[i, k] is
 *                                  the count of isotope k in peak i, in the same
 *                                  element/isotope column order as isotope_names
 *   isotope_names                : list[str], length iso_amount
 *   element_names                : list[str], length element_amount
 *   element_all_iso_calc_amount  : int32[element_amount] -- raw per-element
 *                                  isotope-combination counts (elements[i].all_iso_calc_amount);
 *                                  this is what R's rica=1 branch returns, unfiltered/ungrouped
 */
static PyObject *
pattern_result_to_pyobject(const IsoPatternResult *result)
{
    npy_intp peak_amount = (npy_intp)result->peak_amount;
    npy_intp iso_amount = (npy_intp)result->iso_amount;
    npy_intp element_amount = (npy_intp)result->element_amount;

    PyObject *mass_arr = double_array_from_buffer(result->mass, peak_amount);
    if (mass_arr == NULL) {
        return NULL;
    }
    PyObject *abundance_arr = double_array_from_buffer(result->abundance, peak_amount);
    if (abundance_arr == NULL) {
        Py_DECREF(mass_arr);
        return NULL;
    }

    npy_intp compound_dims[2] = {peak_amount, iso_amount};
    PyObject *compound_arr = PyArray_SimpleNew(2, compound_dims, NPY_INT);
    if (compound_arr == NULL) {
        Py_DECREF(mass_arr);
        Py_DECREF(abundance_arr);
        return NULL;
    }
    if (peak_amount > 0 && iso_amount > 0) {
        memcpy(PyArray_DATA((PyArrayObject *)compound_arr), result->compound,
               (size_t)(peak_amount * iso_amount) * sizeof(int));
    }

    PyObject *isotope_names = PyList_New(iso_amount);
    if (isotope_names == NULL) {
        Py_DECREF(mass_arr);
        Py_DECREF(abundance_arr);
        Py_DECREF(compound_arr);
        return NULL;
    }
    for (npy_intp k = 0; k < iso_amount; k++) {
        /* names[0]="mass", names[1]="abundance", names[2..] = isotopes */
        PyList_SET_ITEM(isotope_names, k, name_from_slot(result->names, (size_t)k + 2));
    }

    PyObject *element_names = PyList_New(element_amount);
    npy_intp element_amount_dim = element_amount;
    PyObject *element_counts = PyArray_SimpleNew(1, &element_amount_dim, NPY_INT);
    if (element_names == NULL || element_counts == NULL) {
        Py_DECREF(mass_arr);
        Py_DECREF(abundance_arr);
        Py_DECREF(compound_arr);
        Py_DECREF(isotope_names);
        Py_XDECREF(element_names);
        Py_XDECREF(element_counts);
        return NULL;
    }
    int *counts_data = (int *)PyArray_DATA((PyArrayObject *)element_counts);
    for (npy_intp e = 0; e < element_amount; e++) {
        char buf[MAX_NAME_SIZE + 1];
        memcpy(buf, result->elements[e].name, MAX_NAME_SIZE);
        buf[MAX_NAME_SIZE] = '\0';
        PyList_SET_ITEM(element_names, e, PyUnicode_FromString(buf));
        counts_data[e] = result->elements[e].all_iso_calc_amount;
    }

    return Py_BuildValue(
        "{s:N,s:N,s:N,s:N,s:N,s:N}",
        "mass", mass_arr,
        "abundance", abundance_arr,
        "compound", compound_arr,
        "isotope_names", isotope_names,
        "element_names", element_names,
        "element_all_iso_calc_amount", element_counts);
}

static PyObject *
call_pattern_algo(PatternAlgoFn algo, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"sum", "peak_limit", "threshold", "rel_to_mono",
                              "iso_list_elem", "iso_list_iso",
                              "iso_list_mass", "iso_list_abu", NULL};
    const char *sum_in;
    int peak_limit;
    double threshold;
    int rel_to_mono;
    PyObject *elem_obj, *iso_obj, *mass_obj, *abu_obj;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "sidiOOOO", kwlist, &sum_in, &peak_limit,
                                      &threshold, &rel_to_mono, &elem_obj, &iso_obj,
                                      &mass_obj, &abu_obj)) {
        return NULL;
    }

    PyObject *elem_seq = PySequence_Fast(elem_obj, "iso_list_elem must be a sequence of str");
    if (elem_seq == NULL) {
        return NULL;
    }
    PyObject *iso_seq = PySequence_Fast(iso_obj, "iso_list_iso must be a sequence of str");
    if (iso_seq == NULL) {
        Py_DECREF(elem_seq);
        return NULL;
    }

    PyArrayObject *mass_arr = (PyArrayObject *)PyArray_FROM_OTF(mass_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *abu_arr = NULL;
    char *sum_copy = NULL;
    char *i_el_ = NULL;
    char *i_il_ = NULL;
    IsoPatternResult result;
    int result_initialized = 0;
    PyObject *ret = NULL;

    if (mass_arr == NULL) {
        goto cleanup;
    }
    abu_arr = (PyArrayObject *)PyArray_FROM_OTF(abu_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (abu_arr == NULL) {
        goto cleanup;
    }
    if (PyArray_NDIM(mass_arr) != 1 || PyArray_NDIM(abu_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "iso_list_mass and iso_list_abu must be 1-D arrays");
        goto cleanup;
    }

    {
        Py_ssize_t n_elem = PySequence_Fast_GET_SIZE(elem_seq);
        Py_ssize_t n_iso = PySequence_Fast_GET_SIZE(iso_seq);
        npy_intp n_mass = PyArray_DIM(mass_arr, 0);
        npy_intp n_abu = PyArray_DIM(abu_arr, 0);

        if (n_elem != n_iso || (npy_intp)n_iso != n_mass || n_mass != n_abu) {
            PyErr_SetString(PyExc_ValueError,
                             "iso_list_elem/iso_list_iso/iso_list_mass/iso_list_abu "
                             "must all have the same length");
            goto cleanup;
        }
        /* iso_amount_global is passed to iso_core as unsigned short. */
        if (n_elem <= 0 || n_elem > 65535) {
            PyErr_SetString(PyExc_ValueError, "iso_list length must be between 1 and 65535");
            goto cleanup;
        }

        size_t n = (size_t)n_elem;
        i_el_ = (char *)malloc(n * MAX_NAME_SIZE);
        i_il_ = (char *)malloc(n * MAX_NAME_SIZE);
        if (i_el_ == NULL || i_il_ == NULL) {
            PyErr_NoMemory();
            goto cleanup;
        }

        for (Py_ssize_t i = 0; i < n_elem; i++) {
            PyObject *e = PySequence_Fast_GET_ITEM(elem_seq, i);
            PyObject *s = PySequence_Fast_GET_ITEM(iso_seq, i);
            const char *e_str = PyUnicode_AsUTF8(e);
            if (e_str == NULL) {
                goto cleanup;
            }
            const char *s_str = PyUnicode_AsUTF8(s);
            if (s_str == NULL) {
                goto cleanup;
            }
            /* MAX_NAME_SIZE-1 chars + NUL, matching the fixed-width packed
               buffer main.c's SEXP wrapper builds via strcpy (see
               src/main.c's iso_pattern/iso_pattern_4). Enforced here
               (unlike the original strcpy, which has no bound check)
               since we control the buffer independently of R's own
               isotope-table invariants. */
            if (strlen(e_str) >= MAX_NAME_SIZE || strlen(s_str) >= MAX_NAME_SIZE) {
                PyErr_Format(PyExc_ValueError,
                             "element/isotope name longer than %d characters at index %zd",
                             MAX_NAME_SIZE - 1, i);
                goto cleanup;
            }
            strcpy(&i_el_[(size_t)i * MAX_NAME_SIZE], e_str);
            strcpy(&i_il_[(size_t)i * MAX_NAME_SIZE], s_str);
        }

        /* Defensive copy of the formula string -- mirrors main.c's own
           R_alloc+strcpy before calling into iso_core (parse_sum_formula_
           vector_R only reads through its sum_formula pointer, but we
           don't rely on that not changing in the future). */
        sum_copy = strdup(sum_in);
        if (sum_copy == NULL) {
            PyErr_NoMemory();
            goto cleanup;
        }

        char errmsg[ISO_CORE_ERRMSG_SIZE] = {0};
        int status = algo(sum_copy, peak_limit, threshold, rel_to_mono,
                           i_el_, i_il_,
                           (double *)PyArray_DATA(mass_arr), (double *)PyArray_DATA(abu_arr),
                           (unsigned short)n_elem, &result, errmsg);
        result_initialized = 1;
        if (status != ISO_CORE_OK) {
            set_core_error(errmsg, "iso_core pattern calculation failed");
            goto cleanup;
        }
    }

    ret = pattern_result_to_pyobject(&result);

cleanup:
    free(i_el_);
    free(i_il_);
    free(sum_copy);
    Py_XDECREF(mass_arr);
    Py_XDECREF(abu_arr);
    Py_DECREF(elem_seq);
    Py_DECREF(iso_seq);
    if (result_initialized) {
        iso_pattern_result_free(&result);
    }
    return ret;
}

static PyObject *
envipat_pattern_algo1(PyObject *self, PyObject *args, PyObject *kwargs)
{
    (void)self;
    return call_pattern_algo(iso_core_pattern_algo1, args, kwargs);
}

static PyObject *
envipat_pattern_algo2(PyObject *self, PyObject *args, PyObject *kwargs)
{
    (void)self;
    return call_pattern_algo(iso_core_pattern_algo4, args, kwargs);
}

/* ------------------------------------------------------------------- */
/* module definition                                                    */
/* ------------------------------------------------------------------- */

static PyMethodDef EnviPatMethods[] = {
    {"ppm_trace", (PyCFunction)(void (*)(void))envipat_ppm_trace, METH_VARARGS | METH_KEYWORDS,
     "ppm_trace(start, end, ppm) -> ndarray[float64]\n\n"
     "Generates a ppm-spaced mass trace between start and end. Mirrors "
     "the R package's iso_ppm_Call .Call entry point."},
    {"centroid", (PyCFunction)(void (*)(void))envipat_centroid, METH_VARARGS | METH_KEYWORDS,
     "centroid(profile_mass, profile_abundance, type) -> (ndarray, ndarray)\n\n"
     "Centroid (type=0), intensoid/local-maxima (type=1) or valley (type=2) "
     "detection on a profile. Mirrors iso_centroid_Call. profile_mass must "
     "be sorted ascending."},
    {"profile_with_trace", (PyCFunction)(void (*)(void))envipat_profile_with_trace, METH_VARARGS | METH_KEYWORDS,
     "profile_with_trace(profile_type, mass, abundance, trace, resolution, "
     "threshold, filter) -> (ndarray, ndarray)\n\n"
     "Convolves an isotope pattern (mass, abundance) onto a mass trace. "
     "profile_type: 0=Gaussian, 1=Cauchy-Lorentz. Mirrors "
     "iso_profile_with_trace_Call."},
    {"pattern_algo1", (PyCFunction)(void (*)(void))envipat_pattern_algo1, METH_VARARGS | METH_KEYWORDS,
     "pattern_algo1(sum, peak_limit, threshold, rel_to_mono, iso_list_elem, "
     "iso_list_iso, iso_list_mass, iso_list_abu) -> dict\n\n"
     "Isotope pattern calculation, algorithm 1. Mirrors iso_pattern."},
    {"pattern_algo2", (PyCFunction)(void (*)(void))envipat_pattern_algo2, METH_VARARGS | METH_KEYWORDS,
     "pattern_algo2(sum, peak_limit, threshold, rel_to_mono, iso_list_elem, "
     "iso_list_iso, iso_list_mass, iso_list_abu) -> dict\n\n"
     "Isotope pattern calculation, algorithm 2. Mirrors iso_pattern_4."},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef enviPat_py_module = {
    PyModuleDef_HEAD_INIT,
    "enviPat_py",
    "CPython/NumPy C-extension wrapping enviPat's pure-C isotope-pattern "
    "core (iso_core.c/h).",
    -1,
    EnviPatMethods,
};

PyMODINIT_FUNC
PyInit_enviPat_py(void)
{
    import_array(); /* on failure, expands to `return NULL;` here */

    PyObject *m = PyModule_Create(&enviPat_py_module);
    if (m == NULL) {
        return NULL;
    }

    EnviPatCoreError = PyErr_NewExceptionWithDoc(
        "envipat._core.enviPat_py.EnviPatCoreError",
        "Raised by this extension's algorithm entry points (pattern_algo1,\n"
        "pattern_algo2, ppm_trace, profile_with_trace, centroid) whenever the\n"
        "underlying iso_core_* call reports a failure (invalid input, or an\n"
        "internal ISO_CORE_ERR_* status -- see iso_core.h). Callers using the\n"
        "pure-Python pipeline (envipat.isopattern, envipat.envelope,\n"
        "envipat.vdetect) will not normally see this directly: those modules\n"
        "catch it and report the same condition via R's own literal \"error\"\n"
        "string sentinel instead of letting it propagate -- see each of\n"
        "their module docstrings.",
        NULL, NULL);
    if (EnviPatCoreError == NULL) {
        Py_DECREF(m);
        return NULL;
    }
    Py_INCREF(EnviPatCoreError);
    if (PyModule_AddObject(m, "EnviPatCoreError", EnviPatCoreError) < 0) {
        Py_DECREF(EnviPatCoreError);
        Py_CLEAR(EnviPatCoreError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
