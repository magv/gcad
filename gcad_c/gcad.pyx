# cython: language_level=3

from sympy import Poly
from fractions import Fraction
from libc.stdlib cimport calloc, free
from cython.cimports.cpython.long cimport PyLong_AsLongAndOverflow, PyLong_FromString
from cython.cimports.cpython.unicode cimport PyUnicode_AsUTF8AndSize

# FLINT <-> Python conversion utils

cdef extern from "<flint/fmpz.h>":
    ctypedef int fmpz
    void fmpz_init(fmpz *z)
    void fmpz_clear(fmpz *z)
    void fmpz_set_si(fmpz *z, long i)
    void fmpz_set_str(fmpz *z, const char *s, int base)
    char *fmpz_get_str(char *s, int base, const fmpz *z)
    long fmpz_get_si(const fmpz *z)
    int fmpz_fits_si(const fmpz *z)

cdef extern from "<flint/flint.h>":
    ctypedef int slong
    ctypedef int ulong

cdef extern from "<flint/fmpq.h>":
    ctypedef int fmpq
    void fmpq_init(fmpq *q)
    void fmpq_clear(fmpq *q)
    fmpz *fmpq_numref(fmpq *q)
    fmpz *fmpq_denref(fmpq *q)

cdef extern from "<flint/fmpq_poly.h>":
    ctypedef int fmpq_poly_struct
    void fmpq_poly_init(fmpq_poly_struct *poly)
    void fmpq_poly_clear(fmpq_poly_struct *poly)
    void fmpq_poly_fit_length(fmpq_poly_struct *poly, long length)
    void fmpq_poly_truncate(fmpq_poly_struct *poly, long n)
    void fmpq_poly_set_coeff_fmpq(fmpq_poly_struct *poly, long n, fmpq *x)
    void fmpq_poly_get_coeff_fmpq(fmpq *x, fmpq_poly_struct *poly, long n)
    long fmpq_poly_length(fmpq_poly_struct *poly)

cdef extern from "<flint/fmpq_mpoly.h>":
    ctypedef int fmpq_mpoly_struct
    ctypedef int fmpq_mpoly_ctx_struct
    cdef enum:
        ORD_LEX
        ORD_DEGLEX
        ORD_DEGREVLEX
    void fmpq_mpoly_ctx_init(fmpq_mpoly_ctx_struct *ctx, long nvars, int ord)
    void fmpq_mpoly_ctx_clear(fmpq_mpoly_ctx_struct *ctx)
    slong fmpq_mpoly_ctx_nvars(fmpq_mpoly_ctx_struct *ctx)
    void fmpq_mpoly_init(fmpq_mpoly_struct *p, fmpq_mpoly_ctx_struct *ctx)
    void fmpq_mpoly_clear(fmpq_mpoly_struct *p, fmpq_mpoly_ctx_struct *ctx)
    void fmpq_mpoly_zero(fmpq_mpoly_struct *p, fmpq_mpoly_ctx_struct *ctx)
    void fmpq_mpoly_reduce(fmpq_mpoly_struct *p, fmpq_mpoly_ctx_struct *ctx)
    long fmpq_mpoly_length(fmpq_mpoly_struct *p, fmpq_mpoly_ctx_struct *ctx)
    void fmpq_mpoly_get_term_coeff_fmpq(fmpq *c, fmpq_mpoly_struct *p, long i, fmpq_mpoly_ctx_struct *ctx)
    void fmpq_mpoly_get_term_exp_si(slong *exps, fmpq_mpoly_struct *p, long i, fmpq_mpoly_ctx_struct *ctx)
    void fmpq_mpoly_push_term_fmpq_ui(fmpq_mpoly_struct *p, fmpq *c, ulong *exp, fmpq_mpoly_ctx_struct *ctx)

cdef void fmpz_set_py_int(fmpz *z, object n):
    """Set an initialised fmpz to a Python int."""
    cdef long val
    cdef int overflow = 0
    val = PyLong_AsLongAndOverflow(n, &overflow)
    if overflow == 0:
        #static_assert(sizeof(slong) >= sizeof(long))
        fmpz_set_si(z, val)
    else:
        s = str(n)
        fmpz_set_str(z, PyUnicode_AsUTF8AndSize(s, NULL), 10)

cdef object fmpz_get_py_int(fmpz *z):
    """Convert an fmpz to a Python int."""
    cdef char *s
    if fmpz_fits_si(z):
        return fmpz_get_si(z)
    else:
        s = fmpz_get_str(NULL, 16, z)
        try:
            return PyLong_FromString(s, NULL, 16)
        finally:
            free(s)

cdef void fmpq_set_py_Fraction(fmpq *q, object frac):
    """Set an initialised fmpq to a Python Fraction."""
    fmpz_set_py_int(fmpq_numref(q), frac.numerator)
    fmpz_set_py_int(fmpq_denref(q), frac.denominator)

cdef object fmpq_get_py_Fraction(fmpq *q):
    """Convert an fmpq to a Python Fraction."""
    num = fmpz_get_py_int(fmpq_numref(q))
    den = fmpz_get_py_int(fmpq_denref(q))
    return Fraction(num, den)

cdef void fmpq_poly_set_sympy_Poly(fmpq_poly_struct *qp, object poly):
    """Set an initialised fmpq_poly to an univariate sympy.Poly."""
    cdef fmpq co
    cdef Py_ssize_t n
    coeffs = poly.all_coeffs()
    n = len(coeffs)
    fmpq_poly_fit_length(qp, n)
    fmpq_init(&co)
    try:
        for i in range(n):
            fmpq_set_py_Fraction(&co, coeffs[n - 1 - i])
            fmpq_poly_set_coeff_fmpq(qp, i, &co)
        fmpq_poly_truncate(qp, n)
    finally:
        fmpq_clear(&co)

cdef object fmpq_poly_get_sympy_Poly(fmpq_poly_struct *qp, object gen):
    """Convert an fmpq_poly to a sympy.Poly with the given generator."""
    cdef fmpq q
    cdef long i
    fmpq_init(&q)
    try:
        coeffs = []
        for i in range(fmpq_poly_length(qp) - 1, -1, -1):
            fmpq_poly_get_coeff_fmpq(&q, qp, i)
            coeffs.append(fmpq_get_py_Fraction(&q))
        return Poly.from_list(coeffs, gen)
    finally:
        fmpq_clear(&q)

cdef void fmpq_mpoly_set_sympy_Poly(fmpq_mpoly_struct *mp, fmpq_mpoly_ctx_struct *ctx, object poly):
    """Set an initialised fmpq_mpoly from a multivariate sympy.Poly."""
    cdef fmpq co
    cdef ulong *exps
    cdef Py_ssize_t nvars
    cdef Py_ssize_t i
    cdef Py_ssize_t n
    gens = poly.gens
    nvars = len(gens)
    fmpq_mpoly_zero(mp, ctx)
    fmpq_init(&co)
    exps = <ulong *>calloc(nvars, sizeof(ulong))
    try:
        for mon, coeff in poly.as_dict().items():
            fmpq_set_py_Fraction(&co, coeff)
            n = len(mon)
            for i in range(n):
                exps[i] = <ulong>mon[i]
            for i in range(n, nvars):
                exps[i] = 0
            fmpq_mpoly_push_term_fmpq_ui(mp, &co, exps, ctx)
        fmpq_mpoly_reduce(mp, ctx)
    finally:
        fmpq_clear(&co)
        free(exps)

cdef object fmpq_mpoly_get_sympy_Poly(fmpq_mpoly_struct *mp, fmpq_mpoly_ctx_struct *ctx, object gens):
    """Convert an fmpq_mpoly to a sympy.Poly with the given generators."""
    cdef fmpq co
    cdef slong *exps
    cdef slong nvars
    cdef slong i
    nvars = fmpq_mpoly_ctx_nvars(ctx)
    fmpq_init(&co)
    exps = <slong *>calloc(nvars, sizeof(slong))
    try:
        terms = {}
        for i in range(fmpq_mpoly_length(mp, ctx)):
            fmpq_mpoly_get_term_coeff_fmpq(&co, mp, i, ctx)
            fmpq_mpoly_get_term_exp_si(exps, mp, i, ctx)
            mon = tuple(exps[j] for j in range(nvars))
            terms[mon] = fmpq_get_py_Fraction(&co)
        return Poly.from_dict(terms, *gens)
    finally:
        fmpq_clear(&co)
        free(exps)

# Round-trip test functions

def _identity_int(i: object) -> int:
    cdef fmpz z
    fmpz_init(&z)
    try:
        fmpz_set_py_int(&z, i)
        return fmpz_get_py_int(&z)
    finally:
        fmpz_clear(&z)

def _identity_Fraction(f: Fraction) -> Fraction:
    cdef fmpq q
    fmpq_init(&q)
    try:
        fmpq_set_py_Fraction(&q, f)
        return fmpq_get_py_Fraction(&q)
    finally:
        fmpq_clear(&q)

def _identity_univariate_Poly(object poly):
    cdef fmpq_poly_struct qp
    fmpq_poly_init(&qp)
    try:
        fmpq_poly_set_sympy_Poly(&qp, poly)
        return fmpq_poly_get_sympy_Poly(&qp, poly.gen)
    finally:
        fmpq_poly_clear(&qp)

def _identity_multivariate_Poly(object poly):
    cdef fmpq_mpoly_struct mp
    cdef fmpq_mpoly_ctx_struct ctx
    cdef long nvars
    nvars = len(poly.gens)
    fmpq_mpoly_ctx_init(&ctx, nvars, ORD_LEX)
    fmpq_mpoly_init(&mp, &ctx)
    try:
        fmpq_mpoly_set_sympy_Poly(&mp, &ctx, poly)
        return fmpq_mpoly_get_sympy_Poly(&mp, &ctx, poly.gens)
    finally:
        fmpq_mpoly_clear(&mp, &ctx)
        fmpq_mpoly_ctx_clear(&ctx)

# API

cdef extern from "shortest_fraction_between.c":
    void _shortest_fraction_between "shortest_fraction_between"(fmpq *res, const fmpq *x, const fmpq *y)

def shortest_fraction_between(a: Fraction, b: Fraction) -> Fraction:
    """Simplest fraction between, or equal to, two rationals."""
    cdef fmpq qa, qb, qres
    fmpq_init(&qa)
    fmpq_init(&qb)
    fmpq_init(&qres)
    try:
        fmpq_set_py_Fraction(&qa, a)
        fmpq_set_py_Fraction(&qb, b)
        _shortest_fraction_between(&qres, &qa, &qb)
        return fmpq_get_py_Fraction(&qres)
    finally:
        fmpq_clear(&qa)
        fmpq_clear(&qb)
        fmpq_clear(&qres)
