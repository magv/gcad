# cython: language_level = 3
# cython: language = c++

from cpython.long cimport PyLong_AsLongAndOverflow, PyLong_FromString
from cpython.unicode cimport PyUnicode_AsUTF8AndSize
from fractions import Fraction
from libc.stdlib cimport calloc, free
from libcpp.vector cimport vector
from sympy import Poly

# FLINT <-> Python conversion utils

cdef extern from "<flint/flint.h>":
    ctypedef int slong
    ctypedef int ulong

cdef extern from "<flint/fmpz.h>":
    ctypedef int fmpz
    void fmpz_init(fmpz *z)
    void fmpz_clear(fmpz *z)
    void fmpz_set_si(fmpz *z, long i)
    void fmpz_set_str(fmpz *z, const char *s, int base)
    char *fmpz_get_str(char *s, int base, const fmpz *z)
    long fmpz_get_si(const fmpz *z)
    int fmpz_fits_si(const fmpz *z)

cdef extern from "<flint/fmpq.h>":
    ctypedef int fmpq
    void fmpq_init(fmpq *q)
    void fmpq_clear(fmpq *q)
    fmpz *fmpq_numref(fmpq *q)
    fmpz *fmpq_denref(fmpq *q)

cdef extern from "<flint/fmpz_poly.h>":
    ctypedef int fmpz_poly_struct
    void fmpz_poly_init(fmpz_poly_struct *poly)
    void fmpz_poly_clear(fmpz_poly_struct *poly)
    void fmpz_poly_fit_length(fmpz_poly_struct *poly, long length)
    void fmpz_poly_truncate(fmpz_poly_struct *poly, long n)
    void fmpz_poly_set_coeff_fmpz(fmpz_poly_struct *poly, long n, fmpz *x)
    void fmpz_poly_get_coeff_fmpz(fmpz *x, fmpz_poly_struct *poly, long n)
    long fmpz_poly_length(fmpz_poly_struct *poly)
    long fmpz_poly_degree(fmpz_poly_struct *poly)
    void fmpz_poly_shift_right(fmpz_poly_struct *res, const fmpz_poly_struct *poly, long n)
    fmpz *fmpz_poly_get_coeff_ptr(const fmpz_poly_struct *poly, long n)

cdef extern from "<flint/fmpz_mpoly.h>":
    ctypedef int fmpz_mpoly_struct
    ctypedef int fmpz_mpoly_ctx_struct
    cdef enum:
        ORD_LEX
        ORD_DEGLEX
        ORD_DEGREVLEX
    void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_struct *ctx, long nvars, int ord)
    void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_struct *ctx)
    slong fmpz_mpoly_ctx_nvars(fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_init(fmpz_mpoly_struct *p, fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_clear(fmpz_mpoly_struct *p, fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_set(fmpz_mpoly_struct *p, fmpz_mpoly_struct *q, fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_zero(fmpz_mpoly_struct *p, fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_reduce(fmpz_mpoly_struct *p, fmpz_mpoly_ctx_struct *ctx)
    long fmpz_mpoly_length(fmpz_mpoly_struct *p, fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_get_term_coeff_fmpz(fmpz *c, fmpz_mpoly_struct *p, long i, fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_get_term_exp_si(slong *exps, fmpz_mpoly_struct *p, long i, fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_push_term_fmpz_ui(fmpz_mpoly_struct *p, const fmpz *c, const ulong *exp, const fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_sort_terms(fmpz_mpoly_struct *p, const fmpz_mpoly_ctx_struct *ctx)
    int fmpz_mpoly_discriminant(fmpz_mpoly_struct *R, const fmpz_mpoly_struct *A, slong var, const fmpz_mpoly_ctx_struct *ctx)
    int fmpz_mpoly_resultant(fmpz_mpoly_struct *R, const fmpz_mpoly_struct *A, const fmpz_mpoly_struct *B, slong var, const fmpz_mpoly_ctx_struct *ctx);

cdef extern from "<flint/fmpz_mpoly_factor.h>":
    ctypedef int fmpz_mpoly_factor_struct
    void fmpz_mpoly_factor_init(fmpz_mpoly_factor_struct *f, fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_factor_clear(fmpz_mpoly_factor_struct *f, fmpz_mpoly_ctx_struct *ctx)
    slong fmpz_mpoly_factor_length(fmpz_mpoly_factor_struct *f, fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_factor_get_constant_fmpz(fmpz *c, fmpz_mpoly_factor_struct *f, fmpz_mpoly_ctx_struct *ctx)
    void fmpz_mpoly_factor_get_base(fmpz_mpoly_struct *B, fmpz_mpoly_factor_struct *f, slong i, fmpz_mpoly_ctx_struct *ctx)
    slong fmpz_mpoly_factor_get_exp_si(fmpz_mpoly_factor_struct *f, slong i, fmpz_mpoly_ctx_struct *ctx)
    int fmpz_mpoly_factor(fmpz_mpoly_factor_struct *f, fmpz_mpoly_struct *A, fmpz_mpoly_ctx_struct *ctx)

cdef void fmpz_set_py_int(fmpz *z, object n):
    """Set an initialised fmpz to a Python int."""
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

cdef void fmpz_poly_set_sympy_Poly(fmpz_poly_struct *qp, object poly):
    """Set an initialised fmpz_poly to an univariate sympy.Poly."""
    cdef fmpz z
    coeffs = poly.all_coeffs()
    n = len(coeffs)
    fmpz_init(&z)
    try:
        for i in range(n):
            fmpz_set_py_int(&z, coeffs[n - 1 - i])
            fmpz_poly_set_coeff_fmpz(qp, i, &z)
    finally:
        fmpz_clear(&z)

cdef object fmpz_poly_get_sympy_Poly(fmpz_poly_struct *qp, object gen):
    """Convert an fmpz_poly to a sympy.Poly with the given generator."""
    cdef fmpz z
    fmpz_init(&z)
    try:
        coeffs = []
        for i in range(fmpz_poly_length(qp) - 1, -1, -1):
            fmpz_poly_get_coeff_fmpz(&z, qp, i)
            coeffs.append(fmpz_get_py_int(&z))
        return Poly.from_list(coeffs, gen)
    finally:
        fmpz_clear(&z)

cdef void fmpz_mpoly_set_sympy_Poly(fmpz_mpoly_struct *mp, fmpz_mpoly_ctx_struct *ctx, object poly):
    """Set an initialised fmpz_mpoly from a multivariate sympy.Poly."""
    cdef fmpz co
    cdef ulong *exps
    gens = poly.gens
    nvars = len(gens)
    fmpz_mpoly_zero(mp, ctx)
    fmpz_init(&co)
    exps = <ulong *>calloc(nvars, sizeof(ulong))
    try:
        for monomial, coeff in poly.as_dict().items():
            fmpz_set_py_int(&co, coeff)
            n = len(monomial)
            assert n == nvars
            for i in range(n):
                exps[i] = <ulong>monomial[i]
            fmpz_mpoly_push_term_fmpz_ui(mp, &co, exps, ctx)
        #fmpz_mpoly_reduce(mp, ctx)
        fmpz_mpoly_sort_terms(mp, ctx)
    finally:
        fmpz_clear(&co)
        free(exps)

cdef object fmpz_mpoly_get_sympy_Poly(fmpz_mpoly_struct *mp, fmpz_mpoly_ctx_struct *ctx, object gens):
    """Convert an fmpz_mpoly to a sympy.Poly with the given generators."""
    cdef fmpz co
    cdef slong *exps
    cdef slong nvars
    cdef slong i
    nvars = fmpz_mpoly_ctx_nvars(ctx)
    fmpz_init(&co)
    exps = <slong *>calloc(nvars, sizeof(slong))
    try:
        terms = {}
        for i in range(fmpz_mpoly_length(mp, ctx)):
            fmpz_mpoly_get_term_coeff_fmpz(&co, mp, i, ctx)
            fmpz_mpoly_get_term_exp_si(exps, mp, i, ctx)
            mon = tuple(exps[j] for j in range(nvars))
            terms[mon] = fmpz_get_py_int(&co)
        return Poly.from_dict(terms, *gens)
    finally:
        fmpz_clear(&co)
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

def _identity_univariate_Poly(poly: object):
    cdef fmpz_poly_struct qp
    fmpz_poly_init(&qp)
    try:
        fmpz_poly_set_sympy_Poly(&qp, poly)
        return fmpz_poly_get_sympy_Poly(&qp, poly.gen)
    finally:
        fmpz_poly_clear(&qp)

def _identity_multivariate_Poly(poly: object):
    cdef fmpz_mpoly_struct mp
    cdef fmpz_mpoly_ctx_struct ctx
    nvars = len(poly.gens)
    fmpz_mpoly_ctx_init(&ctx, nvars, ORD_LEX)
    fmpz_mpoly_init(&mp, &ctx)
    try:
        fmpz_mpoly_set_sympy_Poly(&mp, &ctx, poly)
        return fmpz_mpoly_get_sympy_Poly(&mp, &ctx, poly.gens)
    finally:
        fmpz_mpoly_clear(&mp, &ctx)
        fmpz_mpoly_ctx_clear(&ctx)

# FLINT API

def discriminant(poly: object, varidx: int):
    """
    Find a discriminant of a multivariate poly with integer
    coefficients in the variable with the given index.
    """
    cdef fmpz_mpoly_struct p
    cdef fmpz_mpoly_ctx_struct ctx
    nvars = len(poly.gens)
    fmpz_mpoly_ctx_init(&ctx, nvars, ORD_LEX)
    fmpz_mpoly_init(&p, &ctx)
    try:
        fmpz_mpoly_set_sympy_Poly(&p, &ctx, poly)
        fmpz_mpoly_discriminant(&p, &p, varidx, &ctx)
        return fmpz_mpoly_get_sympy_Poly(&p, &ctx, poly.gens)
    finally:
        fmpz_mpoly_clear(&p, &ctx)
        fmpz_mpoly_ctx_clear(&ctx)

def resultant(poly1: object, poly2: object, varidx: int):
    """
    Find a resultant of two multivariate polys with integer
    coefficients in the variable with the given index.
    """
    cdef fmpz_mpoly_struct p1, p2
    cdef fmpz_mpoly_ctx_struct ctx
    assert poly1.gens == poly2.gens
    nvars = len(poly1.gens)
    fmpz_mpoly_ctx_init(&ctx, nvars, ORD_LEX)
    fmpz_mpoly_init(&p1, &ctx)
    fmpz_mpoly_init(&p2, &ctx)
    try:
        fmpz_mpoly_set_sympy_Poly(&p1, &ctx, poly1)
        fmpz_mpoly_set_sympy_Poly(&p2, &ctx, poly2)
        fmpz_mpoly_resultant(&p1, &p1, &p2, varidx, &ctx)
        return fmpz_mpoly_get_sympy_Poly(&p1, &ctx, poly1.gens)
    finally:
        fmpz_mpoly_clear(&p1, &ctx)
        fmpz_mpoly_clear(&p2, &ctx)
        fmpz_mpoly_ctx_clear(&ctx)

def factor(poly: object):
    """
    Factor a multivariate polynomials with integer coefficients.
    Return (content, [(factor, exp), ...]).
    """
    cdef fmpz_mpoly_ctx_struct ctx
    cdef fmpz_mpoly_struct mp
    cdef fmpz_mpoly_struct base
    cdef fmpz_mpoly_factor_struct fac
    cdef fmpz content
    nvars = len(poly.gens)
    fmpz_mpoly_ctx_init(&ctx, nvars, ORD_LEX)
    fmpz_mpoly_init(&mp, &ctx)
    fmpz_mpoly_init(&base, &ctx)
    fmpz_mpoly_factor_init(&fac, &ctx)
    fmpz_init(&content)
    try:
        fmpz_mpoly_set_sympy_Poly(&mp, &ctx, poly)
        if fmpz_mpoly_factor(&fac, &mp, &ctx) == 0:
            raise RuntimeError("fmpz_mpoly_factor failed")
        fmpz_mpoly_factor_get_constant_fmpz(&content, &fac, &ctx)
        factors = []
        for i in range(fmpz_mpoly_factor_length(&fac, &ctx)):
            fmpz_mpoly_factor_get_base(&base, &fac, i, &ctx)
            exp = int(fmpz_mpoly_factor_get_exp_si(&fac, i, &ctx))
            factors.append((fmpz_mpoly_get_sympy_Poly(&base, &ctx, poly.gens), exp))
        return fmpz_get_py_int(&content), factors
    finally:
        fmpz_clear(&content)
        fmpz_mpoly_factor_clear(&fac, &ctx)
        fmpz_mpoly_clear(&base, &ctx)
        fmpz_mpoly_clear(&mp, &ctx)
        fmpz_mpoly_ctx_clear(&ctx)

# API

cdef extern from "shortest_fraction_between.cpp":
    void _shortest_fraction_between "shortest_fraction_between"(fmpq *res, const fmpq *x, const fmpq *y)

cdef extern from "root_isolation.cpp":
    ctypedef int ZPoly
    ctypedef int Z
    ctypedef int Q
    ctypedef struct RootInterval:
        Q lo, hi
    ZPoly ZPoly_of_fmpz_poly(fmpz_poly_struct *p);
    Z Z_of_fmpz(fmpz *z);
    Q Q_of_fmpq(fmpq *q);
    fmpz *Z_to_fmpz(Z &z);
    fmpq *Q_to_fmpq(Q &z);
    Q _nth_root_ub "nth_root_ub"(ulong n, const Q &c)
    Q _poly_root_ub "poly_root_ub"(ZPoly &p)
    Q _poly_root_lb "poly_root_lb"(ZPoly &p)
    vector[RootInterval] _isolate_positive_roots "isolate_positive_roots"(ZPoly &poly)
    vector[RootInterval] _isolate_roots "isolate_roots"(ZPoly &poly)
    vector[vector[RootInterval]] _isolate_many_roots "isolate_many_roots"(vector[ZPoly] &polys)

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

def isolate_roots(poly: object):
    """
    Isolate all real roots of a square-free polynomial with
    rational coefficients. Return an ordered disjoined list of
    intervals (a, b), each with a single root. The intervals
    should be interpreted as open (i.e. the boundaries are
    excluded), unless a=b.
    """
    cdef fmpz_poly_struct qp
    cdef RootInterval ri
    result = []
    fmpz_poly_init(&qp)
    try:
        fmpz_poly_set_sympy_Poly(&qp, poly)
        intervals = _isolate_roots(ZPoly_of_fmpz_poly(&qp))
        for i in intervals:
            result.append((
                fmpq_get_py_Fraction(Q_to_fmpq(i.lo)),
                fmpq_get_py_Fraction(Q_to_fmpq(i.hi)),
            ))
        return result
    finally:
        fmpz_poly_clear(&qp)

def isolate_many_roots(polys: object):
    """
    Isolate all real roots of a list of square-free co-prime
    polynomials with rational coefficients. For each polynomial,
    return an ordered disjoined list of intervals (a, b), each
    with a single root. The intervals belonging to different
    polynomials are guaranteed to neither overlap, nor touch at
    the boundaries. The intervals are open (i.e. the boundaries
    are excluded), unless a=b.
    """
    cdef vector[ZPoly] cpp_polys
    cdef fmpz_poly_struct qp
    result = []
    for poly in polys:
        fmpz_poly_init(&qp)
        try:
            fmpz_poly_set_sympy_Poly(&qp, poly)
            cpp_polys.push_back(ZPoly_of_fmpz_poly(&qp))
        finally:
            fmpz_poly_clear(&qp)
    result_cpp = _isolate_many_roots(cpp_polys)
    for i in range(result_cpp.size()):
        inner = []
        for ri in result_cpp[i]:
            inner.append((
                fmpq_get_py_Fraction(Q_to_fmpq(ri.lo)),
                fmpq_get_py_Fraction(Q_to_fmpq(ri.hi)),
            ))
        result.append(inner)
    return result
