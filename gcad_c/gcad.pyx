# cython: language_level=3

from fractions import Fraction
from libc.stdlib cimport malloc, free

# FLINT <-> Python conversion utils

cdef extern from "<Python.h>":
    char *PyUnicode_AsUTF8AndSize(object o, long *size)
    long PyLong_AsLong(object o)
    int PyErr_Occurred()
    void PyErr_Clear()

cdef extern from "<flint/fmpz.h>":
    ctypedef int fmpz
    void fmpz_init(fmpz *out)
    void fmpz_clear(fmpz *out)
    void fmpz_set_si(fmpz *out, long x)
    void fmpz_set_str(fmpz *out, const char *s, int base)
    char *fmpz_get_str(char *s, int base, const fmpz *x)
    long fmpz_get_si(const fmpz *f)
    int fmpz_fits_si(const fmpz *f)

cdef extern from "<flint/fmpq.h>":
    ctypedef int fmpq
    void fmpq_init(fmpq *q)
    void fmpq_clear(fmpq *q)
    fmpz *fmpq_numref(fmpq *f)
    fmpz *fmpq_denref(fmpq *f)

cdef extern from "shortest_fraction_between.c":
    void _shortest_fraction_between "shortest_fraction_between"(fmpq *res, const fmpq *x, const fmpq *y)

cdef void fmpz_set_py_int(fmpz *z, object n):
    """Set an initialised fmpz to a Python int."""
    cdef long val
    val = PyLong_AsLong(n)
    if not PyErr_Occurred():
        fmpz_set_si(z, val)
    else:
        PyErr_Clear()
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
            return int(s, 16)
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

# API

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
        result = fmpq_get_py_Fraction(&qres)
    finally:
        fmpq_clear(&qa)
        fmpq_clear(&qb)
        fmpq_clear(&qres)
    return result
