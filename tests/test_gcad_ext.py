from fractions import Fraction
from gcad_ext import *
from gcad_ext import (
    _identity_Fraction,
    _identity_int,
    _identity_multivariate_Poly,
    _identity_univariate_Poly,
)
import sympy as sp

def test_identity_int():
    assert _identity_int(3) == 3
    assert _identity_int(-7) == -7
    assert _identity_int(-2**90) == -2**90
    try:
        _identity_int("zzz")
        assert False
    except Exception:
        pass

def test_identity_fraction():
    assert _identity_Fraction(Fraction(3, 4)) == Fraction(3, 4)
    assert _identity_Fraction(Fraction(-9, 1)) == Fraction(-9, 1)
    assert _identity_Fraction(Fraction(2**150, 3)) == Fraction(2**150, 3)
    assert _identity_Fraction(Fraction(7, 2**160)) == Fraction(7, 2**160)

def test_identity_fraction_conversions():
    assert _identity_Fraction(17) == Fraction(17, 1)
    assert _identity_Fraction(sp.Rational(23, 71)) == Fraction(23, 71)

def test_identity_univariate_Poly():
    x = sp.Symbol("x")
    _, p = sp.Poly(sp.sympify("x**2 - 2*x + 7"), x).clear_denoms(convert=True)
    q = _identity_univariate_Poly(p)
    assert q == p

def test_identity_multivariate_Poly():
    x, y = sp.symbols("x y")
    p = sp.Poly(sp.sympify("x*y**2 - 2*x + 4*y - 7"), x, y)
    q = _identity_multivariate_Poly(p)
    assert q == p
