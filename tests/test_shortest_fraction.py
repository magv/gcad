from gcad_c_ext import *

from fractions import Fraction

def test_shortest_fraction_zero():
    f1 = Fraction(-1, 3)
    f2 = Fraction(2, 7)
    f = shortest_fraction_between(f1, f2)
    assert f == Fraction(0)

def test_shortest_fraction_int():
    f1 = Fraction(5, 3)
    f2 = Fraction(21, 7)
    f = shortest_fraction_between(f1, f2)
    assert f == Fraction(2)

def test_shortest_fraction_negint():
    f1 = Fraction(-50, 3)
    f2 = Fraction(-8, 7)
    f = shortest_fraction_between(f1, f2)
    assert f == Fraction(-2)

def test_shortest_fraction_close():
    f1 = Fraction(31, 107)
    f2 = Fraction(32, 107)
    f = shortest_fraction_between(f1, f2)
    assert f == Fraction(5, 17)

def test_shortest_fraction_closer():
    f1 = Fraction(193, 7121)
    f2 = Fraction(191, 7109)
    f = shortest_fraction_between(f1, f2)
    assert f == Fraction(9, 332)
