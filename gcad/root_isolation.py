"""
Real root isolation via the Vincent-Akritas-Strzeboński
Continued-Fractions (VAS-CF) algorithm, as given in ASV08 and AS05,
used together with root bounds computed via the LMQ algorithm,
as given in AAS08 and ASV08 (based on ASV06).

ASV08:
    A. Akritas, A. Strzeboński, P. Vigklas.
    "Improving the Performance of the Continued Fractions Method Using New Bounds of Positive Roots".
    Nonlinear Anal.: Model. Control, Vol. 13, No. 3, 265--279 (2008).
    https://doi.org/10.15388/NA.2008.13.3.14557

AAS08:
    A. Akritas, A. Argyris, A. Strzeboński.
    "FLQ, the Fastest Quadratic Complexity Bound on the Values of Positive Roots of Polynomials".
    Serdica J. Comput., Vol. 2, No. 2, 145--162 (2008).
    https://doi.org/10.55630/sjc.2008.2.145-162

ASV06:
    A. Akritas, A. Strzeboński, P. Vigklas.
    "Implementations of a New Theorem for Computing Bounds for Positive Roots of Polynomials".
    Comput. 78, 355--367 (2006).
    https://doi.org/10.1007/s00607-006-0186-y

AS05:
    A. Akritas, A. Strzeboński.
    "A Comparative Study of Two Real Root Isolation Methods".
    Nonlinear Anal.: Model. Control, Vol. 10, No. 4, 297--304 (2005).
    https://doi.org/10.15388/NA.2005.10.4.15110
"""

from __future__ import annotations

from sympy import Rational, Poly
import heapq
import math

def sign(r: Rational) -> int:
    return 0 if r.is_zero else +1 if r.is_positive else -1

def sgc(p: Poly) -> int:
    """
    The number of sign changes in the sequence of nonzero
    coefficients of p.
    """
    signs = [sign(c) for c in p.all_coeffs() if c != 0]
    return sum(1 for i in range(1, len(signs)) if signs[i] != signs[i - 1])

def shift(p: Poly, k: int) -> Poly:
    """Find g(x) = p(x + k)."""
    return p.compose(Poly(p.gen + k, p.gen))

def scale(p: Poly, alpha: int) -> Poly:
    """Find g(x) = p(alpha * x)."""
    return p.compose(Poly(alpha * p.gen, p.gen))

def reciprocal_transform(p: Poly) -> Poly:
    """Find g(x) = (x+1)^{deg(p)} * p(1/(x+1))."""
    rev = Poly(list(reversed(p.all_coeffs())), p.gen)
    return rev.compose(Poly(p.gen + 1, p.gen))

def negate(p: Poly) -> Poly:
    """Find g(x) = p(-x)."""
    return p.compose(Poly(-p.gen, p.gen))

def nth_root_ub(n: int, c: Rational) -> Rational:
    """
    Upper bound of the roots of p(x) = x^n - c. I.e. a rational
    that is larger than, or equal to, c^{1/n}.
    """
    assert c > 0
    assert n > 0
    if c == 1 or n == 1:
        return c
    else:
        approx = float(c) ** (1 / n)
        # approx == m * 2**e
        m, e = math.frexp(approx)
        # Making the bound very tight here does not really improve
        # the algorithm, so we won't be trying too hard.
        m = int(math.ceil(1024 * m))
        e = e - 10
        while True:
            y = m * Rational(2) ** e
            if c <= y**n:
                return y
            m += 1

def poly_root_ub(p: Poly) -> Rational:
    """
    A possibly tight upper bound on the real positive roots of
    p, according to the LMQ algorithm from AAS08 (Algorithm 3.1,
    the "Local-Max" Quadratic implementation of Theorem 3). Also
    benchmarked and recommended in ASV08.
    """
    n = p.degree()
    if n <= 0:
        return Rational(0)
    # p = \sum_{i=0}^{n} cl[i] * x^i
    cl = list(reversed(p.all_coeffs()))
    used = [1] * (n + 1)
    tmax = Rational(0)
    for mu in range(n - 1, -1, -1):
        if cl[mu] < 0:
            tmin = None
            for nu in range(n, mu, -1):
                if cl[nu] > 0:
                    t = nth_root_ub(
                        nu - mu, (-cl[mu]) / (cl[nu] / Rational(2) ** used[nu])
                    )
                    used[nu] += 1
                    if tmin is None or t < tmin:
                        tmin = t
            if tmin is not None and tmin > tmax:
                tmax = tmin
    return tmax

def poly_root_lb(p: Poly) -> Rational:
    """
    A possibly tight lower bound on the positive roots of the
    given polynomial, computed as 1/poly_root_ub(x^n p(1/x)).
    """
    assert p.nth(0) != 0
    rev_p = Poly(list(reversed(p.all_coeffs())), p.gen)
    if rev_p.LC().is_negative:
        rev_p = -rev_p
    ub = poly_root_ub(rev_p)
    assert ub > 0
    return 1/ub

def _intrv(a, b, c, d, p: Poly) -> tuple[Rational, Rational]:
    """
    Return the interval (min(b/d, a/c), max(b/d, a/c)). When
    c=0, the interval is (b/d, ∞), but instead of ∞, we use
    an upper bound on the roots of p.
    """
    v1 = Rational(b, d)
    if c == 0:
        # Use poly_root_ub for a finite right endpoint instead of
        # infinity, as the original algorithm would do. This is
        # to help with root refinement later. The only concern
        # is to handle the case of upper bound being exactly the
        # root, as we want to keep our intervals open.
        ub = int(poly_root_ub(p)) + 1
        #assert p.eval(ub) != 0
        v2 = Rational(a * ub + b, d)
        assert v1 < v2
        return (v1, v2)
    else:
        v2 = Rational(a, c)
        return (v1, v2) if v1 <= v2 else (v2, v1)

def isolate_positive_roots(poly: Poly, alpha0: int = 16) -> list[tuple[Rational, Rational]]:
    """
    Isolate all real positive roots of a square-free polynomial
    with rational coefficients. Return an ordered disjoined
    list of intervals (a, b), each with a single root. The
    intervals should be interpreted as open (i.e. the boundaries
    are excluded), unless a=b. We assume no roots at zero, and
    a positive sign of the leading coefficient.
    """
    # Assuming no roots at zero, degree > 1, and positive leading
    # coefficient.
    assert poly.degree() > 1
    assert poly.nth(0) != 0
    assert sign(poly.LC()) > 0
    rootlist: list[tuple[Rational, Rational]] = []
    s = sgc(poly)
    if s == 0:
        return rootlist
    if s == 1:
        # Use poly_root_ub for a finite right endpoint instead of
        # infinity, as the original algorithm would do. This is
        # to help with root refinement later. The only concern
        # is to handle the case of upper bound being exactly the
        # root, as we want to keep our intervals open.
        lo = Rational(0)
        hi = int(poly_root_ub(poly)) + 1
        rootlist.append((lo, hi))
        return rootlist
    intervalstack: list[tuple[int, int, int, int, Poly, int]] = [(1, 0, 0, 1, poly, s)]
    while intervalstack:
        a, b, c, d, p, s = intervalstack.pop()
        # Rounding down the lower bound to an int, as ASV08
        # requests. Keeping the fractional part is also allowed
        # here, but will not improve the performance, and integers
        # are shorter. With this, some exact integer roots will
        # be recognized.
        alpha = int(poly_root_lb(p))
        # Rescale if lower bound is large.
        if alpha > alpha0:
            p = scale(p, alpha)
            a = alpha * a
            c = alpha * c
            alpha = 1
        if alpha >= 1:
            # Shift by the lower bound to be closer to the root.
            p = shift(p, alpha)
            b = alpha * a + b
            d = alpha * c + d
            if p.nth(0) == 0:
                # Exact root at zero. Divide it out.
                rootlist.append((Rational(b, d), Rational(b, d)))
                p = p.quo(Poly(p.gen, p.gen))
                # Dividing only once, because p is supposed to
                # be square-free.
                assert p.nth(0) != 0
            s = sgc(p)
            if s == 0:
                continue
            if s == 1:
                rootlist.append(_intrv(a, b, c, d, p))
                continue
        # The root should be close. Let's split the interval at
        # x=1, and look at both sides:
        # - p1(x) = p(x+1) will cover x ∈ (1; ∞);
        # - p2(x) = (x+1)^n p(1/(x+1)) will cover x ∈ (0; 1).
        p1 = shift(p, 1)
        a1, b1, c1, d1 = a, a + b, c, c + d
        if p1.nth(0) == 0:
            rootlist.append((Rational(b1, d1), Rational(b1, d1)))
            p1 = p1.quo(Poly(p1.gen, p1.gen))
            r = 1
        else:
            r = 0
        s1 = sgc(p1)
        # The whole s2_bound business is tricky, and only valid
        # conjecturally. AS05 mentions the conjecture that
        # sgc(p1)+sgc(p2) <= sgc(p). Here, however, we assume
        # more: we assume that if s2_bound <= 0, then sgc(p2)
        # is 0, and that if s2_bound is 1, then sgc(p2) is 1.
        s2_bound = s - s1 - r
        a2, b2, c2, d2 = b, a + b, d, c + d
        # By the conjecture above, we can skip the computation
        # of p2 in some cases.
        if s2_bound > 1:
            p2 = reciprocal_transform(p)
            if p2.LC().is_negative:
                p2 = -p2
            if p2.nth(0) == 0:
                # p2(0) = p1(0) = p(1) = 0. If we are here,
                # this root was already counted.
                p2 = p2.quo(Poly(p2.gen, p2.gen))
            s2 = sgc(p2)
        else:
            p2 = None # This value will not be used.
            s2 = s2_bound
        # To keep the interval stack small, it's best to push
        # intervals with higher sgc first.
        if s1 < s2:
            a1, b1, c1, d1, p1, s1, a2, b2, c2, d2, p2, s2 = \
                a2, b2, c2, d2, p2, s2, a1, b1, c1, d1, p1, s1
        if s1 == 0:
            continue
        elif s1 == 1:
            rootlist.append(_intrv(a1, b1, c1, d1, p1))
        else:
            intervalstack.append((a1, b1, c1, d1, p1, s1))
        if s2 == 0:
            continue
        elif s2 == 1:
            rootlist.append(_intrv(a2, b2, c2, d2, p2))
        else:
            intervalstack.append((a2, b2, c2, d2, p2, s2))
    rootlist.sort()
    return rootlist

def isolate_roots(p: Poly) -> list[tuple[Rational, Rational]]:
    """
    Isolate all real roots of a square-free polynomial with
    rational coefficients. Return an ordered disjoined list of
    intervals (a, b), each with a single root. The intervals
    should be interpreted as open (i.e. the boundaries are
    excluded), unless a=b.
    """
    roots = []
    if p.degree() <= 0:
        return roots
    if p.nth(0) == 0:
        # Divide out the root at zero.
        p = p.quo(Poly(p.gen, p.gen))
        roots.append((Rational(0), Rational(0)))
    if p.degree() <= 0:
        return roots
    # The input should be square-free.
    assert p.nth(0) != 0
    # Special-case linear polys.
    if p.degree() == 1:
        v = -p.nth(0)/p.nth(1)
        roots.append((v, v))
        return roots
    # Negative roots.
    neg_p = negate(p)
    if neg_p.LC().is_negative:
        neg_p = -neg_p
    for (lo, hi) in isolate_positive_roots(neg_p):
        roots.append((-hi, -lo))
    roots.reverse()
    # Positive roots.
    if p.LC().is_negative:
        p = -p
    for (lo, hi) in isolate_positive_roots(p):
        roots.append((lo, hi))
    return roots

def bisect_root(p: Poly, lo: Rational, hi: Rational) -> tuple[Rational, Rational]:
    """Refine a real root interval for a given polynomial once."""
    sign_lo = sign(p.eval(lo))
    mid = (lo + hi)/2
    sign_mid = sign(p.eval(mid))
    if sign_mid == 0:
        # Exact root at mid.
        return (mid, mid)
    elif sign_mid != sign_lo:
        # Root in (lo, mid) => refine from right.
        return (lo, mid)
    else:
        # Root in (mid, hi) => refine from left.
        return (mid, hi)

def isolate_many_roots(polys: list[Poly]) -> list[list[tuple[Rational, Rational]]]:
    """
    Isolate all real roots of a list of square-free co-prime
    polynomials with rational coefficients. For each polynomial,
    return an ordered disjoined list of intervals (a, b), each
    with a single root. The intervals belonging to different
    polynomials are guaranteed to neither overlap, nor touch at
    the boundaries. The intervals are open (i.e. the boundaries
    are excluded), unless a=b.
    """
    polys = list(polys)
    # Isolate roots of each poly separately.
    intervals = []
    for idx, p in enumerate(polys):
        for lo, hi in isolate_roots(p):
            if lo == hi:
                # We'll divide out exact roots immediately, these
                # don't need bisection. They still need to be
                # accounted for, to make sure other intervals
                # don't overlap with them.
                polys[idx] = p = p.quo(Poly(p.gen - lo, p.gen))
            intervals.append((lo, hi, idx))
    # Sort by the lower boundary.
    heapq.heapify(intervals)
    # Bisect overlapping intervals till none are.
    result = [[] for i in range(len(polys))]
    if not intervals:
        return result
    lo1, hi1, idx1 = heapq.heappop(intervals)
    if not intervals:
        result[idx1].append((lo1, hi1))
        return result
    lo2, hi2, idx2 = heapq.heappop(intervals)
    while True:
        if hi1 < lo2:
            # No intersection. Interval 1 is done.
            if intervals:
                result[idx1].append((lo1, hi1))
                lo1, hi1, idx1 = lo2, hi2, idx2
                lo2, hi2, idx2 = heapq.heappop(intervals)
            else:
                result[idx1].append((lo1, hi1))
                result[idx2].append((lo2, hi2))
                return result
        else:
            # Intersecting intervals. Will bisect and retry.
            while not ((hi1 < lo2) or (hi2 < lo1)):
                if lo1 != hi1:
                    lo1, hi1 = bisect_root(polys[idx1], lo1, hi1)
                if lo2 != hi2:
                    lo2, hi2 = bisect_root(polys[idx2], lo2, hi2)
            if hi2 < lo1:
                # Interval 2 < interval 1.
                lo1, hi1, idx1, lo2, hi2, idx2 = \
                    lo2, hi2, idx2, lo1, hi1, idx1
            lo1, hi1, idx1 = heapq.heappushpop(intervals, (lo1, hi1, idx1))
            lo2, hi2, idx2 = heapq.heappushpop(intervals, (lo2, hi2, idx2))
