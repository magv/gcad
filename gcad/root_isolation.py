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

import heapq
import math
from sympy import Rational, Poly

def sgc(p: Poly) -> int:
    """
    The number of sign changes in the sequence of nonzero
    coefficients of p.
    """
    signs = [1 if c > 0 else -1 for c in p.all_coeffs() if c != 0]
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

def sign(r: Rational) -> int:
    return 0 if r.is_zero else +1 if r.is_positive else -1

def nth_root_ub(n: int, c: Rational) -> Rational:
    """
    Upper bound of the roots of p(x) = x^n + c. I.e. a rational
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
    if n == 0:
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

def lower_bound(p: Poly) -> int:
    """Compute 1 / poly_root_ub(reciprocal_polynomial)."""
    # Strip trailing zeros (roots at x = 0).
    work = p
    while not work.is_zero and work.nth(0) == 0:
        work = work.quo(Poly(work.gen, work.gen))
    if work.is_zero or sgc(work) == 0:
        return 0

    # Reciprocal polynomial: reversed coefficients.
    rev_poly = Poly(list(reversed(work.all_coeffs())), p.gen)
    if rev_poly.LC() < 0:
        rev_poly = -rev_poly

    ub = poly_root_ub(rev_poly)
    assert ub > 0
    # Rounding down to an int here, as ASV08 requests. Keeping
    # the fractional part is also allowed here, but will not
    # improve the performance, and integers are shorter. With
    # this, some exact integer roots will be recognized.
    return int(1 / ub)

def isolate_positive_roots(poly: Poly, alpha0: int = 16) -> list[tuple[Rational, Rational]]:
    """
    Isolate all real positive roots of a square-free polynomial
    with rational coefficients. Return an ordered disjoined list
    of intervals (a, b), each with a single root. The intervals
    should be interpreted as open (i.e. the boundaries are
    excluded), unless a=b.
    """
    # Assuming no root at zero.
    assert poly.nth(0) != 0
    # Make the leading coefficient positive.
    if poly.LC() < 0:
        poly = -poly
    s = sgc(poly)
    def make_interval(b, d, a, c, p: Poly):
        """
        Build (lo, hi) from Möbius endpoints b/d and a/c. When
        c=0 the interval is unbounded above; cap it with an upper
        bound on p. Endpoints are not ordered (per paper).
        """
        v1 = Rational(b, d)
        if c == 0:
            # Use poly_root_ub for a finite right endpoint instead of
            # infinity, as the original algorithm would do. This is
            # to help with root refinement later. The only concern
            # is to handle the case of upper bound being exactly the
            # root, as we want to keep our intervals open.
            ub1 = int(poly_root_ub(p)) + 1
            if p.eval(ub1) == 0:
                # ub is the exact root — return as point interval
                v = Rational(a * ub1 + b, d)
                return (v, v)
            v2 = Rational(a * ub1 + b, d)
            assert v1 <= v2, f"make_interval c=0: v1={v1} > v2={v2}"
            return (v1, v2)
        v2 = Rational(a, c)
        return (v1, v2) if v1 <= v2 else (v2, v1)

    # --- Step 1 ---
    rootlist: list[tuple[Rational, Rational]] = []
    if s == 0:
        return rootlist
    if s == 1:
        # Use poly_root_ub for a finite right endpoint instead of
        # infinity, as the original algorithm would do. This is
        # to help with root refinement later. The only concern
        # is to handle the case of upper bound being exactly the
        # root, as we want to keep our intervals open.
        ub = int(poly_root_ub(poly)) + 1
        if poly.eval(ub) == 0:
            rootlist.append((ub, ub))
        else:
            rootlist.append((Rational(0), ub))
        return rootlist
    # intervalstack stores: (a, b, c, d, p_poly, s)
    intervalstack: list[tuple[int, int, int, int, Poly, int]] = [(1, 0, 0, 1, poly, s)]
    # --- Step 2 (loop) ---
    while intervalstack:
        a, b, c, d, p, s = intervalstack.pop()
        # --- Step 3: compute lower bound α on positive roots of p ---
        alpha = lower_bound(p)
        # --- Step 4: scale if α is large ---
        if alpha > alpha0:
            p = scale(p, alpha)
            a = alpha * a
            c = alpha * c
            alpha = 1
        # --- Step 5: shift by α ---
        step5_extracted = (
            False  # track if we extracted a root (avoids Step 6 duplication)
        )
        if alpha >= 1:
            p = shift(p, alpha)
            b = alpha * a + b
            d = alpha * c + d
            # Check for exact root (alpha equals a root of p)
            if p.nth(0) == 0:  # p(0) == 0
                rootlist.append((Rational(b, d), Rational(b, d)))
                p = p.quo(Poly(p.gen, p.gen))  # divide by x
                step5_extracted = True
            s = sgc(p)
            if s == 0:
                continue  # go to Step 2
            if s == 1:
                rootlist.append(make_interval(b, d, a, c, p))
                continue  # go to Step 2
        # --- Step 6: split — forward shift p(x+1) and reciprocal transform ---
        p1 = shift(p, 1)
        a1, b1, c1, d1 = a, a + b, c, c + d
        # Track whether forward branch will extract the root at x=1 in p's frame.
        # p1(0) = p(1), so p1.nth(0) == 0 means root at x = 1.
        # We need this BEFORE dividing p1 by x, for coordinating with reciprocal branch.
        # When step5_extracted, the root at x=1 in p's frame was already extracted
        # in Step 5. The forward shift p(x+1) detects the same root at p1(0)=0.
        # Skip to avoid duplicate intervals.
        forward_will_extract = p1.nth(0) == 0 and not step5_extracted
        if forward_will_extract:
            rootlist.append((Rational(b1, d1), Rational(b1, d1)))
            p1 = p1.quo(Poly(p1.gen, p1.gen))
        s1 = sgc(p1)
        # --- Step 6 (reciprocal branch): handle root at x = 1 before transform ---
        # p2(0) = p(1) (sum of all coefficients). When p(1) = 0, the reciprocal
        # transform drops degree, creating spurious roots. Fix: if p(1) = 0,
        # extract the root first (if not already extracted), then transform p/(x-1).
        p_sum = p.eval(1)
        if p_sum == 0:
            if not forward_will_extract:
                # Forward branch didn't extract it, so we do.
                rootval = Rational(a + b, c + d)
                rootlist.append((rootval, rootval))
            # Divide p by (x-1) via division, transform the quotient.
            p_for_recip = p.quo(Poly(p.gen - 1, p.gen))
        else:
            p_for_recip = p
        # Reciprocal transformation: p2(x) = (x+1)^m * p(1/(x+1))
        p2 = reciprocal_transform(p_for_recip)
        if p2.LC() < 0:
            p2 = -p2
        a2, b2, c2, d2 = b, a + b, d, c + d
        # p2 should have no leading zeros after the p(1)=0 fix above
        # --- Step 7: recompute s2 for the reciprocal branch ---
        if p2.is_zero:
            s2 = 0
        else:
            s2 = sgc(p2)
        # --- Step 8: prefer the branch with more sign changes on the left ---
        if s1 < s2:
            a1, b1, c1, d1, p1, s1, a2, b2, c2, d2, p2, s2 = \
                a2, b2, c2, d2, p2, s2, a1, b1, c1, d1, p1, s1
        # --- Step 9: process branch 1 (forward) ---
        if s1 == 0:
            pass  # no roots in this sub-interval
        elif s1 == 1:
            rootlist.append(make_interval(b1, d1, a1, c1, p1))
        else:
            intervalstack.append((a1, b1, c1, d1, p1, s1))
        # --- Step 10: process branch 2 (reciprocal) ---
        if s2 == 0:
            pass
        elif s2 == 1:
            rootlist.append(make_interval(b2, d2, a2, c2, p2))
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
    if p.degree() <= 0:
        return []
    if p.nth(0) == 0:
        # Divide out the root at zero.
        p = p.quo(Poly(p.gen, p.gen))
        assert p.nth(0) != 0 # The input should be square-free.
        roots_zero = [(Rational(0), Rational(0))]
        if p.degree() <= 0:
            return roots_zero
    else:
        roots_zero = []
    roots_neg = [(-hi, -lo) for lo, hi in reversed(isolate_positive_roots(negate(p)))]
    roots_pos = isolate_positive_roots(p)
    return roots_neg + roots_zero + roots_pos

def bisect_roots(
    p: Poly, intervals: list[tuple[Rational, Rational]], n: int
) -> list[tuple[Rational, Rational]]:
    """
    Refine real root intervals via bisection, making at least n
    steps on each side of each root interval.
    """
    # Divide out all exact roots from p. This leaves us with
    # only open intervals.
    for lo, hi in intervals:
        if lo == hi:
            p = p.quo(Poly(p.gen - lo, p.gen))
    result = []
    for lo0, hi0 in intervals:
        if lo0 == hi0:
            # We've looked at this root already.
            result.append((lo0, hi0))
        else:
            n_left = n_right = 0
            lo, hi = lo0, hi0
            sign_lo = sign(p.eval(lo))
            assert sign_lo != 0
            while n_left < n or n_right < n:
                mid = (lo + hi) / 2
                sign_mid = sign(p.eval(mid))
                if sign_mid == 0:
                    # Exact root at mid.
                    lo = hi = mid
                    break
                elif sign_mid != sign_lo:
                    # Root in (lo, mid) => refine from right.
                    hi = mid
                    n_right += 1
                else:
                    # Root in (mid, hi) => refine from left.
                    lo = mid
                    sign_lo = sign_mid
                    n_left += 1
            result.append((lo, hi))
    return result

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
