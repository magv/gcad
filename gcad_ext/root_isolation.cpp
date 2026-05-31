/*
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
*/

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <vector>

#include "flint_tools.cpp"

/* The number of sign changes in the sequence of nonzero
 * coefficients of p.
 */
static int
sgc(const ZPoly &p)
{
    slong n = p.length();
    if (n <= 1) return 0;
    int changes = 0;
    int prev = p[0].sign();
    for (slong i = 1; i < n; i++) {
        int s = p[i].sign();
        if ((s != 0) && (s != prev)) {
            changes++;
            prev = s;
        }
    }
    return changes;
}

/* Find g(x) = p(alpha * x).
 */
static void
scale(ZPoly &p, Z &alpha)
{
    slong n = p.length();
    if (n > 1) {
        Z k = alpha;
        for (slong i = 1; i < n-1; i++) {
            p[i] *= k;
            k *= alpha;
        }
        p[n-1] *= k;
    }
}

/* Find g(x) = (x+1)^{deg(p)} * p(1/(x+1)).
 */
static void
reciprocal_transform(ZPoly &p)
{
    slong n = p.length();
    p.set_reverse(n);
    p.set_taylor_shift(Z_of_si(1));
}

/* Find g(x) = p(-x).
 */
static void
negate(ZPoly &p)
{
    slong n = p.length();
    for (slong i = 1; i < n; i += 2) {
        p[i].set_neg();
    }
}

/* Upper bound of the roots of p(x) = x^n - c. I.e. a rational
 * that is larger than, or equal to, c^{1/n}.
 */
static Q
nth_root_ub(ulong n, const Q &c)
{
    if ((n == 1) || c.is_one()) {
        return c;
    } else {
        float approx = powf(c.get_d(), 1.0f / n);
        int e = 0;
        slong m = (slong)ceilf(1024.0 * frexpf(approx, &e));
        e -= 10;
        for (;;) {
            Q y = Q_ldexp_si(m, e);
            if (c <= y.pow_si(n)) return y;
            m++;
        }
    }
}

/* A possibly tight upper bound on the real positive roots of
 * p, according to the LMQ algorithm from AAS08 (Algorithm 3.1,
 * the "Local-Max" Quadratic implementation of Theorem 3). Also
 * benchmarked and recommended in ASV08.
 */
static Q
poly_root_ub(const ZPoly &p)
{
    slong n = p.length();
    if (n <= 1) return Q_of_si(0);
    std::vector<slong> used(n);
    for (slong i = 0; i < n; i++)
        used[i] = 1;
    Q tmax = Q_of_si(0);
    for (slong mu = n - 2; mu >= 0; mu--) {
        const Z &p_mu = p[mu];
        if (p_mu.is_negative()) {
            bool have_tmin = false;
            Q tmin;
            for (slong nu = n - 1; nu > mu; nu--) {
                const Z &p_nu = p[nu];
                if (p_nu.is_positive()) {
                    Q c = p_mu / p_nu;
                    c.set_neg();
                    c.set_mul_2exp(used[nu]);
                    Q t = nth_root_ub(nu - mu, c);
                    used[nu]++;
                    if ((!have_tmin) || (t < tmin)) {
                        tmin = t;
                        have_tmin = true;
                    }
                }
            }
            if (have_tmin && (tmin > tmax)) tmax = tmin;
        }
    }
    return tmax;
}

/* A possibly tight lower bound on the positive roots of the
 * given polynomial, computed as 1/poly_root_ub(x^n p(1/x)).
 */
static Q
poly_root_lb(const ZPoly &p)
{
    slong n = p.length();
    ZPoly rev = p;
    rev.set_reverse(n);
    if (rev.lead().is_negative()) {
        rev.set_neg();
    }
    return poly_root_ub(rev).inv();
}

struct RootInterval {
    Q lo, hi;
};

/* Return the interval (min(b/d, a/c), max(b/d, a/c)). When
 * c=0, the interval is (b/d, ∞), but instead of ∞, we use
 * an upper bound on the roots of p.
 */
static RootInterval
_intrv(const Z &a, const Z &b, const Z &c, const Z &d, const ZPoly &p)
{
    Q v1 = b / d;
    if (c.is_zero()) {
        // Use poly_root_ub for a finite right endpoint instead of
        // infinity, as the original algorithm would do. This is
        // to help with root refinement later. The only concern
        // is to handle the case of upper bound being exactly the
        // root, as we want to keep our intervals open.
        Z ub = poly_root_ub(p).truncate() + Z_of_si(1);
        Q v2 = (a * ub + b) / d;
        assert(v1 < v2);
        return {v1, v2};
    } else {
        Q v2 = a / c;
        return (v1 <= v2) ? (RootInterval){v1, v2} : (RootInterval){v2, v1};
    }
}

/* Isolate all real positive roots of a square-free polynomial
 * with rational coefficients. Return an ordered disjoined list
 * of intervals (a, b), each with a single root. The intervals
 * should be interpreted as open (i.e. the boundaries are
 * excluded), unless a=b. We assume no roots at zero, and
 * a positive sign of the leading coefficient.
 */
static std::vector<RootInterval>
isolate_positive_roots(const ZPoly &poly, slong alpha0 = 16)
{
    // Assuming no roots at zero, degree > 1, and positive leading
    // coefficient.
    assert(poly.length() > 2);
    assert(!poly[0].is_zero());
    assert(poly.lead().is_positive());
    std::vector<RootInterval> rootlist;
    int s = sgc(poly);
    if (s == 0) {
        return rootlist;
    }
    if (s == 1) {
        // Use poly_root_ub for a finite right endpoint instead of
        // infinity, as the original algorithm would do. This is
        // to help with root refinement later. The only concern
        // is to handle the case of upper bound being exactly the
        // root, as we want to keep our intervals open.
        Q lo = Q_of_si(0);
        Q hi = (poly_root_ub(poly).truncate() + Z_of_si(1)) / Z_of_si(1);
        rootlist.push_back({lo, hi});
        return rootlist;
    }
    std::vector<std::tuple<Z, Z, Z, Z, ZPoly, int>> intervalstack;
    intervalstack.push_back(
        {Z_of_si(1), Z_of_si(0), Z_of_si(0), Z_of_si(1), poly, s});
    while (!intervalstack.empty()) {
        Z a, b, c, d;
        ZPoly p;
        int s;
        std::tie(a, b, c, d, p, s) = intervalstack.back();
        intervalstack.pop_back();
        // Rounding down the lower bound to an int, as ASV08
        // requests. Keeping the fractional part is also allowed
        // here, but will not improve the performance, and integers
        // are shorter. With this, some exact integer roots will
        // be recognized.
        Z alpha = poly_root_lb(p).truncate();
        // Rescale if lower bound is large.
        if (alpha > Z_of_si(alpha0)) {
            scale(p, alpha);
            a *= alpha;
            c *= alpha;
            alpha = Z_of_si(1);
        }
        if (alpha >= Z_of_si(1)) {
            // Shift by the lower bound to be closer to the root.
            p.set_taylor_shift(alpha);
            b = a * alpha + b;
            d = c * alpha + d;
            if (p[0].is_zero()) {
                // Exact root at zero. Divide it out.
                rootlist.push_back({b / d, b / d});
                p.set_shift_right(1);
                // Dividing only once, because p is supposed to
                // be square-free.
                assert(!p[0].is_zero());
            }
            s = sgc(p);
            if (s == 0) {
                continue;
            }
            if (s == 1) {
                rootlist.push_back(_intrv(a, b, c, d, p));
                continue;
            }
        }
        // The root should be close. Let's split the interval at
        // x=1, and look at both sides:
        // - p1(x) = p(x+1) will cover x ∈ (1; ∞);
        // - p2(x) = (x+1)^n p(1/(x+1)) will cover x ∈ (0; 1).
        ZPoly p1 = p;
        p1.set_taylor_shift(Z_of_si(1));
        Z a1 = a;
        Z b1 = a + b;
        Z c1 = c;
        Z d1 = c + d;
        int r;
        if (p1[0].is_zero()) {
            rootlist.push_back({b1 / d1, b1 / d1});
            p1.set_shift_right(1);
            r = 1;
        } else {
            r = 0;
        }
        int s1 = sgc(p1);
        // The whole s2_bound business is tricky, and only valid
        // conjecturally. AS05 mentions the conjecture that
        // sgc(p1)+sgc(p2) <= sgc(p). Here, however, we assume
        // more: we assume that if s2_bound <= 0, then sgc(p2)
        // is 0, and that if s2_bound is 1, then sgc(p2) is 1.
        int s2_bound = s - s1 - r;
        Z a2 = b;
        Z b2 = a + b;
        Z c2 = d;
        Z d2 = c + d;
        int s2;
        ZPoly p2;
        // By the conjecture above, we can skip the computation
        // of p2 in some cases.
        if (s2_bound > 1) {
            p2 = p;
            reciprocal_transform(p2);
            if (p2.lead().is_negative()) {
                p2.set_neg();
            }
            if (p2[0].is_zero()) {
                // p2(0) = p1(0) = p(1) = 0. If we are here,
                // this root was already counted.
                p2.set_shift_right(1);
            }
            s2 = sgc(p2);
        } else {
            // The value of p2 will not be used.
            s2 = s2_bound;
        }
        // To keep the interval stack small, it's best to push
        // intervals with higher sgc first.
        if (s1 < s2) {
            std::swap(a1, a2);
            std::swap(b1, b2);
            std::swap(c1, c2);
            std::swap(d1, d2);
            std::swap(s1, s2);
            std::swap(p1, p2);
        }
        if (s1 == 0) {
            continue;
        } else if (s1 == 1) {
            rootlist.push_back(_intrv(a1, b1, c1, d1, p1));
        } else {
            intervalstack.push_back({a1, b1, c1, d1, std::move(p1), s1});
        }
        if (s2 == 0) {
            continue;
        } else if (s2 == 1) {
            rootlist.push_back(_intrv(a2, b2, c2, d2, p2));
        } else {
            intervalstack.push_back({a2, b2, c2, d2, std::move(p2), s2});
        }
    }
    std::sort(rootlist.begin(),
              rootlist.end(),
              [](const RootInterval &a, const RootInterval &b) {
                  return a.lo < b.lo;
              });
    return rootlist;
}

/* Isolate all real roots of a square-free polynomial with
 * rational coefficients. Return an ordered disjoined list of
 * intervals (a, b), each with a single root. The intervals
 * should be interpreted as open (i.e. the boundaries are
 * excluded), unless a=b.
 */
static std::vector<RootInterval>
isolate_roots(ZPoly p)
{
    std::vector<RootInterval> roots;
    if (p.length() <= 1)
        return roots;
    if (p[0].is_zero()) {
        // Divide out the root at zero.
        p.set_shift_right(1);
        roots.push_back({Q_of_si(0), Q_of_si(0)});
    }
    if (p.length() <= 1)
        return roots;
    // The input should be square-free.
    assert(!p[0].is_zero());
    // Special-case linear polys.
    if (p.length() == 2) {
        Q v = p[0] / p[1];
        v.set_neg();
        roots.push_back({v, v});
        return roots;
    }
    // Negative roots.
    ZPoly neg_p = p;
    negate(neg_p);
    if (neg_p.lead().is_negative())
        neg_p.set_neg();
    auto neg_roots = isolate_positive_roots(neg_p);
    roots.insert(roots.end(),
                 std::make_move_iterator(neg_roots.begin()),
                 std::make_move_iterator(neg_roots.end()));
    neg_roots.clear();
    for (auto &r : roots) {
        r.lo.set_neg();
        r.hi.set_neg();
        std::swap(r.lo, r.hi);
    }
    std::reverse(roots.begin(), roots.end());
    // Positive roots.
    if (p.lead().is_negative())
        p.set_neg();
    auto pos_roots = isolate_positive_roots(p);
    roots.insert(roots.end(),
                 std::make_move_iterator(pos_roots.begin()),
                 std::make_move_iterator(pos_roots.end()));
    pos_roots.clear();
    return roots;
}

/* Refine a real root interval for a given polynomial once.
 */
static RootInterval
bisect_root(const ZPoly &p, const Q &lo, const Q &hi)
{
    assert(lo != hi);
    int sign_lo = p.eval(lo).sign();
    Q mid = (lo + hi) / Q_of_si(2);
    int sign_mid = p.eval(mid).sign();
    if (sign_mid == 0) {
        // Exact root at mid.
        return {mid, mid};
    } else if (sign_mid != sign_lo) {
        // Root in (lo, mid) => refine from right.
        return {lo, mid};
    } else {
        // Root in (mid, hi) => refine from left.
        return {mid, hi};
    }
}

struct IndexedRootInterval {
    RootInterval interval;
    int idx;
};

static bool
operator>(const IndexedRootInterval &a, const IndexedRootInterval &b)
{
    return a.interval.lo > b.interval.lo;
}

/* Isolate all real roots of a list of square-free co-prime
 * polynomials with rational coefficients. For each polynomial,
 * return an ordered disjoined list of intervals (a, b), each
 * with a single root. The intervals belonging to different
 * polynomials are guaranteed to neither overlap, nor touch at
 * the boundaries. The intervals are open (i.e. the boundaries
 * are excluded), unless a=b.
 */
static std::vector<std::vector<RootInterval>>
isolate_many_roots(std::vector<ZPoly> &polys)
{
    // Isolate roots of each poly separately.
    std::vector<IndexedRootInterval> intervals;
    for (int idx = 0; idx < (int)polys.size(); idx++) {
        ZPoly &p = polys[idx];
        for (auto &ri : isolate_roots(p)) {
            if (ri.lo == ri.hi) {
                // We'll divide out exact roots immediately, these
                // don't need bisection. They still need to be
                // accounted for, to make sure other intervals
                // don't overlap with them.
                p.set_divexact_root(ri.lo);
            }
            intervals.push_back({ri.lo, ri.hi, idx});
        }
    }
    // Sort by the lower boundary.
    std::make_heap(intervals.begin(), intervals.end(), std::greater<>());
    // Bisect overlapping intervals till none are.
    std::vector<std::vector<RootInterval>> result(polys.size());
    if (intervals.empty())
        return result;
    std::pop_heap(intervals.begin(), intervals.end(), std::greater<>());
    IndexedRootInterval i1 = std::move(intervals.back());
    intervals.pop_back();
    if (intervals.empty()) {
        result[i1.idx].push_back(i1.interval);
        return result;
    }
    std::pop_heap(intervals.begin(), intervals.end(), std::greater<>());
    IndexedRootInterval i2 = std::move(intervals.back());
    intervals.pop_back();
    for (;;) {
        if (i1.interval.hi < i2.interval.lo) {
            // No intersection. Interval 1 is done.
            if (!intervals.empty()) {
                result[i1.idx].push_back(i1.interval);
                i1 = i2;
                std::pop_heap(intervals.begin(), intervals.end(), std::greater<>());
                i2 = std::move(intervals.back());
                intervals.pop_back();
            } else {
                result[i1.idx].push_back(i1.interval);
                result[i2.idx].push_back(i2.interval);
                return result;
            }
        } else {
            // Intersecting intervals. Will bisect and retry.
            while (!((i1.interval.hi < i2.interval.lo) ||
                     (i2.interval.hi < i1.interval.lo))) {
                if (i1.interval.lo != i1.interval.hi) {
                    i1.interval = bisect_root(
                        polys[i1.idx], i1.interval.lo, i1.interval.hi);
                }
                if (i2.interval.lo != i2.interval.hi) {
                    i2.interval = bisect_root(
                        polys[i2.idx], i2.interval.lo, i2.interval.hi);
                }
            }
            if (i2.interval.hi < i1.interval.lo) {
                // Interval 2 < interval 1.
                std::swap(i1, i2);
            }
            intervals.push_back(std::move(i1));
            std::push_heap(intervals.begin(), intervals.end(), std::greater<>());
            std::pop_heap(intervals.begin(), intervals.end(), std::greater<>());
            i1 = std::move(intervals.back());
            intervals.pop_back();
            intervals.push_back(std::move(i2));
            std::push_heap(intervals.begin(), intervals.end(), std::greater<>());
            std::pop_heap(intervals.begin(), intervals.end(), std::greater<>());
            i2 = std::move(intervals.back());
            intervals.pop_back();
        }
    }
}
