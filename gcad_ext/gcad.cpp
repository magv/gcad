#include <flint/fmpq_mpoly.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_poly.h>
#include <set>

using std::set;
using std::vector;

struct PolyRoot {
    // Real root of a polynomial.
    int poly_idx; // Index of the polynomial defining the root.
    int root_idx; // Index of the real root, if roots were ordered by value.
    Q value_lo; // Lower bound on the root's value at the sample point.
    Q value_hi; // Upper bound on the root's value at the sample point.
};

#define PolyRoot_NONE PolyRoot{ -1, 0, Q(), Q() }

struct AxisBound {
    // Cell boundaries along one axis.
    Q point; // A sample var value inside the cell.
    PolyRoot cell_lo; // Lower bounding value, or negative infinity.
    PolyRoot cell_hi; // Higher bounding value, or positive infinity.
};

// The following hack is a testament to the inability of Cython
// to work with nested C structs without converting them to
// Python objects. What a waste of my time.
#define cell_lo_poly_idx cell_lo.poly_idx
#define cell_lo_root_idx cell_lo.root_idx
#define cell_lo_value_lo cell_lo.value_lo
#define cell_lo_value_hi cell_lo.value_hi
#define cell_hi_poly_idx cell_hi.poly_idx
#define cell_hi_root_idx cell_hi.root_idx
#define cell_hi_value_lo cell_hi.value_lo
#define cell_hi_value_hi cell_hi.value_hi

using Cell = vector<AxisBound>;

/* Square-free and relatively prime polynomials multiplicatively
 * generating the product of polys. (Definition 3.2 of [S00])
 */
static vector<fmpz_mpoly_struct>
SFRP(const vector<fmpz_mpoly_struct> &polys, const fmpz_mpoly_ctx_struct &ctx)
{
    trace_scope("SFRP({}p)", polys.size());
    vector<fmpz_mpoly_struct> result;
    auto cmp = [&](const fmpz_mpoly_struct &a,
                   const fmpz_mpoly_struct &b) -> bool {
        return fmpz_mpoly_cmp(&a, &b, &ctx) > 0;
    };
    set<fmpz_mpoly_struct, decltype(cmp)> result_set(cmp);
    fmpz_mpoly_factor_struct fac;
    fmpz_mpoly_factor_init(&fac, &ctx);
    size_t npolys = polys.size();
    double total_C = 0;
    for (size_t i = 0; i < npolys; i++) {
        const auto &poly = polys[npolys-1-i];
        double d = fmpz_mpoly_total_degree_si(&poly, &ctx);
        total_C += d*d*d*d;
    }
    double C = 0;
    trace_progress(C, total_C);
    for (size_t i = 0; i < npolys; i++) {
        const auto &poly = polys[npolys-1-i];
        double d = fmpz_mpoly_total_degree_si(&poly, &ctx);
        slong t = fmpz_mpoly_length(&poly, &ctx);
        {
            trace_scope("factor({}t {}d)", t, (slong)d);
            int ok = fmpz_mpoly_factor(&fac, &poly, &ctx);
            assert(ok);
        }
        C += d*d*d*d;
        trace_progress(C, total_C);
        slong len = fmpz_mpoly_factor_length(&fac, &ctx);
        for (slong j = 0; j < len; j++) {
            fmpz_mpoly_struct &f = fac.poly[j];
            if (result_set.count(f) == 0) {
                fmpz_mpoly_struct r;
                fmpz_mpoly_init(&r, &ctx);
                fmpz_mpoly_swap(&r, &f, &ctx);
                result.push_back(r);
                result_set.insert(r);
            }
        }
    }
    fmpz_mpoly_factor_clear(&fac, &ctx);
    sort(result.begin(), result.end(), cmp);
    return result;
}

/* Obtain the leading coefficient of a polynomial with respect
 * to the given variable.
 */
static void
LC(fmpz_mpoly_struct &res,
   const fmpz_mpoly_struct &poly,
   slong var,
   const fmpz_mpoly_ctx_struct &ctx)
{
    slong deg = fmpz_mpoly_degree_si(&poly, var, &ctx);
    assert(deg >= 0);
    ulong d[1] = {(ulong)deg};
    fmpz_mpoly_get_coeff_vars_ui(&res, &poly, &var, &d[0], 1, &ctx);
}

/* The set of leading coefficients, discriminants, and pairwise
 * resultants of the given list of square-free co-prime polynomials,
 * with respect to the given variable (Definition 3.2 of [S00]).
 */
static vector<fmpz_mpoly_struct>
PR(const vector<fmpz_mpoly_struct> &polys,
   slong var,
   fmpz_mpoly_ctx_struct &ctx)
{
    trace_scope("PR({}p)", polys.size());
    vector<fmpz_mpoly_struct> result;
    auto cmp = [&](const fmpz_mpoly_struct &a,
                   const fmpz_mpoly_struct &b) -> bool {
        return fmpz_mpoly_cmp(&a, &b, &ctx) > 0;
    };
    set<fmpz_mpoly_struct, decltype(cmp)> result_set(cmp);
    auto add = [&](fmpz_mpoly_struct &res) {
        if (!fmpz_mpoly_is_fmpz(&res, &ctx)) {
            if (result_set.count(res) == 0) {
                result.push_back(res);
                result_set.insert(res);
                fmpz_mpoly_init(&res, &ctx);
            }
        }
    };
    size_t npolys = polys.size();
    fmpz_mpoly_struct res;
    fmpz_mpoly_init(&res, &ctx);
    // First pass: complexity analysis.
    double total_C = 0;
    for (size_t i = 0; i < npolys; i++) {
        const auto &p1 = polys[npolys-1-i];
        if (fmpz_mpoly_is_fmpz(&p1, &ctx)) continue;
        double d1 = fmpz_mpoly_degree_si(&p1, var, &ctx);
        double t1 = fmpz_mpoly_length(&p1, &ctx);
        if (d1 >= 1) {
            total_C += 8*d1*d1*d1*t1*t1;
        }
        for (size_t j = i + 1; j < npolys; j++) {
            const auto &p2 = polys[npolys-1-j];
            if (fmpz_mpoly_is_fmpz(&p2, &ctx)) continue;
            double d2 = fmpz_mpoly_degree_si(&p2, var, &ctx);
            double t2 = fmpz_mpoly_length(&p2, &ctx);
            total_C += (d1+d2)*(d1+d2)*(d1+d2)*t1*t2;
        }
    }
    // Second pass: actual calculation.
    double C = 0;
    trace_progress(C, total_C);
    for (size_t i = 0; i < npolys; i++) {
        const auto &p1 = polys[npolys-1-i];
        if (fmpz_mpoly_is_fmpz(&p1, &ctx)) continue;
        double d1 = fmpz_mpoly_degree_si(&p1, var, &ctx);
        double t1 = fmpz_mpoly_length(&p1, &ctx);
        LC(res, p1, var, ctx);
        add(res);
        if (d1 >= 1) {
            {
                trace_scope("discriminant({}t {}d)",
                            fmpz_mpoly_length(&p1, &ctx),
                            fmpz_mpoly_total_degree_si(&p1, &ctx));
                int ok = fmpz_mpoly_discriminant(&res, &p1, var, &ctx);
                assert(ok);
                add(res);
            }
            C += 8*d1*d1*d1*t1*t1;
            trace_progress(C, total_C);
        }
        for (size_t j = i + 1; j < npolys; j++) {
            const auto &p2 = polys[npolys-1-j];
            if (fmpz_mpoly_is_fmpz(&p2, &ctx)) continue;
            double d2 = fmpz_mpoly_degree_si(&p2, var, &ctx);
            double t2 = fmpz_mpoly_length(&p2, &ctx);
            {
                trace_scope("resultant({}t {}d, {}t {}d)",
                            fmpz_mpoly_length(&p1, &ctx),
                            fmpz_mpoly_total_degree_si(&p1, &ctx),
                            fmpz_mpoly_length(&p2, &ctx),
                            fmpz_mpoly_total_degree_si(&p2, &ctx));
                int ok = fmpz_mpoly_resultant(&res, &p1, &p2, var, &ctx);
                assert(ok);
                add(res);
            }
            C += (d1+d2)*(d1+d2)*(d1+d2)*t1*t2;
            trace_progress(C, total_C);
        }
    }
    fmpz_mpoly_clear(&res, &ctx);
    sort(result.begin(), result.end(), cmp);
    return result;
}

/* SFRP(PR(polys, var))
 */
static vector<fmpz_mpoly_struct>
SFRP_PR(const vector<fmpz_mpoly_struct> &polys,
        slong var,
        fmpz_mpoly_ctx_struct &ctx)
{
    vector<fmpz_mpoly_struct> pr = PR(polys, var, ctx);
    vector<fmpz_mpoly_struct> sfrp = SFRP(pr, ctx);
    for (auto &mp : pr) {
        fmpz_mpoly_clear(&mp, &ctx);
    }
    return sfrp;
}

/* The recursive step of RSFC.
 */
static void
_RSFC(const vector<fmpq_mpoly_struct> &positives,
      const vector<vector<fmpq_mpoly_struct>> &pr,
      int k,
      int nvars,
      Cell &cell,
      vector<Cell> &cells,
      size_t max_cells,
      size_t &n_rejected,
      vector<size_t> &n_early_exits,
      fmpq_mpoly_ctx_struct &ctx)
{
    if (k >= nvars) {
        for (const auto &p : positives) {
            assert(fmpq_mpoly_is_fmpq(&p, &ctx));
            Q c;
            fmpq_mpoly_get_fmpq(&c.q, &p, &ctx);
            if (c.sign() < 0) {
                n_rejected++;
                return;
            }
        }
        cells.push_back(cell);
        return;
    }
    // Early exit check.
    for (const auto &p : positives) {
        if (fmpq_mpoly_is_fmpq(&p, &ctx)) {
            Q c;
            fmpq_mpoly_get_fmpq(&c.q, &p, &ctx);
            if (c.sign() < 0) {
                n_early_exits[k]++;
                return;
            }
        }
    }
    // All lower variables are already substituted into pr[k]; we
    // just need to convert it to univariate for root finding.
    vector<ZPoly> pr_k;
    for (const auto &mp0 : pr[k]) {
        ZPoly zp;
        fmpz_mpoly_get_fmpz_poly(
            &zp.p,
            fmpq_mpoly_zpoly_ref((fmpq_mpoly_struct *)&mp0, &ctx),
            k,
            ctx.zctx);
        pr_k.push_back(std::move(zp));
    }
    vector<vector<RootInterval>> intervals = isolate_many_roots(pr_k);
    vector<PolyRoot> roots;
    for (int pidx = 0; pidx < (int)intervals.size(); pidx++) {
        for (int ridx = 0; ridx < (int)intervals[pidx].size(); ridx++) {
            roots.push_back({pidx,
                             ridx,
                             Q_of_fmpq(&intervals[pidx][ridx].lo.q),
                             Q_of_fmpq(&intervals[pidx][ridx].hi.q)});
        }
    }
    sort(roots.begin(), roots.end(), [](const PolyRoot &a, const PolyRoot &b) {
        return a.value_lo < b.value_lo;
    });
    // Now we can recurse into every interval between roots.
    auto recurse_at = [&](const fmpq_t val) {
        vector<vector<fmpq_mpoly_struct>> next_pr(nvars);
        for (int l = k + 1; l < nvars; l++) {
            for (const auto &p : pr[l]) {
                fmpq_mpoly_struct np;
                fmpq_mpoly_init(&np, &ctx);
                fmpq_mpoly_evaluate_one_fmpq(&np, &p, k, val, &ctx);
                next_pr[l].push_back(np);
            }
        }
        vector<fmpq_mpoly_struct> next_positives;
        for (auto &p : positives) {
            fmpq_mpoly_struct np;
            fmpq_mpoly_init(&np, &ctx);
            fmpq_mpoly_evaluate_one_fmpq(&np, &p, k, val, &ctx);
            next_positives.push_back(np);
        }
        _RSFC(next_positives,
              next_pr,
              k + 1,
              nvars,
              cell,
              cells,
              max_cells,
              n_rejected,
              n_early_exits,
              ctx);
        for (auto &p : next_positives) {
            fmpq_mpoly_clear(&p, &ctx);
        }
        for (int l = k + 1; l < nvars; l++) {
            for (auto &p : next_pr[l]) fmpq_mpoly_clear(&p, &ctx);
        }
    };
    if (!roots.empty()) {
        {
            const auto &hi = roots[0];
            Q pt = (hi.value_lo.truncate() - Z_of_si(1)) / Z_of_si(1);
            cell.push_back(AxisBound{pt, PolyRoot_NONE, hi});
            recurse_at(&pt.q);
            cell.pop_back();
            if (max_cells && (cells.size() >= max_cells)) return;
        }
        for (size_t i = 0; i + 1 < roots.size(); i++) {
            const auto &lo = roots[i];
            const auto &hi = roots[i + 1];
            Q a = lo.value_hi;
            Q b = hi.value_lo;
            Q a_adj = (lo.value_lo != lo.value_hi)
                          ? a
                          : (a + (b - a) / Q_of_si(1024));
            Q b_adj = (hi.value_lo != hi.value_hi)
                          ? b
                          : (b - (b - a) / Q_of_si(1024));
            fmpq qa, qb, qres;
            fmpq_init(&qa);
            fmpq_set(&qa, &a_adj.q);
            fmpq_init(&qb);
            fmpq_set(&qb, &b_adj.q);
            fmpq_init(&qres);
            shortest_fraction_between(&qres, &qa, &qb);
            Q mid;
            fmpq_init(&mid.q);
            fmpq_set(&mid.q, &qres);
            fmpq_clear(&qa);
            fmpq_clear(&qb);
            fmpq_clear(&qres);
            cell.push_back(AxisBound{mid, lo, hi});
            recurse_at(&mid.q);
            cell.pop_back();
            if (max_cells && (cells.size() >= max_cells)) return;
        }
        {
            const auto &lo = roots.back();
            Q pt = (lo.value_hi.truncate() + Z_of_si(1)) / Z_of_si(1);
            cell.push_back(AxisBound{pt, lo, PolyRoot_NONE});
            recurse_at(&pt.q);
            cell.pop_back();
        }
    } else {
        Q pt = Q_of_si(0);
        cell.push_back(AxisBound{pt, PolyRoot_NONE, PolyRoot_NONE});
        recurse_at(&pt.q);
        cell.pop_back();
    }
}

/* Recursive Solution Formula Construction (Algorithm 3.5 of [S00]).
 */
vector<Cell>
RSFC(const vector<fmpz_mpoly_struct> &positives,
     const vector<vector<fmpz_mpoly_struct>> &pr_fmpz,
     int nvars,
     size_t max_cells,
     fmpz_mpoly_ctx_struct &zctx)
{
    log_trace_scope("RSFC({}p, {}v)", positives.size(), nvars);
    vector<Cell> cells;
    size_t n_rejected = 0;
    vector<size_t> n_early_exits(nvars, 0);
    fmpq_mpoly_ctx_struct ctx = {{zctx}};
    // Convert to fmpq_mpoly, because we'll need FLINT functions
    // for variable substitution available only on fmpq_mpoly...
    vector<fmpq_mpoly_struct> qpositives;
    for (const auto &p : positives) {
        fmpq_mpoly_struct qp;
        fmpq_mpoly_init(&qp, &ctx);
        fmpq_set_si(fmpq_mpoly_content_ref(&qp, &ctx), 1, 1);
        fmpz_mpoly_set(fmpq_mpoly_zpoly_ref(&qp, &ctx), &p, &zctx);
        qpositives.push_back(qp);
    }
    vector<vector<fmpq_mpoly_struct>> qpr;
    for (int k = 0; k < nvars; k++) {
        vector<fmpq_mpoly_struct> pr_k;
        assert((size_t)k < pr_fmpz.size());
        for (auto &p : pr_fmpz[k]) {
            fmpq_mpoly_struct qp;
            fmpq_mpoly_init(&qp, &ctx);
            fmpq_set_si(fmpq_mpoly_content_ref(&qp, &ctx), 1, 1);
            fmpz_mpoly_set(fmpq_mpoly_zpoly_ref(&qp, &ctx), &p, &zctx);
            pr_k.push_back(qp);
        }
        qpr.push_back(std::move(pr_k));
    }
    Cell cell;
    _RSFC(
        qpositives, qpr, 0, nvars, cell, cells, max_cells, n_rejected, n_early_exits, ctx);
    for (auto &qp : qpositives) {
        fmpq_mpoly_clear(&qp, &ctx);
    }
    for (auto &pr_k : qpr) {
        for (auto &p : pr_k) {
            fmpq_mpoly_clear(&p, &ctx);
        }
    }
    fmpq_mpoly_ctx_clear(&ctx);
    log("Accepted {} cells, rejected {}", cells.size(), n_rejected);
    log("Early exits: {}", n_early_exits);
    return cells;
}
