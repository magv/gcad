#include <set>

struct PolyRoot {
    // Real root of a polynomial.
    int poly_idx; // Index of the polynomial defining the root.
    int root_idx; // Index of the real root, if roots were ordered by value.
    Q value_lo; // Lower bound on the root's value at the sample point.
    Q value_hi; // Upper bound on the root's value at the sample point.
};

struct AxisBound {
    // Cell boundaries along one axis.
    int var; // Variable index corresponding to this axis.
    Q point; // A sample var value inside the cell.
    PolyRoot cell_lo; // Lower bounding value, or negative infinity.
    PolyRoot cell_hi; // Higher bounding value, or positive infinity.
};

using std::set;
using std::vector;
using Cell = vector<AxisBound>;

/* Square-free and relatively prime polynomials multiplicatively
 * generating the product of polys. (Definition 3.2 of [S00])
 */
static vector<fmpz_mpoly_struct>
SFRP(const vector<fmpz_mpoly_struct> &polys, const fmpz_mpoly_ctx_struct &ctx)
{
    LOGBLOCK("SFRP({}p)", polys.size());
    vector<fmpz_mpoly_struct> result;
    auto cmp = [&](const fmpz_mpoly_struct &a,
                   const fmpz_mpoly_struct &b) -> bool {
        return fmpz_mpoly_cmp(&a, &b, &ctx) > 0;
    };
    set<fmpz_mpoly_struct, decltype(cmp)> result_set(cmp);
    fmpz_mpoly_factor_struct fac;
    fmpz_mpoly_factor_init(&fac, &ctx);
    for (const auto &poly : polys) {
        {
            LOGBLOCK("{}/{}: factor({}t {}d)",
                     &poly - &polys[0] + 1,
                     polys.size(),
                     fmpz_mpoly_length(&poly, &ctx),
                     fmpz_mpoly_total_degree_si(&poly, &ctx));
            int ok = fmpz_mpoly_factor(&fac, &poly, &ctx);
            assert(ok);
        }
        slong len = fmpz_mpoly_factor_length(&fac, &ctx);
        for (slong i = 0; i < len; i++) {
            fmpz_mpoly_struct &f = fac.poly[i];
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
    LOGBLOCK("PR({}p)", polys.size());
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
    for (size_t i = 0, idx = 1; i < npolys; i++, idx++) {
        const auto &p1 = polys[i];
        if (fmpz_mpoly_is_fmpz(&p1, &ctx)) continue;
        LC(res, p1, var, ctx);
        add(res);
        if (fmpz_mpoly_degree_si(&p1, var, &ctx) >= 1) {
            LOGBLOCK("{}/{}: discriminant({}t {}d)",
                     idx,
                     npolys*(npolys+1)/2,
                     fmpz_mpoly_length(&p1, &ctx),
                     fmpz_mpoly_total_degree_si(&p1, &ctx));
            int ok = fmpz_mpoly_discriminant(&res, &p1, var, &ctx);
            assert(ok);
            add(res);
        }
        for (size_t j = i + 1; j < npolys; j++) {
            idx++;
            const auto &p2 = polys[j];
            if (fmpz_mpoly_is_fmpz(&p2, &ctx)) continue;
            LOGBLOCK("{}/{}: resultant({}t {}d, {}t {}d)",
                     idx,
                     npolys*(npolys+1)/2,
                     fmpz_mpoly_length(&p1, &ctx),
                     fmpz_mpoly_total_degree_si(&p1, &ctx),
                     fmpz_mpoly_length(&p2, &ctx),
                     fmpz_mpoly_total_degree_si(&p2, &ctx));
            int ok = fmpz_mpoly_resultant(&res, &p1, &p2, var, &ctx);
            assert(ok);
            add(res);
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
