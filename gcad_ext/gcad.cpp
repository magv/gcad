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
    vector<fmpz_mpoly_struct> result;
    auto cmp = [&](const fmpz_mpoly_struct &a,
                   const fmpz_mpoly_struct &b) -> bool {
        return fmpz_mpoly_cmp(&a, &b, &ctx) > 0;
    };
    set<fmpz_mpoly_struct, decltype(cmp)> result_set(cmp);
    fmpz_mpoly_factor_struct fac;
    fmpz_mpoly_factor_init(&fac, &ctx);
    for (const auto &poly : polys) {
        int ok = fmpz_mpoly_factor(&fac, &poly, &ctx);
        assert(ok == 1);
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
