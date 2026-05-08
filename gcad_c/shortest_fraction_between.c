#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <stdlib.h>

// Find the shortest fraction between (or equal to) two rationals.
static void
shortest_fraction_between(fmpq *res, const fmpq *x1, const fmpq *x2)
{
    // SF[x1_, x2_] /; x2 < 0 := -SF[-x2, -x1]
    if (fmpq_sgn(x2) < 0) {
        fmpq_t lo, hi;
        fmpq_init(lo);
        fmpq_init(hi);
        fmpq_neg(lo, x1);
        fmpq_neg(hi, x2);
        shortest_fraction_between(res, hi, lo);
        fmpq_neg(res, res);
        fmpq_clear(lo);
        fmpq_clear(hi);
        return;
    }
    // SF[x1_, x2_] /; x1 <= 0 <= x2 := 0
    if (fmpq_sgn(x1) <= 0) {
        fmpq_zero(res);
        return;
    }
    fmpq_t lo, hi;
    fmpz_t P0, P1, Q0, Q1, c, tmp;
    // lo = x1;
    fmpq_init(lo); fmpq_set(lo, x1);
    // hi = x2;
    fmpq_init(hi); fmpq_set(hi, x2);
    // P0 = 0;
    fmpz_init_set_ui(P0, 0);
    // P1 = 1;
    fmpz_init_set_ui(P1, 1);
    // Q0 = 1;
    fmpz_init_set_ui(Q0, 1);
    // Q1 = 0;
    fmpz_init_set_ui(Q1, 0);
    fmpz_init(c);
    fmpz_init(tmp);
    for (;;) {
        // c = Ceiling[lo];
        fmpz_cdiv_q(c, fmpq_numref(lo), fmpq_denref(lo));
        // If[c <= hi, Break[]];
        fmpz_mul(tmp, c, fmpq_denref(hi));
        if (fmpz_cmp(tmp, fmpq_numref(hi)) <= 0) break;
        // c = c - 1;
        fmpz_sub_ui(c, c, 1);
        // {lo, hi} = {lo - c, hi - c};
        fmpz_submul(fmpq_numref(lo), c, fmpq_denref(lo));
        fmpz_submul(fmpq_numref(hi), c, fmpq_denref(hi));
        // {lo, hi} = {1/hi, 1/lo};
        fmpz_swap(fmpq_numref(lo), fmpq_denref(lo));
        fmpz_swap(fmpq_numref(hi), fmpq_denref(hi));
        fmpq_swap(lo, hi);
        // {P1, P0} = {c * P1 + P0, P1};
        fmpz_swap(P0, P1);
        fmpz_addmul(P1, c, P0);
        // {Q1, Q0} = {c * Q1 + Q0, Q1};
        fmpz_swap(Q0, Q1);
        fmpz_addmul(Q1, c, Q0);
    }
    // {P1, P0} = {c * P1 + P0, P1};
    fmpz_swap(P0, P1);
    fmpz_addmul(P1, c, P0);
    // {Q1, Q0} = {c * Q1 + Q0, Q1};
    fmpz_swap(Q0, Q1);
    fmpz_addmul(Q1, c, Q0);
    // P1/Q1
    fmpz_swap(fmpq_numref(res), P1);
    fmpz_swap(fmpq_denref(res), Q1);
    fmpq_clear(lo);
    fmpq_clear(hi);
    fmpz_clear(P0);
    fmpz_clear(P1);
    fmpz_clear(Q0);
    fmpz_clear(Q1);
    fmpz_clear(c);
    fmpz_clear(tmp);
}
