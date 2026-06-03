"""
Generic Cylindrical Algebraic Decomposition (GCAD),
as described in [S00].

S00:
    A. Strzeboński.
    "Solving Systems of Strict Polynomial Inequalities".
    J. Symb. Comput. 29, 471--480 (2000).
    https://doi.org/10.1006/jsco.1999.0327

DSS04:
    A. Dolzmann, A. Seidl, T. Sturm.
    "Efficient Projection Orders for CAD".
    ISSAC '04.
    https://doi.org/10.1145/1005285.1005303
"""

from __future__ import annotations

from dataclasses import dataclass
from gcad_ext import (
    discriminant,
    factor,
    isolate_many_roots,
    resultant,
    shortest_fraction_between,
)
#from gcad.root_isolation import isolate_many_roots
import sympy as sp
from .log import *

@dataclass(slots=True)
class PolyRoot:
    poly: sp.Poly
    idx: int
    # Isolating interval for the root value at the sample point.
    value_lo: sp.Rational
    value_hi: sp.Rational
    def __repr__(self):
        p = self.poly.as_expr().subs({self.poly.gen: sp.Symbol("#")})
        p = str(p).replace(" ", "").replace("**", "^")
        return f"Root[{p}, {self.idx}]"

@dataclass(slots=True)
class AxisBound:
    var: sp.Symbol
    point: sp.Rational # A sample var value inside the cell.
    cell_lo: PolyRoot | None # None means negative infinity
    cell_hi: PolyRoot | None # None means positive infinity
    def __repr__(self):
        lo = "-∞" if self.cell_lo is None else str(self.cell_lo)
        hi = "∞" if self.cell_hi is None else str(self.cell_hi)
        return f"{lo} < {self.var} < {hi}"

Cell = list[AxisBound]

def uniq(seq: list) -> list:
    """Remove duplicates, preserve order."""
    seen = set()
    result = []
    for item in seq:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result

def SFRP(polys: list[sp.Expr], variables: list[sp.Symbol]) -> list[sp.Poly]:
    """
    Square-free and relatively prime polynomials multiplicatively
    generating the product of polys, as a poly in the given
    variable. (Definition 3.2 of [S00])
    """
    assert isinstance(polys, list)
    for p in polys:
        assert isinstance(p, sp.Poly)
        assert p.gens == tuple(variables)
    result = []
    with logblock("factor"):
        ff = [factor(p) for p in polys]
    for content, factors in ff:
        for poly, exp in factors:
            assert isinstance(poly, sp.Poly)
            assert p.domain == sp.ZZ
            assert poly.gens == tuple(variables)
            assert poly.LC().is_positive
            result.append(poly)
    return uniq(result)

def PR(polys: list[sp.Poly], var: sp.Symbol, rest: list[sp.Symbol]) -> list[sp.Expr]:
    """
    The set of leading coefficients, discriminants, and pairwise
    resultants of the given list of square-free co-prime polynomials,
    with respect to the given variable (Definition 3.2 of [S00]).
    """
    result = []
    for p in polys:
        if len(rest) > 0:
            lc = sp.Poly(sp.LC(p, var), *rest)
            assert lc != 0
            assert isinstance(lc, sp.Poly)
            result.append(lc)
        with logblock(f"discriminant({p.length()}t {p.total_degree()}d)"):
            disc = discriminant(p, p.gens.index(var))
            disc = sp.Poly(disc, *rest)
        if disc != 0:
            assert isinstance(disc, sp.Poly)
            result.append(disc)
    n = len(polys)
    for i in range(n):
        pi = polys[i]
        for j in range(i + 1, n):
            pj = polys[j]
            with logblock(f"resultant({pi.length()}t {pi.total_degree()}d, {pj.length()}t {pj.total_degree()}d)"):
                r = resultant(pi, pj, pi.gens.index(var))
                r = sp.Poly(r, *rest)
            if r != 0:
                assert isinstance(r, sp.Poly)
                result.append(r)
    return uniq(result)

@autolog
def GPROJ(positives: list[sp.Poly], varlist: list[sp.Symbol]) -> list[list[sp.Poly]]:
    """
    Generic projection (Algorithm 3.4 of [S00]). Return the list
    of $pr_i$, as polynomials with integer coefficients.
    """
    n = len(varlist)
    # Note: our lists are 0-indexed, unlike the paper.
    pr: list[list[sp.Poly]] = [[] for k in range(n)]
    with logblock(f"{varlist[n-1]}"):
        pr[n - 1] = SFRP([sp.Poly(p, *varlist) for p in positives], varlist)
    for k in reversed(range(0, n - 1)):
        with logblock(f"{varlist[k]}"):
            pr[k] = SFRP(PR(pr[k + 1], varlist[k + 1], varlist[: k + 1]), varlist[: k + 1])
    # Filter out constants polys from the projection: these never
    # have roots, so we can safely skip them.
    pr = [[p for p in polys if p.degree(varlist[k]) >= 1] for k, polys in enumerate(pr)]
    return pr

@autolog
def greedy_sotd_order(
    relations: sp.Expr | list[sp.Expr], var_groups: list[list[sp.Symbol]]
) -> list[sp.Symbol]:
    """
    Variable order that greedily minimizes the "sum of total
    degrees" metric, as advocated in [DSS04].
    """
    rev_order = []
    varlist = [v for g in var_groups for v in g]
    positives = relations_to_positives(relations, varlist)
    pr = SFRP([sp.Poly(p, *varlist) for p in positives], varlist)
    for group in reversed(var_groups):
        while len(group) > 1:
            n = len(varlist) - len(rev_order)
            with logblock(f"Var #{n}"):
                log(f"Searching among {group}")
                best_sotd = None
                for var in group:
                    rest = [v for v in varlist if v is not var]
                    with logblock(var):
                        new_pr = SFRP(PR(pr, var, rest), rest)
                    new_sotd = sum(sum(m) for p in new_pr for m in p.monoms())
                    if best_sotd is None or new_sotd < best_sotd:
                        best_var = var
                        best_pr = new_pr
                        best_sotd = new_sotd
                log(f"Best var #{len(varlist) - len(rev_order)}: {best_var}")
                rev_order.append(best_var)
                group.remove(best_var)
                pr = best_pr
        if len(group) == 1:
            n = len(varlist) - len(rev_order)
            log(f"Best var #{n}: {group[0]} (the only remaining)")
            rev_order.append(group[0])
    return list(reversed(rev_order))

def isolate_real_roots(pr: list[sp.Poly], subs: dict, var: sp.Symbol) -> list[PolyRoot]:
    """
    Find all real roots of the given list of square-free
    co-prime polynomials (Algorithm 3.5, step 3).
    """
    intervals = isolate_many_roots(
        [p.subs(subs).clear_denoms(convert=True)[1] for p in pr]
    )
    roots = [
        PolyRoot(p, ridx, lo, hi)
        for p, i in zip(pr, intervals)
        for ridx, (lo, hi) in enumerate(i)
    ]
    roots.sort(key=lambda r: r.value_lo)
    return roots

@autolog
def RSFC(
    positives: list[sp.Poly], pr: list[list[sp.Poly]], varlist: list[sp.Symbol]
) -> list[Cell]:
    """Recursive Solution Formula Construction (Algorithm 3.5 of [S00])."""
    def _RSFC(cell: Cell, positives: list[sp.Poly]):
        nonlocal n_rejected_cells, n_early_exits
        k = len(cell)
        if k >= len(varlist):
            if all(sp.sign(p.as_expr()) > 0 for p in positives):
                all_cells.append(cell)
            else:
                n_rejected_cells += 1
        else:
            # Early exit check.
            for p in positives:
                if p.total_degree() <= 0:
                    if sp.sign(p.as_expr()) < 0:
                        n_early_exits[k] += 1
                        return
            var = varlist[k]
            subs = dict(zip(varlist, (ab.point for ab in cell)))
            roots = isolate_real_roots(pr[k], subs, var)
            # Get a sample coordinate strictly inside each region
            # (i.e. between the roots), and recurse to the next var.
            if len(roots) > 0:
                hi = roots[0]
                mid = int(hi.value_lo) - 1
                _RSFC(
                    cell + [AxisBound(var, mid, None, hi)],
                    [p.subs({var: mid}) for p in positives],
                )
                for i in range(len(roots) - 1):
                    lo = roots[i]
                    hi = roots[i + 1]
                    a = lo.value_hi
                    b = hi.value_lo
                    mid = shortest_fraction_between(
                        a if lo.value_lo != lo.value_hi else a + (b - a) / 1024,
                        b if hi.value_lo != hi.value_hi else b - (b - a) / 1024,
                    )
                    _RSFC(
                        cell + [AxisBound(var, mid, lo, hi)],
                        [p.subs({var: mid}) for p in positives],
                    )
                lo = roots[-1]
                mid = int(lo.value_hi) + 1
                _RSFC(
                    cell + [AxisBound(var, mid, lo, None)],
                    [p.subs({var: mid}) for p in positives],
                )
            else:
                mid = sp.Rational(0)
                _RSFC(
                    cell + [AxisBound(var, mid, None, None)],
                    [p.subs({var: mid}) for p in positives],
                )
    positives = [sp.Poly(p, *varlist) for p in positives]
    all_cells = []
    n_rejected_cells = 0
    n_early_exits = [0] * len(varlist)
    _RSFC([], positives)
    log(f"Cells: {len(all_cells)} accepted, {n_rejected_cells} rejected")
    log(f"Early exits: {n_early_exits}")
    return all_cells

@autolog
def relations_to_positives(
    ex: sp.Expr | list[sp.Expr], varlist: list[sp.Symbol]
) -> list[sp.Poly]:
    """
    Turn one or more $a>b$ or $a<b$ relations into a list of positive
    expressions (i.e. turn $a>b$ into $a-b$, and $a<b$ into $b-a$).
    """
    positives = []
    todo = [ex]
    while todo:
        ex = todo.pop()
        if isinstance(ex, sp.StrictGreaterThan):
            positives.append(ex.lhs - ex.rhs)
        elif isinstance(ex, sp.StrictLessThan):
            positives.append(ex.rhs - ex.lhs)
        elif isinstance(ex, list):
            todo.extend(ex)
        elif isinstance(ex, sp.And):
            todo.extend(ex.args)
        elif ex == True:
            pass
        else:
            raise ValueError(f"Not a supported relation: {ex}")
    # We want to consistently have the positives to be multivariate
    # polynomials in all the variables with integer coefficients.
    return [sp.Poly(p, *varlist).clear_denoms(convert=True)[1] for p in positives]

@autolog
def GCAD(relations: sp.Expr | list[sp.Expr], varlist: list[sp.Symbol]) -> list[Cell]:
    """
    Generic Cylindrical Algebraic Decomposition (Algorithm 3.1 of
    [S00]). Take a list of multivariate polynomial inequalities,
    given as expressions that are implied to be positive, and
    decompose the region their conjunction defines into a set of
    cells. (The list of boundaries from the original algorithm
    is not computed here).
    """
    positives = relations_to_positives(relations, varlist)
    pr = GPROJ(positives, varlist)
    cells = RSFC(positives, pr, varlist)
    # At this point users of the API should consider merging
    # cells. We don't do that here.
    return cells

@autolog
def merge(cells: list[Cell]) -> list[Cell]:
    """
    Merge adjacent cells (Remark 3.7 of [S00]) using an aggressive
    merging strategy, merging cells even if the inequalities may
    not hold at the boundaries between the cells.
    """
    if len(cells) == 0:
        return []
    def PolyRoot_eq(r1: PolyRoot | None, r2: PolyRoot | None) -> bool:
        if r1 is None: return False # Infinity is never a separator
        if r2 is None: return False # Infinity is never a separator
        return r1.idx == r2.idx and r1.poly == r2.poly
    def AxisBound_eq(ab1: AxisBound, ab2: AxisBound) -> bool:
        return PolyRoot_eq(ab1.cell_lo, ab2.cell_lo) and \
               PolyRoot_eq(ab1.cell_hi, ab2.cell_hi)
    def Cell_eq(c1: Cell, c2: Cell) -> bool:
        return all(AxisBound_eq(ab1, ab2) for ab1, ab2 in zip(c1, c2))
    # While there's certainly a smarter way to find pairs to
    # merge, no way is more fool-proof than brute force :)
    dim = len(cells[0])
    cells = list(cells)
    n_merges = [0]*dim
    while True:
        merged = False
        for i in range(len(cells)):
            c1 = cells[i]
            if c1 is None: continue
            for j in range(len(cells)):
                if j == i: continue
                c2 = cells[j]
                if c2 is None: continue
                k = 0
                while True:
                    assert k < dim
                    if not AxisBound_eq(c1[k], c2[k]): break
                    k += 1
                if PolyRoot_eq(c1[k].cell_hi, c2[k].cell_lo) and Cell_eq(c1[k+1:], c2[k+1:]):
                    # We've got k matching axis bounds, an adjacent pair,
                    # and matching bounds afterwards. Let's merge.
                    c1 = list(c1)
                    c1[k] = AxisBound(c1[k].var, c1[k].point, c1[k].cell_lo, c2[k].cell_hi)
                    cells[i] = c1
                    # Mark the cell as absent.
                    cells[j] = None
                    merged = True
                    n_merges[k] += 1
        cells = [c for c in cells if c is not None]
        if not merged:
            log(f"Merges: {n_merges}")
            return cells
