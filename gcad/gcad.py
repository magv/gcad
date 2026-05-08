"""
Generic Cylindrical Algebraic Decomposition (GCAD),
as described in S00.

S00:
    A. Strzeboński.
    "Solving Systems of Strict Polynomial Inequalities".
    J. Symb. Comput. 29, 471--480 (2000).
    https://doi.org/10.1006/jsco.1999.0327
"""

from __future__ import annotations

from dataclasses import dataclass
from sympy import Rational
from gcad.root_isolation import isolate_many_roots
from gcad_c_ext import shortest_fraction_between
import sympy as sp

@dataclass(slots=True)
class PolyRoot:
    poly: sp.Poly
    idx: int
    # Isolating interval for the root value at the sample point.
    value_lo: Rational
    value_hi: Rational
    def __repr__(self):
        p = self.poly.as_expr().subs({self.poly.gen: sp.Symbol("#")})
        p = str(p).replace(" ", "").replace("**", "^")
        return f"Root[{p}, {self.idx}]"

@dataclass(slots=True)
class AxisBound:
    var: sp.Symbol
    point: Rational # A sample var value inside the cell.
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

def SFRP(polys: list[sp.Expr], sym: sp.Symbol) -> list[sp.Poly]:
    """
    Square-free and relatively prime polynomials multiplicatively
    generating the product of polys, as a poly in the given
    variable. (Definition 3.2)
    """
    poly_prod = sp.prod(p for p in polys)
    content, factors = sp.factor_list(poly_prod)
    assert len(content.free_symbols) == 0
    polys = uniq(sp.Poly(fac, sym) for fac, _ in factors)
    polys = [-p if p.LC().is_negative == True else p for p in polys]
    return polys

def PR(polys: list[sp.Poly], sym: sp.Symbol) -> list[sp.Expr]:
    """
    The set of leading coefficients, discriminants, and pairwise
    resultants of the given list of square-free co-prime
    polynomials, with respect to the given variable (Definition 3.2).
    """
    for p in polys:
        # Our polys always come from SFRP(..., sym), and are
        # already in this variable.
        assert p.gen == sym
    result = []
    for p in polys:
        lc = p.LC()
        assert lc != 0
        result.append(lc)
        disc = p.discriminant()
        if disc != 0:
            result.append(disc)
    n = len(polys)
    for i in range(n):
        for j in range(i + 1, n):
            r = polys[i].resultant(polys[j])
            if r != 0:
                result.append(r)
    return uniq(result)

def GPROJ(polys: list[sp.Expr], varlist: list[sp.Symbol]) -> list[list[sp.Poly]]:
    """
    Generic projection (Algorithm 3.4).
    Return [pr_1 ... pr_n]
    """
    n = len(varlist)
    # Note: our lists are 0-indexed, unlike the paper.
    F: list[list[sp.Expr]] = [[] for k in range(n)]
    pr: list[list[sp.Poly]] = [[] for k in range(n)]
    F[n - 1] = polys
    pr[n - 1] = SFRP(F[n - 1], varlist[n - 1])
    for k in reversed(range(0, n - 1)):
        F[k] = PR(pr[k + 1], varlist[k + 1])
        pr[k] = SFRP(F[k], varlist[k])
    # Filter out constants polys from the projection. These never
    # have roots, so we can safely skip them.
    pr = [[p for p in polys if p.degree() >= 1] for polys in pr]
    return pr

def isolate_real_roots(
    k: int, pr: list[list[sp.Poly]], subs: dict, var: sp.Symbol
) -> list[PolyRoot]:
    """
    Find all real roots of the given list of square-free
    co-prime polynomials (Algorithm 3.5, step 3).
    """
    intervals = isolate_many_roots(
        [sp.Poly(p.as_expr().subs(subs), var) for p in pr[k]]
    )
    roots = [
        PolyRoot(p, ridx, lo, hi)
        for p, i in zip(pr[k], intervals)
        for ridx, (lo, hi) in enumerate(i)
    ]
    roots.sort(key=lambda r: r.value_lo)
    return roots

def RSFC(
    inequalities: list[sp.Expr], pr: list[list[sp.Poly]], varlist: list[sp.Symbol]
) -> list[Cell]:
    """Recursive Solution Formula Construction (Algorithm 3.5)."""
    def _RSFC(k: int, cell: Cell):
        subs = dict(zip(varlist, (ab.point for ab in cell)))
        if k >= len(varlist):
            if all(sp.sign(eq.subs(subs)) > 0 for eq in inequalities):
                all_cells.append(cell)
        else:
            var = varlist[k]
            roots = isolate_real_roots(k, pr, subs, var)
            # Get a sample coordinate strictly inside each region
            # (i.e. between the roots), and recurse to the next var.
            if len(roots) > 0:
                hi = roots[0]
                _RSFC(k + 1, cell + [AxisBound(var, int(hi.value_lo) - 1, None, hi)])
                for i in range(len(roots) - 1):
                    lo = roots[i]
                    hi = roots[i + 1]
                    a = lo.value_hi
                    b = hi.value_lo
                    mid = shortest_fraction_between(
                        a if lo.value_lo != lo.value_hi else a + (b - a) / 1024,
                        b if hi.value_lo != hi.value_hi else b - (b - a) / 1024,
                    )
                    _RSFC(k + 1, cell + [AxisBound(var, mid, lo, hi)])
                lo = roots[-1]
                _RSFC(k + 1, cell + [AxisBound(var, int(lo.value_hi) + 1, lo, None)])
            else:
                _RSFC(k + 1, cell + [AxisBound(var, sp.Rational(0), None, None)])
    all_cells = []
    _RSFC(0, [])
    return all_cells

def GCAD(inequalities: list[sp.Expr], varlist: list[sp.Symbol]) -> list[Cell]:
    """
    Generic Cylindrical Algebraic Decomposition (Algorithm 3.1).
    Take a list of multivariate polynomial inequalities, given as
    expressions that are implied to be positive, and decompose
    the region their conjunction defines into a set of cells.
    (The list of boundaries from the original algorithm is not
    computed here).
    """
    assert inequalities
    pr = GPROJ(inequalities, varlist)
    cells = RSFC(inequalities, pr, varlist)
    # At this point users of the API should consider merging
    # cells. We don't do that here.
    return cells

def merge(cells: list[Cell]) -> list[Cell]:
    """
    Merge adjacent cells (Remark 3.7) using an aggressive merging
    strategy, merging cells even if the inequalities may not
    hold at the boundaries between the cells.
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
        cells = [c for c in cells if c is not None]
        if not merged:
            return cells
