# Generic Cylindrical Algebraic Decomposition

Given a list of polynomials in $x_1,…,x_n$, GCAD separates $ℝ^n$
into *cells* of constant sign of all polynomials, each given as

$$\bigwedge_{1≤i≤n} l_i(x_1,…,x_{i-1}) < x_i < h_i(x_1,…,x_{i-1}),$$

where the boundaries $l_i$ and $h_i$ are either ±∞, or real roots
of polynomials.

Conventionally, the input polynomials are given as inequalities,
and only the cells that satisfy them are returned.

## How to use

To decompose $2-x^2-y^2-z^2>0 ∧ 1-x^2-(y-1)^2-z^2>0 ∧ 1-z^2>0$ into cells:

```python
import sympy as sp
import gcad

x, y, z = sp.symbols("x y z")
p1 = 2 - x**2 - y**2 - z**2
p2 = 1 - x**2 - (y-1)**2 - z**2
p3 = 1 - z**2
cells = gcad.merge(gcad.GCAD([p1 > 0, p2 > 0, p3 > 0], [x, y, z]))
for i, cell in enumerate(cells):
    print(f"Cell {i}:")
    for axis in cell:
        print(f"  {axis.cell_lo} < {axis.var} < {axis.cell_hi}")
print()
for i, cell in enumerate(cells):
    print(f"Cell {i} example point:")
    for axis in cell:
        print(f"  {axis.var} = {axis.point}")
```

The expected output is:

```
Cell 0:
  Root[#+1, 0] < x < Root[#-1, 0]
  Root[#^2-2*#+x^2, 0] < y < Root[#-1, 0]
  Root[#^2+x^2+y^2-2*y, 0] < z < Root[#^2+x^2+y^2-2*y, 1]
Cell 1:
  Root[#+1, 0] < x < Root[#-1, 0]
  Root[#-1, 0] < y < Root[#^2+x^2-2, 1]
  Root[#^2+x^2+y^2-2, 0] < z < Root[#^2+x^2+y^2-2, 1]

Cell 0 example point:
  x = -7/8
  y = 2/3
  z = 0
Cell 1 example point:
  x = -7/8
  y = 17/16
  z = 0
```

## How to install

You will need a C++ compiler ([GCC] or [Clang]), [FLINT], [MPFR],
and [GMP] development libraries installed on your system. Once
you have those, you can install this package via

```
pip install git+https://github.com/magv/gcad/
```

If you prefer to experiment with the source code, then clone it,
install the dependencies listed in `pyproject.toml`, and build
the library with `make`. Run tests with `make test`.

[GCC]: https://gcc.gnu.org/
[CLANG]: https://clang.llvm.org/
[FLINT]: https://flintlib.org/
[MPFR]: https://www.mpfr.org/
[GMP]: https://gmplib.org/

## References

This implementation is based on the work of A. Strzeboński ([S00]).
The root finding used here is based on the works of A. Akritas,
A. Strzeboński, and others ([ASV08], [AAS08], [ASV06], [AS05]).
The variable ordering heuristics are from [DSS04] and [dRE22].
The handling of rationals and polynomials is done by [FLINT],
[MPFR], and [GMP].

[S00]: https://doi.org/10.1006/jsco.1999.0327
[ASV08]: https://doi.org/10.15388/NA.2008.13.3.14557
[AAS08]: https://doi.org/10.55630/sjc.2008.2.145-162
[ASV06]: https://doi.org/10.1007/s00607-006-0186-y
[AS05]: https://doi.org/10.15388/NA.2005.10.4.15110
[DSS04]: https://doi.org/10.1145/1005285.1005303
[dRE22]: https://doi.org/10.1007/978-3-031-14788-3_17 

## License

This library is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either [version 3
of the License], or (at your option) any later version.

[version 3 of the License]: https://www.gnu.org/licenses/gpl-3.0.txt

## API

### Functions

<a id="gcad.GCAD"></a>
**GCAD** (relations: *Expr | list[Expr]*, varlist: *list[Symbol]*) → *list[list[[AxisBound](#gcad.AxisBound)]]*

Generic Cylindrical Algebraic Decomposition (Algorithm 3.1 of
[S00]). Take one or multiple multivariate polynomial
inequalities, and decompose the region their conjunction
defines into a set of cells. (The list of boundaries from
the original algorithm is not computed here).

<a id="gcad.GPROJ"></a>
**GPROJ** (positives: *list[Poly]*, varlist: *list[Symbol]*) → *list[list[Poly]]*

Generic projection (Algorithm 3.4 of [S00]). Return the list
of $pr_i$, as polynomials with integer coefficients.

<a id="gcad.RSFC"></a>
**RSFC** (positives: *list[Poly]*, pr: *list[list[Poly]]*, varlist: *list[Symbol]*) → *list[list[[AxisBound](#gcad.AxisBound)]]*

Recursive Solution Formula Construction (Algorithm 3.5 of [S00]).

<a id="gcad.greedy_mods_order"></a>
**greedy_mods_order** (relations: *Expr | list[Expr]*, var_groups: *list[list[Symbol]]*) → *list[Symbol]*

Variable order that greedily minimizes the "multiplication
of degree sum" metric, as advocated in [dRE22].

<a id="gcad.greedy_sotd_order"></a>
**greedy_sotd_order** (relations: *Expr | list[Expr]*, var_groups: *list[list[Symbol]]*) → *list[Symbol]*

Variable order that greedily minimizes the "sum of total
degrees" metric, as advocated in [DSS04].

<a id="gcad.merge"></a>
**merge** (cells: *list[list[[AxisBound](#gcad.AxisBound)]]*) → *list[list[[AxisBound](#gcad.AxisBound)]]*

Merge adjacent cells (Remark 3.7 of [S00]) using an aggressive
merging strategy, merging cells even if the inequalities may
not hold at the boundaries between the cells.

<a id="gcad.relations_to_positives"></a>
**relations_to_positives** (relations: *Expr | list[Expr]*, varlist: *list[Symbol]*) → *list[Poly]*

Turn one or more $a>b$ or $a<b$ relations into a list of positive
expressions (i.e. turn $a>b$ into $a-b$, and $a<b$ into $b-a$).

### Classes

<a id="gcad.AxisBound"></a>
**AxisBound**. Cell boundaries along one axis.

* **var**: *Symbol*. Variable name corresponding to this axis.
* **point**: *Rational*. A sample var value inside the cell.
* **cell_lo**: *[PolyRoot](#gcad.PolyRoot) | None*. Lower bounding value. `None` means negative infinity.
* **cell_hi**: *[PolyRoot](#gcad.PolyRoot) | None*. Higher bounding value. `None` means positive infinity.

<a id="gcad.PolyRoot"></a>
**PolyRoot**. Real root of a polynomial.

* **poly**: *Poly*. Polynomial defining the root.
* **idx**: *int*. Index of the real root, if roots were ordered by value.
* **value_lo**: *Rational*. Lower bound on the root's value at the sample point.
* **value_hi**: *Rational*. Upper bound on the root's value at the sample point.

