# Generic Cylindrical Algebraic Decomposition

Given a list of polynomials in $x_1,…,x_n$, GCAD separates $ℝ^n$
into *cells* of constant sign of all polynomials, each given as

$$\bigwedge_{1≤i≤n} l_i(x_1,…,x_{i-1}) < x_i < h_i(x_1,…,x_{i-1}),$$

where the boundaries $l_i$ and $h_i$ are either ±∞, or real roots
of polynomials.

Conventionally, the input polynomials are given as inequalities,
and only the cells that satisfy inequalities are returned.

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
```

## How to install

You will need a C++ compiler (GCC or Clang), [FLINT], [MPFR],
and [GMP] development libraries installed on your system. Once
you have those, you can install this package via

```
pip install git+https://github.com/magv/gcad/
```

If you prefer to experiment with the source code, clone the source
code, install the dependencies listed in `pyproject.toml`, and
build the library with:

```
make
```

[FLINT]: https://flintlib.org/
[MPFR]: https://www.mpfr.org/
[GMP]: https://gmplib.org/

## References

This implementation is based on the work of A. Strzeboński ([S00]).
The root finding used here is based on the works of A. Akritas,
A. Strzeboński, and others ([ASV08], [AAS08], [ASV06], [AS05]).
The handling of rationals and polynomials is done by [FLINT],
[MPFR], and [GMP].

[S00]: https://doi.org/10.1006/jsco.1999.0327
[ASV08]: https://doi.org/10.15388/NA.2008.13.3.14557
[AAS08]: https://doi.org/10.55630/sjc.2008.2.145-162
[ASV06]: https://doi.org/10.1007/s00607-006-0186-y
[AS05]: https://doi.org/10.15388/NA.2005.10.4.15110
