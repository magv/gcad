# Generic Cylindrical Algebraic Decomposition

# How to use

To decompose $2-x^2-y^2-z^2>0 ∧ 1-x^2-(y-1)^2-z^2>0 ∧ 1-z^2>0$ into cells:

```python
import sympy as sp
import gcad

x, y, z = sp.symbols("x y z")
p1 = 2 - x**2 - y**2 - z**2
p2 = 1 - x**2 - (y-1)**2 - z**2
p3 = 1 - z**2
cells = gcad.merge(gcad.GCAD([p1, p2, p3], [x, y, z]))
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

## How to build

In addition to the Python dependencies specified in `pyproject.toml`,
you will need [FLINT], [MPFR], and [GMP] development libraries
installed on your system. Once you have those, build the C
extension this project uses with:

```
make
```

[FLINT]: https://flintlib.org/
[MPFR]: https://www.mpfr.org/
[GMP]: https://gmplib.org/
