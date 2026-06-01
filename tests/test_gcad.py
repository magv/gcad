from gcad.gcad import *
import pytest

x, y, z, t = sp.symbols("x y z t")

def test_unbounded():
    cells = GCAD(x > 0, [x])
    assert len(cells) == 1
    assert cells[0][0].cell_lo.poly == sp.Poly(x, x)
    assert cells[0][0].cell_lo.idx == 0
    assert cells[0][0].cell_hi is None

def test_unbounded_neg():
    cells = GCAD(x < 0, [x])
    assert len(cells) == 1
    assert cells[0][0].cell_lo is None
    assert cells[0][0].cell_hi.poly == sp.Poly(x, x)
    assert cells[0][0].cell_hi.idx == 0

def test_no_real_roots():
    cells = GCAD(x**2 + x + 1 > 0, [x])
    assert len(cells) == 1
    assert cells[0][0].cell_lo is None
    assert cells[0][0].cell_hi is None

def test_quadrant():
    assert len(GCAD([x > 0, y > 0], [x, y])) == 1

def test_cylinder_3d():
    assert len(GCAD([x**2 + y**2 < 1, z**2 < 1], [x, y, z])) == 1

def test_empty_3d():
    assert GCAD([1 + x**2 < 0], [x, y, z]) == []

def test_cube():
    assert len(GCAD([1-x**2 > 0, 1-y**2 > 0, 1-z**2 > 0], [x, y, z])) == 1

def test_boundary():
    assert GCAD([x**2 < 0], [x]) == []

def test_x_y_z_minus_1():
    result = GCAD([0<x, x<1, 0<y, y<1, 0<z, z<1, x + y + z > 1], [z, y, x])
    assert len(result) == 2
    # 0 < z < 1 && 0 < y < 1-z && 1-y-z < x < 1
    assert result[0][0].cell_lo.poly == sp.Poly(z, z)
    assert result[0][0].cell_hi.poly == sp.Poly(z - 1, z)
    assert result[0][1].cell_lo.poly == sp.Poly(y, z, y)
    assert result[0][1].cell_hi.poly == sp.Poly(y + z - 1, z, y)
    assert result[0][2].cell_lo.poly == sp.Poly(x + y + z - 1, z, y, x)
    assert result[0][2].cell_hi.poly == sp.Poly(x - 1, z, y, x)
    # 0 < z < 1 && 1-z < y < 1 && 0 < x <  1
    assert result[1][0].cell_lo.poly == sp.Poly(z, z)
    assert result[1][0].cell_hi.poly == sp.Poly(z - 1, z)
    assert result[1][1].cell_lo.poly == sp.Poly(y + z - 1, z, y)
    assert result[1][1].cell_hi.poly == sp.Poly(y - 1, z, y)
    assert result[1][2].cell_lo.poly == sp.Poly(x, z, y, x)
    assert result[1][2].cell_hi.poly == sp.Poly(x - 1, z, y, x)

def test_hypercube_1():
    # GenericCylindricalDecomposition[{x>0, 1-x>0, y>0, 1-y>0, z>0, 1-z>0}, {x,y,z}]
    cells = GCAD([x>0, 1-x>0, y>0, 1-y>0, z>0, 1-z>0], [x, y, z])
    assert len(cells) == 1
    # 0 < x < 1
    assert cells[0][0].cell_lo.poly == sp.Poly(x, x)
    assert cells[0][0].cell_lo.idx == 0
    assert cells[0][0].cell_hi.poly == sp.Poly(x-1, x)
    assert cells[0][0].cell_hi.idx == 0
    # 0 < y < 1
    assert cells[0][1].cell_lo.poly == sp.Poly(y, x, y)
    assert cells[0][1].cell_lo.idx == 0
    assert cells[0][1].cell_hi.poly == sp.Poly(y-1, x, y)
    assert cells[0][1].cell_hi.idx == 0
    # 0 < z < 1
    assert cells[0][2].cell_lo.poly == sp.Poly(z, x, y, z)
    assert cells[0][2].cell_lo.idx == 0
    assert cells[0][2].cell_hi.poly == sp.Poly(z-1, x, y, z)
    assert cells[0][2].cell_hi.idx == 0

def test_hypercube_2():
    # GenericCylindricalDecomposition[{x>0, 1-x>0, y>0, 1-y>0, z>0, 1-z>0, x-y>0}, {x,y,z}]
    cells = GCAD([x>0, 1-x>0, y>0, 1-y>0, z>0, 1-z>0, x-y>0], [x, y, z])
    assert len(cells) == 1
    # 0 < x < 1
    assert cells[0][0].cell_lo.poly == sp.Poly(x, x)
    assert cells[0][0].cell_lo.idx == 0
    assert cells[0][0].cell_hi.poly == sp.Poly(x-1, x)
    assert cells[0][0].cell_hi.idx == 0
    # 0 < y < x
    assert cells[0][1].cell_lo.poly == sp.Poly(y, x, y)
    assert cells[0][1].cell_lo.idx == 0
    assert cells[0][1].cell_hi.poly == sp.Poly(x-y, x, y)
    assert cells[0][1].cell_hi.idx == 0
    # 0 < z < 1
    assert cells[0][2].cell_lo.poly == sp.Poly(z, x, y, z)
    assert cells[0][2].cell_lo.idx == 0
    assert cells[0][2].cell_hi.poly == sp.Poly(z-1, x, y, z)
    assert cells[0][2].cell_hi.idx == 0

def test_hypercube_3():
    # GenericCylindricalDecomposition[{x>0, 1-x>0, y>0, 1-y>0, z>0, 1-z>0, 1-x-y-z>0}, {x,y,z}]
    cells = GCAD([x>0, 1-x>0, y>0, 1-y>0, z>0, 1-z>0, 1-x-y-z>0], [x, y, z])
    assert len(cells) == 1
    # 0 < x < 1
    assert cells[0][0].cell_lo.poly == sp.Poly(x, x)
    assert cells[0][0].cell_lo.idx == 0
    assert cells[0][0].cell_hi.poly == sp.Poly(x-1, x)
    assert cells[0][0].cell_hi.idx == 0
    # 0 < y < 1 - x
    assert cells[0][1].cell_lo.poly == sp.Poly(y, x, y)
    assert cells[0][1].cell_lo.idx == 0
    assert cells[0][1].cell_hi.poly == sp.Poly(y-(1-x), x, y)
    assert cells[0][1].cell_hi.idx == 0
    # 0 < z < 1 - x - y
    assert cells[0][2].cell_lo.poly == sp.Poly(z, x, y, z)
    assert cells[0][2].cell_lo.idx == 0
    assert cells[0][2].cell_hi.poly == sp.Poly(z-(1-x-y), x, y, z)
    assert cells[0][2].cell_hi.idx == 0

def test_hypercube_4():
    # GenericCylindricalDecomposition[{x>0, 1-x>0, y>0, 1-y>0, z>0, 1-z>0, 2-x^2-y^2-z^2>0}, {x,y,z}]
    cells = GCAD([x>0, 1-x>0, y>0, 1-y>0, z>0, 1-z>0, 1-x**2-y**2-z**2>0], [x, y, z])
    assert len(cells) == 1
    # 0 < x < 1
    assert cells[0][0].cell_lo.poly == sp.Poly(x, x)
    assert cells[0][0].cell_lo.idx == 0
    assert cells[0][0].cell_hi.poly == sp.Poly(x-1, x)
    assert cells[0][0].cell_hi.idx == 0
    # 0 < y < Sqrt[1 - x^2]
    assert cells[0][1].cell_lo.poly == sp.Poly(y, x, y)
    assert cells[0][1].cell_lo.idx == 0
    assert cells[0][1].cell_hi.poly == sp.Poly(y**2-(1-x**2), x, y)
    assert cells[0][1].cell_hi.idx == 1
    # 0 < z < Sqrt[1 - x^2 - y^2]
    assert cells[0][2].cell_lo.poly == sp.Poly(z, x, y, z)
    assert cells[0][2].cell_lo.idx == 0
    assert cells[0][2].cell_hi.poly == sp.Poly(z**2-(1-x**2-y**2), x, y, z)
    assert cells[0][2].cell_hi.idx == 1

def test_hypercube_5():
    # GenericCylindricalDecomposition[{x>0, 1-x>0, y>0, 1-y>0, z>0, 1-z>0, 2-x^2-y^2-z^2>0}, {x,y,z}]
    cells = GCAD([x>0, 1-x>0, y>0, 1-y>0, z>0, 1-z>0, 2-x**2-y**2-z**2>0], [x, y, z])
    assert len(cells) == 2
    cells = merge(cells)
    assert len(cells) == 2
    # 0 < x < 1
    assert cells[0][0].cell_lo.poly == sp.Poly(x, x)
    assert cells[0][0].cell_lo.idx == 0
    assert cells[0][0].cell_hi.poly == sp.Poly(x-1, x)
    assert cells[0][0].cell_hi.idx == 0
    # 0 < y < Sqrt[1 - x^2]
    assert cells[0][1].cell_lo.poly == sp.Poly(y, x, y)
    assert cells[0][1].cell_lo.idx == 0
    assert cells[0][1].cell_hi.poly == sp.Poly(y**2-(1-x**2), x, y)
    assert cells[0][1].cell_hi.idx == 1
    # 0 < z < 1
    assert cells[0][2].cell_lo.poly == sp.Poly(z, x, y, z)
    assert cells[0][2].cell_lo.idx == 0
    assert cells[0][2].cell_hi.poly == sp.Poly(z-1, x, y, z)
    assert cells[0][2].cell_hi.idx == 0
    # 0 < x < 1
    assert cells[1][0].cell_lo.poly == sp.Poly(x, x)
    assert cells[1][0].cell_lo.idx == 0
    assert cells[1][0].cell_hi.poly == sp.Poly(x-1, x)
    assert cells[1][0].cell_hi.idx == 0
    # Sqrt[1 - x^2] < y < 1
    assert cells[1][1].cell_lo.poly == sp.Poly(y**2-(1-x**2), x, y)
    assert cells[1][1].cell_lo.idx == 1
    assert cells[1][1].cell_hi.poly == sp.Poly(y-1, x, y)
    assert cells[1][1].cell_hi.idx == 0
    # 0 < z < Sqrt[2 - x^2 - y^2]
    assert cells[1][2].cell_lo.poly == sp.Poly(z, x, y, z)
    assert cells[1][2].cell_lo.idx == 0
    assert cells[1][2].cell_hi.poly == sp.Poly(z**2-(2-x**2-y**2), x, y, z)
    assert cells[1][2].cell_hi.idx == 1

def test_massive_triangle():
    s, x0, x1 = sp.symbols("s x0 x1")
    F = -s*x0*x1 + x0 + x1
    # Above the threshold
    cells = merge(GCAD([F>0, 0<x0, x0<1, 0<x1, x1<1, 0<1-x0-x1, 0<s, s<10], [s, x0, x1]))
    assert len(cells) == 4
    # Below the threshold
    cells = merge(GCAD([F>0, 0<x0, x0<1, 0<x1, x1<1, 0<1-x0-x1, 0<s, s<4], [s, x0, x1]))
    assert len(cells) == 1

def test_x2_minus_y():
    cells = GCAD([x**2-y>0], [x, y])
    assert len(cells) == 1
    cells = merge(cells)
    assert len(cells) == 1
    # -inf < x < inf
    assert cells[0][0].cell_lo == None
    assert cells[0][0].cell_hi == None
    # -inf < y < x^2
    assert cells[0][1].cell_lo == None
    assert cells[0][1].cell_hi.poly == sp.Poly(x**2-y, x, y)
    assert cells[0][1].cell_hi.idx == 0

def test_y2_minus_x():
    cells = GCAD([y**2-x>0], [x, y])
    assert len(cells) == 3
    cells = merge(cells)
    assert len(cells) == 3

def hanging_test_ex38():
    # Example 3.8 from S00.
    p = [
        1 - x**4 - y**2 - z**2 - t**2 > 0,
        1 - x**2 - y**4 - z**2 - t**2 > 0,
        1 - x**2 - y**2 - z**4 - t**2 > 0,
        1 - x**2 - y**2 - z**2 - t**4 > 0,
        1 - x**2 - y**2 - z**2 - t**2 > 0
    ]
    cells = GCAD(p, [x, y, z, t])
    cells = merge(cells)
    assert len(cells) == 29

def failing_test_ex43():
    # Example 4.3 from S00.
    p = x*t**3 + (x + y + z)*t**2 + (x**2 + y**2 + z**2)*t + x**3 + y**3 + z**3 - 1
    pr = GPROJ(p > 0, [x,y,z,t])
    assert len(pr[3]) == 1
    assert len(pr[2]) == 2
    assert len(pr[1]) == 4
    assert len(pr[0]) == 16
    cells = RSFC([p], pr, [x,y,z,t])
    cells = merge(cells)
    assert len(cells) == 29

def test_ex44_B1_B2():
    # Example 4.4 from S00.
    B1 = [1 - x**2 - y**2 - z**2 > 0]
    B2 = [1 - (x-1)**2 - (y-1)**2 - (z-1)**2 > 0]
    cells = merge(GCAD(B1 + B2, [x, y, z]))
    assert len(cells) == 1

def test_ex44_B1_B4():
    # Example 4.4 from S00.
    B1 = [1 - x**2 - y**2 - z**2 > 0]
    B4 = [1 - (x-sp.Rational(3,2))**2 - (y-2)**2 - z**2 > 0]
    cells = merge(GCAD(B1 + B4, [x, y, z]))
    assert len(cells) == 0

def test_ex44_B1_B2_B3():
    # Example 4.4 from S00.
    B1 = [1 - x**2 - y**2 - z**2 > 0]
    B2 = [1 - (x-1)**2 - (y-1)**2 - (z-1)**2 > 0]
    B3 = [1 - (x-1)**2 - (y-1)**2 - (z+sp.Rational(1,2))**2 > 0]
    cells = merge(GCAD(B1 + B2 + B3, [x, y, z]))
    assert len(cells) == 2

def test_ex44_B1_B2_B4():
    # Example 4.4 from S00.
    B1 = [1 - x**2 - y**2 - z**2 > 0]
    B2 = [1 - (x-1)**2 - (y-1)**2 - (z-1)**2 > 0]
    B4 = [1 - (x-sp.Rational(3,2))**2 - (y-2)**2 - z**2 > 0]
    cells = merge(GCAD(B1 + B2 + B4, [x, y, z]))
    assert len(cells) == 0

def test_ex44_B1_C1():
    # Example 4.4 from S00.
    B1 = [1 - x**2 - y**2 - z**2 > 0]
    C1 = [-x**2 - y**2 - z**2 - 2*y*z + 4*y + 4*z - 3 > 0, z - y + 1 > 0, y + 1 - z > 0]
    cells = merge(GCAD(B1 + C1, [x, y, z]))
    assert len(cells) == 1

def test_ex44_B1_C2():
    # Example 4.4 from S00.
    B1 = [1 - x**2 - y**2 - z**2 > 0]
    C2 = [-x**2 - y**2 - z**2 - 2*y*z + 4*y + 4*z - 3 > 0, z - y - 1 > 0, y + 2 - z > 0]
    cells = merge(GCAD(B1 + C2, [x, y, z]))
    assert len(cells) == 0

def test_ex44_T_C1():
    # Example 4.4 from S00.
    C1 = [-x**2 - y**2 - z**2 - 2*y*z + 4*y + 4*z - 3 > 0, z - y + 1 > 0, y + 1 - z > 0]
    T = [-z**4 - (2*y**2 + 2*x**2 + 6)*z**2 - y**4 - 2*x**2*y**2 + 10*y**2 - x**4 + 10*x**2 - 9 > 0]
    cells = merge(GCAD(T + C1, [x, y, z]))
    assert len(cells) == 9

def test_ex44_T_B2():
    # Example 4.4 from S00.
    B2 = [1 - (x-1)**2 - (y-1)**2 - (z-1)**2 > 0]
    T = [-z**4 - (2*y**2 + 2*x**2 + 6)*z**2 - y**4 - 2*x**2*y**2 + 10*y**2 - x**4 + 10*x**2 - 9 > 0]
    cells = merge(GCAD(T + B2, [x, y, z]))
    assert len(cells) == 1

@pytest.mark.slow
def test_ex44_HB1_HB2_HB3():
    # Example 4.4 from S00.
    B1 = [1 - x**2 - y**2 - z**2 > 0]
    B2 = [1 - (x-1)**2 - (y-1)**2 - (z-1)**2 > 0]
    B3 = [1 - (x-1)**2 - (y-1)**2 - (z+sp.Rational(1,2))**2 > 0]
    HB1 = B1 + [- x - y - z > 0]
    HB2 = B2 + [x + y + z - 3 > 0]
    HB3 = B3 + [sp.Rational(3,2) - x - y - z > 0]
    cells = merge(GCAD(HB1 + HB2 + HB3, [x, y, z]))
    assert len(cells) == 0

@pytest.mark.slow
def test_ex44_HT_HB2_HB3():
    # Example 4.4 from S00.
    B2 = [1 - (x-1)**2 - (y-1)**2 - (z-1)**2 > 0]
    B3 = [1 - (x-1)**2 - (y-1)**2 - (z+sp.Rational(1,2))**2 > 0]
    T = [-z**4 - (2*y**2 + 2*x**2 + 6)*z**2 - y**4 - 2*x**2*y**2 + 10*y**2 - x**4 + 10*x**2 - 9 > 0]
    HB2 = B2 + [x + y + z - 3 > 0]
    HB3 = B3 + [sp.Rational(3,2) - x - y - z > 0]
    HT = T + [- x - y > 0]
    cells = merge(GCAD(HT + HB2 + HB3, [x, y, z]))
    assert len(cells) == 0

def failing_test_ex44_T_C1_B2():
    # Example 4.4 from S00.
    B2 = [1 - (x-1)**2 - (y-1)**2 - (z-1)**2 > 0]
    C1 = [-x**2 - y**2 - z**2 - 2*y*z + 4*y + 4*z - 3 > 0, z - y + 1 > 0, y + 1 - z > 0]
    T = [-z**4 - (2*y**2 + 2*x**2 + 6)*z**2 - y**4 - 2*x**2*y**2 + 10*y**2 - x**4 + 10*x**2 - 9 > 0]
    cells = merge(GCAD(T + C1 + B2, [x, y, z]))
    assert len(cells) == 28 # We get only 19.

@pytest.mark.slow
def test_ex44_HT_C1_HB2():
    # Example 4.4 from S00.
    B2 = [1 - (x-1)**2 - (y-1)**2 - (z-1)**2 > 0]
    C1 = [-x**2 - y**2 - z**2 - 2*y*z + 4*y + 4*z - 3 > 0, z - y + 1 > 0, y + 1 - z > 0]
    T = [-z**4 - (2*y**2 + 2*x**2 + 6)*z**2 - y**4 - 2*x**2*y**2 + 10*y**2 - x**4 + 10*x**2 - 9 > 0]
    HB2 = B2 + [x + y + z - 3 > 0]
    HT = T + [- x - y > 0]
    cells = merge(GCAD(HT + C1 + HB2, [x, y, z]))
    assert len(cells) == 0
