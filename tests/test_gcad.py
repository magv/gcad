from gcad.gcad import *
import pytest

x, y, z, t, u = sp.symbols("x y z t u")

def test_input_unmodified():
    p1 = t*x + u*y - z
    variables = [[t,u],[x,y,z]]
    variables_flat = [t,u,x,y,z]

    greedy_sotd_order([p1>0],variables)
    assert p1 == t*x + u*y - z
    assert variables == [[t,u],[x,y,z]]

    greedy_mods_order([p1>0],variables)
    assert p1 == t*x + u*y - z
    assert variables == [[t,u],[x,y,z]]

    greedy_t1_order([p1>0],variables)
    assert p1 == t*x + u*y - z
    assert variables == [[t,u],[x,y,z]]

    GCAD([p1>0],variables_flat)
    assert p1 == t*x + u*y - z
    assert variables_flat == [t,u,x,y,z]

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
    
def test_ex422_bubble():
    # Example 4.2.2 from JOS25.
    b, x1, x2 = sp.symbols("b x1 x2")
    F = x1**2 + x2**2 - 2*b*x1*x2
    # Above the threshold
    cells = merge(GCAD([F>0, 0<x1, 0<x2, 1<b], [b, x1, x2]))
    assert len(cells) == 2
    cells = merge(GCAD([F<0, 0<x1, 0<x2, 1<b], [b, x1, x2]))
    assert len(cells) == 1

def test_ex411_box():
    # Example 4.1.1 from JOS25.
    s12, s13, x1, x2, x3, x4 = sp.symbols("s12 s13 x1 x2 x3 x4")
    F = -s12*x1*x2 -s13*x3*x4
    # Above the threshold
    cells = merge(GCAD([F<0, 0<x1, 0<x2, 0<x3, 0<x4, 0<s12, s13<0], [s12, s13, x1, x2, x3, x4]))
    assert len(cells) == 1
    # Below the threshold
    cells = merge(GCAD([F>0, 0<x1, 0<x2, 0<x3, 0<x4, s12<0, s13<0], [s12, s13, x1, x2, x3, x4]))
    assert len(cells) == 1

def test_ex412_box():
    # Example 4.1.2 from JOS25.
    s12, s13, p1s, x1, x2, x3, x4 = sp.symbols("s12 s13 p1s x1 x2 x3 x4")
    F = -s12*x1*x2 -s13*x3*x4 -p1s*x1*x3
    # Above the s12 threshold
    cells = merge(GCAD([F<0, 0<x1, 0<x2, 0<x3, 0<x4, 0<s12, s13<0, 0<p1s], [s12, s13, p1s, x1, x2, x3, x4]))
    assert len(cells) == 1
    cells = merge(GCAD([F>0, 0<x1, 0<x2, 0<x3, 0<x4, 0<s12, s13<0, 0<p1s], [s12, s13, p1s, x1, x2, x3, x4]))
    assert len(cells) == 1
    
def test_ex413_pentagon():
    # Example 4.1.3 from JOS25.
    s12, s23, s34, s45, s51, x1, x2, x3, x4, x5 = sp.symbols("s12 s23 s34 s45 s51 x1 x2 x3 x4 x5")
    F = -s12*x2*x5 -s23*x1*x3 -s34*x2*x4 -s45*x3*x5 -s51*x1*x4
    # Note: gmods order is [s51, s45, s34, s23, s12, x5, x4, x3, x2, x1] which gives 3 + 3 cells
    # Above the s12, s34, s51 threshold
    cells = merge(GCAD([F<0, 0<x1, 0<x2, 0<x3, 0<x4, 0<x5, 0<s12, 0<s34, 0<s51, s23<0, s45<0], [s12, s23, s34, s45, s51, x1, x2, x4, x5, x3]))
    assert len(cells) == 1
    cells = merge(GCAD([F>0, 0<x1, 0<x2, 0<x3, 0<x4, 0<x5, 0<s12, 0<s34, 0<s51, s23<0, s45<0], [s12, s23, s34, s45, s51, x1, x2, x4, x5, x3]))
    assert len(cells) == 1
    # GenericCylindricalDecomposition[f < 0 &&
    # {x1, x2, x3, x4, x5} \[Element] PositiveReals &&
    # 0 < s12 && 0 < s34 && 0 < s51 && s23 < 0 && s45 < 0,
    # {s12, s23, s34, s45, s51, x1, x2, x4, x5, x3}]]
    
def test_ex414_bnp6():
    # Example 4.1.4 from JOS25.
    s12, s23, x1, x2, x3, x4, x5, x6 = sp.symbols("s12 s23 x1 x2 x3 x4 x5 x6")
    F = -s12*x2*x3*x6 -s23*x1*x2*x4 +(s12+s23)*x1*x3*x5
    # Above the s12 threshold
    cells = merge(GCAD([F<0, 0<x1, 0<x2, 0<x3, 0<x4, 0<x5, 0<x6, 0<s12, s23<0, -s12<s23], [s23, s12, x5, x4, x3, x2, x6, x1]))
    assert len(cells) == 1
    cells = merge(GCAD([F>0, 0<x1, 0<x2, 0<x3, 0<x4, 0<x5, 0<x6, 0<s12, s23<0, -s12<s23], [s23, s12, x5, x4, x3, x2, x6, x1]))
    assert len(cells) == 1
    
def test_ex415_bnp7():
    # Example 4.1.5 from JOS25.
    s12, s23, x1, x2, x3, x4, x5, x6, x7 = sp.symbols("s12 s23 x1 x2 x3 x4 x5 x6 x7")
    F = -s12*(x3*x4*x6+x2*x5*x7+x2*x3*x7+x2*x3*x6+x2*x3*x5+x2*x3*x4)-s23*x1*x5*x6+(s12+s23)*x1*x4*x7
    # Above the s12 threshold
    cells = merge(GCAD([F<0, 0<x1, 0<x2, 0<x3, 0<x4, 0<x5, 0<x6, 0<x7, 0<s12, s23<0, -s12<s23], [s23, s12, x7, x6, x5, x4, x3, x2, x1]))
    assert len(cells) == 1
    cells = merge(GCAD([F>0, 0<x1, 0<x2, 0<x3, 0<x4, 0<x5, 0<x6, 0<x7, 0<s12, s23<0, -s12<s23], [s23, s12, x7, x6, x5, x4, x3, x2, x1]))
    assert len(cells) == 1
    # GenericCylindricalDecomposition[
    # f < 0 && {x1, x2, x3, x4, x5, x6, x7} \[Element] PositiveReals &&
    # 0 < s12 && s23 < 0 && -s12 < s23,
    # {s23, s12, x7, x6, x5, x4, x3, x2, x1}]]

@pytest.mark.slow
def test_massless_planar_elliptic():
    # Example of a planar massless elliptic integral
    x0, x1, x2, x3, x4, x5, x6 = sp.symbols("x0 x1 x2 x3 x4 x5 x6")
    p1 = 7*x0*x1*x2 + 7*x0*x1*x3 - 23*x0*x2*x3 - 23*x1*x2*x3 + 7*x0*x1*x4 + 29*x0*x2*x4 + 29*x1*x2*x4 + 31*x0*x3*x4 + 31*x1*x3*x4 - 2*x0*x2*x5 - 3*x1*x2*x5 - 2*x0*x3*x5 - 3*x1*x3*x5 - 23*x2*x3*x5 - 2*x0*x4*x5 - 3*x1*x4*x5 + 29*x2*x4*x5 + 31*x3*x4*x5 + 7*x0*x1*x6 + 11*x0*x2*x6 - 13*x0*x3*x6 - 17*x1*x3*x6 - 23*x2*x3*x6 + 41*x0*x4*x6 + 19*x1*x4*x6 + 29*x2*x4*x6 + 31*x3*x4*x6 - 2*x0*x5*x6 - 3*x1*x5*x6 + 5*x2*x5*x6 + 37*x3*x5*x6
    problem = [ p1 > 0, x0 > 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, x6 > 0]
    variables = [[x0, x1, x2, x3, x4, x5, x6]]
    greedy_variables = greedy_mods_order(problem, variables)
    assert greedy_variables == [x4, x3, x2, x6, x5, x1, x0]
    cells = GCAD(problem, greedy_variables)
    assert len(cells) == 138678
    # TODO: too slow to compute merged cells
    # mathematica: 793 cells in 561 seconds
    #merged_cells = merge(cells)
    #assert len(cells) == ???

@pytest.mark.slow
def test_elliptic2l_physical():
    # pySecDec example: elliptic2L_physical
    x0, x1, x2, x3, x4, x5, x6 = sp.symbols("x0 x1 x2 x3 x4 x5 x6")
    s, t, pp4, msq = sp.symbols("s t pp4 msq")
    ff =  + (msq)*x5*x6**2 + (msq)*x4*x6**2 + (2*msq - t)*x4*x5*x6 + (msq)*x4**2*x6 + (msq)*x4**2*x5 + (2*msq - pp4)*x3*x5*x6 + (2*msq - pp4)*x3*x4*x6 + (2*msq)*x3*x4*x5 + (msq)*x3*x4**2 + (msq)*x3**2*x5 + (msq)*x3**2*x4 + (msq)*x1*x6**2 + (2*msq - pp4)*x1*x5*x6 + (2*msq)*x1*x4*x6 + (2*msq)*x1*x4*x5 + (2*msq - pp4)*x1*x3*x6 + (2*msq)*x1*x3*x5 + (2*msq)*x1*x3*x4 + (msq)*x1*x3**2 + (msq)*x1**2*x6 + (msq)*x1**2*x5 + (msq)*x1**2*x3 + (msq)*x0*x6**2 + (2*msq)*x0*x5*x6 + (2*msq)*x0*x4*x6 + (2*msq)*x0*x4*x5 + (2*msq - pp4)*x0*x3*x6 + (2*msq - s)*x0*x3*x5 + (2*msq)*x0*x3*x4 + (msq)*x0*x3**2 + (2*msq - s)*x0*x1*x6 + (2*msq - s)*x0*x1*x5 + (2*msq - s)*x0*x1*x3 + (msq)*x0**2*x6 + (msq)*x0**2*x5 + (msq)*x0**2*x3
    ff = ff.subs({msq:1})
    problem = [ ff > 0, x0 > 0, x1 > 0, x2 > 0, x3 > 0, x5 > 0, x6 > 0, pp4 < 4, pp4 > 0, s > 4, t < 0]
    variables = [[s,t,pp4],[x0, x1, x2, x3, x4, x5, x6]] # interesting to switch x4 <> x5 then gmods takes a long time, t1 is still fast
    greedy_variables = greedy_t1_order(problem, variables)
    assert greedy_variables == [t, pp4, s, x6, x3, x5, x1, x0, x4, x2]
    cells = GCAD(problem, greedy_variables)
    assert len(cells) == 28794
    # TODO: too slow to compute merged cells
    
@pytest.mark.slow
def test_ex_e5():
    # Integral E5 from LLSSV25.
    x0, x1, x2, x3, x4, x5 = sp.symbols("x0 x1 x2 x3 x4 x5")
    p1s, ms = sp.symbols("p1s ms")
    ff =  + (ms)*x2*x3*x5**2 + (2*ms - p1s)*x2*x3*x4*x5 + (ms)*x2*x3*x4**2 + (ms)*x2*x3**2*x5 + (ms)*x2*x3**2*x4 + (ms)*x2**2*x3*x5 + (ms)*x2**2*x3*x4 + (ms)*x1*x3*x5**2 + (2*ms - p1s)*x1*x3*x4*x5 + (ms)*x1*x3*x4**2 + (ms)*x1*x3**2*x5 + (ms)*x1*x3**2*x4 + (ms)*x1*x2*x5**2 + (2*ms - p1s)*x1*x2*x4*x5 + (ms)*x1*x2*x4**2 + (4*ms)*x1*x2*x3*x5 + (4*ms - p1s)*x1*x2*x3*x4 + (ms)*x1*x2*x3**2 + (ms)*x1*x2**2*x5 + (ms)*x1*x2**2*x4 + (ms)*x1*x2**2*x3 + (ms)*x1**2*x3*x5 + (ms)*x1**2*x3*x4 + (ms)*x1**2*x2*x5 + (ms)*x1**2*x2*x4 + (ms)*x1**2*x2*x3 + (ms)*x0*x3*x5**2 + (2*ms - p1s)*x0*x3*x4*x5 + (ms)*x0*x3*x4**2 + (ms)*x0*x3**2*x5 + (ms)*x0*x3**2*x4 + (ms)*x0*x2*x5**2 + (2*ms - p1s)*x0*x2*x4*x5 + (ms)*x0*x2*x4**2 + (4*ms - p1s)*x0*x2*x3*x5 + (4*ms)*x0*x2*x3*x4 + (ms)*x0*x2*x3**2 + (ms)*x0*x2**2*x5 + (ms)*x0*x2**2*x4 + (ms)*x0*x2**2*x3 + (2*ms - p1s)*x0*x1*x3*x5 + (2*ms - p1s)*x0*x1*x3*x4 + (2*ms - p1s)*x0*x1*x2*x5 + (2*ms - p1s)*x0*x1*x2*x4 + (2*ms - p1s)*x0*x1*x2*x3 + (ms)*x0**2*x3*x5 + (ms)*x0**2*x3*x4 + (ms)*x0**2*x2*x5 + (ms)*x0**2*x2*x4 + (ms)*x0**2*x2*x3
    ff = ff.subs({ms:1})
    # Below threshold
    problem1 = [ ff < 0, x0 > 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, p1s < 4]
    problem2 = [ ff > 0, x0 > 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, p1s < 4]
    variables = [[p1s],[x0, x1, x2, x3, x4, x5]]
    greedy_variables = greedy_t1_order(problem1, variables)
    assert greedy_variables == [p1s, x4, x5, x3, x2, x1, x0]
    cells = GCAD(problem1, greedy_variables)
    assert len(cells) == 0
    mcells = merge(cells)
    assert len(mcells) == 0
    cells = GCAD(problem2, greedy_variables)
    assert len(cells) == 25160
    # TODO: too slow to compute merged cells
    # mathematica: 1 cell in 17 seconds
    #mcells = merge(cells)
    #assert len(mcells) == 1
    # Above 2-particle on-shell threshold in p1s
    problem3 = [ ff < 0, x0 > 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, p1s > 4, p1s < 16]
    problem4 = [ ff > 0, x0 > 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, p1s > 4, p1s < 16]
    cells = GCAD(problem3, greedy_variables)
    assert len(cells) == 31457
    # TODO: too slow to compute merged cells
    # mathematica: X cells in Y seconds (too slow to determine)
    #mcells = merge(cells)
    #assert len(mcells) == 0
    cells = GCAD(problem4, greedy_variables)
    assert len(cells) == 105322
    # TODO: too slow to compute merged cells
    # mathematica: X cells in Y seconds (too slow to determine)
    #mcells = merge(cells)
    #assert len(mcells) == 0

def test_melih_box():
    # Double box integral provided by Melih Ozcelik
    x1, x2, x3, x4, x5, x6, x7 = sp.symbols("x1 x2 x3 x4 x5 x6 x7")
    ff = -4*x1*x4*x5 + 2*x1*x3*x6 + 2*x2*x3*x6 + 2*x2*x5*x6 + 2*x3*x5*x6 + x1*x6**2 + x2*x6**2 + x5*x6**2 + 2*x2*x3*x7 + 2*x2*x4*x7 + 2*x2*x5*x7 + 2*x3*x5*x7 + 2*x2*x6*x7 + 2*x3*x6*x7 + 2*x5*x6*x7 + x6**2*x7 + x3*x7**2 + x4*x7**2 + x5*x7**2 + x6*x7**2
    problem1 = [ ff < 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, x6 >0, x7 > 0]
    problem2 = [ ff > 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, x6 >0, x7 > 0]
    variables = [[x1, x2, x3, x4, x5, x6, x7]]
    greedy_variables = greedy_t1_order(problem1, variables)
    assert greedy_variables == [x7, x5, x6, x3, x4, x2, x1]
    mcells = merge(GCAD(problem1, greedy_variables))
    assert len(mcells) == 1
    mcells = merge(GCAD(problem2, greedy_variables))
    assert len(mcells) == 2
    
def test_triangle2l_split():
    # pySecDec example: triangle2L_split
    x0, x1, x2, x3, x4, x5 = sp.symbols("x0 x1 x2 x3 x4 x5")
    ff = + (1)*x3**2*x5 + (1)*x3**2*x4 + (1)*x2*x3*x5 + (1)*x2*x3*x4 + (1)*x1*x3*x5 + (1)*x1*x3*x4 + (1)*x1*x3**2 + (-1)*x1*x2*x4 + (1)*x1*x2*x3 + (1)*x0*x3*x4 + (1)*x0*x3**2 + (1)*x0*x2*x3 + (-1)*x0*x1*x5 + (-1)*x0*x1*x4 + (-1)*x0*x1*x3 + (-1)*x0*x1*x2
    problem1 = [ ff < 0, x0 > 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0]
    problem2 = [ ff > 0, x0 > 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0]
    variables = [[x0, x1, x2, x3, x4, x5]]
    greedy_variables = greedy_t1_order(problem1, variables)
    assert greedy_variables == [x3, x1, x0, x2, x4, x5]
    mcells = merge(GCAD(problem1, greedy_variables))
    assert len(mcells) == 8
    mcells = merge(GCAD(problem2, greedy_variables))
    assert len(mcells) == 10

@pytest.mark.slow
@pytest.mark.intractable
def test_gluza_ex30():
    # Example lh_np30 from [2201.02576] provided by Krzysztof Grzanka/Janus Gluza on 12.08.2025
    x0, x1, x2, x3, x4, x5, x6, x7 = sp.symbols("x0 x1 x2 x3 x4 x5 x6 x7")
    ff = + (-1)*x4*x5*x6*x7 + (-1)*x3*x4*x6*x7 + (-1)*x3*x4*x5*x7 + (1)*x3*x4*x5*x6 + (1)*x3**2*x6*x7 + (1)*x3**2*x5*x7 + (1)*x3**2*x5*x6 + (1)*x3**2*x4*x6 + (1)*x3**2*x4*x5 + (-1)*x2*x4*x6*x7 + (-1)*x2*x4*x5*x6 + (-1)*x2*x3*x4*x7 + (-1)*x2*x3*x4*x5 + (1)*x2*x3**2*x7 + (1)*x2*x3**2*x5 + (1)*x2*x3**2*x4 + (-1)*x1*x4*x6*x7 + (-1)*x1*x4*x5*x7 + (1)*x1*x3*x4*x6 + (1)*x1*x3*x4*x5 + (1)*x1*x3**2*x6 + (1)*x1*x3**2*x5 + (-1)*x1*x2*x4*x7 + (-1)*x1*x2*x4*x6 + (-1)*x1*x2*x4*x5 + (1)*x1*x2*x3*x4 + (1)*x1*x2*x3**2 + (-1)*x0*x5*x6*x7 + (-1)*x0*x4*x5*x7 + (-1)*x0*x4*x5*x6 + (-1)*x0*x3*x6*x7 + (1)*x0*x3*x5*x7 + (-1)*x0*x3*x5*x6 + (-1)*x0*x3*x4*x7 + (-1)*x0*x3*x4*x6 + (1)*x0*x3**2*x7 + (1)*x0*x3**2*x5 + (1)*x0*x3**2*x4 + (-1)*x0*x2*x6*x7 + (-1)*x0*x2*x5*x6 + (-1)*x0*x2*x4*x7 + (-1)*x0*x2*x4*x6 + (-1)*x0*x2*x4*x5 + (1)*x0*x2*x3*x7 + (1)*x0*x2*x3*x5 + (1)*x0*x2*x3*x4 + (-1)*x0*x1*x6*x7 + (-1)*x0*x1*x5*x7 + (-1)*x0*x1*x4*x7 + (-1)*x0*x1*x4*x6 + (-1)*x0*x1*x4*x5 + (-1)*x0*x1*x3*x6 + (-1)*x0*x1*x3*x5 + (1)*x0*x1*x3*x4 + (1)*x0*x1*x3**2 + (-1)*x0*x1*x2*x7 + (-1)*x0*x1*x2*x6 + (-1)*x0*x1*x2*x5 + (1)*x0*x1*x2*x3
    ff = ff.subs({x3:1}) # t1/mods is very slow with x3 active
    problem1 = [ ff < 0, x0 > 0, x1 > 0, x2 > 0, x4 > 0, x5 > 0, x6 > 0, x7 > 0]
    problem2 = [ ff > 0, x0 > 0, x1 > 0, x2 > 0, x4 > 0, x5 > 0, x6 > 0, x7 > 0]
    variables = [[x0, x1, x2, x4, x5, x6, x7]]
    greedy_variables = greedy_t1_order(problem1, variables)
    assert greedy_variables == [x0, x6, x4, x2, x5, x7, x1]
    mcells = merge(GCAD(problem1, greedy_variables))
    #assert len(mcells) == 7 # TODO: example too slow to complete
    mcells = merge(GCAD(problem2, greedy_variables))
    #assert len(mcells) == 10 # TODO: example too slow to complete

@pytest.mark.slow
@pytest.mark.intractable
def test_gluza_ex33():
    # Example 33-35 from [2201.02576] provided by Krzysztof Grzanka/Janus Gluza on 12.08.2025
    MT, x0, x1, x2, x3, x4, x5, x6, x7 = sp.symbols("MT x0 x1 x2 x3 x4 x5 x6 x7")
    ff = + (MT**2)*x5*x6*x7**2 + (MT**2)*x5*x6**2*x7 + (2*MT**2 - 1)*x4*x5*x6*x7 + (MT**2)*x4*x5*x6**2 + (MT**2)*x4**2*x5*x6 + (MT**2)*x3*x6*x7**2 + (MT**2)*x3*x6**2*x7 + (MT**2)*x3*x5*x7**2 + (2*MT**2 - 1)*x3*x5*x6*x7 + (MT**2)*x3*x5*x6**2 + (2*MT**2 - 1)*x3*x4*x6*x7 + (MT**2)*x3*x4*x6**2 + (2*MT**2 - 1)*x3*x4*x5*x7 + (2*MT**2)*x3*x4*x5*x6 + (MT**2)*x3*x4**2*x6 + (MT**2)*x3*x4**2*x5 + (MT**2)*x2*x6*x7**2 + (MT**2)*x2*x6**2*x7 + (2*MT**2)*x2*x5*x6*x7 + (MT**2)*x2*x5*x6**2 + (2*MT**2 - 1)*x2*x4*x6*x7 + (MT**2)*x2*x4*x6**2 + (2*MT**2 - 1)*x2*x4*x5*x6 + (MT**2)*x2*x4**2*x6 + (MT**2)*x2*x3*x7**2 + (2*MT**2 - 1)*x2*x3*x6*x7 + (2*MT**2)*x2*x3*x5*x7 + (2*MT**2 - 1)*x2*x3*x5*x6 + (2*MT**2 - 1)*x2*x3*x4*x7 + (2*MT**2 - 1)*x2*x3*x4*x6 + (2*MT**2 - 1)*x2*x3*x4*x5 + (MT**2)*x2*x3*x4**2 + (MT**2)*x2**2*x6*x7 + (MT**2)*x2**2*x5*x6 + (MT**2)*x2**2*x4*x6 + (MT**2)*x2**2*x3*x7 + (MT**2)*x2**2*x3*x5 + (MT**2)*x2**2*x3*x4 + (MT**2)*x1*x6*x7**2 + (MT**2)*x1*x6**2*x7 + (MT**2)*x1*x5*x7**2 + (2*MT**2)*x1*x5*x6*x7 + (2*MT**2 - 1)*x1*x4*x6*x7 + (MT**2)*x1*x4*x6**2 + (2*MT**2 - 1)*x1*x4*x5*x7 + (2*MT**2)*x1*x4*x5*x6 + (MT**2)*x1*x4**2*x6 + (MT**2)*x1*x4**2*x5 + (2*MT**2 - 1)*x1*x3*x6*x7 + (MT**2)*x1*x3*x6**2 + (2*MT**2 - 1)*x1*x3*x5*x7 + (2*MT**2)*x1*x3*x5*x6 + (2*MT**2)*x1*x3*x4*x6 + (2*MT**2)*x1*x3*x4*x5 + (MT**2)*x1*x2*x7**2 + (4*MT**2)*x1*x2*x6*x7 + (MT**2)*x1*x2*x6**2 + (2*MT**2)*x1*x2*x5*x7 + (2*MT**2)*x1*x2*x5*x6 + (2*MT**2 - 1)*x1*x2*x4*x7 + (4*MT**2 - 1)*x1*x2*x4*x6 + (2*MT**2 - 1)*x1*x2*x4*x5 + (MT**2)*x1*x2*x4**2 + (2*MT**2 - 1)*x1*x2*x3*x7 + (2*MT**2 - 1)*x1*x2*x3*x6 + (2*MT**2 - 1)*x1*x2*x3*x5 + (2*MT**2)*x1*x2*x3*x4 + (MT**2)*x1*x2**2*x7 + (MT**2)*x1*x2**2*x6 + (MT**2)*x1*x2**2*x5 + (MT**2)*x1*x2**2*x4 + (MT**2)*x1*x2**2*x3 + (MT**2)*x1**2*x6*x7 + (MT**2)*x1**2*x5*x7 + (MT**2)*x1**2*x4*x6 + (MT**2)*x1**2*x4*x5 + (MT**2)*x1**2*x3*x6 + (MT**2)*x1**2*x3*x5 + (MT**2)*x1**2*x2*x7 + (MT**2)*x1**2*x2*x6 + (MT**2)*x1**2*x2*x5 + (MT**2)*x1**2*x2*x4 + (MT**2)*x1**2*x2*x3 + (MT**2)*x0*x5*x7**2 + (2*MT**2 - 1)*x0*x5*x6*x7 + (2*MT**2 - 1)*x0*x4*x5*x7 + (2*MT**2 - 1)*x0*x4*x5*x6 + (MT**2)*x0*x4**2*x5 + (MT**2)*x0*x3*x7**2 + (2*MT**2 - 1)*x0*x3*x6*x7 + (2*MT**2)*x0*x3*x5*x7 + (2*MT**2 - 1)*x0*x3*x5*x6 + (2*MT**2 - 1)*x0*x3*x4*x7 + (2*MT**2 - 1)*x0*x3*x4*x6 + (2*MT**2 - 1)*x0*x3*x4*x5 + (MT**2)*x0*x3*x4**2 + (MT**2)*x0*x2*x7**2 + (2*MT**2 - 1)*x0*x2*x6*x7 + (2*MT**2)*x0*x2*x5*x7 + (2*MT**2 - 1)*x0*x2*x5*x6 + (2*MT**2 - 1)*x0*x2*x4*x7 + (2*MT**2 - 1)*x0*x2*x4*x6 + (2*MT**2 - 1)*x0*x2*x4*x5 + (MT**2)*x0*x2*x4**2 + (2*MT**2)*x0*x2*x3*x7 + (2*MT**2)*x0*x2*x3*x5 + (2*MT**2)*x0*x2*x3*x4 + (MT**2)*x0*x2**2*x7 + (MT**2)*x0*x2**2*x5 + (MT**2)*x0*x2**2*x4 + (MT**2)*x0*x1*x7**2 + (2*MT**2 - 1)*x0*x1*x6*x7 + (2*MT**2 - 1)*x0*x1*x5*x7 + (2*MT**2 - 1)*x0*x1*x4*x7 + (2*MT**2 - 1)*x0*x1*x4*x6 + (2*MT**2 - 1)*x0*x1*x4*x5 + (MT**2)*x0*x1*x4**2 + (2*MT**2 - 1)*x0*x1*x3*x7 + (2*MT**2 - 1)*x0*x1*x3*x6 + (2*MT**2 - 1)*x0*x1*x3*x5 + (2*MT**2)*x0*x1*x3*x4 + (4*MT**2 - 1)*x0*x1*x2*x7 + (2*MT**2 - 1)*x0*x1*x2*x6 + (2*MT**2 - 1)*x0*x1*x2*x5 + (4*MT**2)*x0*x1*x2*x4 + (2*MT**2)*x0*x1*x2*x3 + (MT**2)*x0*x1*x2**2 + (MT**2)*x0*x1**2*x7 + (MT**2)*x0*x1**2*x4 + (MT**2)*x0*x1**2*x3 + (MT**2)*x0*x1**2*x2 + (MT**2)*x0**2*x5*x7 + (MT**2)*x0**2*x4*x5 + (MT**2)*x0**2*x3*x7 + (MT**2)*x0**2*x3*x5 + (MT**2)*x0**2*x3*x4 + (MT**2)*x0**2*x2*x7 + (MT**2)*x0**2*x2*x5 + (MT**2)*x0**2*x2*x4 + (MT**2)*x0**2*x1*x7 + (MT**2)*x0**2*x1*x4 + (MT**2)*x0**2*x1*x3 + (MT**2)*x0**2*x1*x2
    ff = ff.subs({MT:1})
    problem1 = [ ff < 0, x0 > 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, x6 > 0, x7 > 0]
    problem2 = [ ff > 0, x0 > 0, x1 > 0, x2 > 0, x3 > 0, x4 > 0, x5 > 0, x6 > 0, x7 > 0]
    variables = [[x0, x1, x2, x3, x4, x5, x6, x7]]
    greedy_variables = greedy_t1_order(problem1, variables)
    #assert greedy_variables == [x0, x6, x4, x2, x3, x5, x7, x1] # TODO: example too slow to complete
    mcells = merge(GCAD(problem1, greedy_variables))
    #assert len(mcells) == 7 # TODO: example too slow to complete
    mcells = merge(GCAD(problem2, greedy_variables))
    #assert len(mcells) == 10 # TODO: example too slow to complete

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

def test_ex44_HT_C1_HB2():
    # Example 4.4 from S00.
    B2 = [1 - (x-1)**2 - (y-1)**2 - (z-1)**2 > 0]
    C1 = [-x**2 - y**2 - z**2 - 2*y*z + 4*y + 4*z - 3 > 0, z - y + 1 > 0, y + 1 - z > 0]
    T = [-z**4 - (2*y**2 + 2*x**2 + 6)*z**2 - y**4 - 2*x**2*y**2 + 10*y**2 - x**4 + 10*x**2 - 9 > 0]
    HB2 = B2 + [x + y + z - 3 > 0]
    HT = T + [- x - y > 0]
    cells = merge(GCAD(HT + C1 + HB2, [x, y, z]))
    assert len(cells) == 0

def test_greedy_sotd_order():
    # Figure 1 and Figure 2 from DSS04.
    c1 = (x+3)**2 + (y+1)**2 - 4
    c2 = (x-3)**2 + (y-1)**2 - 4
    assert len(GCAD([c1>0,c2>0], [x,y])) == 7
    assert len(GCAD([c1>0,c2>0], [y,x])) == 9
    assert greedy_sotd_order([c1>0,c2>0], [[x,y]]) == [x,y]
    assert greedy_sotd_order([c1>0,c2>0], [[x],[y]]) == [x,y]
    assert greedy_sotd_order([c1>0,c2>0], [[y],[x]]) == [y,x]

def test_greedy_mods_order():
    p = x**4*y**2*z + x*y*z + 1
    assert greedy_mods_order([p>0], [[x,y,z]]) == [x,y,z]
