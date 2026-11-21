#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: the cut-cell quadrature rule of the plasma integrals (setting `lcfs_cut_cells`)
#
# test_gs_cut_cell.F90 writes the rule's points and weights for vertex values swept through every change of the
# clipped region; the test checks them against independent clipping and for continuity in the swept parameter.
import os
import sys

import numpy as np
import pytest

test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT

oft_in = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/
"""


def clip(poly, vals):
    '''Keep the part of a convex polygon (barycentric vertices) where a linear function (vertex values) is positive'''
    out = []
    n = len(poly)
    for i in range(n):
        j = (i + 1) % n
        if vals[i] > 0.0:
            out.append(poly[i])
        if (vals[i] > 0.0) != (vals[j] > 0.0):
            s = vals[i]/(vals[i] - vals[j])
            out.append(poly[i] + s*(poly[j] - poly[i]))
    return out


def triangle_rule(n=10):
    '''Collapsed Gauss-Legendre rule on the reference triangle (exact to degree 2n-2), weights summing to one'''
    x, w = np.polynomial.legendre.leggauss(n)
    x = 0.5*(x + 1.0)
    w = 0.5*w
    pts, wts = [], []
    for xi, wi in zip(x, w):
        for eta, we in zip(x, w):
            pts.append([xi, (1.0 - xi)*eta])
            wts.append(2.0*wi*we*(1.0 - xi))
    return np.array(pts), np.array(wts)


def poly_f(l1, l2):
    return 1.0 + 2.0*l1 - l2 + 3.0*l1**2 + l1*l2 - 2.0*l2**3 + (l1*l2)**2 + l1**5


def reference(dv, hv):
    '''Area fraction and integral of poly_f over {dv > 0} and every {hv_k > 0}, by independent clipping'''
    poly = [np.eye(3)[i] for i in range(3)]
    poly = clip(poly, list(dv))
    for h in hv:
        if len(poly) < 3:
            break
        poly = clip(poly, [float(np.dot(p, h)) for p in poly])
    if len(poly) < 3:
        return 0.0, 0.0
    qp, qw = triangle_rule()
    area, integral = 0.0, 0.0
    a = poly[0]
    for k in range(1, len(poly) - 1):
        b, c = poly[k], poly[k + 1]
        frac = (b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0])
        l1 = a[0] + (b[0] - a[0])*qp[:, 0] + (c[0] - a[0])*qp[:, 1]
        l2 = a[1] + (b[1] - a[1])*qp[:, 0] + (c[1] - a[1])*qp[:, 1]
        area += frac
        integral += frac*np.sum(qw*poly_f(l1, l2))
    return area, integral


@pytest.mark.coverage
def test_cut_cell_rule():
    '''Cut-cell rule (`gs_cut_cell_rule`) on one triangle with OFT's order 5 and 9 rules.

    A wholly-inside cell must get the cell's own rule bit for bit; the kept area and the integral of a degree-5
    polynomial must match an independent clipping with a collapsed Gauss rule; and the integral of a smooth
    function must change continuously as the clipping function and half-planes sweep through every transition
    (a tenfold finer sweep shrinks the largest step change about tenfold, while a jump would not shrink; and at
    each transition the rule just before and just after it gives the same integrals).
    '''
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in)
    assert run_OFT('./test_gs_cut_cell', 1, 120)
    # A wholly-inside cell gets exactly the cell's rule
    whole = np.loadtxt('gs_cut_cell_whole.results')
    assert np.all(whole[:, 2] == 0.0) and np.all(whole[:, 3] == 0.0)
    data = np.loadtxt('gs_cut_cell.results')
    for order in (5, 9):
        rows = data[data[:, 0] == order]
        # Exactness on the whole cell and the 1000-step sweeps, which pass through every region shape of the
        # 10000-step sweeps (those are used for the continuity check only, which keeps the test fast)
        for row in rows[rows[:, 2] < 2]:
            nplanes = int(row[4])
            dv, hv = row[9:12], [row[12:15], row[15:18]][:nplanes]
            area, integral = reference(dv, hv)
            # The rule integrates polynomials up to its order exactly over the clipped region
            assert row[6] == pytest.approx(area, abs=1.0E-13)
            assert row[7] == pytest.approx(integral, abs=1.0E-12)
        # Continuity in the swept parameter: refining the sweep tenfold shrinks the largest step change tenfold
        for sweep in range(1, 7):
            steps = []
            for ires in (1, 2):
                sel = rows[(rows[:, 1] == sweep) & (rows[:, 2] == ires)]
                sel = sel[np.argsort(sel[:, 5])]
                steps.append(np.max(np.abs(np.diff(sel[:, 8]))))
            assert steps[0] > 0.0
            assert steps[1] < 0.2*steps[0]
    # Continuity directly at each transition (a change of the point count or of a vertex value's sign, located by
    # bisection): the sum of the weights and both integrals just before and just after it (parameter -+ 1e-12)
    jumps = np.loadtxt('gs_cut_cell_jumps.results')
    for order in (5, 9):
        for sweep in range(1, 7):
            sel = jumps[(jumps[:, 0] == order) & (jumps[:, 1] == sweep)]
            assert len(sel) > 0
            assert np.max(np.abs(sel[:, 4:7] - sel[:, 8:11])) < 1.0E-10
