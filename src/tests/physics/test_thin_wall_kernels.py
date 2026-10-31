import os
import sys

import numpy as np
import pytest
from scipy.integrate import dblquad, quad
from scipy.optimize import minimize_scalar
from scipy.special import ellipe, ellipk

test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

oft_in = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/
"""


def triangle_integrals(pt, scale):
    '''Potential int 1/|r-r'| dA' over the triangle (0,0,0), (1,0,0), (0,1,0) (times scale) and its gradient with
    respect to the observation point, by adaptive quadrature (smooth integrands for points off the triangle)'''
    def integrate(f):
        return dblquad(lambda v, u: f(scale*np.array([u, v, 0.0])), 0.0, 1.0, 0.0, lambda u: 1.0 - u,
                       epsabs=0.0, epsrel=1.E-12)[0]*scale**2
    phi = integrate(lambda rp: 1.0/np.linalg.norm(pt - rp))
    grad = [integrate(lambda rp, k=k: -(pt[k] - rp[k])/np.linalg.norm(pt - rp)**3) for k in range(3)]
    return phi, np.array(grad)


@pytest.mark.coverage
def test_triangle_potential():
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in)
    assert run_OFT('./test_thin_wall_kernels', 1, 60)
    data = np.loadtxt('thin_wall_kernels.results')
    for row in data:
        scale, pt, phi, grad = row[0], row[1:4], row[4], row[5:8]
        phi_ref, grad_ref = triangle_integrals(pt, scale)
        assert np.all(np.isfinite(row))
        assert phi == pytest.approx(phi_ref, rel=1.E-9)
        assert np.allclose(grad, grad_ref, rtol=0.0, atol=1.E-9*np.linalg.norm(grad_ref))


def test_lowrank_qrcp():
    """Truncated pivoted QR + SVD compression against the rank of the SVD with the same truncation rule"""
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in)
    assert run_OFT('./test_thin_wall_kernels', 1, 60)
    with open('thin_wall_qrcp.results') as fid:
        lines = fid.readlines()
    i = 0
    while i < len(lines):
        m, n, tol, rank, err = lines[i].split()
        m, n, tol, rank, err = int(m), int(n), float(tol), int(rank), float(err)
        A = np.array([[float(v) for v in lines[i + 1 + k].split()] for k in range(m)])
        i += 1 + m
        s = np.linalg.svd(A, compute_uv=False)
        tail = np.sqrt(np.cumsum(s[::-1]**2))[::-1]  # tail[k]: norm of the discarded part at rank k
        svd_rank = max(1, int(np.sum(tail > tol)))
        if svd_rank >= m*n//(m + n):
            assert rank == -1
        else:
            assert svd_rank <= rank <= svd_rank + 1
            assert err <= tol


def test_aca_plus():
    """ACA+ of a smooth kernel block (absolute and relative tolerance) against the dense block; an incompressible
    block does not converge within a few terms"""
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in)
    assert run_OFT('./test_thin_wall_kernels', 1, 60)
    data = np.loadtxt('thin_wall_aca.results')
    assert data[0, 1] > 0 and data[0, 2] < 1.E-8
    assert data[1, 1] > 0 and data[1, 2] < 1.E-6
    assert data[2, 1] == -1


def seg_seg_reference(a0, a1, b0, b1):
    """Double integral of 1/|x-y| over two segments: exact along B, adaptive quadrature along A with breakpoints at
    the closest points to B's ends"""
    la, lb = np.linalg.norm(a1 - a0), np.linalg.norm(b1 - b0)
    ua, ub = (a1 - a0)/la, (b1 - b0)/lb
    def inner(s):
        w = a0 + s*ua - b0
        s0 = w @ ub
        rho = np.sqrt(max(w @ w - s0*s0, 0.0))
        return np.arcsinh((lb - s0)/rho) + np.arcsinh(s0/rho)
    pts = sorted(set(min(la, max(0.0, (p - a0) @ ua)) for p in (b0, b1)) - {0.0, la})
    return quad(inner, 0.0, la, points=pts or None, epsabs=0.0, epsrel=1.E-13, limit=2000)[0]


def test_seg_seg_int():
    """Straight segment pair integrals for close parallel, skew and collinear pairs against quadrature"""
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in)
    assert run_OFT('./test_thin_wall_kernels', 1, 60)
    for row in np.loadtxt('thin_wall_segseg.results'):
        p = row[:12].reshape(4, 3)
        assert row[12] == pytest.approx(seg_seg_reference(*p), rel=1.E-9)


def test_ring_poly_softened():
    """Softened ring kernel integrated along a polygon of the same ring against the softened ring-ring mutual"""
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in)
    assert run_OFT('./test_thin_wall_kernels', 1, 60)
    ring_poly, ring_ring = np.loadtxt('thin_wall_ringpoly.results')
    assert ring_poly == pytest.approx(ring_ring, rel=1.E-4)


def ring_seg_reference(a, b):
    """Integral along a segment of u . (closed integral of dl'/|r-r'| around the ring R=1, Z=0), with the ring
    integral from elliptic integrals and adaptive quadrature along the segment with a breakpoint at the closest point"""
    lb = np.linalg.norm(b - a)
    u = (b - a)/lb
    def kernel(s):
        y = a + s*u
        rho = np.hypot(y[0], y[1])
        q = (rho + 1.0)**2 + y[2]**2
        m = 4.0*rho/q
        aphi = 4.0/np.sqrt(q)*((2.0 - m)*ellipk(m) - 2.0*ellipe(m))/m
        return aphi*(u @ np.array([-y[1], y[0], 0.0]))/rho
    def dist(s):
        y = a + s*u
        return np.hypot(np.hypot(y[0], y[1]) - 1.0, y[2])
    s0 = minimize_scalar(dist, bounds=(0.0, lb), method='bounded', options={'xatol': 1.E-12*lb}).x
    return quad(kernel, 0.0, lb, points=[s0], epsabs=0.0, epsrel=1.E-13, limit=2000)[0]


def test_ring_seg_close():
    """Exact ring kernel along straight segments passing close to the ring against quadrature"""
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in)
    assert run_OFT('./test_thin_wall_kernels', 1, 60)
    for row in np.loadtxt('thin_wall_ringseg.results'):
        assert row[6] == pytest.approx(ring_seg_reference(row[:3], row[3:6]), rel=1.E-7)
