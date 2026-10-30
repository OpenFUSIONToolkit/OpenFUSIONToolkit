import os
import sys

import numpy as np
import pytest
from scipy.integrate import dblquad

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
