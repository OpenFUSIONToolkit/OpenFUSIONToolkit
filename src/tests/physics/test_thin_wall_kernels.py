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
