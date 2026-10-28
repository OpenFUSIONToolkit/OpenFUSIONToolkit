import os
import sys

import numpy as np
import pytest
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


def green_integral(r, z, rc, zc, nphi=4000):
    '''Defining integral of the filament's Green's function, -(r/4 pi) int rc cos(phi)/|x-x'| dphi, and its R and Z
    derivatives (trapezoidal rule, exponentially convergent for these smooth periodic integrands)'''
    phi = 2.0*np.pi*np.arange(nphi)/nphi
    dist = np.sqrt(r*r + rc*rc - 2.0*r*rc*np.cos(phi) + (z - zc)**2)
    w = 2.0*np.pi/nphi
    I = w*np.sum(rc*np.cos(phi)/dist)
    dIdr = -w*np.sum(rc*np.cos(phi)*(r - rc*np.cos(phi))/dist**3)
    dIdz = -w*np.sum(rc*np.cos(phi)*(z - zc)/dist**3)
    return -r*I/(4.0*np.pi), -(I + r*dIdr)/(4.0*np.pi), -r*dIdz/(4.0*np.pi)


def green_closed_form(r, z, rc, zc):
    '''Closed form with SciPy's complete elliptic integrals (parameter m = k^2)'''
    m = 4.0*r*rc/((r + rc)**2 + (z - zc)**2)
    return -np.sqrt(r*rc/m)*((2.0 - m)*ellipk(m) - 2.0*ellipe(m))/(2.0*np.pi)


@pytest.mark.coverage
def test_green():
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in)
    assert run_OFT('./test_axi_green', 1, 60)
    data = np.loadtxt('axi_green.results', skiprows=1)
    rc, zc = np.loadtxt('axi_green.results', max_rows=1)
    for r, z, g, fg, dgdr, dgdz in data:
        g_ref, dgdr_ref, dgdz_ref = green_integral(r, z, rc, zc)
        assert g == pytest.approx(g_ref, rel=1.0E-9)
        assert fg == pytest.approx(g_ref, rel=1.0E-9)
        assert g == pytest.approx(green_closed_form(r, z, rc, zc), rel=1.0E-9)
        assert dgdr == pytest.approx(dgdr_ref, rel=1.0E-8)
        assert dgdz == pytest.approx(dgdz_ref, rel=1.0E-8)
        # Far field: magnetic dipole, -r^2 rc^2/(4 rho^3) with rho the distance from the filament's center
        rho = np.hypot(r, z - zc)
        if rho > 10.0*rc:
            eps = (rc/rho)**2
            assert g == pytest.approx(-(r*rc)**2/(4.0*rho**3), rel=5.0*eps)
            assert dgdr == pytest.approx(-r*rc*rc*(2.0*rho**2 - 3.0*r*r)/(4.0*rho**5), rel=5.0*eps)
            assert dgdz == pytest.approx(3.0*(r*rc)**2*(z - zc)/(4.0*rho**5), rel=5.0*eps)
