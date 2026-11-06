#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: X-point isoflux boundary helpers in TokaMaker.util
import pytest
import numpy as np
from tokamaker_testing import load_fixture
from OpenFUSIONToolkit.TokaMaker.util import create_isoflux_xpts, xpoints_from_moments


isoflux_xpts_expected = load_fixture('isoflux_xpts_expected.json')


def test_xpoints_from_moments():
    '''X-point locations from Miller shape moments are the exact analytic coordinates, for symmetric and up/down
    asymmetric shapes.'''
    r0, z0, a = 6.2, 0.3, 2.0
    kU, dU, kL, dL = 1.7, 0.33, 2.0, 0.5
    # Up-down symmetric (lower defaults to upper)
    xpts = xpoints_from_moments(r0, z0, a, kU, dU)
    assert xpts.shape == (2, 2)
    assert np.allclose(xpts[0], [r0 - a * dU, z0 + a * kU])
    assert np.allclose(xpts[1], [r0 - a * dU, z0 - a * kU])
    # Asymmetric
    xpts = xpoints_from_moments(r0, z0, a, kU, dU, kappa_lower=kL, delta_lower=dL)
    assert np.allclose(xpts[0], [r0 - a * dU, z0 + a * kU])
    assert np.allclose(xpts[1], [r0 - a * dL, z0 - a * kL])


@pytest.mark.parametrize('case', ['symmetric', 'asymmetric'])
def test_create_isoflux_xpts(case):
    '''The generated double-null boundary matches the reference contour, has exactly npts distinct points and passes
    through both X-points.'''
    ref = isoflux_xpts_expected[case]
    params = ref['params']
    pts = create_isoflux_xpts(**params)
    # Shape: exactly npts points, no duplicate closing point
    assert pts.shape == (params['npts'], 2)
    closed = np.vstack([pts, pts[0]])
    gaps = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    assert np.all(gaps > 1e-9), "boundary contains a duplicate point"
    # Both X-points must appear on the contour
    xpts = np.asarray(ref['xpoints'])
    for xpt in xpts:
        assert np.min(np.linalg.norm(pts - xpt, axis=1)) < 1e-8
    # Regression against the stored reference contour
    assert np.allclose(pts, np.asarray(ref['points']), rtol=1e-10, atol=1e-12)


def test_create_isoflux_xpts_requires_min_npts():
    '''Fewer than three boundary points is rejected.'''
    with pytest.raises(ValueError):
        create_isoflux_xpts(2, 6.2, 0.0, 2.0, 1.7, 0.33)
