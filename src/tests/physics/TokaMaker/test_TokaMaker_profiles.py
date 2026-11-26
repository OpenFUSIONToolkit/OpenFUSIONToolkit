#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: flux-function profile helpers in TokaMaker.util (pure Python)
import os
import sys

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'python')))
from OpenFUSIONToolkit.TokaMaker.util import create_spline_flux_fun


@pytest.mark.coverage
def test_spline_flux_fun_normalized_to_axis():
    '''A spline profile with a nonzero axis value is scaled to 1 on axis'''
    x = np.array([0.0, 0.5, 1.0])
    y = np.array([2.0, 1.0, 0.0])
    prof = create_spline_flux_fun(40, x, y)
    expected = CubicSpline(x, y, bc_type=[(1, 0.0), (1, 0.0)])(np.linspace(0.0, 1.0, 40))
    assert prof['type'] == 'linterp'
    assert prof['y'][0] == pytest.approx(1.0)
    assert np.allclose(prof['y'], expected/2.0)


@pytest.mark.coverage
@pytest.mark.parametrize('peak', [3.0, -2.0])
def test_spline_flux_fun_zero_axis(peak):
    '''A spline profile that is zero on axis is scaled by its maximum absolute value instead of dividing by zero'''
    x = np.array([0.0, 0.5, 1.0])
    y = np.array([0.0, peak, 0.0])
    prof = create_spline_flux_fun(40, x, y)
    expected = CubicSpline(x, y, bc_type=[(1, 0.0), (1, 0.0)])(np.linspace(0.0, 1.0, 40))
    assert np.all(np.isfinite(prof['y']))
    assert np.allclose(prof['x'], np.linspace(0.0, 1.0, 40))
    assert np.max(np.abs(prof['y'])) == pytest.approx(1.0)
    assert np.allclose(prof['y'], expected/np.max(np.abs(expected)))
    assert np.allclose(create_spline_flux_fun(40, x, y, normalize=False)['y'], expected)


@pytest.mark.coverage
@pytest.mark.parametrize('axis,peak', [(1.e-12, 3.0), (-1.e-12, 3.0), (1.e-12, -2.0)])
def test_spline_flux_fun_tiny_axis(axis, peak):
    '''A spline profile whose axis value is tiny relative to its peak (eg. round-off) is
    scaled by its maximum absolute value, keeping its sign, rather than blown up by the axis value'''
    x = np.array([0.0, 0.5, 1.0])
    y = np.array([axis, peak, 0.0])
    prof = create_spline_flux_fun(40, x, y)
    expected = CubicSpline(x, y, bc_type=[(1, 0.0), (1, 0.0)])(np.linspace(0.0, 1.0, 40))
    assert np.max(np.abs(prof['y'])) == pytest.approx(1.0)
    assert np.allclose(prof['y'], expected/np.max(np.abs(expected)))


@pytest.mark.coverage
def test_spline_flux_fun_all_zero_warns():
    '''A spline profile that is zero everywhere is left unnormalized with a warning'''
    x = np.array([0.0, 0.5, 1.0])
    y = np.zeros(3)
    with pytest.warns(RuntimeWarning, match='not normalized'):
        prof = create_spline_flux_fun(40, x, y)
    assert np.all(np.isfinite(prof['y']))
    assert np.all(prof['y'] == 0.0)
