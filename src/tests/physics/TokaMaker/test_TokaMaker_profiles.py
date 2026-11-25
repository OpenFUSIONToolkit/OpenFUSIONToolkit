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


def test_spline_flux_fun_normalized_to_axis():
    '''A spline profile with a nonzero axis value is scaled to 1 on axis'''
    x = np.array([0.0, 0.5, 1.0])
    y = np.array([2.0, 1.0, 0.0])
    prof = create_spline_flux_fun(40, x, y)
    expected = CubicSpline(x, y, bc_type=[(1, 0.0), (1, 0.0)])(np.linspace(0.0, 1.0, 40))
    assert prof['type'] == 'linterp'
    assert prof['y'][0] == pytest.approx(1.0)
    assert np.allclose(prof['y'], expected/2.0)


def test_spline_flux_fun_zero_axis():
    '''A spline profile that is zero on axis keeps its finite sampled shape instead of dividing by zero'''
    x = np.array([0.0, 0.5, 1.0])
    y = np.array([0.0, 1.0, 0.0])
    prof = create_spline_flux_fun(40, x, y)
    expected = CubicSpline(x, y, bc_type=[(1, 0.0), (1, 0.0)])(np.linspace(0.0, 1.0, 40))
    assert np.all(np.isfinite(prof['y']))
    assert np.allclose(prof['x'], np.linspace(0.0, 1.0, 40))
    assert np.allclose(prof['y'], expected)
    assert np.allclose(prof['y'], create_spline_flux_fun(40, x, y, normalize=False)['y'])
