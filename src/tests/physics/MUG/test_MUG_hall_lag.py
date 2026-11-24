#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: rotation of the polarization of a shear Alfven wave in a triply periodic box (xMHD with all fields on
# Lagrange elements) by the Hall term in Ohm's law, with and without significant electron inertia. Same cases as
# test_MUG_hall.py, using the test_MUG_gyro_lag driver (the combined Hall and gyroviscosity case is in
# test_MUG_gyro_lag.py)
import pytest
from test_MUG_gyro_lag import gyro_case
from test_MUG_gyro import validate_result

#============================================================================
# Each case runs with the nonlinear and the linearized model, which have separate regression values
#============================================================================
# Hall term
@pytest.mark.parametrize("linear,err_exp", ((False, 1.381234E-03), pytest.param(True, 1.712936E-03, marks=pytest.mark.linear)))
def test_hall_rotation(linear, err_exp):
    '''Shear Alfven wave with the Hall term in Ohm's law, which splits the whistler and ion cyclotron branches and
    rotates the polarization at half their split, k^2 v_A d_i/2 (k d_i = 0.32). The electron mass factor is 1, so
    electron inertia is negligible (k^2 d_e^2 = 3E-5).'''
    assert gyro_case(-1., hall=True, me_factor=1., linear=linear)
    assert validate_result(0.0, err_exp, hall=True, me_factor=1.)

# Electron inertia
@pytest.mark.parametrize("linear,err_exp", ((False, 1.376975E-03), pytest.param(True, 1.726795E-03, marks=pytest.mark.linear)))
def test_hall_inertia(linear, err_exp):
    '''Shear Alfven wave with the Hall term and an electron mass factor of 1000 (k^2 d_e^2 = 0.028), where electron
    inertia shifts the phase of V by ~0.035 relative to the inertia-free result.'''
    assert gyro_case(-1., hall=True, me_factor=1000., linear=linear)
    assert validate_result(0.0, err_exp, hall=True, me_factor=1000.)
