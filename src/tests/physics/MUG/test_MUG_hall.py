#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: rotation of the polarization of a shear Alfven wave in a triply periodic box (xMHD, H(curl) magnetic
# field) by the Hall term in Ohm's law, with and without significant electron inertia. Uses the test_MUG_gyro
# driver, input template, and linear theory (the combined Hall and gyroviscosity case is in test_MUG_gyro.py)
import pytest
from test_MUG_gyro import gyro_case, validate_result

#============================================================================
# Each case runs with the nonlinear and the linearized model, which have separate regression values
#============================================================================
# Hall term
@pytest.mark.parametrize("linear,petsc_flag,err_exp", ((False, False, 7.690471E-04),
                                                  pytest.param(True, False, 1.100484E-03, marks=pytest.mark.linear),
                                                  (False, True, 7.690471E-04)))
def test_hall_rotation(linear, petsc_flag, err_exp):
    '''Shear Alfven wave with the Hall term in Ohm's law, which splits the whistler and ion cyclotron branches and
    rotates the polarization at half their split, k^2 v_A d_i/2 (k d_i = 0.32). The opposite sign of the Hall term
    would rotate it the other way. The electron mass factor is 1, so electron inertia is negligible
    (k^2 d_e^2 = 3E-5).'''
    assert gyro_case(-1., hall=True, me_factor=1., linear=linear, petsc=petsc_flag)
    assert validate_result(0.0, err_exp, hall=True, me_factor=1.)

# Electron inertia
@pytest.mark.parametrize("linear,petsc_flag,err_exp", ((False, False, 9.079000E-04),
                                                  pytest.param(True, False, 1.258472E-03, marks=pytest.mark.linear),
                                                  (False, True, 9.079000E-04)))
def test_hall_inertia(linear, petsc_flag, err_exp):
    '''Shear Alfven wave with the Hall term and an electron mass factor of 1000 (k^2 d_e^2 = 0.028), where electron
    inertia shifts the phase of V by ~0.035 relative to the inertia-free result.'''
    assert gyro_case(-1., hall=True, me_factor=1000., linear=linear, petsc=petsc_flag)
    assert validate_result(0.0, err_exp, hall=True, me_factor=1000.)
