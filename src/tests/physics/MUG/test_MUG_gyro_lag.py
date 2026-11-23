#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: rotation of the polarization of a shear Alfven wave in a triply periodic box (xMHD with all fields on
# Lagrange elements) by Braginskii gyroviscosity and by the Hall term in Ohm's law. Same cases, input template, and
# linear theory as test_MUG_gyro.py, run with the test_MUG_gyro_lag driver. Shared with test_MUG_hall_lag.py, which
# holds the Hall-only cases
import pytest
from test_MUG_gyro import gyro_case as _gyro_case, validate_result

def gyro_case(nu_gyro, hall=False, me_factor=100., linear=False):
    return _gyro_case(nu_gyro, hall=hall, me_factor=me_factor, linear=linear, exe='./test_MUG_gyro_lag')

#============================================================================
# Each case runs with the nonlinear and the linearized model, which share the linear theory and the deviation
# tolerance but have separate regression values
#============================================================================
# Reference case without gyroviscosity
@pytest.mark.parametrize("linear,err_exp", ((False, 4.266953E-04), pytest.param(True, 9.117855E-04, marks=pytest.mark.linear)))
def test_alfven_ref(linear, err_exp):
    '''Linearly polarized shear Alfven wave without gyroviscosity, which keeps its polarization.'''
    assert gyro_case(-1., linear=linear)
    assert validate_result(0.0, err_exp)

# Gyroviscous rotation of the polarization
@pytest.mark.parametrize("linear,err_exp", ((False, 1.767624E-03), pytest.param(True, 2.074131E-03, marks=pytest.mark.linear)))
def test_gyro_rotation(linear, err_exp):
    '''Shear Alfven wave with Braginskii's gyroviscous stress (S. I. Braginskii, Rev. Plasma Phys. 1, 205 (1965)),
    whose polarization rotates from x toward y (B along +z) at half the split between the circular polarizations,
    eta_4 k^2/rho = k^2 T_i/B0. The opposite sign of the stress would rotate it the other way.'''
    assert gyro_case(1., linear=linear)
    assert validate_result(1.0, err_exp)

#============================================================================
# Gyroviscosity with the Hall term (Hall-only cases are in test_MUG_hall_lag.py)
@pytest.mark.parametrize("linear,err_exp", ((False, 3.098049E-03), pytest.param(True, 3.009142E-03, marks=pytest.mark.linear)))
def test_hall_gyro(linear, err_exp):
    '''Shear Alfven wave with both gyroviscosity and the Hall term (default electron mass factor of 100). Their
    rotation rates are in the ratio beta_i/2 ~ 0.96 and, with the conventions above, add, which checks the relative
    sign of the two terms.'''
    assert gyro_case(1., hall=True, linear=linear)
    assert validate_result(1.0, err_exp, hall=True)
