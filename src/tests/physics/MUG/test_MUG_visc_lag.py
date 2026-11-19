#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: viscous damping of a sound wave in a triply periodic box (xMHD with all fields on Lagrange elements),
# with a uniform field that sets the direction of the parallel stress but is not advanced. Same cases, input template,
# and linear theory as test_MUG_visc.py, run with the test_MUG_visc_lag driver
from test_MUG_visc import visc_case as _visc_case, validate_result, nu_visc, b_magic

def visc_case(visc_type,nu_visc,nu_par,b_dir):
    return _visc_case(visc_type,nu_visc,nu_par,b_dir,exe='./test_MUG_visc_lag')

#============================================================================
# Isotropic viscosity
def test_kin():
    '''Sound wave damped by the `kin` viscous stress (nu grad V), which acts on the longitudinal velocity with
    coefficient nu.'''
    assert visc_case('kin',nu_visc,-1.,'0.,0.,1.')
    assert validate_result(nu_visc,6.246030554160897E-003)

def test_iso():
    '''Sound wave damped by the `iso` viscous stress (nu W, with W the rate of strain tensor), which acts on the
    longitudinal velocity with coefficient 4/3 nu.'''
    assert visc_case('iso',nu_visc,-1.,'0.,0.,1.')
    assert validate_result(4.0*nu_visc/3.0,7.594990997046444E-003)

#============================================================================
# Parallel viscosity
def test_par_aligned():
    '''Sound wave along B damped by Braginskii's parallel stress (S. I. Braginskii, Rev. Plasma Phys. 1, 205
    (1965)) with eta_0 = m_i n0 nu_par, which acts on the longitudinal velocity with coefficient 4/3 nu_par, the same
    as `iso`.'''
    assert visc_case('kin',0.,nu_visc,'0.,0.,1.')
    assert validate_result(4.0*nu_visc/3.0,7.728128491080409E-003)

def test_par_perp():
    '''Sound wave across B damped by the parallel stress, which acts on the longitudinal velocity with coefficient
    1/3 nu_par, a quarter of the field-aligned case.'''
    assert visc_case('kin',0.,nu_visc,'1.,0.,0.')
    assert validate_result(nu_visc/3.0,3.337417272638117E-003)

def test_par_magic():
    '''Sound wave at the angle to B where b.W.b = 0 (cos^2 = 1/3), so the parallel stress vanishes and the wave
    is undamped.'''
    assert visc_case('kin',0.,nu_visc,b_magic)
    assert validate_result(0.0,1.452759337776555E-003)
