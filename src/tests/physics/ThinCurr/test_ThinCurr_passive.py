#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThinCurr tests: passive filament (vcoil-only) model eigenmodes, time-dependent, and frequency-response cases
import pytest
from thincurr_testing import ThinCurr_setup, validate_eigs, validate_fr, validate_td
from oft_testing import approx_list


@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_eig_passive(direct_flag):
    '''L/R eigenvalues of passive filaments (voltage coils with no applied voltage) with a small dummy mesh placed far
    away, so the filaments carry all the current'''
    eigs = approx_list((1.489156E-1, 6.274424E-2, 3.041702E-2, 2.793962E-2), rel=1.E-5)
    assert ThinCurr_setup(None,'eig',direct_flag,eta=1.E4,
                           vcoils=((0.5, 0.1), (0.5, 0.05),
                                   (0.5, -0.05), (0.5, -0.1)))
    assert validate_eigs(eigs)


@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_eig_passive_circular(direct_flag):
    '''Filament-only eigenvalues with exact circular filaments instead of polygons'''
    eigs = approx_list((1.489295E-1, 6.274780E-2, 3.041867E-2, 2.794114E-2), rel=1.E-5)
    assert ThinCurr_setup(None,'eig',direct_flag,eta=1.E4,
                           vcoils=((0.5, 0.1), (0.5, 0.05),
                                   (0.5, -0.05), (0.5, -0.1)),circular=True)
    assert validate_eigs(eigs)


@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_passive(direct_flag):
    '''Time-dependent response of flux loops to a ramped coil coupled to a passive filament'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((7.938330E-4, 8.069765E-4), rel=1.E-3)
    assert ThinCurr_setup(None,'td',direct_flag,eta=1.E4,
                          icoils=((0.5, 0.1),),
                          vcoils=((0.5, 0.0),),
                          floops=((0.5, -0.05), (0.5, -0.1)),
                          curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)))
    assert validate_td(sigs_final)


@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_fr_passive(direct_flag):
    '''Frequency response at 5 kHz of flux loops to a coil coupled to a passive filament'''
    fr_real = approx_list((1.844588E-1, 1.917030E-1), rel=1.E-4)
    fr_imag = approx_list((-2.272322E-4, -1.629831E-4), rel=1.E-4)
    assert ThinCurr_setup(None,'fr',direct_flag,eta=1.E4,freq=5.E3,fr_limit=0,
                           icoils=((0.5, 0.1),),
                           vcoils=((0.5, 0.0),),
                           floops=((0.5, -0.05), (0.5, -0.1)))
    assert validate_fr(fr_real, fr_imag)


@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_passive_volt(direct_flag):
    '''Time-dependent response to one voltage-driven filament coupled to a passive one, including both currents'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((4.296507E-4, 4.366419E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([-673.8428, 1721.1650], rel=1.E-3)
    assert ThinCurr_setup(None,'td',direct_flag,eta=1.E4,
                          vcoils=((0.5, 0.0), (0.5, 0.1)),
                          floops=((0.5, -0.05), (0.5, -0.1)),
                          volt_waveform=((0.0, 0.0, 1.0), (1.0, 0.0, 1.0)))
    assert validate_td(sigs_final, jumpers_final)
