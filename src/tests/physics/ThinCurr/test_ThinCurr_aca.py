#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThinCurr tests: high-resolution cylinder, dense and with ACA+ (HODLR) matrix compression
import pytest
from thincurr_testing import ThinCurr_setup, validate_eigs, validate_fr, validate_td
from oft_testing import approx_list


# HODLR compression of the inductance matrix with ACA+ (Grasedyck, Computing 74, 205 (2005)) on a cylinder fine
# enough (>2000 elements) for the compressed and the block-Jacobi preconditioned dense solvers to be used
@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_eig_aca(convert_xml):
    '''L/R eigenmodes of the high-resolution cylinder with a HODLR-compressed inductance matrix.'''
    eigs = approx_list((2.659575E-2, 1.254552E-2, 1.254536E-2, 1.208636E-2), rel=1.E-5)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'eig','F',use_aca=True,jumper_start=2,convert_xml=convert_xml)
    assert validate_eigs(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_aca(convert_xml):
    '''Time-dependent eddy currents in the high-resolution cylinder with HODLR compression, including a reduced model
    built from 30 eigenmodes that must reproduce the slowest eigenvalues and run on its own.'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((7.205529E-4, 6.100309E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([5.453372E3, 5.453372E3], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'td','F',use_aca=True,convert_xml=convert_xml,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           jumper_start=2,run_reduced=True)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_fr_cyl_hr(direct_flag):
    '''Frequency response of the dense high-resolution cylinder, where the iterative solve uses the block-Jacobi
    preconditioner.'''
    # Dense model with >2000 elements, so the iterative solve uses the block-Jacobi preconditioner
    fr_real = approx_list((5.884322E-2, 4.882137E-3), rel=1.E-4)
    fr_imag = approx_list((-2.016855E-3, -2.314362E-3), rel=1.E-4)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'fr',direct_flag,freq=5.E3,fr_limit=0,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           jumper_start=2)
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_fr_aca(convert_xml):
    '''Frequency response of the high-resolution cylinder with HODLR compression.'''
    fr_real = approx_list((5.888736E-2, 4.881440E-3), rel=1.E-3)
    fr_imag = approx_list((-2.017045E-3, -2.313881E-3), rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'fr','F',use_aca=True,freq=5.E3,fr_limit=0,convert_xml=convert_xml,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           jumper_start=2)
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.coverage
def test_td_aca_single():
    '''The HODLR time-dependent case with single-precision block-Jacobi preconditioners.'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((7.205529E-4, 6.100309E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([5.453372E3, 5.453372E3], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'td','F',use_aca=True,hodlr_opts={'pre_single': True},
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.coverage
def test_fr_aca_single():
    '''The HODLR frequency-response case with single-precision block-Jacobi preconditioners.'''
    fr_real = approx_list((5.888736E-2, 4.881440E-3), rel=1.E-3)
    fr_imag = approx_list((-2.017045E-3, -2.313881E-3), rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'fr','F',use_aca=True,hodlr_opts={'pre_single': True},freq=5.E3,fr_limit=0,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           jumper_start=2)
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_volt_aca(convert_xml):
    '''Voltage-driven time-dependent response of the high-resolution cylinder with HODLR compression.'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((1.823946E-4, 1.560471E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([1.348327E3, 1348.327, 788.7134, 788.7133], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'td','F',use_aca=True,convert_xml=convert_xml,
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.parametrize("pre_single,pre_neigs", ((False, 0), (True, 0), (False, 2)))
def test_td_volt_aca_coarse(pre_single,pre_neigs):
    '''The voltage-driven HODLR case with the two-level (coarse-corrected) preconditioner, optionally in single
    precision or with the slowest eigenmodes in the coarse space.'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((1.823946E-4, 1.560471E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([1.348327E3, 1348.327, 788.7134, 788.7133], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'td','F',use_aca=True,hodlr_opts={'pre_coarse': True, 'pre_single': pre_single,
                                                                          'pre_neigs': pre_neigs},
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)
