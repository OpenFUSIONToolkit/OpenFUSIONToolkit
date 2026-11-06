#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThinCurr tests: torus eigenmodes, time-dependent, and frequency-response cases
import pytest
from thincurr_testing import ThinCurr_setup, circular_coil_defs, validate_eigs, validate_fr, validate_td
from oft_testing import approx_list


@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_eig_torus(direct_flag,convert_xml):
    '''Four slowest L/R eigenvalues of a closed torus, with its two topological holes, using direct and iterative
    solvers and XML or mesh-file inputs'''
    eigs = approx_list((4.751344E-2, 2.564491E-2, 2.555695E-2, 2.285850E-2), rel=1.E-5)
    assert ThinCurr_setup("tw_test-torus.h5",'eig',direct_flag,convert_xml=convert_xml)
    assert validate_eigs(eigs)


@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_torus(direct_flag,convert_xml):
    '''Time-dependent response of flux loops inside and outside the torus to a pair of ramped coils'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((4.775732E-4, 3.408042E-5), rel=1.E-3)
    assert ThinCurr_setup("tw_test-torus.h5",'td',direct_flag,convert_xml=convert_xml,
                           icoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)),
                           lin_tol=1.E-11)
    assert validate_td(sigs_final)


@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("coil_def", ('xml', 'mesh', 'converted'))
def test_td_torus_circular(direct_flag,coil_def):
    '''Torus time-dependent case with exact circular coils defined in XML, in the mesh file, or converted from point
    lists'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((4.778715E-4, 3.407226E-5), rel=1.E-3)
    assert ThinCurr_setup("tw_test-torus.h5",'td',direct_flag,
                           icoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)),
                           lin_tol=1.E-11,**circular_coil_defs[coil_def])
    assert validate_td(sigs_final)


@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_fr_torus(direct_flag,convert_xml):
    '''Frequency response at 5 kHz of flux loops to a pair of coils around the torus'''
    fr_real = approx_list((-2.807268E-3, -1.196091E-4), rel=1.E-4)
    fr_imag = approx_list((-1.870323E-3, -1.248355E-4), rel=1.E-4)
    assert ThinCurr_setup("tw_test-torus.h5",'fr',direct_flag,freq=5.E3,fr_limit=0,convert_xml=convert_xml,
                           icoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)))
    assert validate_fr(fr_real, fr_imag)


@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_torus_volt(direct_flag,convert_xml):
    '''Time-dependent response to a pair of voltage-driven coils around the torus, including the recorded currents'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((7.423015E-5, 5.296107E-6), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([None, -780.7178, 485.56802, 485.56820], rel=1.E-3)
    assert ThinCurr_setup("tw_test-torus.h5",'td',direct_flag,convert_xml=convert_xml,
                           vcoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           lin_tol=1.E-11)
    assert validate_td(sigs_final,jumpers_final)
