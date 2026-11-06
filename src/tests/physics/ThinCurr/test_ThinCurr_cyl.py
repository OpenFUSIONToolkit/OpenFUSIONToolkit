#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThinCurr tests: cylinder eigenmodes, time-dependent, and frequency-response cases
import pytest
from thincurr_testing import ThinCurr_setup, circular_coil_defs, validate_eigs, validate_fr, validate_td
from oft_testing import approx_list


@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_eig_cyl(direct_flag,convert_xml):
    '''The four slowest L/R eigenmodes of a cylindrical shell with jumper regions, with the iterative and the direct
    eigensolver, and the resistivity read from XML or converted into the mesh file.'''
    eigs = approx_list((2.657195E-2, 1.248071E-2, 1.247103E-2, 1.200566E-2), rel=1.E-5)
    assert ThinCurr_setup("tw_test-cyl.h5",'eig',direct_flag,jumper_start=2,convert_xml=convert_xml)
    assert validate_eigs(eigs)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_cyl(direct_flag,convert_xml):
    '''Time-dependent eddy currents in a cylindrical shell with jumper regions driven by a ramped current in two
    circular drive coils: the final flux-loop signals and jumper currents.'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((7.178595E-4, 6.040177E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([5.451791E3, 5.451791E3], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl.h5",'td',direct_flag,convert_xml=convert_xml,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_fr_cyl(direct_flag,convert_xml):
    '''Frequency response of a cylindrical shell with jumper regions at 5 kHz to two circular drive coils: real and
    imaginary parts of the flux-loop signals.'''
    fr_real = approx_list((6.113541E-2, 4.356824E-3), rel=1.E-4)
    fr_imag = approx_list((-1.911861E-3, -2.283967E-3), rel=1.E-4)
    assert ThinCurr_setup("tw_test-cyl.h5",'fr',direct_flag,freq=5.E3,fr_limit=0,convert_xml=convert_xml,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           jumper_start=2)
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_cyl_volt(direct_flag,convert_xml):
    '''Time-dependent response of a cylindrical shell with jumper regions to voltage-driven passive coils: the final
    flux-loop signals and jumper currents and coil currents.'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((1.813054E-4, 1.542203E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([1.3453041E3, 1345.304, 787.155, 787.150], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl.h5",'td',direct_flag,convert_xml=convert_xml,
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("coil_def", ('xml', 'mesh', 'converted'))
def test_td_cyl_volt_circular(direct_flag,coil_def):
    '''The voltage-driven cylinder case with the coils given as circular filaments, defined in the XML file, in the
    mesh file, or converted from point lists by the mesh tool.'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((1.813188E-4, 1.541818E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([1.344134E3, 1344.134, 786.6749, 786.6708], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl.h5",'td',direct_flag,
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           jumper_start=2,**circular_coil_defs[coil_def])
    assert validate_td(sigs_final,jumpers_final)
