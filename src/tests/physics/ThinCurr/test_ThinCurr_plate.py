#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThinCurr tests: flat plate eigenmodes, time-dependent, and frequency-response cases
import pytest
from thincurr_testing import ThinCurr_setup, validate_eigs, validate_fr, validate_td
from oft_testing import approx_list


@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_eig_plate(direct_flag,convert_xml):
    '''The four slowest L/R eigenmodes of a conducting square plate, with the iterative and the direct eigensolver, and
    the resistivity read from XML or converted into the mesh file.'''
    eigs = approx_list((9.735667E-3, 6.532314E-3, 6.532201E-3, 5.251598E-3), rel=1.E-5)
    assert ThinCurr_setup("tw_test-plate.h5",'eig',direct_flag,convert_xml=convert_xml)
    assert validate_eigs(eigs)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_plate(direct_flag,convert_xml):
    '''Time-dependent eddy currents in a conducting square plate driven by a ramped current in a circular drive coil
    above it: the final flux-loop signals.'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((8.442110E-4, 7.117118E-4), rel=1.E-3)
    assert ThinCurr_setup("tw_test-plate.h5",'td',direct_flag,convert_xml=convert_xml,
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)))
    assert validate_td(sigs_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_fr_plate(direct_flag,convert_xml):
    '''Frequency response of a conducting square plate at 5 kHz to a circular drive coil above it: real and imaginary
    parts of the flux-loop signals.'''
    fr_real = approx_list((6.801291E-2, 7.203606E-2), rel=1.E-4)
    fr_imag = approx_list((-3.012114E-3, -2.177358E-3), rel=1.E-4)
    assert ThinCurr_setup("tw_test-plate.h5",'fr',direct_flag,freq=5.E3,fr_limit=0,convert_xml=convert_xml,
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)))
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_plate_volt(direct_flag,convert_xml):
    '''Time-dependent response of a conducting square plate to voltage-driven passive coils: the final flux-loop
    signals and coil currents.'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((4.699522E-4, 3.954280E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([1740.513], rel=1.E-3)
    assert ThinCurr_setup("tw_test-plate.h5",'td',direct_flag,convert_xml=convert_xml,
                           vcoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           volt_waveform=((0.0, 1.0), (1.0, 1.0)))
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.coverage
def test_td_output_dir():
    '''Restart and history files of a time-dependent run go to the requested output directory rather than the working
    directory.'''
    assert ThinCurr_setup("tw_test-plate.h5", 'td', 'F',
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           basepath='td_output_dir_test')
