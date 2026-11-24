import os
import sys
import pytest
import numpy as np
import h5py
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT
from OpenFUSIONToolkit.io import XDMF_plot_file

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={0}
 nbase={1}
 grid_order=1
/

&cube_options
 mesh_type=1
/

&test_io_options
 test_id={2}
/
"""

# Common setup function and process handling
def io_setup(test_id, parallel = False):
    nbase = 3
    if parallel:
        nlevels = 4
        nproc = 2
    else:
        nlevels = 3
        nproc = 1
    #
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nlevels,nbase,test_id))
    return run_OFT("./test_io", nproc, 60)

# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("test_id", (1,2,3,4,5))
@pytest.mark.parametrize("parallel", (False,True))
def test_base(test_id,parallel):
    assert io_setup(test_id,parallel)


# Test for time interpolation of plot file fields
@pytest.mark.coverage
def test_xdmf_get_field_time(tmp_path):
    '''Check that `XDMF_plot_mesh.get_field(name, time=t)` linearly interpolates between saved timesteps
    and returns the stored field exactly at saved times (issue #460).'''
    # One-triangle plot file with field "F" saved at uneven times, so interior intervals are exercised
    times = (0.0, 1.0, 3.0)
    values = (0.0, 1.0, 5.0)
    filename = str(tmp_path / 'interp_test.h5')
    with h5py.File(filename, 'w') as fid:
        mesh = fid.create_group('test/mesh')
        mesh['TYPE'] = 31
        mesh['R'] = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        mesh['LC'] = np.array([[0, 1, 2]], dtype=np.int32)
        for i, (time, value) in enumerate(zip(times, values)):
            step = mesh.create_group('{0:04d}'.format(i+1))
            step['TIME'] = np.array([time])
            step['F'] = np.array([value])
    plot_mesh = XDMF_plot_file(filename)['test']['mesh']
    # Saved times must return the stored data exactly
    for time, value in zip(times, values):
        assert plot_mesh.get_field('F', time=time)[0] == pytest.approx(value)
    # Interior points are linear interpolants of the bracketing timesteps
    for time in (0.25, 0.5, 2.0, 2.5):
        assert plot_mesh.get_field('F', time=time)[0] == pytest.approx(np.interp(time, times, values))