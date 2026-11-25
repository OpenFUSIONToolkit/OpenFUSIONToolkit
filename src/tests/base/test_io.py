import os
import sys
import pytest
import numpy as np
import h5py
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

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


def _write_plot_file(filename, times, values):
    # Write a one-triangle plot file with field "F" saved at the given times (no timesteps if times is empty)
    with h5py.File(filename, 'w') as fid:
        mesh = fid.create_group('test/mesh')
        mesh['TYPE'] = 31
        mesh['R'] = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        mesh['LC'] = np.array([[0, 1, 2]], dtype=np.int32)
        static = mesh.create_group('0000')
        static['S'] = np.array([7.0])
        for i, (time, value) in enumerate(zip(times, values)):
            step = mesh.create_group('{0:04d}'.format(i+1))
            step['TIME'] = np.array([time])
            step['F'] = np.array([value])


# Test for time interpolation of plot file fields
@pytest.mark.coverage
def test_xdmf_get_field_time(tmp_path):
    '''Check that `XDMF_plot_mesh.get_field(name, time=t)` linearly interpolates between saved timesteps
    and returns a copy of the stored field at saved times (issue #460).'''
    # Import here, not at module level, so collecting the Fortran-driven tests above does not load liboftpy
    from OpenFUSIONToolkit.io import XDMF_plot_file
    # One-triangle plot file with field "F" saved at uneven times, so interior intervals are exercised
    times = (0.0, 1.0, 3.0)
    values = (0.0, 1.0, 5.0)
    filename = str(tmp_path / 'interp_test.h5')
    _write_plot_file(filename, times, values)
    plot_mesh = XDMF_plot_file(filename)['test']['mesh']
    # Saved times must return the stored data exactly
    for time, value in zip(times, values):
        assert plot_mesh.get_field('F', time=time)[0] == pytest.approx(value)
    # Interior points are linear interpolants of the bracketing timesteps
    for time in (0.25, 0.5, 2.0, 2.5):
        assert plot_mesh.get_field('F', time=time)[0] == pytest.approx(np.interp(time, times, values))
    # Results at saved times are copies, so modifying them does not change the stored data
    result = plot_mesh.get_field('F', time=1.0)
    result *= 2.0
    assert plot_mesh.get_field('F', time=1.0)[0] == pytest.approx(1.0)
    assert plot_mesh.get_field('F', time=2.0)[0] == pytest.approx(3.0)


@pytest.mark.coverage
def test_xdmf_get_field_timestep_range(tmp_path):
    '''`timestep` accepts Python-style negative indices; values outside [-N, N-1] raise ValueError.'''
    from OpenFUSIONToolkit.io import XDMF_plot_file as plot_file
    filename = str(tmp_path / 'timestep_range.h5')
    _write_plot_file(filename, (0.0, 1.0, 2.0), (3.0, 4.0, 5.0))
    plot_mesh = plot_file(filename)['test']['mesh']
    for timestep, value in enumerate((3.0, 4.0, 5.0)):
        assert plot_mesh.get_field('F', timestep=timestep)[0] == value
        assert plot_mesh.get_field('F', timestep=timestep-3)[0] == value
    for timestep in (-4, -10, 3, 10):
        with pytest.raises(ValueError):
            plot_mesh.get_field('F', timestep=timestep)


@pytest.mark.coverage
def test_xdmf_get_field_duplicate_times(tmp_path):
    '''Repeated TIME values do not produce NaN and return the first stored timestep with that time.'''
    from OpenFUSIONToolkit.io import XDMF_plot_file as plot_file
    filename = str(tmp_path / 'duplicate_times.h5')
    times = (0.0, 0.0, 1.0, 1.0, 2.0)
    values = (1.0, 2.0, 3.0, 4.0, 6.0)
    _write_plot_file(filename, times, values)
    plot_mesh = plot_file(filename)['test']['mesh']
    for time, expected in ((0.0, 1.0), (1.0, 3.0), (2.0, 6.0), (0.5, 2.5), (1.5, 5.0)):
        result = plot_mesh.get_field('F', time=time)
        assert np.all(np.isfinite(result))
        assert result[0] == pytest.approx(expected)
    # Single timestep: the stored time is in range and returns the stored data
    filename = str(tmp_path / 'single_time.h5')
    _write_plot_file(filename, (0.5,), (9.0,))
    plot_mesh = plot_file(filename)['test']['mesh']
    assert plot_mesh.get_field('F', time=0.5)[0] == 9.0
    with pytest.raises(ValueError):
        plot_mesh.get_field('F', time=0.6)


@pytest.mark.coverage
def test_xdmf_get_field_static_only(tmp_path):
    '''A mesh with no timesteps raises a clear ValueError for `time` and `timestep` requests.'''
    from OpenFUSIONToolkit.io import XDMF_plot_file as plot_file
    filename = str(tmp_path / 'static_only.h5')
    _write_plot_file(filename, (), ())
    plot_mesh = plot_file(filename)['test']['mesh']
    assert plot_mesh.get_field('S')[0] == 7.0
    for kwargs in ({'time': 0.0}, {'timestep': 0}):
        for name in ('S', 'F'):
            with pytest.raises(ValueError, match='no time-dependent data'):
                plot_mesh.get_field(name, **kwargs)


@pytest.mark.coverage
def test_xdmf_get_field_decreasing_times(tmp_path):
    '''Decreasing TIME values load and allow `timestep` access, but `time` interpolation raises ValueError.'''
    from OpenFUSIONToolkit.io import XDMF_plot_file as plot_file
    filename = str(tmp_path / 'decreasing_times.h5')
    _write_plot_file(filename, (0.0, 2.0, 1.0), (0.0, 2.0, 1.0))
    plot_mesh = plot_file(filename)['test']['mesh']
    assert plot_mesh.get_field('F', timestep=2)[0] == 1.0
    for time in (0.5, 1.0, 1.5):
        with pytest.raises(ValueError, match='decrease between timesteps 1'):
            plot_mesh.get_field('F', time=time)
