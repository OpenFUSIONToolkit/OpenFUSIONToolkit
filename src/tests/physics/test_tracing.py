#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Field line tracing tests (3D tracers on an analytic helical field)
import os
import sys
import h5py
import numpy as np
import pytest
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
 mesh_type={3}
 rscale=2.d0,2.d0,2.d0
 shift=-1.d0,-1.d0,-1.d0
 ref_per={7}
/

&test_tracing_options
 tracer_type={2}
 test_poincare={4}
 field_jump={5}
 tol={6}
/
"""

# Common setup function and process handling
def tracing_setup(tracer_type, parallel=False, mesh_type=1, poincare=False, field_jump=0.0, tol=1.E-10, periodic=False):
    # Both give the same final mesh, as the level where the mesh is distributed is not refined
    if parallel:
        nlevels = 4
        nbase = 3
        nproc = 2
    else:
        nlevels = 3
        nbase = 3
        nproc = 1
    poincare_flag = 'F'
    if poincare:
        poincare_flag = 'T'
    ref_per = 'F,F,F'
    if periodic:
        ref_per = 'F,F,T'
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nlevels, nbase, tracer_type, mesh_type, poincare_flag, field_jump, tol, ref_per))
    return run_OFT("./test_tracing", nproc, 120)


def jump_trace(parallel, mesh_type, tol):
    '''Trace a field line with jumps in the field strength between cells and return its points'''
    assert tracing_setup(4, parallel, mesh_type, field_jump=0.5, tol=tol)
    with h5py.File('trace_test', 'r') as fid:
        return np.asarray(fid['trace'])

# Tracer type 1 is LSODE (Adams method), type 4 is the cell-limited
# Bogacki-Shampine 3(2) Runge-Kutta tracer
# Mesh type 1 is tetrahedral, 2 is hexahedral (different cell exit tests for type 4)
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_tracing(tracer_type, parallel, mesh_type):
    '''Single field line traced with tracing_line, checked against the analytic helical field line'''
    assert tracing_setup(tracer_type, parallel, mesh_type)

# Poincare sections of several field lines traced concurrently (tracing_poincare),
# which spawns per-line tracers from the template using the tracer's copy method
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_poincare(tracer_type, parallel, mesh_type):
    '''Poincare section of several field lines (tracing_poincare), checked against the analytic crossings'''
    assert tracing_setup(tracer_type, parallel, mesh_type, poincare=True)

# Mesh that is periodic in z: field lines are not continued across periodic faces, so the
# line stops at z=1 as at the mesh boundary
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_tracing_periodic(tracer_type, parallel, mesh_type):
    '''Single field line on a mesh that is periodic in z, which stops at the periodic face'''
    assert tracing_setup(tracer_type, parallel, mesh_type, periodic=True)


# Field whose strength jumps at every cell face, as for finite element fields. The field line
# is still the exact helix (checked by the test program), and the time along it (fourth
# column) is a piecewise linear function of the angle, which is compared to a serial trace with
# a tight tolerance on the same mesh. This covers the retries of steps that overshoot a face
# too far and, in parallel, the steps that end on another domain.
@pytest.mark.coverage
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_tracing_jump(parallel, mesh_type):
    '''Field line through a field that jumps between cells, checked against a tight-tolerance reference'''
    ref = jump_trace(False, mesh_type, 1.E-11)
    trace = jump_trace(parallel, mesh_type, 1.E-8)
    phase_ref = np.unwrap(np.arctan2(ref[:,1], ref[:,0]))
    phase = np.unwrap(np.arctan2(trace[:,1], trace[:,0]))
    inside = phase <= phase_ref[-1]
    err = np.abs(trace[inside,3] - np.interp(phase[inside], phase_ref, ref[:,3]))
    print('Max time error = {0:.3E}, points = {1}'.format(err.max(), inside.sum()))
    assert err.max() < 2.E-6
