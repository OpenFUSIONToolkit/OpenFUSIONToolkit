import os
import sys
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
/

&test_tracing_options
 tracer_type={2}
 test_poincare={4}
/
"""

# Common setup function and process handling
def tracing_setup(tracer_type, parallel=False, mesh_type=1, poincare=False):
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
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nlevels, nbase, tracer_type, mesh_type, poincare_flag))
    return run_OFT("./test_tracing", nproc, 120)

# Tracer type 1 is LSODE (Adams method), type 4 is the cell-limited
# Bogacki-Shampine 3(2) Runge-Kutta tracer
# Mesh type 1 is tetrahedral, 2 is hexahedral (different cell exit tests for type 4)
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_tracing(tracer_type, parallel, mesh_type):
    assert tracing_setup(tracer_type, parallel, mesh_type)

# Poincare sections of several field lines traced concurrently (tracing_poincare),
# which spawns per-line tracers from the template using the tracer's copy method
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_poincare(tracer_type, parallel, mesh_type):
    assert tracing_setup(tracer_type, parallel, mesh_type, poincare=True)
