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
 mesh_type=1
 rscale=2.d0,2.d0,2.d0
 shift=-1.d0,-1.d0,-1.d0
/

&test_tracing_options
 tracer_type={2}
/
"""

# Common setup function and process handling
def tracing_setup(tracer_type, parallel=False):
    if parallel:
        nlevels = 4
        nbase = 3
        nproc = 2
    else:
        nlevels = 3
        nbase = 3
        nproc = 1
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nlevels, nbase, tracer_type))
    return run_OFT("./test_tracing", nproc, 120)

# Tracer type 1 uses CVODE when available, type 4 is always LSODE
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
def test_tracing(tracer_type, parallel):
    assert tracing_setup(tracer_type, parallel)
