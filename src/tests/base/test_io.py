from __future__ import print_function
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