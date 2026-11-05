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

&test_cvode_options
 test_id={0}
/
"""

# Common setup function and process handling
def cvode_setup(test_id):
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(test_id))
    return run_OFT("./test_cvode", 1, 60)

# Test runner for CVODE interface
@pytest.mark.coverage
@pytest.mark.parametrize("test_id", (1,2,3))
def test_cvode(test_id):
    assert cvode_setup(test_id)
