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
 debug={1}
 test_run=T
/

&test_replicated_options
 n={0}
/
"""

# Common setup function and process handling
def replicated_setup(nproc, n=1000, debug=1):
    os.chdir(test_dir)
    if os.path.isfile('replicated.results'):
        os.remove('replicated.results')
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(n, debug))
    return run_OFT("./test_replicated_seam", nproc, 60)

def validate_result():
    if not os.path.isfile('replicated.results'):
        print("FAILED: No results file")
        return False
    with open('replicated.results','r') as fid:
        result = fid.readline().strip()
    if result != 'PASS':
        print("FAILED: {0}".format(result))
        return False
    return True

#============================================================================
# Serial test (replicated layout with a single owner)
@pytest.mark.coverage
def test_serial():
    assert replicated_setup(1)
    assert validate_result()

# Parallel tests (power-of-two and non-power-of-two processor counts)
@pytest.mark.coverage
@pytest.mark.parametrize("nproc", (2, 3))
def test_mpi(nproc):
    assert replicated_setup(nproc)
    assert validate_result()

# Parallel test with more processors than some owners have entries
@pytest.mark.coverage
def test_mpi_small():
    assert replicated_setup(3, n=5)
    assert validate_result()
