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

&test_strumpack_options
 nsize={0}
 comp_type={1}
 leaf_size={2}
 min_sep_size={3}
 rel_tol={4}
 iter_refine={5}
/
"""

# Common setup function and process handling
def strumpack_setup(nsize, comp_type, leaf_size=-1, min_sep_size=-1, rel_tol=-1.0, iter_refine='F'):
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nsize, comp_type, leaf_size, min_sep_size, rel_tol, iter_refine))
    return run_OFT("./test_strumpack", 1, 60)

def validate_result(max_error):
    with open('strumpack.results','r') as fid:
        errors = [float(fid.readline()) for _ in range(3)]
    for i, error in enumerate(errors):
        if error > max_error:
            print("FAILED: Solution error too large (solve {0})!".format(i+1))
            print("  Expected < {0}".format(max_error))
            print("  Actual =   {0}".format(error))
            return False
    return True

#============================================================================
# Test runner for sparse solves (0 -> none, 1 -> HSS, 2 -> BLR)
@pytest.mark.coverage
@pytest.mark.parametrize("comp_type", (0, 1, 2))
def test_sparse(comp_type):
    assert strumpack_setup(60, comp_type, leaf_size=16, min_sep_size=16, rel_tol=1.E-10)
    assert validate_result(1.E-8)

# Test runner for sparse solves with loose compression and iterative refinement (GMRES)
@pytest.mark.parametrize("comp_type", (1, 2))
def test_sparse_refine(comp_type):
    assert strumpack_setup(60, comp_type, leaf_size=16, min_sep_size=16, rel_tol=1.E-2, iter_refine='T')
    assert validate_result(1.E-5)
