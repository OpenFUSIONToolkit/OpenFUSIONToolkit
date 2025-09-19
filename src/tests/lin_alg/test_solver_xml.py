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
 use_petsc={0}
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels=3
 nbase=3
 grid_order=1
/

&cube_options
 mesh_type=1
/
"""
#
xml_template = """
<oft>
  {0}
</oft>
"""

# Common setup function and process handling
def xml_setup(xml_file, petsc=False):
    os.chdir(test_dir)
    petsc_flag=('T' if petsc else 'F')
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(petsc_flag))
    with open('oft_in.xml','w+') as fid:
        fid.write(xml_template.format(xml_file))
    return run_OFT("./test_solver_xml oft.in oft_in.xml", 1, 60)

def validate_result(iteration_count,converged_error):
    #
    fid = open('xml.results','r')
    its_test = int(fid.readline())
    if iteration_count is not None:
        if abs(iteration_count-its_test) >= max(3,0.05*iteration_count):
            print("FAILED: Iteration count incorrect!")
            print("  Expected = {0}".format(iteration_count))
            print("  Actual =   {0}".format(its_test))
            return False
    error_test = float(fid.readline())
    if abs(converged_error-error_test) > 1.E-6:
        print("FAILED: Residual error incorrect!")
        print("  Expected = {0}".format(converged_error))
        print("  Actual =   {0}".format(error_test))
        return False
    return True

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_gmres(petsc_flag):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
    </solver>
    """
    iteration_count = 530
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_gmres_jacobi(petsc_flag):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="jacobi"></pre>
    </solver>
    """
    iteration_count = 492
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_cg(petsc_flag):
    xml_file = """
    <solver type="cg">
        <its>-3</its>
    </solver>
    """
    iteration_count = 145
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_cg_jacobi(petsc_flag):
    xml_file = """
    <solver type="cg">
        <its>-3</its>
        <pre type="jacobi"></pre>
    </solver>
    """
    iteration_count = 154
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)
