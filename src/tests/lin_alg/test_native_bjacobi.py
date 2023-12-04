from __future__ import print_function
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

nSolvers = 5

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
 nbase=3
 grid_order=1
/

&cube_options
 mesh_type=1
/

&test_bj_options
 nlocal={1}
 sol_type={2}
/
"""

# Common setup function and process handling
def native_bjacobi_setup(nprocs=1, nlocal=1, sol_type=1):
    nlevels = 3
    if nprocs > 1:
        nlevels = 4
    #
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nlevels, nlocal, sol_type))
    return run_OFT("./test_native_bjacobi", nprocs, 60)

def validate_result(iteration_count,converged_error):
    #
    fid = open('bjacobi.results','r')
    its_test = int(fid.readline())
    if iteration_count != None:
        if its_test>iteration_count+1:
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
#
def check_solver(nprocs,nlocal,sol_type,iteration_count,converged_error):
    assert native_bjacobi_setup(nprocs,nlocal,sol_type)
    assert validate_result(iteration_count,converged_error)

# Test runner for base test case
@pytest.mark.parametrize("solver_type", range(nSolvers))
def test_base(solver_type):
    iteration_count = 1
    converged_error = 1.4688646289584524
    check_solver(1, 1, solver_type+1, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", range(nSolvers))
def test_base_nlocal2(solver_type):
    iteration_count = 42
    converged_error = 1.4688646289584524
    check_solver(1, 2, solver_type+1, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", range(nSolvers))
def test_base_nlocal2_part(solver_type):
    iteration_count = 47
    converged_error = 1.4688646289584524
    check_solver(1, -2, solver_type+1, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", range(nSolvers))
def test_mpi2(solver_type):
    iteration_count = 50
    converged_error = 1.4688646289584524
    check_solver(2, 1, solver_type+1, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", range(nSolvers))
def test_mpi2_nlocal2(solver_type):
    iteration_count = 57
    converged_error = 1.4688646289584524
    check_solver(2, 2, solver_type+1, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", range(nSolvers))
def test_mpi2_nlocal2_part(solver_type):
    iteration_count = 57
    converged_error = 1.4688646289584524
    check_solver(2, -2, solver_type+1, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", range(nSolvers))
def test_mpi4(solver_type):
    iteration_count = 59
    converged_error = 1.4688646289584524
    check_solver(4, 1, solver_type+1, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", range(nSolvers))
def test_mpi4_nlocal2(solver_type):
    iteration_count = 65
    converged_error = 1.4688646289584524
    check_solver(4, 2, solver_type+1, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", range(nSolvers))
def test_mpi4_nlocal2_part(solver_type):
    iteration_count = 65
    converged_error = 1.4688646289584524
    check_solver(4, -2, solver_type+1, iteration_count, converged_error)
