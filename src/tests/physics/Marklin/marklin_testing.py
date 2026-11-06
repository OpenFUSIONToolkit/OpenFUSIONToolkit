#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Shared helpers for Marklin Python interface tests: process runners and result validation
import os
import sys
import numpy as np
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
import oft_testing
from OpenFUSIONToolkit._interface import oftpy_dump_cov


def mp_run(target,args,timeout=180):
    '''Run `target` in its own process from the Marklin test folder (see `oft_testing.mp_run`)'''
    return oft_testing.mp_run(target,args,test_dir,timeout=timeout)


def clear_cache():
    '''Remove cached eigenmode/vacuum field files (`*.rst`) from previous runs'''
    oft_testing.run_command("rm -f *.rst", cwd=test_dir)


def run_marklin(meshfile,nmodes,order,grid_order,meshfile2,mp_q):
    '''Compute Taylor eigenmodes on `meshfile`, then on `meshfile2` (if given) with a second Marklin instance
    in the same OFT environment, returning the eigenvalues of the last instance (None on failure)'''
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.Marklin import Marklin
        myOFT = OFT_env(nthreads=-1)
        taylor_solver = Marklin(myOFT)
        taylor_solver.setup_mesh(mesh_file=meshfile,grid_order=grid_order)
        taylor_solver.setup(order,minlev=1)
        taylor_solver.compute_eig(nmodes,cache_file='Marklin_{0}.rst'.format(meshfile.split(".")[0]))
        result = True
    except BaseException as e:
        print(e)
        result = False
    if (meshfile2 is not None) and result:
        try:
            # Keep the first instance alive so both coexist in the same environment
            taylor_solver_old = taylor_solver
            taylor_solver = Marklin(myOFT)
            taylor_solver.setup_mesh(mesh_file=meshfile2,grid_order=grid_order)
            taylor_solver.setup(order,minlev=1)
            taylor_solver.compute_eig(nmodes,cache_file='Marklin_{0}.rst'.format(meshfile2.split(".")[0]))
            result = True
        except BaseException as e:
            print(e)
            result = False
    oftpy_dump_cov()
    if result:
        mp_q.put(taylor_solver.eig_vals)
    else:
        mp_q.put(None)


def validate_eigs(results,eigs):
    '''Check the lowest eigenvalue against `eigs` to a relative tolerance of 1E-4'''
    if results is None:
        print("FAILED: error in solve!")
        return False
    test_result = True
    if abs((results[0]-eigs)/eigs) > 1.E-4:
        print("FAILED: Eigenvalue error too high!")
        print("  Expected = {0}".format(eigs))
        print("  Actual =   {0}".format(results[0]))
        test_result = False
    return test_result


def run_marklin_vac(meshfile,order,grid_order,mp_q):
    '''Compute the vacuum field for one flux jump plane and return |B| at (1, 0, 0) (None on failure)'''
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.Marklin import Marklin
        myOFT = OFT_env(nthreads=-1)
        taylor_solver = Marklin(myOFT)
        taylor_solver.setup_mesh(mesh_file=meshfile,grid_order=grid_order)
        taylor_solver.setup(order,minlev=1)
        # One flux jump plane through (1, 0, 0) with normal along the toroidal direction there
        nh = 1
        hcpc = np.array([[1.0,0.0,0.0],])
        hcpv = np.array([[0.0,1.5,0.0],])
        taylor_solver.compute_vac(nh,hcpc,hcpv)
        binterp_obj = taylor_solver.get_binterp(vac_facs=np.r_[1.0,])
        b = binterp_obj.eval(np.r_[1.0,0.0,0.0])
        bnorm = np.linalg.norm(b)
        result = True
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    if result:
        mp_q.put(bnorm)
    else:
        mp_q.put(None)
