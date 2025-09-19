from __future__ import print_function
import os
import sys
import time
import multiprocessing
import pytest
import numpy as np
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..','..','python')))
from oft_testing import run_command
from OpenFUSIONToolkit._interface import oftpy_dump_cov


def mp_run(target,args,timeout=180):
    if os.environ.get('OFT_DEBUG_TEST', 0):
        timeout *= 4
    os.chdir(test_dir)
    mp_q = multiprocessing.Queue()
    p = multiprocessing.Process(target=target, args=args + (mp_q,))
    p.start()
    start = time.time()
    while time.time() - start <= timeout:
        if not p.is_alive():
            break
        time.sleep(.5)
    else: # Reached timeout
        print("Timeout reached")
        p.terminate()
        p.join()
        return None
    # Completed successfully?
    try:
        test_result = mp_q.get(timeout=5)
    except:
        print("Failed to get output")
        return None
    p.join()
    return test_result


def run_marklin(meshfile,nmodes,order,grid_order,meshfile2,mp_q):
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
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.Marklin import Marklin
        myOFT = OFT_env(nthreads=-1)
        taylor_solver = Marklin(myOFT)
        taylor_solver.setup_mesh(mesh_file=meshfile,grid_order=grid_order)
        taylor_solver.setup(order,minlev=1)
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


# Test runners for cylinder taylor states (straight mesh)
def marklin_lin_cell(order):
    eigs = [5.029946,4.958005,4.957746,4.957744]
    if order == 1:
        run_command("rm -f *.rst", cwd=test_dir)
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,order,1,None))
    assert validate_eigs(results,eigs[order-1])

@pytest.mark.coverage
def test_marklin_g1_p1():
    marklin_lin_cell(1)

@pytest.mark.coverage
def test_marklin_g1_p2():
    marklin_lin_cell(2)

def test_marklin_g1_p3():
    marklin_lin_cell(3)

def test_marklin_g1_p4():
    marklin_lin_cell(4)

# Test runners for cylinder taylor states (curved mesh)
def marklin_quad_cell(order):
    eigs = [5.027066,4.955207,4.954955,4.954955]
    if order == 1:
        run_command("rm -f *.rst", cwd=test_dir)
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,order,2,None))
    assert validate_eigs(results,eigs[order-1])

@pytest.mark.coverage
def test_marklin_g2_p1():
    marklin_quad_cell(1)

@pytest.mark.coverage
def test_marklin_g2_p2():
    marklin_quad_cell(2)

def test_marklin_g2_p3():
    marklin_quad_cell(3)

def test_marklin_g2_p4():
    marklin_quad_cell(4)

# Test runner for multiple concurrent Marklin instances
def test_marklin_multiple():
    cyl_eig = 4.955207
    tor_eig = 6.317001
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,2,2,'torus_test.h5'))
    assert validate_eigs(results,tor_eig)
    results = mp_run(run_marklin,('torus_test.h5',1,2,2,'cyl_Marklin.h5'))
    assert validate_eigs(results,cyl_eig)

# Test runners for torus vacuum field (curved mesh)
def marklin_torus_vac(order,mag_val):
    if order == 1:
        run_command("rm -f *.rst", cwd=test_dir)
    results = mp_run(run_marklin_vac,('torus_test.h5',order,2))
    print(results)
    assert abs((results-mag_val)/mag_val) < 1.E-4

@pytest.mark.coverage
def test_marklin_vac_g2_p1():
    field_mag = 1.145501
    marklin_torus_vac(1,field_mag)

@pytest.mark.coverage
def test_marklin_vac_g2_p2():
    field_mag = 1.188539
    marklin_torus_vac(2,field_mag)

def test_marklin_vac_g2_p3():
    field_mag = 1.188352
    marklin_torus_vac(3,field_mag)

def test_marklin_vac_g2_p4():
    field_mag = 1.188126
    marklin_torus_vac(4,field_mag)

# # Example of how to run single test without pytest
# if __name__ == '__main__':
#     multiprocessing.freeze_support()
#     mp_q = multiprocessing.Queue()
#     run_command("rm -f *.rst")
#     run_marklin('cyl_Marklin.h5',1,2,2,'torus_test.h5',mp_q)