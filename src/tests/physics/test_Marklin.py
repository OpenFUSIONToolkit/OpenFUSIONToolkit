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
from OpenFUSIONToolkit.util import oftpy_dump_cov


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


def run_marklin(meshfile,nmodes,order,grid_order,mp_q):
    try:
        from OpenFUSIONToolkit.Marklin import Marklin
        taylor_solver = Marklin(grid_order=grid_order)
        taylor_solver.setup_mesh(mesh_file=meshfile)
        taylor_solver.compute(nmodes,order,minlev=1)
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
        print("FAILED: psi error too high!")
        print("  Expected = {0}".format(eigs))
        print("  Actual =   {0}".format(results[0]))
        test_result = False
    return test_result

# Test runners for cylinder taylor states (straight mesh)
def marklin_lin_cell(order):
    eigs = [5.029946,4.958005,4.957746,4.957744]
    if order == 1:
        run_command("rm -f *.rst")
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,order,1))
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
        run_command("rm -f *.rst")
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,order,2))
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
