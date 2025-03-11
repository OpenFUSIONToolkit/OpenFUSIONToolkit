from __future__ import print_function
import os
import sys
import time
import multiprocessing
import pytest
import numpy as np
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..','..','python')))


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
    test_results = []
    try:
        for _ in range(3):
            result = mp_q.get(timeout=5)
            test_results.append(result)
    except:
        print("Failed to get output")
        return None
    p.join()
    return test_results


def run_all(dummy,mp_q):
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.Marklin import Marklin
    from OpenFUSIONToolkit.ThinCurr import ThinCurr
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
    myOFT = OFT_env(nthreads=-1)
    # Run Marklin
    try:
        taylor_solver = Marklin(myOFT)
        taylor_solver.setup_mesh(mesh_file='cyl_Marklin.h5',grid_order=1)
        taylor_solver.setup(2,minlev=1)
        taylor_solver.compute(1)
        mp_q.put(taylor_solver.eig_vals)
    except:
        mp_q.put(None)
        mp_q.put(None)
        mp_q.put(None)
        return
    # Run ThinCurr
    try:
        tw_model = ThinCurr(myOFT)
        tw_model.setup_model(mesh_file="tw_test-plate.h5",xml_filename='oft_in.xml')
        tw_model.setup_io()
        tw_model.compute_Mcoil()
        tw_model.compute_Lmat()
        tw_model.compute_Rmat()
        eig_vals, _ = tw_model.get_eigs(4,direct=True)
        mp_q.put(eig_vals[0])
    except:
        mp_q.put(None)
        mp_q.put(None)
        return
    # Run TokaMaker
    try:
        gs_mesh = gs_Domain()
        gs_mesh.define_region('plasma',0.05,'plasma')
        gs_mesh.add_rectangle(0.5,0.5,1.0,1.0,'plasma')
        mesh_pts, mesh_lc, _ = gs_mesh.build_mesh()
        mygs = TokaMaker(myOFT)
        mygs.setup_mesh(mesh_pts,mesh_lc)
        mygs.settings.free_boundary = False
        mygs.setup(order=2)
        mygs.pnorm=0.0
        ffp_prof={
            'type': 'linterp',
            'x': [0.0,1.0],
            'y': [1.0,0.0],
        }
        mygs.set_profiles(ffp_prof=ffp_prof,pp_prof={'type': 'flat'})
        mygs.settings.nl_tol = 1.E-12
        mygs.settings.maxits = 100
        mygs.urf = 0.0
        mygs.update_settings()
        mygs.init_psi()
        mygs.solve()
        mp_q.put(mygs.alam)
    except:
        mp_q.put(None)
        return


def validate_all(results, expected, rtol=1.E-3):
    if len(results) != 3:
        print("FAILED: Number of results does not match expected")
        return False
    retval = True
    for i in range(3):
        if abs((results[i]-expected[i])/expected[i]) > rtol:
            print("FAILED: Result {0} outside tolerance!".format(i))
            print("  Expected = {0}".format(expected[i]))
            print("  Actual =   {0}".format(results[i]))
            retval = False
    return retval

@pytest.mark.coverage
def test_packages_serial():
    with open("oft_in.xml", "w+") as fid:
        fid.write("""<oft>
  <thincurr>
    <eta>1.256637E-5</eta>
  </thincurr>
</oft>
""")
    expected = [4.95800457, 9.735667E-3, 49.131116]
    results = mp_run(run_all,(None,))
    assert validate_all(results, expected)