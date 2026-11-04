import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import mp_run


def run_all(dummy,mp_q):
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.Marklin import Marklin
    from OpenFUSIONToolkit.ThinCurr import ThinCurr
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
    myOFT = OFT_env(nthreads=-1)
    # Results from each package, left as None if that package (or an earlier one) fails
    results = [None, None, None]
    try:
        # Run Marklin
        taylor_solver = Marklin(myOFT)
        taylor_solver.setup_mesh(mesh_file='cyl_Marklin.h5',grid_order=1)
        taylor_solver.setup(2,minlev=1)
        taylor_solver.compute_eig(1)
        results[0] = taylor_solver.eig_vals
        # Run ThinCurr
        tw_model = ThinCurr(myOFT)
        tw_model.setup_model(mesh_file="tw_test-plate.h5",xml_filename='oft_in.xml')
        tw_model.setup_io()
        tw_model.compute_Mcoil()
        tw_model.compute_Lmat()
        tw_model.compute_Rmat()
        eig_vals, _ = tw_model.get_eigs(4,direct=True)
        results[1] = eig_vals[0]
        # Run TokaMaker
        gs_mesh = gs_Domain()
        gs_mesh.define_region('plasma',0.05,'plasma')
        gs_mesh.add_rectangle(0.5,0.5,1.0,1.0,'plasma')
        mesh_pts, mesh_lc, _ = gs_mesh.build_mesh()
        mygs = TokaMaker(myOFT)
        mygs.setup_mesh(mesh_pts,mesh_lc)
        mygs.settings.free_boundary = False
        mygs.setup(order=2)
        mygs.p_scale=0.0
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
        results[2] = mygs.ffp_scale
    except BaseException as e:
        print(e)
    mp_q.put(results)


def validate_all(results, expected, rtol=1.E-3):
    if results is None:
        print("FAILED: error in run!")
        return False
    if len(results) != 3:
        print("FAILED: Number of results does not match expected")
        return False
    retval = True
    for i in range(3):
        if results[i] is None:
            print("FAILED: No result {0}!".format(i))
            retval = False
            continue
        if abs((results[i]-expected[i])/expected[i]) > rtol:
            print("FAILED: Result {0} outside tolerance!".format(i))
            print("  Expected = {0}".format(expected[i]))
            print("  Actual =   {0}".format(results[i]))
            retval = False
    return retval

@pytest.mark.coverage
def test_packages_serial():
    with open(os.path.join(test_dir,"oft_in.xml"), "w+") as fid:
        fid.write("""<oft>
  <thincurr>
    <eta>1.256637E-5</eta>
  </thincurr>
</oft>
""")
    expected = [4.95800457, 9.735667E-3, 49.131116]
    results = mp_run(run_all,(None,),test_dir)
    assert validate_all(results, expected)