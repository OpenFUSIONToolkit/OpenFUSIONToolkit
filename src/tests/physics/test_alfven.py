from __future__ import print_function
import os
import sys
import h5py
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
 use_petsc={10}
/

&mesh_options
 meshname='cube'
 cad_type={11}
 nlevels={1}
 nbase={0}
 grid_order=1
/

&t3d_options
 filename='cube.t3d'
 inpname='cube.inp'
 reflect='xyz'
 ref_per=T,T,T
/

&cube_options
 mesh_type=2
 ni=1,2,4
 rscale=1.,1.,2.
 shift=-0.5,-0.5,-1.
 ref_per=T,T,T
/

&test_alfven_options
 order={2}
 minlev={3}
 delta=1.e-3
 linear={5}
/

&xmhd_options
 xmhd_jcb=T
 xmhd_advec=F
 xmhd_hall=F
 xmhd_adv_temp=F
 xmhd_adv_den=F
 bbc='bc'
 vbc='all'
 dt={6}
 eta=1.e-3
 nu_par=1.e-3
 nsteps={7}
 rst_freq=125
 lin_tol={8}
 nl_tol=1.E-5
 nu_xmhd={4}
 nclean=5000
 rst_ind=0
 maxextrap=2
 ittarget=100
 xmhd_mfnk={9}
/

&xmhd_plot_options
 t0=.9e-5
 dt=1.e-5
 rst_start=0
 rst_end={7}
/
"""

# Common setup function and process handling
def alfven_setup(nbase,nlevels,order,minlev,nu_mhd='-1',linear=False,
                 mf=False,petsc='F',hex_mesh=False):
    mesh_type=1
    if hex_mesh:
        mesh_type=92
    dt='4.E-7'
    its='250'
    tol='1.E-12'
    mf_flag='F'
    if linear:
        lin_flag='T'
    else:
        lin_flag='F'
        if order==4:
            tol='1.E-13'
        if mf:
            mf_flag='T'
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, minlev, nu_mhd,
                                       lin_flag, dt, its, tol, mf_flag,
                                       petsc, mesh_type))
    return run_OFT("./test_alfven", nproc, 4000)

def validate_result(berr_exp,verr_exp,steps_exp=11,linear=False):
    retval = True
    with open('alfven.results','r') as fid:
        berr_test = float(fid.readline())
        if berr_test > 1.01*berr_exp:
            print("FAILED: Magnetic error too high!")
            print("  Expected = {0}".format(berr_exp))
            print("  Actual =   {0}".format(berr_test))
            retval = False
        verr_test = float(fid.readline())
        if verr_test > 1.01*verr_exp:
            print("FAILED: Velocity error too high!")
            print("  Expected = {0}".format(verr_exp))
            print("  Actual =   {0}".format(verr_test))
            retval = False
    #
    step_count=0
    B0_found = False
    with h5py.File("{0}.{1}.h5".format('oft_xdmf',str(1).zfill(4)),'r') as h5_file:
        if 'MUG' in h5_file:
            for _, mesh_obj in h5_file['MUG'].items():
                if mesh_obj['TYPE'][()] > 30:
                    for i in range(9999):
                        timestep = mesh_obj.get('{0:04d}'.format(i),None)
                        if timestep is None:
                            break
                        step_count += 1
                        if 'B0' in timestep:
                            B0_found = True
        else:
            print("FAILED: MUG plot group not found in output file")
            retval = False
    if step_count != steps_exp:
        print("FAILED: Incorrect number of time steps!")
        print("  Expected = {0}".format(steps_exp))
        print("  Actual =   {0}".format(step_count))
        retval = False
    if linear and (not B0_found):
        print("FAILED: Linearization plot fields missing for linear run!")
        retval = False
    return retval

#============================================================================
# Non-Linear test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p2(petsc_flag):
    berr_exp = 3.6337154826251369E-002
    verr_exp = 3.2114524936411738E-002
    assert alfven_setup(1,1,2,2,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p2_mpi(mf,petsc_flag):
    berr_exp = 3.6337154826251369E-002
    verr_exp = 3.2114524936411738E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p3(petsc_flag):
    berr_exp = 3.1344661129302444E-003
    verr_exp = 1.8333585622266604E-003
    assert alfven_setup(1,1,3,2,'0,2,2',petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p3_mpi(mf,petsc_flag):
    berr_exp = 3.1344661129302444E-003
    verr_exp = 1.8333585622266604E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p4(petsc_flag):
    berr_exp = 3.4446909270430422E-004
    verr_exp = 2.5485805438810947E-004
    assert alfven_setup(1,1,4,2,'0,10,4,2',petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p4_mpi(mf,petsc_flag):
    berr_exp = 3.4446909270430422E-004
    verr_exp = 2.5485805438810947E-004
    assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p2(petsc_flag):
    berr_exp = 3.6339799347819605E-002
    verr_exp = 3.2104584635227973E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p2_mpi(petsc_flag):
    berr_exp = 3.6339799347819605E-002
    verr_exp = 3.2104584635227973E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p3(petsc_flag):
    berr_exp = 3.1424512164583028E-003
    verr_exp = 1.8563760659227942E-003
    assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p3_mpi(petsc_flag):
    berr_exp = 3.1424512164583028E-003
    verr_exp = 1.8563760659227942E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p4(petsc_flag):
    berr_exp = 6.5622744829070675E-004
    verr_exp = 3.8358125487992524E-004
    assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p4_mpi(petsc_flag):
    berr_exp = 6.5622744829070675E-004
    verr_exp = 3.8358125487992524E-004
    assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Non-Linear test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p2(petsc_flag):
    berr_exp = 1.8217184894941060E-002
    verr_exp = 1.8207419422505190E-002
    assert alfven_setup(1,1,2,2,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p2_mpi(mf,petsc_flag):
    berr_exp = 1.8217184894941060E-002
    verr_exp = 1.8207419422505190E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p3(petsc_flag):
    berr_exp = 2.2785287899224150E-003
    verr_exp = 2.2636641474315156E-003
    #assert alfven_setup(1,1,3,2,'0,2,2',petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,3,3,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p3_mpi(mf,petsc_flag):
    berr_exp = 2.2785287899224150E-003
    verr_exp = 2.2636641474315156E-003
    #assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,3,4,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p4(petsc_flag):
    berr_exp = 1.8358574384029231E-004
    verr_exp = 1.8591003585664032E-004
    #assert alfven_setup(1,1,4,2,'0,10,4,2',petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,4,3,'0,0,0,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p4_mpi(mf,petsc_flag):
    berr_exp = 1.8358574384029231E-004
    verr_exp = 1.8591003585664032E-004
    #assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p2(petsc_flag):
    berr_exp = 1.8214844434637032E-002
    verr_exp = 1.8184066089285304E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p2_mpi(petsc_flag):
    berr_exp = 1.8214844434637032E-002
    verr_exp = 1.8184066089285304E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p3(petsc_flag):
    berr_exp = 2.2981509470689828E-003
    verr_exp = 2.2815413576994501E-003
    #assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,3,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p3_mpi(petsc_flag):
    berr_exp = 2.2981509470689828E-003
    verr_exp = 2.2815413576994501E-003
    #assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,3,4,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p4(petsc_flag):
    berr_exp = 1.5943398846056970E-003
    verr_exp = 1.3887929641779406E-003
    #assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,4,3,'0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p4_mpi(petsc_flag):
    berr_exp = 1.5943398846056970E-003
    verr_exp = 1.3887929641779406E-003
    #assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
