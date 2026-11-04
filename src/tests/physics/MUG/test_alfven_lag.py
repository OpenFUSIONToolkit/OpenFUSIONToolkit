import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT
from mug_testing import validate_errors, validate_mug_output

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
 cad_type=92
 nlevels={1}
 nbase={0}
 grid_order=1
/

&cube_options
 mesh_type={11}
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
                 mf=False,petsc=False,hex_mesh=False):
    os.chdir(test_dir)
    mesh_type=1
    if hex_mesh:
        mesh_type=2
    dt='4.E-7'
    its='250'
    tol='1.E-12'
    mf_flag='F'
    petsc_flag=('T' if petsc else 'F')
    lin_flag=('T' if linear else 'F')
    if not linear:
        if order==4:
            tol='1.E-13'
        if mf:
            mf_flag='T'
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, minlev, nu_mhd,
                                       lin_flag, dt, its, tol, mf_flag,
                                       petsc_flag, mesh_type))
    return run_OFT("./test_alfven_lag", nproc, 2000)

def validate_result(berr_exp,verr_exp,steps_exp=11,linear=False):
    errors_ok = validate_errors('alfven.results', (('Magnetic', berr_exp), ('Velocity', verr_exp)), 1.01)
    output_ok = validate_mug_output(steps_exp, linear=linear)
    return errors_ok and output_ok

#============================================================================
# Non-Linear test runners for NP=2
def test_nl_r1_p2(petsc_flag=False):
    berr_exp = 4.2378931381112289E-002
    verr_exp = 4.1830283541690064E-002
    assert alfven_setup(1,1,2,2,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p2_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 4.2378931381112289E-002
    verr_exp = 4.1830283541690064E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
def test_nl_r1_p3(petsc_flag=False):
    berr_exp = 2.1613939097864013E-003
    verr_exp = 2.0399660410128605E-003
    assert alfven_setup(1,1,3,2,'0,2,2',petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p3_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 2.1613939097864013E-003
    verr_exp = 2.0399660410128605E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
def test_nl_r1_p4(petsc_flag=False):
    berr_exp = 2.6353342214151861E-004
    verr_exp = 2.4711224138100239E-004
    assert alfven_setup(1,1,4,2,'0,10,4,2',petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p4_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 2.6353342214151861E-004
    verr_exp = 2.4711224138100239E-004
    assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
def test_lin_r1_p2(petsc_flag=False):
    berr_exp = 4.2445515350299369E-002
    verr_exp = 4.1907987520465033E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(petsc_flag):
    berr_exp = 4.2445515350299369E-002
    verr_exp = 4.1907987520465033E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p3(petsc_flag=False):
    berr_exp = 2.1743513983662793E-003
    verr_exp = 2.0589634472042181E-003
    assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(petsc_flag):
    berr_exp = 2.1743513983662793E-003
    verr_exp = 2.0589634472042181E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p4(petsc_flag=False):
    berr_exp = 3.3133689642593500E-004
    verr_exp = 3.1807433060756504E-004
    assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(petsc_flag):
    berr_exp = 3.3133689642593500E-004
    verr_exp = 3.1807433060756504E-004
    assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Non-Linear test runners for NP=2
def test_hex_nl_r1_p2(petsc_flag=False):
    berr_exp = 1.8204622112726550E-002
    verr_exp = 1.8204690721409598E-002
    assert alfven_setup(1,1,2,2,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p2_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 1.8204622112726550E-002
    verr_exp = 1.8204690721409598E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
def test_hex_nl_r1_p3(petsc_flag=False):
    berr_exp = 2.2675346818877107E-003
    verr_exp = 2.2675252946567641E-003
    #assert alfven_setup(1,1,3,2,'0,2,2',petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,3,3,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p3_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 2.2675346818877107E-003
    verr_exp = 2.2675252946567641E-003
    #assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,3,4,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
def test_hex_nl_r1_p4(petsc_flag=False):
    berr_exp = 1.7790763796394384E-004
    verr_exp = 1.7781738147450092E-004
    #assert alfven_setup(1,1,4,2,'0,10,4,2',petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,4,3,'0,0,0,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p4_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 1.7790763796394384E-004
    verr_exp = 1.7781738147450092E-004
    #assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
def test_hex_lin_r1_p2(petsc_flag=False):
    berr_exp = 1.8185749172894013E-002
    verr_exp = 1.8185749172889051E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p2_mpi(petsc_flag):
    berr_exp = 1.8185749172894013E-002
    verr_exp = 1.8185749172889051E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p3(petsc_flag=False):
    berr_exp = 2.2830097215737651E-003
    verr_exp = 2.2830075204362107E-003
    #assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,3,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p3_mpi(petsc_flag):
    berr_exp = 2.2830097215737651E-003
    verr_exp = 2.2830075204362107E-003
    #assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,3,4,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p4(petsc_flag=False):
    berr_exp = 2.6757952245217109E-004
    verr_exp = 2.6757961586721258E-004
    #assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,4,3,'0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p4_mpi(petsc_flag):
    berr_exp = 2.6757952245217109E-004
    verr_exp = 2.6757961586721258E-004
    #assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
