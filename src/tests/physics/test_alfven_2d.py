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
 use_petsc={9}
/

&mesh_options
 meshname='test'
 cad_type=92
 nlevels={1}
 nbase={0}
 grid_order=2
 fix_boundary=T
/

&cube_options
 mesh_type={8}
 ni=4,1,0
 rscale=2.0,0.5,0.
 shift=-1.0,-0.25,0.
 ref_per=T,T,F
/

&xmhd_options
  linear={3}
  nl_tol=1.E-8
  cyl=F
  order={2}
  n0=1E19,
  psi0=0.d0,
  velx0=0,
  vely0=0,
  velz0=0,
  t0=1E3,
  by0=0.d0,
  bx0=0.d0,
  bz0=0.d0,
  chi=1.E0
  eta=0.00795774715
  nu=0.59880239
  gamma=1.67
  D_diff=1.E0
  dt={4}
  nsteps={5}
  rst_freq=5
  use_mfnk={6}
  lin_tol={7}
  pm=F
  den_scale=1.E19
/

&xmhd_plot_options
 rst_start=0
 rst_end={5}
/
"""

# Common setup function and process handling
def alfven_2d_setup(nbase,nlevels,order,
                    linear=False,mf=False,petsc=False,hex_mesh=False):
    mesh_type=1
    if hex_mesh:
        mesh_type=2
    dt='4.0E-7'
    its='250'
    tol='1.E-10'
    mf_flag='F'
    petsc_flag=('T' if petsc else 'F')
    lin_flag=('T' if linear else 'F')
    if (not linear) and mf:
        tol='1.E-9'
        mf_flag='T'
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, lin_flag,
                                         dt, its, mf_flag, tol, mesh_type, petsc_flag))
    return run_OFT("./test_alfven_2d", nproc, 4000)

def validate_result(verr_exp,berr_exp,steps_exp=11,linear=False):
    retval = True
    with open('alfven_2d.results','r') as fid:
        verr_test = float(fid.readline())
        if verr_test > 1.05*verr_exp:
            print("FAILED: Velocity error too high!")
            print("  Expected = {0:.8E}".format(verr_exp))
            print("  Actual =   {0:.8E}".format(verr_test))
            retval = False
        berr_test = float(fid.readline())
        if berr_test > 1.05*berr_exp:
            print("FAILED: Magnetic error too high!")
            print("  Expected = {0:.8E}".format(berr_exp))
            print("  Actual =   {0:.8E}".format(berr_test))
            retval = False
    return retval

#============================================================================
# Non-Linear test runners for NP=2
@pytest.mark.parametrize("mf", (False, True))
def test_nl_r1_p2(mf,petsc_flag=False):
    verr_exp = 3.627799135938326E-2
    berr_exp = 2.7401604608693758E-2
    assert alfven_2d_setup(1,1,2,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p2_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    verr_exp = 3.627799135938326E-2
    berr_exp = 2.7401604608693758E-2
    assert alfven_2d_setup(1,2,2,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.parametrize("mf", (False, True))
def test_nl_r1_p3(mf,petsc_flag=False):
    verr_exp = 3.1531577649676595E-2
    berr_exp = 2.7201129077752486E-2
    assert alfven_2d_setup(1,1,3,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p3_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    verr_exp = 3.1531577649676595E-2
    berr_exp = 2.7201129077752486E-2
    assert alfven_2d_setup(1,2,3,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.parametrize("mf", (False, True))
def test_nl_r1_p4(mf,petsc_flag=False):
    verr_exp = 3.13872526257328E-2
    berr_exp = 2.719581233696034E-2
    assert alfven_2d_setup(1,1,4,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p4_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    verr_exp = 3.13872526257328E-2
    berr_exp = 2.719581233696034E-2
    assert alfven_2d_setup(1,2,4,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

#============================================================================
#Linear test runners for NP=2
def test_lin_r1_p2(petsc_flag=False):
    verr_exp = 1.811020939389565E-2
    berr_exp = 4.923999966698998E-3
    assert alfven_2d_setup(1,1,2, linear = True, petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(petsc_flag):
    verr_exp = 1.811020939389565E-2
    berr_exp = 4.923999966698998E-3
    assert alfven_2d_setup(1,2,2, petsc=petsc_flag, linear = True)
    assert validate_result(verr_exp, berr_exp)

#============================================================================
#Linear test runners for NP=3
def test_lin_r1_p3(petsc_flag=False):
    verr_exp = 2.5746106154166067E-3
    berr_exp = 5.148998451109316E-4
    assert alfven_2d_setup(1,1,3, linear = True, petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(petsc_flag):
    verr_exp = 2.5746106154166067E-3
    berr_exp = 5.148998451109316E-4
    assert alfven_2d_setup(1,2,3, petsc=petsc_flag, linear = True)
    assert validate_result(verr_exp, berr_exp) 

#============================================================================
#Linear test runners for NP=4
def test_lin_r1_p4(petsc_flag=False):
    verr_exp = 4.997848969169074E-4
    berr_exp = 4.504851399342237E-4
    assert alfven_2d_setup(1,1,4, linear = True, petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(petsc_flag):
    verr_exp = 4.997848969169074E-4
    berr_exp = 4.504851399342237E-4
    assert alfven_2d_setup(1,2,4, petsc=petsc_flag, linear = True)
    assert validate_result(verr_exp, berr_exp) 