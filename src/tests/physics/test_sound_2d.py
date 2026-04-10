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
 rscale=2.,0.5,0.
 shift=-1.0,-0.25,0.
 ref_per=T,T,F
/

&xmhd_options
  linear={7}
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
  chi=1.E-16
  eta=1.E-16
  nu=1.E-16
  gamma=1.67
  D_diff=1.E-16
  dt={3}
  nsteps={4}
  rst_freq=5
  use_mfnk={5}
  lin_tol={6}
  nl_tol = 1.E-8
  pm=F
  den_scale=1.E19
/

&xmhd_plot_options
 rst_start=0
 rst_end={4}
/
"""

# Common setup function and process handling
def sound_2d_setup(nbase,nlevels,order,
                   linear=False,mf=False,petsc=False,hex_mesh=False):
    mesh_type=1
    if hex_mesh:
        mesh_type=2
    dt='1.E-6'
    its='100'
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
        fid.write(oft_in_template.format(nbase, nlevels, order,
                                         dt, its, mf_flag, tol, lin_flag, mesh_type, petsc_flag))
    return run_OFT("./test_sound_2d", nproc, 4000)

def validate_result(verr_exp,nerr_exp,terr_exp, steps_exp=11,linear=False):
    retval = True
    with open('sound_2d.results','r') as fid:
        verr_test = float(fid.readline())
        if verr_test > 1.05*verr_exp:
            print("FAILED: Velocity error too high!")
            print("  Expected = {0:.8E}".format(verr_exp))
            print("  Actual =   {0:.8E}".format(verr_test))
            retval = False
        nerr_test = float(fid.readline())
        if nerr_test > 1.05*nerr_exp:
            print("FAILED: Density error too high!")
            print("  Expected = {0:.8E}".format(nerr_exp))
            print("  Actual =   {0:.8E}".format(nerr_test))
            retval = False
        terr_test = float(fid.readline())
        if terr_test > 1.05*terr_exp:
            print("FAILED: Temperature error too high!")
            print("  Expected = {0:.8E}".format(terr_exp))
            print("  Actual =   {0:.8E}".format(terr_test))
            retval = False
    return retval
#============================================================================
# Non-Linear test runners for NP=2
@pytest.mark.parametrize("mf", (False, True))
def test_nl_r1_p2(mf,petsc_flag=False):
    verr_exp = 3.4757492919345006E-3
    nerr_exp = 1.0042257109873824E-2
    terr_exp = 1.0042257501380656E-2
    assert sound_2d_setup(1,1,2,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p2_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    verr_exp = 3.4757492919345006E-3
    nerr_exp = 1.0042257109873824E-2
    terr_exp = 1.0042257501380656E-2
    assert sound_2d_setup(1,2,2,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.parametrize("mf", (False, True))
def test_nl_r1_p3(mf,petsc_flag=False):
    verr_exp = 2.1428547309174486E-3
    nerr_exp = 2.1814302475110835E-3
    terr_exp = 2.1800283394622167E-3
    assert sound_2d_setup(1,1,3,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p3_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    verr_exp = 2.1428547309174486E-3
    nerr_exp = 2.1814302475110835E-3
    terr_exp = 2.1800283394622167E-3
    assert sound_2d_setup(1,2,3,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.parametrize("mf", (False, True))
def test_nl_r1_p4(mf,petsc_flag=False):
    verr_exp = 2.081834197322753E-3
    nerr_exp = 2.083623647676870E-3
    terr_exp = 2.083635171836626E-3
    assert sound_2d_setup(1,1,4,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p4_mpi(mf,petsc_flag):
    if mf and petsc_flag:
        pytest.skip()
    verr_exp = 2.081834197322753E-3
    nerr_exp = 2.083623647676870E-3
    terr_exp = 2.083635171836626E-3
    assert sound_2d_setup(1,2,4,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

#============================================================================
# Linear test runners for NP=2
def test_lin_r1_p2(petsc_flag=False):
    verr_exp = 8.54387175E-03
    nerr_exp = 1.18926533E-02
    terr_exp = 1.18926532E-02
    assert sound_2d_setup(1,1,2,linear=True,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(petsc_flag):
    verr_exp = 8.54387175E-03
    nerr_exp = 1.18926533E-02
    terr_exp = 1.18926532E-02
    assert sound_2d_setup(1,2,2,linear=True,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

#============================================================================
# Linear test runners for NP=3
def test_lin_r1_p3(petsc_flag=False):
    verr_exp = 8.74577449E-03
    nerr_exp = 8.75411608E-03
    terr_exp = 8.75411607E-03
    assert sound_2d_setup(1,1,3,linear=True,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(petsc_flag):
    verr_exp = 8.74577449E-03
    nerr_exp = 8.75411608E-03
    terr_exp = 8.75411607E-03
    assert sound_2d_setup(1,2,3,linear=True,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

#============================================================================
# Linear test runners for NP=4
def test_lin_r1_p4(petsc_flag=False):
    verr_exp = 8.74382139E-03
    nerr_exp = 8.74369833E-03
    terr_exp = 8.74369838E-03
    assert sound_2d_setup(1,1,4,linear=True,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(petsc_flag):
    verr_exp = 8.74382139E-03
    nerr_exp = 8.74369833E-03
    terr_exp = 8.74369838E-03
    assert sound_2d_setup(1,2,4,linear=True,petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)
