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
 use_petsc={8}
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
 mesh_type={9}
 ni=1,4,0
 rscale=2.,2.,0.
 shift=1.0,-1.0,0.
 ref_per=T,T,F
/

&xmhd_options
  linear={3}
  lin_tol ={6}
  nl_tol = 1.E-8
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
  dt={4}
  nsteps={5}
  rst_freq=5
  use_mfnk={7}
  pm=F
  den_scale=1.d19
  
/

&xmhd_plot_options
 t0=.9e-5
 dt=1.e-5
 rst_start=0
 rst_end={4}
/
"""

# Common setup function and process handling
def sound_2d_setup(nbase,nlevels,order,linear=False, 
                 mf=False,petsc=False, hex_mesh=False):
    mesh_type=1
    if hex_mesh:
        mesh_type=2
    dt='2.500E-8'
    its='4000'
    tol='1.E-12'
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
        fid.write(oft_in_template.format(nbase, nlevels, order,lin_flag,
                                       dt, its, tol, mf_flag,petsc_flag,
                                    mesh_type))
    return run_OFT("./test_sound_2d", nproc, 4000)

def validate_result(verr_exp,nerr_exp,terr_exp, steps_exp=11,linear=False):
    retval = True
    with open('sound_2d.results','r') as fid:
        verr_test = float(fid.readline())
        if verr_test > 1.01*verr_exp:
            print("FAILED: Velocity error too high!")
            print("  Expected = {0}".format(verr_exp))
            print("  Actual =   {0}".format(verr_test))
            retval = False
        nerr_test = float(fid.readline())
        if nerr_test > 1.01*nerr_exp:
            print("FAILED: Density error too high!")
            print("  Expected = {0}".format(nerr_exp))
            print("  Actual =   {0}".format(nerr_test))
            retval = False
        terr_test = float(fid.readline())
        if terr_test > 1.01*nerr_exp:
            print("FAILED: Density error too high!")
            print("  Expected = {0}".format(terr_exp))
            print("  Actual =   {0}".format(terr_test))
            retval = False
    return retval
#============================================================================
#Non-Linear test runners for NP=2

def test_nl_p2():
    verr_exp = 1.4477427222428568E-002
    nerr_exp = 1.8098066719625182E-002
    terr_exp = 1.8188555555408339E-002
    assert sound_2d_setup(1,1,2)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_p2_mpi(mf, petsc_flag):
    verr_exp = 1.4477427222428568E-002
    nerr_exp = 1.8098066719625182E-002
    terr_exp = 1.8188555555408339E-002
    assert sound_2d_setup(1,1,2, mf=mf, petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

#============================================================================
# Non-Linear test runners for NP=3

def test_nl_p3():
    verr_exp = 8.1141504705730411E-003
    nerr_exp = 8.1134814444124882E-003
    terr_exp = 8.1545550097106480E-003
    assert sound_2d_setup(1,1,3)
    assert validate_result(verr_exp, nerr_exp, terr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_p3_mpi(mf, petsc_flag):
    verr_exp = 8.1141504705730411E-003
    nerr_exp = 8.1134814444124882E-003
    terr_exp = 8.1545550097106480E-003
    assert sound_2d_setup(1,1,3, mf=mf, petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

#============================================================================
# Non-Linear test runners for NP=4

def test_nl_p4():
    verr_exp = 7.9690160287047319E-003
    nerr_exp = 7.9629023996221435E-003
    terr_exp = 8.0027129244056965E-003
    assert sound_2d_setup(1,1,4)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_p4_mpi(mf, petsc_flag):
    verr_exp = 7.9690160287047319E-003
    nerr_exp = 7.9629023996221435E-003
    terr_exp = 8.0027129244056965E-003
    assert sound_2d_setup(1,1,4, mf=mf, petsc=petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)
#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
def test_lin_p2():
    verr_exp = 1.4475321807125909E-002
    nerr_exp = 1.8096366550025649E-002
    terr_exp =1.8186848072669490E-002
    assert sound_2d_setup(1,1,2, linear = True)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_p2_mpi(mf, petsc_flag):
    verr_exp = 1.4475321807125909E-002
    nerr_exp = 1.8096366550025649E-002
    terr_exp =1.8186848072669490E-002
    assert sound_2d_setup(1,1,2, mf=mf, petsc=petsc_flag, linear = True)
    assert validate_result(verr_exp, nerr_exp, terr_exp)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
def test_lin_p3():
    verr_exp = 8.1092620709466330E-003
    nerr_exp = 8.1067291952836486E-003
    terr_exp = 8.1472470392866935E-003
    assert sound_2d_setup(1,1,3, linear = True)

@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_p3_mpi(mf, petsc_flag):
    verr_exp = 8.1092620709466330E-003
    nerr_exp = 8.1067291952836486E-003
    terr_exp = 8.1472470392866935E-003
    assert sound_2d_setup(1,1,3, linear = True, petsc= petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)
#============================================================================
#Linear test runners for NP=4
@pytest.mark.linear
def test_lin_p4():
    verr_exp = 7.9648052134242055E-003
    nerr_exp = 7.9586787434810218E-003
    terr_exp = 7.9984691203724943E-003
    assert sound_2d_setup(1,1,4, linear = True)
    assert validate_result(verr_exp, nerr_exp, terr_exp)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_p4_mpi(mf, petsc_flag):
    verr_exp = 7.9648052134242055E-003
    nerr_exp = 7.9586787434810218E-003
    terr_exp = 7.9984691203724943E-003
    assert sound_2d_setup(1,1,4, linear = True,petsc = petsc_flag)
    assert validate_result(verr_exp, nerr_exp, terr_exp)
