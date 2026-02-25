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
 ni=2,4,0
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
  chi=1.E0
  eta=1.E-8
  nu=1.E-8
  gamma=1.67
  D_diff=1.E0
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
def alfven_2d_setup(nbase,nlevels,order,linear=False, 
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
    return run_OFT("./test_alfven_2d", nproc, 4000)

def validate_result(verr_exp,berr_exp,steps_exp=11,linear=False):
    retval = True
    with open('alfven_2d.results','r') as fid:
        verr_test = float(fid.readline())
        if verr_test > 1.01*verr_exp:
            print("FAILED: Velocity error too high!")
            print("  Expected = {0}".format(verr_exp))
            print("  Actual =   {0}".format(verr_test))
            retval = False
        berr_test = float(fid.readline())
        if berr_test > 1.01*berr_exp:
            print("FAILED: Magnetic error too high!")
            print("  Expected = {0}".format(berr_exp))
            print("  Actual =   {0}".format(berr_test))
            retval = False
    return retval

# #============================================================================
# # Non-Linear test runners for NP=2
# def test_nl_p2(petsc_flag=False):
#     berr_exp = 5.3027098135759567E-003
#     verr_exp = 1.1873721610050637E-002
#     assert alfven_2d_setup(1,1,2)
#     assert validate_result(verr_exp, berr_exp)
#============================================================================
# Non-Linear test runners for NP=3
def test_nl_p3(petsc_flag=False):
    berr_exp = 1.6560625574643158E-003
    verr_exp = 2.2167700412833043E-003
    assert alfven_2d_setup(1,1,3)
    assert validate_result(verr_exp, berr_exp)
