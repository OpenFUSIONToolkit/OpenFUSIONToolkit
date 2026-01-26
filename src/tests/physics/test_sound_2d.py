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
 mesh_type={6}
 ni=16,16,0
 rscale=2.,2.,0.
 shift=-1.0,-1.0,0.
 ref_per=T,T,F
/

&xmhd_options
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
def sound_2d_setup(nbase,nlevels,order,
                 mf=False,hex_mesh=False):
    mesh_type=1
    if hex_mesh:
        mesh_type=2
    dt='2.0E-7'
    its='500'
    tol='1.E-8'
    mf_flag='F'
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order,
                                       dt, its, mf_flag,
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

# Non-Linear test runners for NP=2
def test_nl(petsc_flag=False):
    verr_exp = 3.9219635487227018E-002
    nerr_exp = 3.9144202292077303E-002
    terr_exp = 3.9339922643013781E-002
    assert sound_2d_setup(1,1,2)
    assert validate_result(verr_exp, nerr_exp, terr_exp)
