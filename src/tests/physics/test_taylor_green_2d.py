from __future__ import print_function
import math
import os
import sys
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
 mesh_type=1
 ni={3},{3},0
 rscale=1.,1.,0.
 shift=0.,0.,0.
 ref_per=T,T,F
/

&xmhd_options
  order={2}
  n0=3.E24
  chi=1.E0
  eta=1.E0
  nu={6}
  gamma=1.67
  D_diff=1.E0
  m_i=205
  dt={4}
  nsteps={5}
  rst_freq=1000
  use_mfnk={7}
  timestep_cn=T
  ittarget=1000
  lin_tol=1.E-9
  nl_tol=1.E-7
  pm=F
  den_scale=1.E24
/
"""

def taylor_green_setup(nbase, nlevels, order, mf=False, ni=8, dt='2.5E-2', nsteps=20, nu='1.E-2'):
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    xml_file = 'oft_in_taylor_green.xml'
    if order == 2 and mf and nproc == 1:
        xml_file = 'oft_in_taylor_green_p2_mf.xml'
    mf_flag = ('T' if mf else 'F')
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, ni, dt, nsteps, nu, mf_flag))
    # Block-Jacobi makes the periodic Jacobian suitable for local LU solves.
    # Use two blocks only for serial P2 matrix-free to avoid a near-exact
    # inner GMRES solve.
    try:
        os.remove('taylor_green_2d.results')
    except FileNotFoundError:
        pass
    return run_OFT("./test_taylor_green_2d oft.in {0}".format(xml_file), nproc, 1000)

def validate_result(vxerr_exp, vzerr_exp, perr_exp):
    retval = True
    with open('taylor_green_2d.results', 'r') as fid:
        vxerr_test = float(fid.readline())
        if not math.isfinite(vxerr_test) or vxerr_test > 1.05*vxerr_exp:
            print("FAILED: v_x error too high!")
            print("  Expected = {0:.8E}".format(vxerr_exp))
            print("  Actual =   {0:.8E}".format(vxerr_test))
            retval = False
        vzerr_test = float(fid.readline())
        if not math.isfinite(vzerr_test) or vzerr_test > 1.05*vzerr_exp:
            print("FAILED: v_z error too high!")
            print("  Expected = {0:.8E}".format(vzerr_exp))
            print("  Actual =   {0:.8E}".format(vzerr_test))
            retval = False
        perr_test = float(fid.readline())
        if not math.isfinite(perr_test) or perr_test > 1.05*perr_exp:
            print("FAILED: Pressure error too high!")
            print("  Expected = {0:.8E}".format(perr_exp))
            print("  Actual =   {0:.8E}".format(perr_test))
            retval = False
    return retval

@pytest.mark.coverage
@pytest.mark.parametrize("error_index", (0, 1, 2), ids=("vx", "vz", "pressure"))
@pytest.mark.parametrize("nonfinite", (float('nan'), float('inf'), float('-inf')),
                         ids=("nan", "positive-inf", "negative-inf"))
def test_validate_result_rejects_nonfinite(tmp_path, monkeypatch, error_index, nonfinite):
    """@brief Reject nonfinite values in each reported physics error.
    @param tmp_path Isolated directory for the results file.
    @param monkeypatch Fixture that restores the working directory.
    @param error_index Position of the nonfinite error in the results file.
    @param nonfinite Nonfinite value to reject.
    """
    errors = [0.1, 0.2, 0.3]
    errors[error_index] = nonfinite
    (tmp_path / 'taylor_green_2d.results').write_text(
        ''.join('{0}\n'.format(error) for error in errors))
    monkeypatch.chdir(tmp_path)
    assert not validate_result(0.1, 0.2, 0.3)

@pytest.mark.coverage
def test_validate_result_accepts_finite(tmp_path, monkeypatch):
    """@brief Accept finite physics errors within the existing bounds.
    @param tmp_path Isolated directory for the results file.
    @param monkeypatch Fixture that restores the working directory.
    """
    (tmp_path / 'taylor_green_2d.results').write_text('0.1\n0.2\n0.3\n')
    monkeypatch.chdir(tmp_path)
    assert validate_result(0.1, 0.2, 0.3)

#============================================================================
# Taylor-Green vortex, incompressible, NP=2
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p2(mf):
    vxerr_exp = 2.9214986850512860E-2
    vzerr_exp = 2.9214986892598150E-2
    perr_exp = 8.2921136337263515E-2
    assert taylor_green_setup(1, 1, 2, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p2_mpi(mf):
    vxerr_exp = 2.9214986850512860E-2
    vzerr_exp = 2.9214986892598150E-2
    perr_exp = 8.2921136337263515E-2
    assert taylor_green_setup(1, 2, 2, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)

#============================================================================
# Taylor-Green vortex, incompressible, NP=3
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p3(mf):
    vxerr_exp = 1.4264674522612337E-3
    vzerr_exp = 1.4264676548348997E-3
    perr_exp = 5.8785206705912092E-3
    assert taylor_green_setup(1, 1, 3, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p3_mpi(mf):
    vxerr_exp = 1.4264674522612337E-3
    vzerr_exp = 1.4264676548348997E-3
    perr_exp = 5.8785206705912092E-3
    assert taylor_green_setup(1, 2, 3, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)

#============================================================================
# Taylor-Green vortex, incompressible, NP=4
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p4(mf):
    vxerr_exp = 5.6846822406755744E-5
    vzerr_exp = 5.6845321594545741E-5
    perr_exp = 2.4035705160280522E-4
    assert taylor_green_setup(1, 1, 4, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p4_mpi(mf):
    vxerr_exp = 5.6846822406755744E-5
    vzerr_exp = 5.6845321594545741E-5
    perr_exp = 2.4035705160280522E-4
    assert taylor_green_setup(1, 2, 4, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)
