#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: decaying Taylor-Green vortex with the incompressible 2D xMHD solver, compared with the exact
# viscous decay of velocity and pressure (G. I. Taylor and A. E. Green, Proc. R. Soc. Lond. A 158, 499
# (1937))
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT
from mug_testing import validate_errors

# Input file template, filled in by the setup function
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
    # The level above nbase is the base mesh distributed over the processes (a transfer level, not a refinement),
    # so MPI runs reproduce the serial errors
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    # Block-Jacobi makes the periodic Jacobian suitable for local LU solves; serial quadratic matrix-free runs use
    # two blocks to avoid a near-exact inner GMRES solve
    xml_file = 'oft_in_taylor_green.xml'
    if order == 2 and mf and nproc == 1:
        xml_file = 'oft_in_taylor_green_p2_mf.xml'
    mf_flag = ('T' if mf else 'F')
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, ni, dt, nsteps, nu, mf_flag))
    try:
        os.remove('taylor_green_2d.results')
    except FileNotFoundError:
        pass
    return run_OFT("./test_taylor_green_2d oft.in {0}".format(xml_file), nproc, 1000)

# Errors are regression values from a reference run; each may exceed its stored value by up to 5%
def validate_result(vxerr_exp, vzerr_exp, perr_exp):
    return validate_errors('taylor_green_2d.results', (('v_x', vxerr_exp), ('v_z', vzerr_exp), ('Pressure', perr_exp)), 1.05)

@pytest.mark.coverage
@pytest.mark.parametrize("error_index", (0, 1, 2), ids=("vx", "vz", "pressure"))
@pytest.mark.parametrize("nonfinite", (float('nan'), float('inf'), float('-inf')),
                         ids=("nan", "positive-inf", "negative-inf"))
def test_validate_result_rejects_nonfinite(tmp_path, monkeypatch, error_index, nonfinite):
    '''A NaN or infinite error in any field fails validation, so a diverged run cannot pass.'''
    errors = [0.1, 0.2, 0.3]
    errors[error_index] = nonfinite
    (tmp_path / 'taylor_green_2d.results').write_text(
        ''.join('{0}\n'.format(error) for error in errors))
    monkeypatch.chdir(tmp_path)
    assert not validate_result(0.1, 0.2, 0.3)

@pytest.mark.coverage
def test_validate_result_accepts_finite(tmp_path, monkeypatch):
    '''Finite errors at their expected values pass validation.'''
    (tmp_path / 'taylor_green_2d.results').write_text('0.1\n0.2\n0.3\n')
    monkeypatch.chdir(tmp_path)
    assert validate_result(0.1, 0.2, 0.3)

#============================================================================
# Quadratic elements
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p2(mf):
    '''Taylor-Green vortex with quadratic elements, serial.
    Also covers matrix-free Newton-Krylov.'''
    vxerr_exp = 2.9214986850512860E-2
    vzerr_exp = 2.9214986892598150E-2
    perr_exp = 8.2921136337263515E-2
    assert taylor_green_setup(1, 1, 2, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p2_mpi(mf):
    '''Taylor-Green vortex with quadratic elements, on two MPI processes.
    Also covers matrix-free Newton-Krylov.'''
    vxerr_exp = 2.9214986850512860E-2
    vzerr_exp = 2.9214986892598150E-2
    perr_exp = 8.2921136337263515E-2
    assert taylor_green_setup(1, 2, 2, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)

#============================================================================
# Cubic elements
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p3(mf):
    '''Taylor-Green vortex with cubic elements, serial.
    Also covers matrix-free Newton-Krylov.'''
    vxerr_exp = 1.4264674522612337E-3
    vzerr_exp = 1.4264676548348997E-3
    perr_exp = 5.8785206705912092E-3
    assert taylor_green_setup(1, 1, 3, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p3_mpi(mf):
    '''Taylor-Green vortex with cubic elements, on two MPI processes.
    Also covers matrix-free Newton-Krylov.'''
    vxerr_exp = 1.4264674522612337E-3
    vzerr_exp = 1.4264676548348997E-3
    perr_exp = 5.8785206705912092E-3
    assert taylor_green_setup(1, 2, 3, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)

#============================================================================
# Quartic elements
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p4(mf):
    '''Taylor-Green vortex with quartic elements, serial.
    Also covers matrix-free Newton-Krylov.'''
    vxerr_exp = 5.6846822406755744E-5
    vzerr_exp = 5.6845321594545741E-5
    perr_exp = 2.4035705160280522E-4
    assert taylor_green_setup(1, 1, 4, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
def test_taylor_green_p4_mpi(mf):
    '''Taylor-Green vortex with quartic elements, on two MPI processes.
    Also covers matrix-free Newton-Krylov.'''
    vxerr_exp = 5.6846822406755744E-5
    vzerr_exp = 5.6845321594545741E-5
    perr_exp = 2.4035705160280522E-4
    assert taylor_green_setup(1, 2, 4, mf=mf)
    assert validate_result(vxerr_exp, vzerr_exp, perr_exp)
