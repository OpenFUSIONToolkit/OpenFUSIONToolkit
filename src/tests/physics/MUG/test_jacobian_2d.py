#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: consistency of the 2D xMHD approximate Jacobian with the nonlinear residual (planar and
# cylindrical geometry, compressible and incompressible models), and input checks of the 2D driver
import os
import sys
import re
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT

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
 nlevels=1
 nbase=1
 grid_order=2
 fix_boundary=T
/

&cube_options
 mesh_type=1
 ni=4,4,0
 rscale=1.,1.,0.
 shift=1.,-0.5,0.
 ref_per=F,F,F
/

&xmhd_jac_options
 order={0}
 cyl={1}
 incomp={2}
 b0={3}
 check_abort={4}
/
"""

# Largest allowed relative difference between Jacobian and finite-difference products
JAC_TOL = 1.E-6


def jacobian_2d_setup(order, cyl, incomp, b0, check_abort=False, expect_abort=None):
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(order, 'T' if cyl else 'F', 'T' if incomp else 'F',
                                         ','.join('{0:.4E}'.format(val) for val in b0),
                                         'T' if check_abort else 'F'))
    return run_OFT("./test_jacobian_2d", 1, 600, return_stdout=True, expect_abort=expect_abort)


def validate_result(result):
    status, std_out = result
    if not status:
        return False
    match = re.search(r'MAX_REL_DIFF =\s*(\S+)', std_out)
    if match is None:
        print('FAILED: maximum relative difference not reported')
        return False
    max_diff = float(match.group(1))
    if max_diff > JAC_TOL:
        print('FAILED: Jacobian differs from finite-difference derivative ({0:.4E} > {1:.4E})'.format(max_diff, JAC_TOL))
        return False
    return True


#============================================================================
@pytest.mark.parametrize("order", (2, 3))
@pytest.mark.parametrize("incomp", (False, True))
def test_jacobian_planar(order, incomp):
    '''Jacobian consistency in planar geometry, with a background field in all three directions.'''
    assert validate_result(jacobian_2d_setup(order, False, incomp, (1.E-2, 2.E-2, 1.5E-2)))


@pytest.mark.parametrize("order", (2, 3))
@pytest.mark.parametrize("incomp", (False, True))
def test_jacobian_cyl(order, incomp):
    '''Jacobian consistency in cylindrical geometry, with a Z-directed background field.'''
    assert validate_result(jacobian_2d_setup(order, True, incomp, (0.0, 0.0, 1.E-2)))


def test_cyl_b0_abort():
    '''A background field with R or phi components is rejected in cylindrical geometry.'''
    status, _ = jacobian_2d_setup(2, True, False, (1.E-2, 0.0, 0.0), check_abort=True,
                                  expect_abort='Only vertical B_0')
    assert status
