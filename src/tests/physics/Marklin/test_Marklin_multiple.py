#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Marklin tests: multiple concurrent Marklin instances (cylinder and torus) in one OFT environment
from marklin_testing import mp_run, run_marklin, validate_eigs


def test_marklin_multiple():
    '''Two Marklin objects on different meshes (cylinder and torus) in one OFT environment, created in both
    orders: the eigenvalue of the second instance must match its single-instance value.'''
    cyl_eig = 4.955207
    tor_eig = 6.317001
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,2,2,'torus_test.h5'))
    assert validate_eigs(results,tor_eig)
    results = mp_run(run_marklin,('torus_test.h5',1,2,2,'cyl_Marklin.h5'))
    assert validate_eigs(results,cyl_eig)
