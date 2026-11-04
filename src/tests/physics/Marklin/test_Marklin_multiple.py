#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Marklin tests: multiple concurrent Marklin instances (cylinder and torus) in one OFT environment
from marklin_testing import mp_run, run_marklin, validate_eigs


# Test runner for multiple concurrent Marklin instances
def test_marklin_multiple():
    cyl_eig = 4.955207
    tor_eig = 6.317001
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,2,2,'torus_test.h5'))
    assert validate_eigs(results,tor_eig)
    results = mp_run(run_marklin,('torus_test.h5',1,2,2,'cyl_Marklin.h5'))
    assert validate_eigs(results,cyl_eig)
