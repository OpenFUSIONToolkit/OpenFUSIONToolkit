#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Marklin tests: vacuum field in a torus with a toroidal flux jump plane (curved mesh)
import pytest
from marklin_testing import mp_run, clear_cache, run_marklin_vac


# Test runners for torus vacuum field (curved mesh)
def marklin_torus_vac(order,mag_val):
    if order == 1:
        clear_cache()
    results = mp_run(run_marklin_vac,('torus_test.h5',order,2))
    print(results)
    assert abs((results-mag_val)/mag_val) < 1.E-4

@pytest.mark.coverage
def test_marklin_vac_g2_p1():
    field_mag = 1.145501
    marklin_torus_vac(1,field_mag)

@pytest.mark.coverage
def test_marklin_vac_g2_p2():
    field_mag = 1.188539
    marklin_torus_vac(2,field_mag)

def test_marklin_vac_g2_p3():
    field_mag = 1.188352
    marklin_torus_vac(3,field_mag)

def test_marklin_vac_g2_p4():
    field_mag = 1.188126
    marklin_torus_vac(4,field_mag)
