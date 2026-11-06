#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: time-dependent setup for a mesh with coils but no conducting regions
import pytest
from tokamaker_testing import mp_run
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain


def run_td_noconductor_case(use_vcoil,mp_q):
    gs_mesh = gs_Domain(rextent=1.0,zextents=[0.0,1.0])
    gs_mesh.define_region('air',0.1,'boundary')
    gs_mesh.define_region('plasma',0.1,'plasma')
    gs_mesh.define_region('coil1',0.1,'coil')
    gs_mesh.define_region('coil2',0.1,'coil')
    gs_mesh.add_rectangle(0.4,0.4,0.2,0.2,'plasma')
    gs_mesh.add_rectangle(0.8,0.8,0.1,0.1,'coil1')
    gs_mesh.add_rectangle(0.8,0.4,0.1,0.1,'coil2')
    mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
    coil_dict = gs_mesh.get_coils()
    cond_dict = gs_mesh.get_conductors()
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    mygs.setup(order=2)
    mygs.set_coil_currents({'COIL1': 1.E3})
    mygs.vac_solve()
    if use_vcoil:
        mygs.set_vcoils({'COIL2': 1.E-3})
    try:
        mygs.setup_td(1.E-3,1.E-13,1.E-11)
        if use_vcoil:
            mygs.step_td(0.0,1.E-3)
    except Exception as e:
        mp_q.put(str(e))
        oftpy_dump_cov()
        return
    mp_q.put('')
    oftpy_dump_cov()


@pytest.mark.coverage
def test_td_no_conductors():
    '''Time-dependent setup on a mesh with coils but no conducting regions fails with a clear error.'''
    result = mp_run(run_td_noconductor_case,(False,))
    assert result is not None
    assert 'requires at least one conducting region' in result


@pytest.mark.coverage
def test_td_vcoil_only():
    '''A voltage-driven coil alone provides the conducting structure a time-dependent step needs.'''
    result = mp_run(run_td_noconductor_case,(True,))
    assert result == ''
