#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: cut-cell quadrature of the plasma integrals at the LCFS (setting `lcfs_cut_cells`)
import pytest
import numpy as np
from tokamaker_testing import (mp_run, validate_dict, ITER_eq_dict, ITER_isoflux_pts, ITER_x_point,
                               create_ITER_mesh, setup_ITER_gs, set_ITER_constraints, solve_ITER_eq)
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun


# Reference values for the edge case of `run_ITER_edge_case`, with cut cells, on the ITER mesh refined by two
# (mesh_resolution=2.0, order 2); checked against a live solve by `test_ITER_edge_cut_cells_reference`
ITER_edge_ref_dict = {
    'Ip': 15599994.426932985,
    'W_MHD': 261817359.70739892,
    'vol': 817.1927484319564,
    'q_0': 0.8527581619021124,
    'q_95': 2.7848873842917548,
    'l_i': 0.8554081360691612,
    'beta_pol': 45.80556660545955,
}


def run_ITER_eq_cut_cells(fe_order,mp_q):
    try:
        create_ITER_mesh()
    except Exception as e:
        print(e)
        mp_q.put(None)
        return
    myOFT = OFT_env(nthreads=-1)
    mygs = setup_ITER_gs(myOFT,fe_order)
    set_ITER_constraints(mygs)
    mygs.settings.lcfs_cut_cells = True
    mygs.update_settings()
    try:
        EQ_obj, nl_its = solve_ITER_eq(mygs)
    except ValueError:
        mp_q.put(None)
        return
    eq_info = EQ_obj.get_stats(li_normalization='ITER')
    Lmat = mygs.get_coil_Lmat()
    eq_info['LCS1'] = Lmat[mygs.coil_sets['CS1U']['id'],mygs.coil_sets['CS1U']['id']]
    eq_info['MCS1_plasma'] = Lmat[mygs.coil_sets['CS1U']['id'],-1]
    eq_info['Lplasma'] = Lmat[-1,-1]
    eq_info['nl_its'] = nl_its
    mp_q.put([eq_info])
    oftpy_dump_cov()


def run_ITER_edge_case(mesh_resolution,mp_q):
    # Mesh refined by `mesh_resolution` in its own file, so the shared ITER_mesh.h5 keeps the default resolution
    mesh_file = 'ITER_mesh.h5' if mesh_resolution == 1.0 else 'ITER_mesh_r{0:g}.h5'.format(mesh_resolution)
    try:
        create_ITER_mesh(mesh_resolution,mesh_file)
    except Exception as e:
        print(e)
        mp_q.put(None)
        return
    myOFT = OFT_env(nthreads=-1)
    mygs = setup_ITER_gs(myOFT,2,mesh_file=mesh_file)
    set_ITER_constraints(mygs)
    # Shape and X-point targets moved up by 6 cm
    shift = np.array([0.0, 0.06])
    mygs.set_isoflux(np.vstack((ITER_isoflux_pts,ITER_x_point))+shift)
    mygs.set_saddles(ITER_x_point+shift)
    ffp_prof = create_power_flux_fun(40,1.5,2.0)
    pp_prof = create_power_flux_fun(40,4.0,1.0)
    # P' = 1 - psi_N^4 + 0.3: a finite pressure gradient at the LCFS, as in pedestal profiles
    pp_prof['y'] += 0.3
    mygs.set_profiles(ffp_prof=ffp_prof,pp_prof=pp_prof)
    mygs.settings.lcfs_cut_cells = True
    mygs.update_settings()
    mygs.init_psi(6.3, 0.56, 2.0, 1.4, 0.0)
    try:
        _, nl_its = mygs.solve(return_its=True)
    except ValueError:
        mp_q.put(None)
        return
    eq_info = mygs.get_stats(li_normalization='ITER')
    eq_info['nl_its'] = nl_its
    # Tilt psi vertically in 400 small steps (the LCFS moves by about 1 mm in total, bounds recomputed at each
    # step) and record the plasma volume, an integral over the plasma part of the cells. The volume is used
    # because `set_psi` keeps the flux functions' normalization from the solve, so current and stored energy
    # integrands would themselves change where psi falls below the old LCFS value.
    psi0 = mygs.get_psi(False)
    tilt = (mygs.r[:,1]-mygs.o_point[1])*abs(mygs.psi_bounds[1]-mygs.psi_bounds[0])
    vol_vals = []
    for eps in np.linspace(-2.E-4,2.E-4,401):
        mygs.set_psi(psi0+eps*tilt,update_bounds=True)
        vol_vals.append(mygs.get_globals()[2])
    steps = np.abs(np.diff(vol_vals))
    eq_info['vol_step_ratio'] = np.max(steps)/np.median(steps)
    mp_q.put([eq_info])
    oftpy_dump_cov()


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))
def test_ITER_eq_cut_cells(order):
    '''The ITER equilibrium of test_ITER_eq solved with `lcfs_cut_cells`: the cut-cell quadrature changes the
    discrete answer only by the boundary discretization error, so global quantities and the plasma-coil
    inductances (which use the same quadrature) match the reference values of the default quadrature'''
    eq_dict = {key: val for key, val in ITER_eq_dict.items() if not key.startswith('fsa_')}
    results = mp_run(run_ITER_eq_cut_cells,(order,))
    assert validate_dict(results,eq_dict)


@pytest.mark.coverage
def test_ITER_edge_cut_cells():
    '''Free-boundary ITER case with a finite pressure gradient at the LCFS (P' = 1 - psi_N^4 + 0.3, shape
    targets moved up by 6 cm). With the default point test at the LCFS the source changes by a finite amount
    whenever a quadrature point crosses the boundary, and this case can end in a period-2 cycle of the Picard
    iteration instead of converging (which inputs do so depends on the platform). With `lcfs_cut_cells` it must
    converge, agree with the solution on a mesh refined by two, and give plasma integrals that change smoothly
    with psi: in a small vertical tilt of psi no step of the plasma volume stands out from the trend'''
    results = mp_run(run_ITER_edge_case,(1.0,))
    assert results is not None
    assert results[0]['nl_its'] > 0
    assert validate_dict(results,ITER_edge_ref_dict)
    # A jump from quadrature points crossing the LCFS would make one step many times the median
    assert results[0]['vol_step_ratio'] < 3.0


@pytest.mark.slow
def test_ITER_edge_cut_cells_reference():
    '''The refined-mesh reference values of test_ITER_edge_cut_cells, from a live solve'''
    results = mp_run(run_ITER_edge_case,(2.0,),timeout=180)
    assert results is not None
    assert results[0]['nl_its'] > 0
    assert validate_dict(results,ITER_edge_ref_dict)
    assert results[0]['vol_step_ratio'] < 3.0
