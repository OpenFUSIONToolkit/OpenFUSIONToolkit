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
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
from OpenFUSIONToolkit.TokaMaker.reconstruction import reconstruction


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


def lcfs_area_volume(mygs):
    '''Area and volume inside the traced LCFS contour (independent of the plasma-region quadrature)'''
    lcfs = mygs.trace_surf(1.0)
    R, Z = lcfs[:,0], lcfs[:,1]
    R2, Z2 = np.roll(R,-1), np.roll(Z,-1)
    cross = R*Z2-R2*Z
    return abs(np.sum(cross))/2.0, abs(np.pi*np.sum((R+R2)*cross)/3.0)


def plasma_integral_checks(mygs,Ip_target,lcfs=True):
    '''Relative differences between plasma integrals that must agree when they all use the cut-cell quadrature

    - `Ip_fit`, `dflux_fit`: the plasma current and diamagnetic flux fit constraints (`gs_itor_nl`, `gs_dflux`)
      against `get_globals` (`gs_comp_globals`); the same integrands over the same points, so equal to round-off
    - `Ip_target`: `get_globals` current against the target that the solve imposes on its source (`gs_source`)
    - `vol_lcfs`, `area_lcfs` (if `lcfs`): `get_globals` volume and the area from `compute_flux_integral`
      (`gs_flux_int`) against the region inside the traced LCFS contour, which differs only by the boundary
      discretization
    '''
    Ip, _, vol, _, dflux, _, _ = mygs.get_globals()
    myrecon = reconstruction(mygs)
    myrecon.set_Ip(Ip,1.0)
    myrecon.set_DFlux(dflux,1.0)
    myrecon.setup_constraints()
    con_vals = myrecon.eval_error()[:,1]
    myrecon.destroy_constraints()
    checks = {
        'Ip_fit': con_vals[0]/Ip-1.0,
        'dflux_fit': con_vals[1]/dflux-1.0,
        'Ip_target': Ip/Ip_target-1.0
    }
    if lcfs:
        area = mygs.compute_flux_integral(np.array([0.0,1.0]),np.array([1.0,1.0]))
        area_lcfs, vol_lcfs = lcfs_area_volume(mygs)
        checks['vol_lcfs'] = vol/vol_lcfs-1.0
        checks['area_lcfs'] = area/area_lcfs-1.0
    return checks


# Bounds for `plasma_integral_checks`, measured on the edge case at mesh resolutions 1 and 2: `Ip_fit` and
# `dflux_fit` 0 (bit for bit), `Ip_target` -1.2e-7 and -3.6e-7, `vol_lcfs` -1.9e-5 and -1.9e-6, `area_lcfs`
# 2.0e-5 and 1.0e-5. The point test in any one of these integrals moves its value by 6e-9 to 1.3e-4.
plasma_integral_bounds = {'Ip_fit': 1.E-10, 'dflux_fit': 1.E-10, 'Ip_target': 3.E-6, 'vol_lcfs': 5.E-5,
                          'area_lcfs': 5.E-5}


def check_plasma_integrals(results,keys=None):
    '''Check the entries `keys` (default: all) of `plasma_integral_checks` against `plasma_integral_bounds`'''
    test_result = True
    for key in (plasma_integral_bounds if keys is None else keys):
        bound = plasma_integral_bounds[key]
        if abs(results[key]) >= bound:
            print('FAILED: {0} = {1:.3E}, bound {2:.1E}'.format(key,results[key],bound))
            test_result = False
    return test_result


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
    eq_info.update(plasma_integral_checks(mygs,15.6E6))
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


def run_DIIID_xpoint_pair_cut_cells(mp_q):
    # The case of test_TokaMaker_xpoint_cut.py (DIII-D example mesh, a second lower X-point 9 cm from the active one,
    # which skews the active X-point's core axis), solved with cut cells
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('DIIID_mesh.h5')
    mygs.setup_model(mesh_pts, mesh_lc, reg=mesh_reg, cond_dict=cond_dict, coil_dict=coil_dict, order=2, F0=2.96)
    mygs.set_coil_vsc({'F9A': 1.0, 'F9B': -1.0})
    mygs.set_targets(Ip=1.19E6,pax=6.0E4)
    x_points = np.array([[1.183,-1.227],[1.11,-1.17]])
    lcfs = np.array([[2.287,0.067],[2.275,0.205],[2.241,0.343],[2.187,0.487],[2.101,0.636],[1.972,0.78],
                     [1.793,0.908],[1.568,0.968],[1.346,0.904],[1.195,0.725],[1.128,0.512],[1.094,0.323],
                     [1.077,0.151],[1.076,-0.017],[1.083,-0.193],[1.09,-0.405],[1.117,-0.686],[1.542,-0.97],
                     [1.792,-0.768],[1.959,-0.607],[2.077,-0.465],[2.164,-0.331],[2.227,-0.201],[2.271,-0.069]])
    mygs.set_isoflux_constraints(np.vstack((lcfs,x_points)))
    mygs.set_saddle_constraints(x_points)
    reg_terms = [mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E-2) for name in mygs.coil_sets]
    reg_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E-2))
    mygs.set_coil_reg(reg_terms=reg_terms)
    psi_prof = np.array([0.0,0.1,0.3,0.5,0.7,0.85,0.925,0.95,0.975,1.0])
    ffp_prof = np.array([1.0,0.571,0.355,0.233,0.155,0.112,0.071,-0.002,-0.205,0.0])
    pp_prof = np.array([1.0,1.868,1.758,1.235,0.881,0.75,0.91,1.642,3.084,0.0])
    mygs.set_profiles(ffp_prof={'type': 'linterp', 'x': psi_prof, 'y': ffp_prof},
                      pp_prof={'type': 'linterp', 'x': psi_prof, 'y': pp_prof})
    mygs.settings.lcfs_cut_cells = True
    mygs.update_settings()
    mygs.init_psi(1.7,0.0,0.6,1.7,0.4)
    try:
        _, nl_its = mygs.solve(return_its=True)
    except ValueError:
        mp_q.put(None)
        return
    eq_info = {'nl_its': nl_its, 'diverted': mygs.diverted, 'nx_points': len(mygs.get_xpoints()[0])}
    eq_info.update(plasma_integral_checks(mygs,1.19E6,lcfs=False))
    # Plasma volume and area on the same psi with the point test, which applies both X-point cuts
    # (`gs_test_bounds`), against the cut-cell quadrature
    area_fun = (np.array([0.0,1.0]),np.array([1.0,1.0]))
    vol_cut, area_cut = mygs.get_globals()[2], mygs.compute_flux_integral(*area_fun)
    mygs.settings.lcfs_cut_cells = False
    mygs.update_settings()
    vol_pt, area_pt = mygs.get_globals()[2], mygs.compute_flux_integral(*area_fun)
    eq_info['vol_point'] = vol_cut/vol_pt-1.0
    eq_info['area_point'] = area_cut/area_pt-1.0
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
    with psi: in a small vertical tilt of psi no step of the plasma volume stands out from the trend. The fit
    constraint integrals, the reported globals, the flux-function integral and the solved source must also agree
    with each other and with the traced LCFS (`plasma_integral_checks`), as they do when all use the same quadrature'''
    results = mp_run(run_ITER_edge_case,(1.0,))
    assert results is not None
    assert results[0]['nl_its'] > 0
    assert validate_dict(results,ITER_edge_ref_dict)
    # A jump from quadrature points crossing the LCFS would make one step many times the median
    assert results[0]['vol_step_ratio'] < 3.0
    assert check_plasma_integrals(results[0])


@pytest.mark.slow
def test_ITER_edge_cut_cells_reference():
    '''The refined-mesh reference values of test_ITER_edge_cut_cells, from a live solve'''
    results = mp_run(run_ITER_edge_case,(2.0,),timeout=180)
    assert results is not None
    assert results[0]['nl_its'] > 0
    assert validate_dict(results,ITER_edge_ref_dict)
    assert results[0]['vol_step_ratio'] < 3.0
    assert check_plasma_integrals(results[0])


@pytest.mark.coverage
def test_DIIID_xpoint_pair_cut_cells():
    '''The X-point pair case of test_TokaMaker_xpoint_cut.py with `lcfs_cut_cells`. The second cut at the active
    X-point (across its core axis) removes a patch of private flux that the cut from the O-point alone lets into
    the plasma region; the cut-cell quadrature must apply it too. With both cuts the cut-cell plasma volume and
    area differ from those of the point test on the same psi by 8e-6 and 1.1e-5 (the boundary discretization on
    this mesh); integrating the patch adds about 8e-5 and 1.4e-4, and dropping the cut from the O-point alone
    adds 4e-5 and 8e-5.'''
    results = mp_run(run_DIIID_xpoint_pair_cut_cells,(),timeout=60)
    assert results is not None
    results = results[0]
    assert results['nl_its'] > 0
    assert results['diverted'] and (results['nx_points'] > 1)
    assert check_plasma_integrals(results,('Ip_fit','dflux_fit','Ip_target'))
    assert abs(results['vol_point']) < 4.E-5
    assert abs(results['area_point']) < 4.E-5
