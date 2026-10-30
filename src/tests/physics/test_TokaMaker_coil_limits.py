#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Regression tests for TokaMaker post-solve coil limit checks

All limits used here are illustrative values chosen relative to the computed
equilibrium; they are not engineering limits for any machine.
'''
import os
import sys
import time
import json
import multiprocessing
import numpy as np
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..', 'python')))
from OpenFUSIONToolkit.TokaMaker.coil_limits import (CoilLimitReport, check_coil_limits, check_keepout,
                                                     compute_coil_forces, evaluate_limits)


def mp_run(target, args, timeout=60):
    '''! Run a test case in a separate process (only one OFT_env is allowed per process)

    @param target Function to run, the last argument is a result queue
    @param args Arguments to pass to `target`
    @param timeout Timeout [s]
    @result Value put in the queue by `target`, or `None` on failure
    '''
    if os.environ.get('OFT_DEBUG_TEST', 0):
        timeout *= 4
    os.chdir(test_dir)
    mp_q = multiprocessing.Queue()
    p = multiprocessing.Process(target=target, args=args + (mp_q,))
    p.start()
    start = time.time()
    while time.time() - start <= timeout:
        if not p.is_alive():
            break
        time.sleep(.5)
    else:
        print("Timeout reached")
        p.terminate()
        p.join()
        return None
    try:
        test_result = mp_q.get(timeout=5)
    except:
        print("Failed to get output")
        return None
    p.join()
    return test_result


#============================================================================
# Limit evaluation (no equilibrium needed)
#============================================================================
def synthetic_quantities():
    '''! Build a small set of coil quantities by hand

    @result Quantities dictionary in the format of `compute_coil_quantities`
    '''
    def coil(current, nturns, area, force_r, force_z, coil_set):
        return {'coil_set': coil_set, 'current': current, 'nturns': nturns, 'ampere_turns': current*nturns,
                'current_density': current*nturns/area, 'force_r': force_r, 'force_z': force_z, 'area': area}
    coils = {
        'A': coil(10.0E3, 100.0, 0.5, 2.0E6, -3.0E6, 'A'),
        'BU': coil(-4.0E3, 50.0, 0.1, 1.0E6, 5.0E6, 'B'),
        'BL': coil(-4.0E3, -50.0, 0.1, 1.0E6, -1.0E6, 'B'),
    }
    coil_sets = {
        'A': dict(coils['A']),
        'B': {'current': -4.0E3, 'ampere_turns': 0.0, 'current_density': -2.0E6,
              'force_r': 2.0E6, 'force_z': 4.0E6, 'area': 0.2},
    }
    return {'coils': coils, 'coil_sets': coil_sets, 'energy_coils': 1.0E8, 'energy_total': 1.2E8}


def find_check(checks, name, quantity):
    '''! Get the check record for a name and quantity

    @param checks List of check records
    @param name Coil, group, or "system"
    @param quantity Quantity name
    @result Matching check record
    '''
    for check in checks:
        if (check['name'] == name) and (check['quantity'] == quantity):
            return check
    raise KeyError((name, quantity))


def test_limit_fractions():
    '''! Fraction of limit used for scalar, asymmetric, and one-sided limits'''
    quantities = synthetic_quantities()
    checks = evaluate_limits(quantities, coil_limits={
        'A': {'current': 20.0E3, 'force_z': [-6.0E6, 1.0E6]},
        'B': {'current': [-5.0E3, 0.0], 'force_z': [None, 8.0E6]},
        'BU': {'current_density': 1.5E6},
    }, global_limits={'energy_coils': 4.0E8})
    assert find_check(checks, 'A', 'current')['fraction_used'] == pytest.approx(0.5)
    assert find_check(checks, 'A', 'force_z')['fraction_used'] == pytest.approx(0.5)
    assert find_check(checks, 'A', 'force_z')['limit'] == pytest.approx(-6.0E6)
    assert find_check(checks, 'B', 'current')['fraction_used'] == pytest.approx(0.8)
    assert find_check(checks, 'B', 'force_z')['fraction_used'] == pytest.approx(0.5)
    assert find_check(checks, 'BU', 'current_density')['fraction_used'] == pytest.approx(4.0/3.0)
    assert find_check(checks, 'system', 'energy_coils')['fraction_used'] == pytest.approx(0.25)
    report = CoilLimitReport(checks)
    assert not report.ok
    assert [(c['name'], c['quantity']) for c in report.violations] == [('BU', 'current_density')]
    assert report.binding['name'] == 'BU'
    assert 'Binding limit: BU current_density' in report.summary()


def test_limit_wrong_sign():
    '''! A value on the side of a zero limit (eg. unipolar supply) is flagged'''
    checks = evaluate_limits(synthetic_quantities(), coil_limits={'B': {'current': [0.0, 5.0E3]}})
    assert np.isinf(checks[0]['fraction_used'])
    assert not checks[0]['ok']


def test_group_limits():
    '''! Group sums with default and explicit weights'''
    quantities = synthetic_quantities()
    checks = evaluate_limits(quantities, group_limits={
        'all_z': {'coils': ['A', 'BU', 'BL'], 'force_z': 10.0E6},
        'split_z': {'coils': {'BU': 1.0, 'BL': -1.0}, 'force_z': 10.0E6},
    })
    assert find_check(checks, 'all_z', 'force_z')['value'] == pytest.approx(1.0E6)
    assert find_check(checks, 'split_z', 'force_z')['value'] == pytest.approx(6.0E6)
    assert find_check(checks, 'split_z', 'force_z')['fraction_used'] == pytest.approx(0.6)


@pytest.mark.parametrize('kwargs', [
    {'coil_limits': {'missing': {'current': 1.0}}},
    {'coil_limits': {'A': {'voltage': 1.0}}},
    {'coil_limits': {'A': {'current': -1.0}}},
    {'coil_limits': {'A': {'current': [1.0, 2.0]}}},
    {'group_limits': {'g': {'force_z': 1.0}}},
    {'group_limits': {'g': {'coils': ['A'], 'current': 1.0}}},
    {'global_limits': {'energy': 1.0}},
])
def test_limit_bad_specs(kwargs):
    '''! Invalid limit specifications raise errors

    @param kwargs Arguments to `evaluate_limits`
    '''
    with pytest.raises((KeyError, ValueError)):
        evaluate_limits(synthetic_quantities(), **kwargs)


def test_keepout_geometry():
    '''! Signed clearance for points inside and outside a polygon'''
    square = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    outside = np.array([[1.5, 0.5], [2.0, 2.0]])
    records = check_keepout(outside, {'sq': square, 'sq_margin': {'polygon': square, 'min_clearance': 0.6}})
    assert records[0]['n_inside'] == 0
    assert records[0]['min_clearance'] == pytest.approx(0.5)
    assert records[0]['ok']
    assert not records[1]['ok']
    inside = np.vstack((outside, [[0.5, 0.2]]))
    records = check_keepout(inside, {'sq': np.vstack((square, square[:1, :]))})
    assert records[0]['n_inside'] == 1
    assert records[0]['min_clearance'] == pytest.approx(-0.2)
    assert not records[0]['ok']


def test_plot_margins():
    '''! Margin bar plot runs'''
    matplotlib = pytest.importorskip('matplotlib')
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    checks = evaluate_limits(synthetic_quantities(), coil_limits={'A': {'current': 20.0E3}, 'B': {'force_z': 1.0E6}})
    fig, ax = plt.subplots()
    CoilLimitReport(checks).plot_margins(fig, ax)
    assert len(ax.patches) == 2
    plt.close(fig)


#============================================================================
# ITER-like equilibrium (public geometry from ITER_geom.json)
#============================================================================
def create_ITER_mesh():
    '''! Build the ITER test mesh (same as `test_TokaMaker.py`) if it does not exist'''
    from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain, save_gs_mesh
    with open('ITER_geom.json', 'r') as fid:
        ITER_geom = json.load(fid)
    gs_mesh = gs_Domain()
    gs_mesh.define_region('air', 0.6, 'boundary')
    gs_mesh.define_region('plasma', 0.15, 'plasma')
    gs_mesh.define_region('vacuum1', 0.3, 'vacuum')
    gs_mesh.define_region('vacuum2', 0.3, 'vacuum')
    gs_mesh.define_region('vv1', 0.3, 'conductor', eta=6.9E-7)
    gs_mesh.define_region('vv2', 0.3, 'conductor', eta=6.9E-7)
    for key in ITER_geom['coils']:
        if not key.startswith('VS'):
            gs_mesh.define_region(key, 0.2, 'coil')
    gs_mesh.define_region('VSU', 0.2, 'coil', coil_set='VS', nTurns=1.0)
    gs_mesh.define_region('VSL', 0.2, 'coil', coil_set='VS', nTurns=-1.0)
    gs_mesh.add_polygon(ITER_geom['limiter'], 'plasma', parent_name='vacuum1')
    gs_mesh.add_annulus(ITER_geom['inner_vv'][0], 'vacuum1', ITER_geom['inner_vv'][1], 'vv1', parent_name='vacuum2')
    gs_mesh.add_annulus(ITER_geom['outer_vv'][0], 'vacuum2', ITER_geom['outer_vv'][1], 'vv2', parent_name='air')
    for key, coil in ITER_geom['coils'].items():
        parent = 'vacuum1' if key.startswith('VS') else 'air'
        gs_mesh.add_rectangle(coil['rc'], coil['zc'], coil['w'], coil['h'], key, parent_name=parent)
    mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
    save_gs_mesh(mesh_pts, mesh_lc, mesh_reg, gs_mesh.get_coils(), gs_mesh.get_conductors(), 'ITER_mesh.h5')


def solve_ITER():
    '''! Solve the ITER-like test equilibrium

    @result TokaMaker object with solved equilibrium
    '''
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
    if not os.path.exists('ITER_mesh.h5'):
        create_ITER_mesh()
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh('ITER_mesh.h5')
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.setup(order=2, F0=5.3*6.2)
    mygs.set_coil_vsc({'VS': 1.0})
    mygs.set_coil_bounds({key: [-50.E6, 50.E6] for key in mygs.coil_sets})
    mygs.set_targets(Ip=15.6E6, pax=6.2E5)
    isoflux_pts = np.array([
        [8.20, 0.41], [8.06, 1.46], [7.51, 2.62], [6.14, 3.78], [4.51, 3.02],
        [4.26, 1.33], [4.28, 0.08], [4.49, -1.34], [7.28, -1.89], [8.00, -0.68]
    ])
    x_point = np.array([[5.125, -3.4],])
    mygs.set_isoflux(np.vstack((isoflux_pts, x_point)))
    mygs.set_saddles(x_point)
    reg_terms = []
    for name in mygs.coil_sets:
        weight = 2.E-2 if name.startswith('CS1') else 1.E-2
        reg_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=weight))
    reg_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E2))
    mygs.set_coil_reg(reg_terms=reg_terms)
    mygs.set_profiles(ffp_prof=create_power_flux_fun(40, 1.5, 2.0), pp_prof=create_power_flux_fun(40, 4.0, 1.0))
    mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
    mygs.solve()
    return mygs


def centroid_forces(mygs):
    r'''! Independent coil force calculation using one-point (centroid) quadrature per cell

    Uses the force density convention from `compute_forces_components` usage
    (\f$ f_r = J_\phi B_z \f$, \f$ f_z = -J_\phi B_r \f$).

    @param mygs TokaMaker object with solved equilibrium
    @result `{coil_name: [F_r, F_z]}` [N]
    '''
    _, region_currents = mygs.get_coil_currents()
    field_eval = mygs.get_field_eval('B')
    forces = {}
    for name, coil in mygs._coil_dict.items():
        cells = mygs.lc[mygs.reg == coil['reg_id'], :]
        p0, p1, p2 = mygs.r[cells[:, 0], :2], mygs.r[cells[:, 1], :2], mygs.r[cells[:, 2], :2]
        area = 0.5*np.abs((p1[:, 0]-p0[:, 0])*(p2[:, 1]-p0[:, 1]) - (p2[:, 0]-p0[:, 0])*(p1[:, 1]-p0[:, 1]))
        centroids = (p0 + p1 + p2)/3.0
        B = field_eval.eval(centroids)
        J_phi = region_currents[coil['reg_id']-1]/area.sum()
        forces[name] = [
            2.0*np.pi*J_phi*np.sum(B[:, 2]*centroids[:, 0]*area),
            -2.0*np.pi*J_phi*np.sum(B[:, 0]*centroids[:, 0]*area),
        ]
    return forces


def run_ITER_limits(mp_q):
    '''! Run coil limit checks on the ITER-like equilibrium

    @param mp_q Result queue
    '''
    try:
        mygs = solve_ITER()
        report = check_coil_limits(mygs)
        quantities = report.quantities
        # Illustrative limits relative to the solution: PF6 current at 80% (binding),
        # CS2L vertical force at 110% (violated), other currents at 50%
        coil_limits = {}
        for name, values in quantities['coil_sets'].items():
            coil_limits[name] = {'current': 2.0*abs(values['current'])}
        coil_limits['PF6']['current'] = abs(quantities['coil_sets']['PF6']['current'])/0.8
        coil_limits['CS2L']['force_z'] = abs(quantities['coil_sets']['CS2L']['force_z'])/1.1
        coil_limits['VSU'] = {'current_density': 2.0*abs(quantities['coils']['VSU']['current_density'])}
        upper_cs = ['CS1U', 'CS2U', 'CS3U']
        group_limits = {'CS_upper': {'coils': upper_cs, 'force_z': 1.0E9}}
        global_limits = {'energy_total': 2.0*quantities['energy_total']}
        keepout = {
            'around_xpoint': np.array([[4.8, -3.7], [5.4, -3.7], [5.4, -3.1], [4.8, -3.1]]),
            'outboard': {'polygon': np.array([[10.0, -1.0], [11.0, -1.0], [11.0, 1.0], [10.0, 1.0]]), 'min_clearance': 0.5},
        }
        report = check_coil_limits(mygs, coil_limits=coil_limits, group_limits=group_limits,
                                   global_limits=global_limits, keepout=keepout)
        # Independent checks
        psi = mygs.get_psi(normalized=False)
        J_plasma = mygs.calc_jtor_plasma()
        energy_direct = mygs.compute_area_integral(psi*J_plasma, reg_mask=1)
        for name, values in quantities['coils'].items():
            energy_direct += values['current_density']*mygs.compute_area_integral(psi, reg_mask=mygs._coil_dict[name]['reg_id'])
        energy_direct *= np.pi  # psi is in Wb/rad
        J_delstar = mygs.calc_delstar_curr(psi)
        delstar_avg = {}
        for name, values in quantities['coils'].items():
            reg_id = mygs._coil_dict[name]['reg_id']
            delstar_avg[name] = mygs.compute_area_integral(J_delstar, reg_mask=reg_id)/values['area']
        result = {
            'checks': report.checks,
            'binding': report.binding,
            'ok': report.ok,
            'summary': report.summary(),
            'keepout': report.keepout,
            'quantities': report.quantities,
            'forces': compute_coil_forces(mygs),
            'centroid_forces': centroid_forces(mygs),
            'energy_direct': energy_direct,
            'delstar_avg': delstar_avg,
            'limits': coil_limits,
        }
    except Exception as e:
        print(e)
        result = None
    mp_q.put(result)


@pytest.mark.coverage
def test_ITER_coil_limits():
    '''! Coil limit checks on the ITER-like equilibrium'''
    results = mp_run(run_ITER_limits, ())
    assert results is not None
    checks = results['checks']
    quantities = results['quantities']
    # Margins: fraction used is |value|/limit for symmetric limits
    for check in checks:
        if check['kind'] == 'coil' and check['quantity'] == 'current':
            expected = abs(quantities['coil_sets'][check['name']]['current'])/results['limits'][check['name']]['current']
            assert check['fraction_used'] == pytest.approx(expected, rel=1.E-10)
    assert find_check(checks, 'PF6', 'current')['fraction_used'] == pytest.approx(0.8, rel=1.E-10)
    # Over-limit check is flagged and is the binding limit
    over = find_check(checks, 'CS2L', 'force_z')
    assert over['fraction_used'] == pytest.approx(1.1, rel=1.E-10)
    assert not over['ok']
    assert not results['ok']
    assert results['binding']['name'] == 'CS2L'
    assert results['binding']['quantity'] == 'force_z'
    assert 'Binding limit: CS2L force_z' in results['summary']
    # Group sum
    group = find_check(checks, 'CS_upper', 'force_z')
    assert group['value'] == pytest.approx(sum(quantities['coil_sets'][name]['force_z'] for name in ['CS1U', 'CS2U', 'CS3U']))
    # Current density sign and magnitude agree with Del*(psi) current in coil regions
    for name in ['PF1', 'PF6', 'CS1U']:
        assert results['delstar_avg'][name] == pytest.approx(quantities['coils'][name]['current_density'], rel=0.1)
    # Forces agree with an independent centroid quadrature of J x B
    force_scale = max(abs(np.array(list(results['forces'][name].values()))).max() for name in results['forces'])
    for name, forces in results['forces'].items():
        assert abs(forces['force_r'] - results['centroid_forces'][name][0]) < 1.E-2*force_scale
        assert abs(forces['force_z'] - results['centroid_forces'][name][1]) < 1.E-2*force_scale
    # Net vertical force on all coils vanishes (plasma is in vertical force balance)
    total_fz = sum(values['force_z'] for values in quantities['coils'].values())
    assert abs(total_fz) < 1.E-3*force_scale
    # Stored energy from inductance matrix matches pi*int(psi*J dA)
    assert quantities['energy_total'] == pytest.approx(results['energy_direct'], rel=1.E-3)
    assert quantities['energy_coils'] < quantities['energy_total']
    # Keep-out: separatrix legs pass through the box around the X-point, outboard box is clear
    keepout = {rec['name']: rec for rec in results['keepout']}
    assert keepout['around_xpoint']['n_inside'] > 0
    assert keepout['around_xpoint']['min_clearance'] < 0.0
    assert not keepout['around_xpoint']['ok']
    assert keepout['outboard']['n_inside'] == 0
    assert keepout['outboard']['min_clearance'] > 0.5
    assert keepout['outboard']['ok']
