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
import warnings
import multiprocessing
import numpy as np
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..', 'python')))
from OpenFUSIONToolkit.TokaMaker.coil_limits import CoilLimitReport, check_coil_limits, evaluate_limits

mu0 = 4.E-7*np.pi


def mp_run(target, args, timeout=60):
    '''Run a test case in a separate process (only one OFT_env is allowed per process)'''
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
    '''Small set of hand-built coil quantities, in the format of `compute_coil_quantities`'''
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
    '''Get the check record for a name and quantity'''
    for check in checks:
        if (check['name'] == name) and (check['quantity'] == quantity):
            return check
    raise KeyError((name, quantity))


def test_limit_fractions():
    '''Fraction of limit used for scalar, asymmetric, and one-sided limits'''
    quantities = synthetic_quantities()
    checks = evaluate_limits(quantities, coil_limits={
        'A': {'current': 20.0E3, 'force_z': [-6.0E6, 1.0E6], 'force_r': None},
        'B': {'current': [-5.0E3, 0.0], 'force_z': [None, 8.0E6]},
        'BU': {'current_density': np.array(1.5E6)},
    }, global_limits={'energy_coils': 4.0E8, 'energy_total': None})
    assert len(checks) == 6  # `None` limits are skipped
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
    '''A value on the side of a zero limit (eg. unipolar supply) is flagged'''
    checks = evaluate_limits(synthetic_quantities(), coil_limits={'B': {'current': [0.0, 5.0E3]}})
    assert np.isinf(checks[0]['fraction_used'])
    assert not checks[0]['ok']


def test_nonfinite_value():
    '''Non-finite values fail loudly instead of passing as 0% used'''
    quantities = synthetic_quantities()
    quantities['energy_total'] = np.nan
    with pytest.raises(ValueError, match='not finite'):
        evaluate_limits(quantities, global_limits={'energy_total': 1.0E9})


def test_group_limits():
    '''Group sums with default and explicit weights'''
    quantities = synthetic_quantities()
    checks = evaluate_limits(quantities, group_limits={
        'all_z': {'coils': ['A', 'BU', 'BL'], 'force_z': 10.0E6},
        'split_z': {'coils': {'BU': 1.0, 'BL': -1.0}, 'force_z': 10.0E6},
    })
    assert find_check(checks, 'all_z', 'force_z')['value'] == pytest.approx(1.0E6)
    assert find_check(checks, 'split_z', 'force_z')['value'] == pytest.approx(6.0E6)
    assert find_check(checks, 'split_z', 'force_z')['fraction_used'] == pytest.approx(0.6)


def test_group_double_count_warning():
    '''Listing a coil set and one of its regions in a group warns'''
    with pytest.warns(UserWarning, match='twice'):
        evaluate_limits(synthetic_quantities(), group_limits={'g': {'coils': ['B', 'BU'], 'force_z': 1.0E7}})
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        evaluate_limits(synthetic_quantities(), group_limits={'g': {'coils': ['A', 'BU', 'BL'], 'force_z': 1.0E7}})


@pytest.mark.parametrize('kwargs', [
    {'coil_limits': {'missing': {'current': 1.0}}},
    {'coil_limits': {'A': {'voltage': 1.0}}},
    {'coil_limits': {'A': {'current': -1.0}}},
    {'coil_limits': {'A': {'current': np.nan}}},
    {'coil_limits': {'A': {'current': 'big'}}},
    {'coil_limits': {'A': {'current': [1.0, 2.0]}}},
    {'coil_limits': {'A': {'current': [1.0, 2.0, 3.0]}}},
    {'coil_limits': {'A': {'current': ['low', 2.0]}}},
    {'group_limits': {'g': {'force_z': 1.0}}},
    {'group_limits': {'g': {'coils': ['A'], 'current': 1.0}}},
    {'global_limits': {'energy': 1.0}},
])
def test_limit_bad_specs(kwargs):
    '''Invalid limit specifications raise errors'''
    with pytest.raises((KeyError, ValueError)):
        evaluate_limits(synthetic_quantities(), **kwargs)


def test_plot_margins():
    '''Margin bar plot runs'''
    matplotlib = pytest.importorskip('matplotlib')
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    checks = evaluate_limits(synthetic_quantities(), coil_limits={'A': {'current': 20.0E3}, 'B': {'force_z': 1.0E6}})
    fig, ax = plt.subplots()
    CoilLimitReport(checks).plot_margins(fig, ax)
    assert len(ax.patches) == 2
    plt.close(fig)


#============================================================================
# Vacuum fields from simple coils
#============================================================================
## Coil geometry for the vacuum test mesh: name: (rc, zc, w, h, nturns, coil_set)
VAC_COILS = {
    'C5': (3.0, 3.5, 0.6, 0.6, 5.0, 'C5'),
    'PU': (1.5, 1.0, 0.4, 0.8, 1.0, 'PU'),
    'PL': (1.5, -1.0, 0.4, 0.8, 1.0, 'PL'),
    'MIXU': (6.0, 2.0, 0.5, 0.5, 3.0, 'MIX'),
    'MIXL': (6.0, -2.0, 0.5, 0.5, -2.0, 'MIX'),
}


def hoop_force(rc, w, h, ampere_turns):
    '''Analytic total radial self-force of a thin ring with rectangular cross-section

    F_r = mu0 I^2/2 (ln(8R/g) - 1), with geometric mean distance g ~ 0.2235(w+h)
    '''
    gmd = 0.2235*(w + h)
    return 0.5*mu0*ampere_turns**2*(np.log(8.0*rc/gmd) - 1.0)


def run_vacuum_coils(mp_q):
    '''Vacuum solves with single coils, coil pairs, and a mixed-turn coil set'''
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit._interface import oftpy_dump_cov
        from OpenFUSIONToolkit.TokaMaker import TokaMaker
        from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
        from OpenFUSIONToolkit.TokaMaker.coil_limits import compute_coil_quantities
        gs_mesh = gs_Domain(rextent=8.0, zextents=[-5.0, 5.0])
        gs_mesh.define_region('air', 0.4, 'boundary')
        gs_mesh.define_region('plasma', 0.2, 'plasma')
        for name, (_, _, _, _, nturns, coil_set) in VAC_COILS.items():
            gs_mesh.define_region(name, 0.05, 'coil', coil_set=coil_set, nTurns=nturns)
        gs_mesh.add_rectangle(4.0, 0.0, 1.5, 2.0, 'plasma')
        for name, (rc, zc, w, h, _, _) in VAC_COILS.items():
            gs_mesh.add_rectangle(rc, zc, w, h, name)
        mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
        myOFT = OFT_env(nthreads=-1)
        mygs = TokaMaker(myOFT)
        mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
        mygs.setup_regions(cond_dict=gs_mesh.get_conductors(), coil_dict=gs_mesh.get_coils())
        mygs.setup(order=2)
        results = {}
        # After `vac_solve` get_globals reports a profile Ip that is not in psi;
        # after `solve(vacuum=True)` it reports Ip = 0. Both must give coil-only energy.
        cases = {
            'single': ({'C5': 2.0E5}, 'vac_solve'),
            'single_solve': ({'C5': 2.0E5}, 'solve'),
            'same_sign': ({'PU': 1.0E6, 'PL': 1.0E6}, 'vac_solve'),
            'opposite_sign': ({'PU': 1.0E6, 'PL': -1.0E6}, 'vac_solve'),
            'mixed': ({'MIX': 1.0E5}, 'vac_solve'),
        }
        for case, (currents, method) in cases.items():
            mygs.set_coil_currents(currents)
            if method == 'solve':
                mygs.solve(vacuum=True)
            else:
                mygs.replace_eq(mygs.vac_solve())
            quantities = compute_coil_quantities(mygs)
            psi = mygs.get_psi(normalized=False)
            J_delstar = mygs.calc_delstar_curr(psi)
            energy_direct = 0.0
            delstar_turns = {}
            for name, values in quantities['coils'].items():
                reg_id = mygs._coil_dict[name]['reg_id']
                energy_direct += values['current_density']*mygs.compute_area_integral(psi, reg_mask=reg_id)
                delstar_turns[name] = mygs.compute_area_integral(J_delstar, reg_mask=reg_id)
            report = check_coil_limits(mygs, global_limits={'energy_total': 2.0*quantities['energy_total']})
            results[case] = {
                'quantities': quantities,
                'energy_direct': np.pi*energy_direct,  # psi is in Wb/rad
                'delstar_turns': delstar_turns,
                'Ip': mygs.get_globals()[0],
                'energy_check': report.checks[0],
            }
        # Non-uniform current distributions are not supported
        mygs.set_coil_current_dist('C5', 2.0*np.ones((mygs.np,)))
        try:
            check_coil_limits(mygs)
            results['dist_error'] = None
        except NotImplementedError as e:
            results['dist_error'] = str(e)
        oftpy_dump_cov()
    except Exception as e:
        print(e)
        results = None
    mp_q.put(results)


@pytest.mark.coverage
def test_vacuum_coils():
    '''Coil bookkeeping, forces, and energy for vacuum fields from simple coils'''
    results = mp_run(run_vacuum_coils, ())
    assert results is not None
    assert results['single']['Ip'] != 0.0
    assert results['single_solve']['Ip'] == 0.0
    for case in ('single', 'single_solve', 'same_sign', 'opposite_sign', 'mixed'):
        res = results[case]
        # Vacuum solve: no plasma current in psi, energy_total is finite and equals energy_coils
        assert np.isfinite(res['quantities']['energy_total'])
        assert res['quantities']['energy_total'] == res['quantities']['energy_coils']
        assert res['quantities']['energy_coils'] == pytest.approx(res['energy_direct'], rel=1.E-3)
        assert res['energy_check']['fraction_used'] == pytest.approx(0.5)
        # Region A-turns match the integral of the Del*(psi) current (smeared over region edges)
        for name, values in res['quantities']['coils'].items():
            if values['current'] != 0.0:
                assert res['delstar_turns'][name] == pytest.approx(values['ampere_turns'], rel=0.1)
    # Single 5-turn coil at 200 kA: A-turns and hoop force for 1 MA-turn
    c5 = results['single']['quantities']['coils']['C5']
    rc, _, w, h, _, _ = VAC_COILS['C5']
    assert c5['ampere_turns'] == pytest.approx(1.0E6)
    assert c5['current_density'] == pytest.approx(1.0E6/c5['area'])
    assert c5['force_r'] == pytest.approx(hoop_force(rc, w, h, 1.0E6), rel=1.E-2)
    c5_solve = results['single_solve']['quantities']['coils']['C5']
    assert c5_solve['force_r'] == pytest.approx(c5['force_r'], rel=1.E-6)
    # Coaxial pair: attract with same-sign currents, repel with opposite signs, equal and opposite
    for case, sign in (('same_sign', -1.0), ('opposite_sign', 1.0)):
        coils = results[case]['quantities']['coils']
        assert np.sign(coils['PU']['force_z']) == sign
        assert np.sign(coils['PL']['force_z']) == -sign
        assert coils['PL']['force_z'] == pytest.approx(-coils['PU']['force_z'], rel=1.E-2)
    # Mixed-turn set (3 and -2 turns) at 100 kA per turn
    mixed = results['mixed']['quantities']
    upper, lower, mix_set = mixed['coils']['MIXU'], mixed['coils']['MIXL'], mixed['coil_sets']['MIX']
    assert upper['ampere_turns'] == pytest.approx(3.0E5)
    assert lower['ampere_turns'] == pytest.approx(-2.0E5)
    assert upper['current_density'] == pytest.approx(3.0E5/upper['area'])
    assert lower['current_density'] == pytest.approx(-2.0E5/lower['area'])
    assert mix_set['current'] == pytest.approx(1.0E5)
    assert mix_set['ampere_turns'] == pytest.approx(1.0E5)
    assert mix_set['current_density'] == pytest.approx(upper['current_density'])
    # Opposite A-turns repel: upper pushed up, lower pushed down, equal and opposite
    assert upper['force_z'] > 0.0
    assert lower['force_z'] == pytest.approx(-upper['force_z'], rel=1.E-2)
    # Current distributions raise
    assert results['dist_error'] is not None
    assert 'C5' in results['dist_error']


#============================================================================
# ITER-like equilibrium (public geometry from ITER_geom.json)
#============================================================================
def create_ITER_mesh():
    '''Build the ITER test mesh (same as `test_TokaMaker.py`) if it does not exist'''
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
    '''Solve the ITER-like test equilibrium and return the TokaMaker object'''
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
    mygs.set_isoflux_constraints(np.vstack((isoflux_pts, x_point)))
    mygs.set_saddle_constraints(x_point)
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


def run_ITER_limits(mp_q):
    '''Run coil limit checks on the ITER-like equilibrium'''
    try:
        from OpenFUSIONToolkit._interface import oftpy_dump_cov
        mygs = solve_ITER()
        quantities = check_coil_limits(mygs).quantities
        # Illustrative limits relative to the solution: PF6 current at 80%,
        # CS2L vertical force at 110% (violated, binding), other currents at 50%
        coil_limits = {}
        for name, values in quantities['coil_sets'].items():
            coil_limits[name] = {'current': 2.0*abs(values['current'])}
        coil_limits['PF6']['current'] = abs(quantities['coil_sets']['PF6']['current'])/0.8
        coil_limits['CS2L']['force_z'] = abs(quantities['coil_sets']['CS2L']['force_z'])/1.1
        coil_limits['VSU'] = {'current_density': 2.0*abs(quantities['coils']['VSU']['current_density'])}
        group_limits = {'CS_upper': {'coils': ['CS1U', 'CS2U', 'CS3U'], 'force_z': 1.0E9}}
        global_limits = {'energy_total': 2.0*quantities['energy_total']}
        report = check_coil_limits(mygs, coil_limits=coil_limits, group_limits=group_limits,
                                   global_limits=global_limits)
        # Independent checks: energy from pi*int(psi*J dA) and Del*(psi) current in coils
        psi = mygs.get_psi(normalized=False)
        energy_direct = mygs.compute_area_integral(psi*mygs.calc_jtor_plasma(), reg_mask=1)
        J_delstar = mygs.calc_delstar_curr(psi)
        delstar_avg = {}
        for name, values in quantities['coils'].items():
            reg_id = mygs._coil_dict[name]['reg_id']
            energy_direct += values['current_density']*mygs.compute_area_integral(psi, reg_mask=reg_id)
            delstar_avg[name] = mygs.compute_area_integral(J_delstar, reg_mask=reg_id)/values['area']
        result = {
            'checks': report.checks,
            'binding': report.binding,
            'ok': report.ok,
            'summary': report.summary(),
            'quantities': report.quantities,
            'energy_direct': np.pi*energy_direct,  # psi is in Wb/rad
            'delstar_avg': delstar_avg,
            'limits': coil_limits,
        }
        oftpy_dump_cov()
    except Exception as e:
        print(e)
        result = None
    mp_q.put(result)


@pytest.mark.coverage
def test_ITER_coil_limits():
    '''Coil limit checks on the ITER-like equilibrium'''
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
    # Net vertical force on all coils vanishes (plasma is in vertical force balance),
    # relative to the largest coil force
    largest_force = max(max(abs(values['force_r']), abs(values['force_z'])) for values in quantities['coils'].values())
    total_fz = sum(values['force_z'] for values in quantities['coils'].values())
    assert abs(total_fz) < 1.E-3*largest_force
    # Stored energy from inductance matrix matches pi*int(psi*J dA)
    assert quantities['energy_total'] == pytest.approx(results['energy_direct'], rel=1.E-3)
