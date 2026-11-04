#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Tests for TokaMaker divertor geometry metrics (`TokaMaker.divertor`)'''
import os
import sys
import pytest
import numpy as np
from contourpy import contour_generator
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..', 'python')))
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
from OpenFUSIONToolkit.TokaMaker.divertor import compute_divertor_metrics, _divertor_metrics
from test_TokaMaker import mp_run, create_ITER_mesh


#------------------------------------------------------------------------------
# Independent checks shared by the analytic and ITER tests
#------------------------------------------------------------------------------
def contour_L_par(psi_fun, B_fun, psi_level, start, end, R_range, Z_range, dx):
    '''L_par from `start` to `end` by integrating |B|/B_p along a gridded contour of `psi_level`

    The contour comes from contourpy on a uniform (R, Z) grid, so it is independent of the
    field-line tracer in `divertor.py`.
    '''
    R = np.arange(R_range[0], R_range[1] + 0.5*dx, dx)
    Z = np.arange(Z_range[0], Z_range[1] + 0.5*dx, dx)
    RR, ZZ = np.meshgrid(R, Z)
    psi = psi_fun(np.column_stack((RR.ravel(), ZZ.ravel()))).reshape(RR.shape)
    lines = contour_generator(R, Z, np.ma.masked_invalid(psi)).lines(psi_level)
    line = min(lines, key=lambda pts: np.hypot(pts[:, 0]-start[0], pts[:, 1]-start[1]).min())
    i0 = np.hypot(line[:, 0]-start[0], line[:, 1]-start[1]).argmin()
    i1 = np.hypot(line[:, 0]-end[0], line[:, 1]-end[1]).argmin()
    if i0 > i1:
        i0, i1 = i1, i0
        start, end = end, start
    pts = np.vstack((start, line[i0+1:i1, :], end))
    mid = 0.5*(pts[1:, :] + pts[:-1, :])
    B = B_fun(mid)
    return (np.linalg.norm(B, axis=1)/np.hypot(B[:, 0], B[:, 2])*np.linalg.norm(np.diff(pts, axis=0), axis=1)).sum()


def strike_separation_ratio(metrics, end):
    '''Along-wall separation of a SOL strike point from its separatrix strike point, mapped to the
    flux-surface spacing (separation*sin(angle_poloidal)) and divided by the prediction sol_width*f_x'''
    target = metrics['targets'][end['target']]
    separation = np.hypot(end['R_target'] - target['R_target'], end['Z_target'] - target['Z_target'])
    spacing = separation*np.sin(np.radians(target['angle_poloidal_target']))
    return spacing/(metrics['sol']['width']*target['poloidal_flux_expansion'])


def on_wall(R, Z, wall, tol=1.E-8):
    '''Check whether (R, Z) lies on a closed wall contour'''
    pt = np.r_[R, Z]
    starts = wall[:-1, :]
    vecs = np.diff(wall, axis=0)
    t = np.clip(((pt - starts)*vecs).sum(axis=1)/(vecs**2).sum(axis=1), 0.0, 1.0)
    return np.sqrt((((starts + t[:, None]*vecs) - pt)**2).sum(axis=1)).min() < tol


#------------------------------------------------------------------------------
# Analytic field: straight current filaments (no OFT runtime required)
#------------------------------------------------------------------------------
class FilamentField():
    '''Poloidal field of straight current filaments, psi = -sum_k I_k ln(d_k), with B_phi = F0/R'''
    def __init__(self, filaments, F0=5.0):
        self.filaments = np.asarray(filaments, dtype=np.float64)
        self.F0 = F0

    def psi(self, pts):
        pts = np.atleast_2d(pts)
        out = np.zeros((pts.shape[0],))
        for R, Z, I in self.filaments:
            out -= 0.5*I*np.log((pts[:, 0]-R)**2 + (pts[:, 1]-Z)**2)
        return out

    def grad_psi(self, pts):
        pts = np.atleast_2d(pts)
        out = np.zeros(pts.shape)
        for R, Z, I in self.filaments:
            d2 = (pts[:, 0]-R)**2 + (pts[:, 1]-Z)**2
            out[:, 0] -= I*(pts[:, 0]-R)/d2
            out[:, 1] -= I*(pts[:, 1]-Z)/d2
        return out

    def B(self, pts):
        pts = np.atleast_2d(pts)
        grad = self.grad_psi(pts)
        return np.column_stack((-grad[:, 1]/pts[:, 0], self.F0/pts[:, 0], grad[:, 0]/pts[:, 0]))

    def find_xpoint(self, guess):
        '''Locate an X-point with Newton's method'''
        x = np.array(guess, dtype=np.float64)
        h = 1.E-6
        for _ in range(50):
            hess = np.column_stack((
                (self.grad_psi(x + [h, 0.0])[0] - self.grad_psi(x - [h, 0.0])[0])/(2.0*h),
                (self.grad_psi(x + [0.0, h])[0] - self.grad_psi(x - [0.0, h])[0])/(2.0*h)
            ))
            x = x - np.linalg.solve(hess, self.grad_psi(x)[0])
        return x


ANALYTIC_WALL = np.array([[1.2, -1.05], [2.8, -1.05], [2.8, 1.05], [1.2, 1.05], [1.2, -1.05]])


def analytic_metrics(I_upper, I_lower, sol_width=3.E-3, wall=ANALYTIC_WALL):
    '''Plasma filament at (2, 0) with divertor filaments at (2, +/-1.2); returns field, metrics'''
    field = FilamentField([[2.0, 0.0, 1.0], [2.0, 1.2, I_upper], [2.0, -1.2, I_lower]])
    x_points = np.array([field.find_xpoint([2.0, 0.85]), field.find_xpoint([2.0, -0.85])])
    return field, _divertor_metrics(field, [2.0, 0.0], x_points, [wall], sol_width=sol_width)


def test_divertor_analytic_symmetric_dn():
    '''Up/down symmetric double null: dR_sep ~ 0, symmetric strike points, independent L_par and f_x checks'''
    field, metrics = analytic_metrics(0.35, 0.35, sol_width=1.E-3)
    assert abs(metrics['dR_sep']) < 1.E-6
    assert set(metrics['targets'].keys()) == {'upper_inner', 'upper_outer', 'lower_inner', 'lower_outer'}
    for target in metrics['targets'].values():
        assert target['status'] == 'wall'
        assert on_wall(target['R_target'], target['Z_target'], ANALYTIC_WALL)
        psi_x = metrics['xpoints'][target['xpoint']]['psi']
        assert abs(field.psi(np.r_[target['R_target'], target['Z_target']])[0] - psi_x) < 1.E-4
        assert target['poloidal_flux_expansion'] > 1.0
        assert 0.0 < target['angle_of_incidence_target'] < target['angle_poloidal_target'] <= 90.0
    for side in ('inner', 'outer'):
        upper = metrics['targets']['upper_' + side]
        lower = metrics['targets']['lower_' + side]
        assert upper['R_target'] == pytest.approx(lower['R_target'], abs=1.E-5)
        assert upper['Z_target'] == pytest.approx(-lower['Z_target'], abs=1.E-5)
        assert upper['poloidal_flux_expansion'] == pytest.approx(lower['poloidal_flux_expansion'], rel=1.E-4)
    # Outboard SOL connects only to the two outer targets
    assert sorted(end['target'] for end in metrics['sol']['ends']) == ['lower_outer', 'upper_outer']
    start = np.r_[metrics['sol']['R_OMP'], 0.0]
    for end in metrics['sol']['ends']:
        target = metrics['targets'][end['target']]
        # L_par against a direct integral along a gridded contour
        L_contour = contour_L_par(field.psi, field.B, metrics['sol']['psi'], start, np.r_[end['R_target'], end['Z_target']],
                                  (1.2, 2.8), (-1.05, 1.05), 1.E-3)
        assert target['connection_length_target'] == pytest.approx(L_contour, rel=1.E-3)
        assert 0.0 < target['connection_length_divertor'] < target['connection_length_target']
        # Strike-point separation on the wall against the flux expansion
        assert strike_separation_ratio(metrics, end) == pytest.approx(1.0, rel=2.E-2)


def test_divertor_analytic_unbalanced_dn():
    '''Upper-biased double null: dR_sep > 0 with the upper X-point primary; L_par increases toward separatrix'''
    _, metrics = analytic_metrics(0.40, 0.35)
    primary = metrics['xpoints'][metrics['primary_xpoint']]
    assert primary['Z'] > 0.0
    assert metrics['dR_sep'] > 0.0
    _, metrics_flip = analytic_metrics(0.35, 0.40)
    assert metrics_flip['dR_sep'] == pytest.approx(-metrics['dR_sep'], rel=1.E-4)
    L_par = []
    for width in (1.E-3, 3.E-3, 1.E-2):
        _, metrics = analytic_metrics(0.40, 0.35, sol_width=width)
        L_par.append([metrics['targets'][name]['connection_length_target'] for name in ('upper_inner', 'upper_outer')])
    L_par = np.array(L_par)
    assert np.all(np.isfinite(L_par))
    assert np.all(np.diff(L_par, axis=0) < 0.0)


def test_divertor_analytic_errors():
    '''Invalid SOL widths and a wall that excludes the innermost (primary) X-point raise ValueError'''
    for width in (0.0, -1.E-3, np.nan):
        with pytest.raises(ValueError, match='sol_width'):
            analytic_metrics(0.35, 0.35, sol_width=width)
    with pytest.raises(ValueError, match='outside the wall'):
        analytic_metrics(0.35, 0.35, sol_width=0.5)
    # Upper X-point is primary but the wall cuts below it: the plasma would be limited
    wall = ANALYTIC_WALL.copy()
    wall[2:4, 1] = 0.85
    with pytest.raises(ValueError, match='Primary X-point'):
        analytic_metrics(0.40, 0.35, wall=wall)


#------------------------------------------------------------------------------
# ITER lower single null equilibrium (same setup as test_TokaMaker.py)
#------------------------------------------------------------------------------
def setup_ITER(myOFT):
    '''Load the ITER test mesh and set up coils, targets and profiles as in test_TokaMaker.py'''
    if not os.path.exists('ITER_mesh.h5'):
        create_ITER_mesh()
    mygs = TokaMaker(myOFT)
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh('ITER_mesh.h5')
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.setup(order=2, F0=5.3*6.2)
    mygs.set_coil_vsc({'VS': 1.0})
    mygs.set_coil_bounds({key: [-50.E6, 50.E6] for key in mygs.coil_sets})
    regularization_terms = []
    for name in mygs.coil_sets:
        weight = 2.E-2 if name.startswith('CS1') else 1.E-2
        regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=weight))
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E2))
    mygs.set_coil_reg(reg_terms=regularization_terms)
    mygs.set_profiles(ffp_prof=create_power_flux_fun(40, 1.5, 2.0), pp_prof=create_power_flux_fun(40, 4.0, 1.0))
    return mygs


def solve_ITER_diverted(mygs):
    '''Lower single null ITER equilibrium from test_TokaMaker.py'''
    mygs.set_targets(Ip=15.6E6, pax=6.2E5)
    isoflux_pts = np.array([
        [8.20, 0.41], [8.06, 1.46], [7.51, 2.62], [6.14, 3.78], [4.51, 3.02],
        [4.26, 1.33], [4.28, 0.08], [4.49, -1.34], [7.28, -1.89], [8.00, -0.68]
    ])
    x_point = np.array([[5.125, -3.4], ])
    mygs.set_isoflux_constraints(np.vstack((isoflux_pts, x_point)))
    mygs.set_saddle_constraints(x_point)
    mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
    return mygs.solve()


def solve_ITER_limited(mygs):
    '''Shift the ITER shape 0.2 m inward: the plasma limits on the inner wall but the lower X-point stays inside the wall'''
    mygs.set_targets(Ip=15.6E6, pax=6.2E5)
    isoflux_pts = np.array([
        [8.20, 0.41], [8.06, 1.46], [7.51, 2.62], [6.14, 3.78], [4.51, 3.02],
        [4.26, 1.33], [4.28, 0.08], [4.49, -1.34], [7.28, -1.89], [8.00, -0.68]
    ]) - [0.2, 0.0]
    x_point = np.array([[4.925, -3.4], ])
    mygs.set_isoflux_constraints(np.vstack((isoflux_pts, x_point)))
    mygs.set_saddle_constraints(x_point)
    mygs.init_psi(6.1, 0.5, 2.0, 1.4, 0.0)
    mygs.solve()


def summarize(metrics, mygs, psi_eval):
    '''Picklable summary of the metrics for the primary targets'''
    summary = {key: metrics[key] for key in ('diverted', 'R_OMP', 'B_pol_OMP', 'dR_sep', 'psi_sep')}
    summary['primary_Z'] = metrics['xpoints'][metrics['primary_xpoint']]['Z']
    summary['psi_range'] = abs(mygs.psi_bounds[1] - mygs.psi_bounds[0])
    summary['sol'] = {key: metrics['sol'][key] for key in ('width', 'R_OMP', 'psi')}
    summary['ends'] = [{key: end[key] for key in ('target', 'R_target', 'Z_target', 'connection_length')}
                       for end in metrics['sol']['ends']]
    summary['targets'] = {}
    for name, target in metrics['targets'].items():
        if not target['primary']:
            continue
        entry = {key: target[key] for key in (
            'status', 'R_target', 'Z_target', 'R_OMP', 'B_target', 'poloidal_flux_expansion',
            'toroidal_flux_expansion', 'total_flux_expansion', 'angle_poloidal_target',
            'angle_of_incidence_target', 'connection_length_target', 'connection_length_divertor')}
        entry['on_wall'] = bool(on_wall(target['R_target'], target['Z_target'], mygs.lim_contour, tol=1.E-6))
        entry['dpsi'] = psi_eval.eval(np.r_[target['R_target'], target['Z_target']])[0] - metrics['psi_sep']
        summary['targets'][name] = entry
    summary['strike_ratio'] = {end['target']: strike_separation_ratio(metrics, end) for end in metrics['sol']['ends']}
    return summary


def run_ITER_divertor_case(mp_q):
    '''Diverted ITER metrics at several SOL widths, independent checks, equilibrium-object call and error paths'''
    results = {}
    try:
        myOFT = OFT_env(nthreads=-1)
        mygs = setup_ITER(myOFT)
        EQ_obj = solve_ITER_diverted(mygs)
        psi_eval = mygs.get_field_eval('PSI')
        B_eval = mygs.get_field_eval('B')
        results['widths'] = []
        for sol_width in (1.E-3, 3.E-3, 1.E-2):
            metrics = compute_divertor_metrics(mygs, sol_width=sol_width)
            results['widths'].append(summarize(metrics, mygs, psi_eval))
        # Independent L_par: integrate |B|/B_p along a gridded contour from the OMP to each target (3 mm SOL)
        sol = results['widths'][1]
        grid_ranges = {'lower_outer': ((4.9, 8.35), (-4.45, 0.8)), 'lower_inner': ((4.0, 8.35), (-4.0, 4.7))}
        results['L_par_contour'] = {}
        for end in sol['ends']:
            R_range, Z_range = grid_ranges[end['target']]
            results['L_par_contour'][end['target']] = contour_L_par(
                lambda pts: psi_eval.eval(pts)[:, 0], B_eval.eval, sol['sol']['psi'],
                np.r_[sol['sol']['R_OMP'], mygs.o_point[1]], np.r_[end['R_target'], end['Z_target']],
                R_range, Z_range, 5.E-3)
        # Same result through a TokaMaker_equilibrium object
        metrics_eq = compute_divertor_metrics(EQ_obj, sol_width=3.E-3)
        results['eq_obj'] = summarize(metrics_eq, mygs, psi_eval)
        # Error paths
        errors = {}
        for name, kwargs in (('zero', {'sol_width': 0.0}), ('negative', {'sol_width': -1.E-3}),
                             ('outside', {'sol_width': 0.5})):
            try:
                compute_divertor_metrics(mygs, **kwargs)
                errors[name] = None
            except ValueError as err:
                errors[name] = str(err)
        solve_ITER_limited(mygs)
        x_points, diverted = mygs.get_xpoints()
        results['limited_has_xpoint'] = (x_points is not None) and (not diverted)
        try:
            compute_divertor_metrics(mygs)
            errors['limited'] = None
        except ValueError as err:
            errors['limited'] = str(err)
        results['errors'] = errors
    except Exception as err:
        print(err)
        mp_q.put(None)
        return
    mp_q.put(results)
    oftpy_dump_cov()


# Values from the default 5 mm step and 3 mm SOL; they change by <5E-4 when the step is reduced to 0.625 mm
REGRESSION = {
    'dR_sep': -0.042130,
    'targets': {
        'lower_inner': {
            'R_target': 4.20842, 'Z_target': -3.82620, 'poloidal_flux_expansion': 5.86393,
            'total_flux_expansion': 0.54454, 'angle_poloidal_target': 39.69971, 'angle_of_incidence_target': 2.26071,
            'connection_length_target': 139.31518, 'connection_length_divertor': 30.6143
        },
        'lower_outer': {
            'R_target': 5.56500, 'Z_target': -4.23657, 'poloidal_flux_expansion': 5.02459,
            'total_flux_expansion': 0.71958, 'angle_poloidal_target': 26.43035, 'angle_of_incidence_target': 1.83705,
            'connection_length_target': 71.69449, 'connection_length_divertor': 26.3181
        }
    }
}


def test_divertor_ITER():
    '''Lower single null ITER equilibrium: strike points, flux expansion, angles and connection lengths'''
    results = mp_run(run_ITER_divertor_case, (), timeout=180)
    assert results is not None
    for summary in results['widths']:
        assert summary['diverted']
        assert summary['primary_Z'] < 0.0
        # Upper X-point separatrix lies outside the active (lower) separatrix
        assert summary['dR_sep'] < 0.0
        assert set(summary['targets'].keys()) == {'lower_inner', 'lower_outer'}
        assert sorted(end['target'] for end in summary['ends']) == ['lower_inner', 'lower_outer']
        for target in summary['targets'].values():
            assert target['status'] == 'wall'
            assert target['on_wall']
            assert abs(target['dpsi']) < 1.E-4*summary['psi_range']
            assert target['poloidal_flux_expansion'] > 1.0
            assert 0.0 < target['angle_of_incidence_target'] < target['angle_poloidal_target'] <= 90.0
            assert 0.0 < target['connection_length_divertor'] < target['connection_length_target']
        inner = summary['targets']['lower_inner']
        outer = summary['targets']['lower_outer']
        assert inner['R_target'] < outer['R_target']
        assert outer['connection_length_target'] < inner['connection_length_target']
    # Connection length increases as the SOL surface approaches the separatrix
    for name in ('lower_inner', 'lower_outer'):
        L_par = [summary['targets'][name]['connection_length_target'] for summary in results['widths']]
        assert L_par[0] > L_par[1] > L_par[2]
    # Strike-point separation on the wall matches sol_width*f_x (1 mm SOL)
    for ratio in results['widths'][0]['strike_ratio'].values():
        assert ratio == pytest.approx(1.0, rel=2.E-2)
    # Independent L_par along a gridded contour (3 mm SOL)
    summary = results['widths'][1]
    for name, L_contour in results['L_par_contour'].items():
        assert summary['targets'][name]['connection_length_target'] == pytest.approx(L_contour, rel=5.E-4)
    # TokaMaker_equilibrium input gives the same answer
    for name in ('lower_inner', 'lower_outer'):
        for key in ('R_target', 'Z_target', 'poloidal_flux_expansion', 'connection_length_target'):
            assert results['eq_obj']['targets'][name][key] == pytest.approx(summary['targets'][name][key], rel=1.E-10)
    # Regression (converged to ~1E-4 in step size)
    assert summary['dR_sep'] == pytest.approx(REGRESSION['dR_sep'], rel=1.E-2)
    for name, expected in REGRESSION['targets'].items():
        target = summary['targets'][name]
        assert target['R_target'] == pytest.approx(expected['R_target'], abs=5.E-3)
        assert target['Z_target'] == pytest.approx(expected['Z_target'], abs=5.E-3)
        for key in ('poloidal_flux_expansion', 'total_flux_expansion', 'angle_poloidal_target',
                    'angle_of_incidence_target', 'connection_length_target', 'connection_length_divertor'):
            assert target[key] == pytest.approx(expected[key], rel=1.E-2)
    # Error paths
    errors = results['errors']
    assert 'sol_width' in errors['zero']
    assert 'sol_width' in errors['negative']
    assert 'outside the wall' in errors['outside']
    assert results['limited_has_xpoint']
    assert 'limited' in errors['limited']

