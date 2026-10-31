#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Tests for TokaMaker divertor geometry metrics (`TokaMaker.divertor`)'''
import os
import sys
import time
import json
import multiprocessing
import pytest
import numpy as np
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..', 'python')))
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain, save_gs_mesh, load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
from OpenFUSIONToolkit.TokaMaker.divertor import compute_divertor_metrics, _divertor_metrics


def mp_run(target, args, timeout=60):
    '''! Run `target` in a separate process (only one OFT_env may exist per process)

    @param target Function to run, called with `args + (queue,)`
    @param args Arguments for `target`
    @param timeout Timeout [s]
    @result Object placed on the queue by `target` or `None` on failure
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


#------------------------------------------------------------------------------
# Analytic field: straight current filaments (no OFT runtime required)
#------------------------------------------------------------------------------
class _FilamentField():
    r'''! Poloidal field of straight current filaments, \f$ \psi = -\sum_k I_k \ln d_k \f$, with \f$ B_\phi = F_0/R \f$'''
    def __init__(self, filaments, F0=5.0):
        r'''! Create field

        @param filaments Filament (R, Z, I) [n,3]
        @param F0 Toroidal field function \f$ R B_\phi \f$
        '''
        self.filaments = np.asarray(filaments, dtype=np.float64)
        self.F0 = F0

    def psi(self, pts):
        '''! Evaluate flux

        @param pts Locations [n,2]
        @result Flux [n]
        '''
        pts = np.atleast_2d(pts)
        out = np.zeros((pts.shape[0],))
        for R, Z, I in self.filaments:
            out -= 0.5*I*np.log((pts[:, 0]-R)**2 + (pts[:, 1]-Z)**2)
        return out

    def grad_psi(self, pts):
        '''! Evaluate flux gradient

        @param pts Locations [n,2]
        @result Gradient [n,2]
        '''
        pts = np.atleast_2d(pts)
        out = np.zeros(pts.shape)
        for R, Z, I in self.filaments:
            d2 = (pts[:, 0]-R)**2 + (pts[:, 1]-Z)**2
            out[:, 0] -= I*(pts[:, 0]-R)/d2
            out[:, 1] -= I*(pts[:, 1]-Z)/d2
        return out

    def B(self, pts):
        '''! Evaluate field

        @param pts Locations [n,2]
        @result (B_R, B_phi, B_Z) [n,3]
        '''
        pts = np.atleast_2d(pts)
        grad = self.grad_psi(pts)
        return np.column_stack((-grad[:, 1]/pts[:, 0], self.F0/pts[:, 0], grad[:, 0]/pts[:, 0]))

    def find_xpoint(self, guess):
        '''! Locate an X-point with Newton's method

        @param guess Initial guess [2]
        @result X-point [2]
        '''
        x = np.array(guess, dtype=np.float64)
        h = 1.E-6
        for _ in range(50):
            hess = np.column_stack((
                (self.grad_psi(x + [h, 0.0])[0] - self.grad_psi(x - [h, 0.0])[0])/(2.0*h),
                (self.grad_psi(x + [0.0, h])[0] - self.grad_psi(x - [0.0, h])[0])/(2.0*h)
            ))
            x = x - np.linalg.solve(hess, self.grad_psi(x)[0])
        return x


_ANALYTIC_WALL = np.array([[1.2, -1.05], [2.8, -1.05], [2.8, 1.05], [1.2, 1.05], [1.2, -1.05]])


def _analytic_metrics(I_upper, I_lower, sol_width=3.E-3):
    '''! Compute metrics for a plasma filament with upper/lower divertor filaments

    @param I_upper Upper divertor filament current
    @param I_lower Lower divertor filament current
    @param sol_width SOL width at the OMP [m]
    @result Field, metrics dictionary
    '''
    field = _FilamentField([[2.0, 0.0, 1.0], [2.0, 1.2, I_upper], [2.0, -1.2, I_lower]])
    x_points = np.array([field.find_xpoint([2.0, 0.85]), field.find_xpoint([2.0, -0.85])])
    return field, _divertor_metrics(field, [2.0, 0.0], x_points, [_ANALYTIC_WALL], sol_width=sol_width)


def _on_wall(R, Z, wall, tol=1.E-8):
    '''! Check whether a point lies on a closed wall contour

    @param R Radial location [m]
    @param Z Vertical location [m]
    @param wall Closed contour [n,2]
    @param tol Distance tolerance [m]
    @result True if on the wall
    '''
    pt = np.r_[R, Z]
    starts = wall[:-1, :]
    vecs = np.diff(wall, axis=0)
    t = np.clip(((pt - starts)*vecs).sum(axis=1)/(vecs**2).sum(axis=1), 0.0, 1.0)
    return np.sqrt((((starts + t[:, None]*vecs) - pt)**2).sum(axis=1)).min() < tol


def test_divertor_analytic_symmetric_dn():
    '''! Up/down symmetric double null: dR_sep ~ 0, four symmetric strike points on wall and separatrix'''
    field, metrics = _analytic_metrics(0.35, 0.35)
    assert abs(metrics['dR_sep']) < 1.E-6
    assert set(metrics['targets'].keys()) == {'upper_inner', 'upper_outer', 'lower_inner', 'lower_outer'}
    for name, target in metrics['targets'].items():
        assert target['status'] == 'wall'
        assert _on_wall(target['R_target'], target['Z_target'], _ANALYTIC_WALL)
        psi_x = metrics['xpoints'][target['xpoint']]['psi']
        assert abs(field.psi(np.r_[target['R_target'], target['Z_target']])[0] - psi_x) < 1.E-4
        assert target['poloidal_flux_expansion'] > 1.0
        assert 0.0 < target['angle_of_incidence_target'] < target['angle_poloidal_target'] <= 90.0
        assert target['toroidal_flux_expansion'] == pytest.approx(target['R_target']/target['R_OMP'])
        assert target['total_flux_expansion'] == pytest.approx(
            target['poloidal_flux_expansion']*target['toroidal_flux_expansion'])
    for side in ('inner', 'outer'):
        upper = metrics['targets']['upper_' + side]
        lower = metrics['targets']['lower_' + side]
        assert upper['R_target'] == pytest.approx(lower['R_target'], abs=1.E-5)
        assert upper['Z_target'] == pytest.approx(-lower['Z_target'], abs=1.E-5)
        assert upper['poloidal_flux_expansion'] == pytest.approx(lower['poloidal_flux_expansion'], rel=1.E-4)
    # Outboard SOL connects only to the two outer targets, with equal lengths
    ends = sorted(end['target'] for end in metrics['sol']['ends'])
    assert ends == ['lower_outer', 'upper_outer']
    L_upper = metrics['targets']['upper_outer']['connection_length_target']
    L_lower = metrics['targets']['lower_outer']['connection_length_target']
    assert L_upper == pytest.approx(L_lower, rel=1.E-3)
    assert 0.0 < metrics['targets']['upper_outer']['connection_length_divertor'] < L_upper


def test_divertor_analytic_unbalanced_dn():
    '''! Upper-biased double null: dR_sep > 0, primary X-point is the upper one; L_par increases toward separatrix'''
    _, metrics = _analytic_metrics(0.40, 0.35)
    primary = metrics['xpoints'][metrics['primary_xpoint']]
    assert primary['Z'] > 0.0
    assert metrics['dR_sep'] > 0.0
    _, metrics_flip = _analytic_metrics(0.35, 0.40)
    assert metrics_flip['dR_sep'] == pytest.approx(-metrics['dR_sep'], rel=1.E-4)
    L_par = []
    for width in (1.E-3, 3.E-3, 1.E-2):
        _, metrics = _analytic_metrics(0.40, 0.35, sol_width=width)
        L_par.append([metrics['targets'][name]['connection_length_target'] for name in ('upper_inner', 'upper_outer')])
    L_par = np.array(L_par)
    assert np.all(np.isfinite(L_par))
    assert np.all(np.diff(L_par, axis=0) < 0.0)


#------------------------------------------------------------------------------
# ITER lower single null equilibrium (same setup as test_TokaMaker.py)
#------------------------------------------------------------------------------
def _create_ITER_mesh():
    '''! Create ITER mesh used by TokaMaker tests'''
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


def run_ITER_divertor_case(sol_widths, mp_q):
    '''! Solve ITER lower single null equilibrium and compute divertor metrics

    @param sol_widths SOL widths to evaluate [m]
    @param mp_q Queue for results
    '''
    try:
        if not os.path.exists('ITER_mesh.h5'):
            _create_ITER_mesh()
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
        x_point = np.array([[5.125, -3.4], ])
        mygs.set_isoflux_constraints(np.vstack((isoflux_pts, x_point)))
        mygs.set_saddle_constraints(x_point)
        regularization_terms = []
        for name in mygs.coil_sets:
            weight = 2.E-2 if name.startswith('CS1') else 1.E-2
            regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=weight))
        regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E2))
        mygs.set_coil_reg(reg_terms=regularization_terms)
        mygs.set_profiles(ffp_prof=create_power_flux_fun(40, 1.5, 2.0), pp_prof=create_power_flux_fun(40, 4.0, 1.0))
        mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
        mygs.solve()
        psi_eval = mygs.get_field_eval('PSI')
        wall = mygs.lim_contour
        results = []
        for sol_width in sol_widths:
            metrics = compute_divertor_metrics(mygs, sol_width=sol_width)
            summary = {
                'diverted': metrics['diverted'],
                'R_OMP': metrics['R_OMP'],
                'B_pol_OMP': metrics['B_pol_OMP'],
                'dR_sep': metrics['dR_sep'],
                'primary_Z': metrics['xpoints'][metrics['primary_xpoint']]['Z'],
                'psi_range': abs(mygs.psi_bounds[1] - mygs.psi_bounds[0]),
                'targets': {}
            }
            for name, target in metrics['targets'].items():
                if not target['primary']:
                    continue
                entry = {key: target[key] for key in (
                    'status', 'R_target', 'Z_target', 'R_OMP', 'poloidal_flux_expansion', 'toroidal_flux_expansion',
                    'total_flux_expansion', 'angle_poloidal_target', 'angle_of_incidence_target',
                    'connection_length_target', 'connection_length_divertor')}
                entry['on_wall'] = bool(_on_wall(target['R_target'], target['Z_target'], wall, tol=1.E-6))
                entry['dpsi'] = psi_eval.eval(np.r_[target['R_target'], target['Z_target']])[0] - metrics['psi_sep']
                summary['targets'][name] = entry
            results.append(summary)
    except Exception as err:
        print(err)
        mp_q.put(None)
        return
    mp_q.put(results)


def test_divertor_ITER():
    '''! Lower single null ITER equilibrium: strike points, flux expansion, angles and connection lengths'''
    sol_widths = (1.E-3, 3.E-3, 1.E-2)
    results = mp_run(run_ITER_divertor_case, (sol_widths,))
    assert results is not None
    for summary in results:
        assert summary['diverted']
        assert summary['primary_Z'] < 0.0
        # Upper X-point separatrix lies outside the active (lower) separatrix
        assert summary['dR_sep'] < 0.0
        assert set(summary['targets'].keys()) == {'lower_inner', 'lower_outer'}
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
        L_par = [summary['targets'][name]['connection_length_target'] for summary in results]
        assert L_par[0] > L_par[1] > L_par[2]
    # Loose regression on separatrix geometry
    targets = results[1]['targets']
    assert targets['lower_inner']['R_target'] == pytest.approx(4.208, abs=2.E-2)
    assert targets['lower_inner']['Z_target'] == pytest.approx(-3.826, abs=2.E-2)
    assert targets['lower_outer']['R_target'] == pytest.approx(5.565, abs=2.E-2)
    assert targets['lower_outer']['Z_target'] == pytest.approx(-4.237, abs=2.E-2)
    assert targets['lower_inner']['poloidal_flux_expansion'] == pytest.approx(5.86, rel=0.1)
    assert targets['lower_outer']['poloidal_flux_expansion'] == pytest.approx(5.03, rel=0.1)
