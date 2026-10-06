"""Analysis suite for TokaMaker flux functions on normalized toroidal flux (Phi_N).

Not part of the automated tests. It regenerates the data and figures used to validate the
toroidal-flux map and the cut-cell q/F backend on the ITER (diverted) and LTX (limited) test cases:
  - cut-cell map against a map built from densely traced q surfaces
  - round trips: psi_N profiles rewritten on Phi_N nodes and re-solved
  - solve_bootstrap with all profiles on Phi_N
  - comparison between the cut-cell q/F map against the field-line tracer in accuracy, cost, and in-solve cost

Usage (run from this directory, ~15 min on 8 cores):
    python tokamaker_torflux_analysis.py [--out DIR] [--steps bench iter nodes ltx boot plots]
                                         [--orders 2 3] [--threads 1 4]
Each case runs in its own subprocess (one OFT environment per process). Timings are the minimum
over repeated runs; use a dedicated node for meaningful numbers.
"""
import os, sys, json, time, argparse, subprocess
import numpy as np
from scipy.integrate import cumulative_trapezoid
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(HERE, '..', '..', 'python')))


# ============================================================ setup
def _oft():
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain, save_gs_mesh, load_gs_mesh
    return OFT_env, TokaMaker, gs_Domain, save_gs_mesh, load_gs_mesh


def _iter_mesh():
    _, _, gs_Domain, save_gs_mesh, _ = _oft()
    with open(os.path.join(HERE, 'ITER_geom.json'), 'r') as fid:
        geom = json.load(fid)
    m = gs_Domain()
    m.define_region('air', 0.6, 'boundary'); m.define_region('plasma', 0.15, 'plasma')
    m.define_region('vacuum1', 0.3, 'vacuum'); m.define_region('vacuum2', 0.3, 'vacuum')
    m.define_region('vv1', 0.3, 'conductor', eta=6.9E-7); m.define_region('vv2', 0.3, 'conductor', eta=6.9E-7)
    for key in geom['coils']:
        if not key.startswith('VS'):
            m.define_region(key, 0.2, 'coil')
    m.define_region('VSU', 0.2, 'coil', coil_set='VS', nTurns=1.0)
    m.define_region('VSL', 0.2, 'coil', coil_set='VS', nTurns=-1.0)
    m.add_polygon(geom['limiter'], 'plasma', parent_name='vacuum1')
    m.add_annulus(geom['inner_vv'][0], 'vacuum1', geom['inner_vv'][1], 'vv1', parent_name='vacuum2')
    m.add_annulus(geom['outer_vv'][0], 'vacuum2', geom['outer_vv'][1], 'vv2', parent_name='air')
    for key, c in geom['coils'].items():
        m.add_rectangle(c['rc'], c['zc'], c['w'], c['h'], key, parent_name='vacuum1' if key.startswith('VS') else 'air')
    pts, lc, reg = m.build_mesh()
    save_gs_mesh(pts, lc, reg, m.get_coils(), m.get_conductors(), os.path.join(HERE, 'ITER_mesh.h5'))


def _ltx_mesh():
    _, _, gs_Domain, save_gs_mesh, _ = _oft()
    with open(os.path.join(HERE, 'LTX_geom.json'), 'r') as fid:
        g = json.load(fid)
    m = gs_Domain()
    m.define_region('air', 0.05, 'boundary'); m.define_region('plasma', 0.02, 'plasma')
    m.define_region('shellU', 0.015, 'conductor', eta=4.E-7, noncontinuous=True)
    m.define_region('shellL', 0.015, 'conductor', eta=4.E-7, noncontinuous=True)
    for i, seg in enumerate(g['vv']):
        m.define_region('vv{0}'.format(i), 0.015, 'conductor', eta=seg[1])
    for key, c in g['coils'].items():
        if key.startswith('OH'):
            m.define_region(key, 0.02, 'coil', nTurns=c['nturns'], coil_set='OH')
        else:
            m.define_region(key, 0.02, 'coil', nTurns=c['nturns'])
    m.add_polygon(g['limiter'], 'plasma', parent_name='air')
    m.add_polygon(g['shell'], 'shellU', parent_name='air')
    sl = np.array(g['shell'].copy()); sl[:, 1] *= -1.0
    m.add_polygon(sl, 'shellL', parent_name='air')
    for i, seg in enumerate(g['vv']):
        m.add_polygon(seg[0], 'vv{0}'.format(i), parent_name='air')
    for key, c in g['coils'].items():
        m.add_rectangle(c['rc'], c['zc'], c['w'], c['h'], key, parent_name='air')
    pts, lc, reg = m.build_mesh()
    save_gs_mesh(pts, lc, reg, m.get_coils(), m.get_conductors(), os.path.join(HERE, 'LTX_mesh.h5'))


def setup_iter(nthreads=4, order=2, Ip=15.6E6, pax=6.2E5):
    OFT_env, TokaMaker, _, _, load_gs_mesh = _oft()
    if not os.path.exists(os.path.join(HERE, 'ITER_mesh.h5')):
        _iter_mesh()
    mygs = TokaMaker(OFT_env(nthreads=nthreads))
    mesh = load_gs_mesh(os.path.join(HERE, 'ITER_mesh.h5'))
    mygs.setup_mesh(*mesh[:3])
    mygs.setup_regions(cond_dict=mesh[4], coil_dict=mesh[3])
    mygs.settings.maxits = 100
    mygs.setup(order=order, F0=5.3*6.2)
    mygs.set_coil_vsc({'VS': 1.0})
    mygs.set_coil_bounds({key: [-50.E6, 50.E6] for key in mygs.coil_sets})
    mygs.set_targets(Ip=Ip, pax=pax)
    iso = np.array([[8.20, 0.41], [8.06, 1.46], [7.51, 2.62], [6.14, 3.78], [4.51, 3.02],
                    [4.26, 1.33], [4.28, 0.08], [4.49, -1.34], [7.28, -1.89], [8.00, -0.68]])
    xpt = np.array([[5.125, -3.4]])
    mygs.set_isoflux(np.vstack((iso, xpt)))
    mygs.set_saddles(xpt)
    reg = [mygs.coil_reg_term({n: 1.0}, target=0.0, weight=2.E-2 if n.startswith('CS1') else 1.E-2)
           for n in mygs.coil_sets]
    reg.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E2))
    mygs.set_coil_reg(reg_terms=reg)
    return mygs


def init_iter(mygs):
    mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)


def setup_ltx(nthreads=4):
    from OpenFUSIONToolkit.TokaMaker.util import create_isoflux
    OFT_env, TokaMaker, _, _, load_gs_mesh = _oft()
    if not os.path.exists(os.path.join(HERE, 'LTX_mesh.h5')):
        _ltx_mesh()
    mygs = TokaMaker(OFT_env(nthreads=nthreads))
    mesh = load_gs_mesh(os.path.join(HERE, 'LTX_mesh.h5'))
    mygs.setup_mesh(*mesh[:3]); mygs.setup_regions(cond_dict=mesh[4], coil_dict=mesh[3])
    mygs.setup(order=2, F0=0.10752)
    mygs.set_coil_vsc({'INTERNALU': 1.0, 'INTERNALL': -1.0})
    mygs.set_targets(Ip=8.0E4, Ip_ratio=2.0)
    mygs.set_isoflux(create_isoflux(20, 0.40, 0.0, 0.22, 1.5, 0.1))
    reg = []
    for n in mygs.coil_sets:
        if n[:-1] == 'YELLOW':
            reg.append(mygs.coil_reg_term({n: 1.0}, target=0.0, weight=1.E4)); continue
        if n == 'OH':
            reg.append(mygs.coil_reg_term({n: 1.0}, target=0.0, weight=1.E-1)); continue
        elif n[-1] == 'L':
            continue
        reg.append(mygs.coil_reg_term({n: 1.0}, target=0.0, weight=1.E-1))
        reg.append(mygs.coil_reg_term({n: 1.0, n[:-1]+'L': -1.0}, target=0.0, weight=1.E2))
    reg.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E-4))
    mygs.set_coil_reg(reg_terms=reg)
    return mygs


def init_ltx(mygs):
    mygs.init_psi(0.42, 0.0, 0.15, 1.5, 0.6)


# ============================================================ testing-only Fortran hooks
_hooks = {}


def _load_hooks():
    if not _hooks:
        from ctypes import c_int, c_double, c_void_p, c_char_p
        from OpenFUSIONToolkit._interface import oftpy_lib, ctypes_subroutine, ctypes_numpy_array
        _hooks['backend'] = ctypes_subroutine(oftpy_lib.tokamaker_torflux_backend, [c_int, c_double])
        _hooks['qgeom'] = ctypes_subroutine(oftpy_lib.tokamaker_torflux_qgeom,
            [c_void_p, c_int, ctypes_numpy_array(np.float64, 1), ctypes_numpy_array(np.float64, 1), c_int, c_char_p])
    return _hooks


def set_backend(backend=0, tol=-1.0):
    '''Map q/F backend inside the solve (0 cut-cell, 1 tracer) and tracer tolerance (<0 default)'''
    _load_hooks()['backend'](backend, tol)


def qgeom(mygs, psi_n, backend):
    '''q/F at psi_N (0 at axis) from cut-cell (0) or tracer (1)'''
    x = np.ascontiguousarray(1.0-np.asarray(psi_n, dtype=np.float64))
    g = np.zeros_like(x)
    err = mygs._oft_env.get_c_errorbuff()
    _load_hooks()['qgeom'](mygs._tMaker_equil._equil_ptr, x.shape[0], x, g, backend, err)
    if err.value != b'':
        raise Exception(err.value)
    return np.abs(g)


# ============================================================ helpers
def traced_map(mygs):
    '''Phi_N(psi_N), dPhi_N/dpsi_N and Q from q traced on ~1700 surfaces (trapezoid, analytic end corrections)'''
    psi = np.unique(np.concatenate((np.linspace(1.E-4, 0.99, 1500), 1.0-np.logspace(-6, -2, 200))))
    psi, q, _, _, _, _ = mygs.get_q(psi=psi)
    phi = cumulative_trapezoid(q, psi, initial=0.0)
    s0 = psi[0]; q0 = q[0]-(q[1]-q[0])*s0/(psi[1]-psi[0])
    phi = phi+0.5*(q0+q[0])*s0
    s = 1.0-psi[-2:]
    if mygs.diverted:
        b = (q[-1]-q[-2])/np.log(s[1]/s[0]); a = q[-1]-b*np.log(s[1])
        tail = a*s[1]+b*(s[1]*np.log(s[1])-s[1])
    else:
        tail = q[-1]*s[1]
    Q = phi[-1]+tail
    return np.concatenate(([0.0], psi, [1.0])), np.concatenate(([0.0], phi/Q, [1.0])), \
        np.concatenate(([q0], q, [np.nan]))/Q, Q


def conv_map(mygs):
    '''psi_N -> Phi_N conversion (traced q) for building Phi_N profiles; J = inf at a diverted LCFS'''
    ps, phi, jac, _ = traced_map(mygs)
    jac[-1] = np.inf if mygs.diverted else jac[-2]
    return ps, phi, jac


def to_tor(prof, pmap, coord):
    '''psi_N linterp profile -> Phi_N nodes ('phi_n': y/J, 'phi_n_relabel': y unchanged)'''
    x = np.interp(prof['x'], pmap[0], pmap[1]); x[0] = 0.0; x[-1] = 1.0
    y = prof['y']/np.interp(prof['x'], pmap[0], pmap[2]) if coord == 'phi_n' else prof['y']
    return {'type': prof.get('type', 'linterp'), 'x': x, 'y': y, 'coord': coord}


def snapshot(mygs, npsi=200):
    psi_n, q, _, _, _, _ = mygs.get_q(npsi=npsi, psi_pad=1.E-3)
    pr = mygs.get_profiles(npsi=npsi, psi_pad=1.E-3)
    st = mygs.get_stats()
    return {'psi': mygs.get_psi(False), 'psi_n': psi_n, 'q': q, 'prof_psi': pr[0], 'F': pr[1], 'Fp': pr[2],
            'P': pr[3], 'Pp': pr[4], 'Ip': st['Ip']}


def diff(a, b):
    rq = np.abs(a['q']-b['q'])/np.abs(b['q']); x = np.asarray(b['psi_n'])
    return {'psi_err': float(np.linalg.norm(a['psi']-b['psi'])/np.linalg.norm(b['psi'])),
            'q_err': float(np.max(rq)), 'q_err_mid': float(np.max(rq[(x >= 0.01) & (x <= 0.95)])),
            'Ip_err': float(abs(a['Ip']-b['Ip'])/b['Ip'])}


def field(mygs):
    '''Normalized flux on mesh nodes (0 at axis) and the LCFS for contour plots'''
    try:
        lcfs = mygs.trace_surf(0.9999)
    except Exception:
        lcfs = None
    return {'r': mygs.r, 'lc': mygs.lc, 'psi_n': mygs.get_psi(True), 'lim': mygs.lim_contour, 'lcfs': lcfs}


def solve(mygs, init, ffp, pp, nrep=5):
    '''Solve nrep times from the same initial state; minimum wall time'''
    mygs.set_profiles(ffp_prof=ffp, pp_prof=pp)
    ts = []
    for _ in range(nrep):
        init(mygs)
        t0 = time.perf_counter(); _, its = mygs.solve(return_its=True); ts.append(time.perf_counter()-t0)
    return its, float(np.min(ts))


def timed(fun, nrep):
    ts = []
    for _ in range(nrep):
        t0 = time.perf_counter(); out = fun(); ts.append(time.perf_counter()-t0)
    return out, float(np.min(ts))


class _Enc(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.generic):
            return o.item()
        return super().default(o)


def save_json(out_dir, name, obj):
    with open(os.path.join(out_dir, 'data', name), 'w') as fid:
        json.dump(obj, fid, cls=_Enc)


def load_json(out_dir, name):
    with open(os.path.join(out_dir, 'data', name)) as fid:
        return json.load(fid)


# ============================================================ cases
def case_bench(out_dir, order, nthreads):
    '''Cut-cell vs tracer q/F: accuracy, cost vs tolerance, cost vs surfaces, full solves with each map backend'''
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
    tols = [1.E-5, 1.E-6, 1.E-7, 1.E-8, 1.E-9, 1.E-10, 1.E-11, 1.E-12]  # 1e-8 = default, 1e-10 = high res
    mygs = setup_iter(nthreads=nthreads, order=order)
    ffp = create_power_flux_fun(40, 1.5, 2.0); pp = create_power_flux_fun(40, 4.0, 1.0)
    its_psi, t_psi = solve(mygs, init_iter, ffp, pp)
    out = {'order': order, 'nthreads': nthreads}
    # Accuracy: g = q/F on each backend; reference = tracer at tol 1e-12 (F cancels in g/g_ref)
    psi_n = np.unique(np.concatenate((np.linspace(0.005, 0.98, 196), 1.0-np.logspace(-4, np.log10(0.02), 40))))
    _, q_trace, _, _, _, _ = mygs.get_q(psi=psi_n)
    g_cut, t_cut = timed(lambda: qgeom(mygs, psi_n, 0), 5)
    g_tr, t_tr = {}, {}
    for tol in tols:
        set_backend(0, tol)
        g_tr[tol], t_tr[tol] = timed(lambda: qgeom(mygs, psi_n, 1), 3 if tol > 1.E-11 else 1)
    set_backend(0, -1.0)
    g_ref = g_tr[1.E-12]
    out['acc'] = {'psi_n': psi_n, 'q_ref': q_trace/qgeom(mygs, psi_n, 1)*g_ref, 'rel_cut': g_cut/g_ref-1.0, 't_cut': t_cut,
                  'rel_tr': {str(t): g_tr[t]/g_ref-1.0 for t in tols}, 't_tr': {str(t): t_tr[t] for t in tols}}
    # Cost vs number of surfaces
    tim = {'nr': [10, 25, 50, 100, 200, 400], 'cut': [], 'tr_def': [], 'tr_high': []}
    for nr in tim['nr']:
        x = np.linspace(0.01, 0.99, nr)
        tim['cut'].append(timed(lambda: qgeom(mygs, x, 0), 7)[1])
        for key, tol in (('tr_def', 1.E-8), ('tr_high', 1.E-10)):
            set_backend(0, tol); tim[key].append(timed(lambda: qgeom(mygs, x, 1), 3)[1])
    set_backend(0, -1.0)
    out['timing'] = tim
    # Full 'phi_n' solves with each map backend
    ref = snapshot(mygs); pmap = conv_map(mygs)
    out['solve'] = {'psi_n': {'t': t_psi, 'its': its_psi}}
    snaps = {}
    for label, backend, tol in (('cutcell', 0, -1.0), ('tracer_default', 1, 1.E-8), ('tracer_high', 1, 1.E-10)):
        set_backend(backend, tol)
        its, t = solve(mygs, init_iter, to_tor(ffp, pmap, 'phi_n'), to_tor(pp, pmap, 'phi_n'))
        snaps[label] = snapshot(mygs)
        out['solve'][label] = dict(t=t, its=its, **diff(snaps[label], ref))
        print(label, out['solve'][label], flush=True)
    set_backend(0, -1.0)
    out['solve']['cut_vs_trhigh'] = diff(snaps['cutcell'], snaps['tracer_high'])
    out['solve']['trdef_vs_trhigh'] = diff(snaps['tracer_default'], snaps['tracer_high'])
    save_json(out_dir, 'bench_o{0}_t{1}.json'.format(order, nthreads), out)


def map_check(mygs):
    '''Cut-cell map (get_torflux_map) against the traced-q map'''
    x = np.linspace(0.0, 1.0, 401)
    phi_m, jac_m = mygs.get_torflux_map(x)
    psi_back, _ = mygs.get_torflux_map(phi_m, inverse=True)
    ps, phi_r, jac_r, Q = traced_map(mygs)
    xm = np.linspace(0.02, 0.98, 49)
    _, jm = mygs.get_torflux_map(xm)
    _, qm, _, _, _, _ = mygs.get_q(psi=xm)
    Q_map = float(np.median(qm/jm))
    b = mygs.psi_bounds
    return {'x': x, 'phi_map': phi_m, 'jac_map': jac_m, 'psi_back': psi_back,
            'psi_ref': ps, 'phi_ref': phi_r, 'jac_ref': jac_r, 'Q_ref': Q, 'Q_map': Q_map,
            'phi_edge_map': 2.0*np.pi*Q_map*abs(b[1]-b[0]), 'phi_edge_ref': 2.0*np.pi*Q*abs(b[1]-b[0]),
            'tflux': mygs.get_stats()['tflux']}


def roundtrip(mygs, init, ffp, pp, with_fields=False):
    out = {}
    its, t = solve(mygs, init, ffp, pp)
    ref = snapshot(mygs); pmap = conv_map(mygs)
    out['psi_n'] = {'its': its, 't': t, 'snap': ref}
    if with_fields:
        out['psi_n']['field'] = field(mygs)
    for coord in ('phi_n', 'phi_n_relabel'):
        its, t = solve(mygs, init, to_tor(ffp, pmap, coord), to_tor(pp, pmap, coord))
        s = snapshot(mygs)
        out[coord] = dict(its=its, t=t, snap=s, **diff(s, ref))
        if coord == 'phi_n':
            out[coord]['map'] = map_check(mygs)
        if with_fields:
            out[coord]['field'] = field(mygs)
        print(coord, {k: v for k, v in out[coord].items() if k not in ('snap', 'map', 'field')}, flush=True)
    return out


def case_iter(out_dir):
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
    mygs = setup_iter()
    save_json(out_dir, 'map_iter.json', roundtrip(mygs, init_iter, create_power_flux_fun(40, 1.5, 2.0),
                                                  create_power_flux_fun(40, 4.0, 1.0), with_fields=True))


def case_nodes(out_dir):
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
    mygs = setup_iter()
    out = {}
    for n in (20, 40, 80, 160, 320, 640):
        r = roundtrip(mygs, init_iter, create_power_flux_fun(n, 1.5, 2.0), create_power_flux_fun(n, 4.0, 1.0))
        out[n] = {c: {k: r[c][k] for k in ('its', 't', 'psi_err', 'q_err', 'q_err_mid', 'Ip_err')}
                  for c in ('phi_n', 'phi_n_relabel')}
    save_json(out_dir, 'map_nodes.json', out)


def case_ltx(out_dir):
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
    mygs = setup_ltx()
    save_json(out_dir, 'map_ltx.json', roundtrip(mygs, init_ltx, create_power_flux_fun(50, 1.5, 2.0),
                                                 create_power_flux_fun(50, 4.0, 1.0), with_fields=True))


def case_boot(out_dir):
    '''solve_bootstrap on psi_N; Phi_N nodes from that equilibrium's traced q; solve_bootstrap on Phi_N'''
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
    from OpenFUSIONToolkit.TokaMaker.bootstrap import Hmode_profiles
    mygs = setup_iter(Ip=13.0E6)
    n = 257; x = np.linspace(0.0, 1.0, n)
    ne = Hmode_profiles(edge=0.35, ped=0.6, core=1.1, rgrid=n, expin=1.6, expout=1.6, widthp=0.35, xphalf=0.965)*1.E20
    Te = Hmode_profiles(edge=1500., ped=5000., core=21000., rgrid=n, expin=1.3, expout=1.7, widthp=0.1, xphalf=0.965)
    jind = create_power_flux_fun(n, 2.25, 2.5)['y']
    solve(mygs, init_iter, create_power_flux_fun(40, 1.5, 2.0), create_power_flux_fun(40, 4.0, 1.0), nrep=1)

    def run(x_in, coord):
        t0 = time.perf_counter()
        r = mygs.solve_bootstrap(
            ffp_prof={'type': 'jphi-split-bootstrap', 'x': x_in, 'y': jind},
            te_prof={'type': 'linterp', 'x': x_in, 'y': Te/1.E3}, ne_prof={'type': 'linterp', 'x': x_in, 'y': ne},
            ti_prof={'type': 'linterp', 'x': x_in, 'y': Te/1.E3}, ni_prof={'type': 'linterp', 'x': x_in, 'y': ne},
            Zeff=1.5, Ip_target=13.0E6, coord=coord)
        r = {k: r[k] for k in ('psi_n', 'j_bs_final', 'total_j_phi')}
        r['t'] = time.perf_counter()-t0
        return r
    out = {'x': x, 'Te': Te}
    out['psi'] = run(x, 'psi_n'); ref = snapshot(mygs)
    pmap = conv_map(mygs)
    x_phi = np.interp(x, pmap[0], pmap[1]); x_phi[0] = 0.0; x_phi[-1] = 1.0
    out['x_phi'] = x_phi
    out['phi'] = run(x_phi, 'phi_n'); s = snapshot(mygs)
    out['diff'] = diff(s, ref); out['snap_psi'] = ref; out['snap_phi'] = s
    save_json(out_dir, 'map_boot.json', out)


# ============================================================ figures
def make_plots(out_dir, orders, threads):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    from matplotlib.path import Path
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.3, 'figure.dpi': 130})
    A = lambda v: np.asarray(v, dtype=float)
    CC = {'psi_n': 'k', 'phi_n': 'C0', 'phi_n_relabel': 'C3'}
    LBL = {'psi_n': r'$\hat\psi$ profiles (reference)', 'phi_n': r"$\hat\Phi$, 'phi_n'", 'phi_n_relabel': r"$\hat\Phi$, 'phi_n_relabel'"}
    # Colour = method and thread count, fill = FE order (hollow lowest, solid highest)
    col_cut = {4: 'k', 1: 'C3'}; col_tr = {4: 'C0', 1: 'C1'}; col_psi = {4: '0.45', 1: '0.75'}
    fill = lambda o: o != min(orders)

    def mk(color, o):
        return {'color': color, 'mec': color, 'mfc': color if fill(o) else 'none'}

    def save(fig, name):
        fig.savefig(os.path.join(out_dir, 'figures', name), bbox_inches='tight'); plt.close(fig)
        print('wrote', name, flush=True)

    def contours(ax, fld, levels, key='psi_n', **kw):
        '''Contours of a nodal field, restricted to triangles inside the LCFS'''
        r = A(fld['r']); lc = np.asarray(fld['lc'], dtype=int)
        tri = Triangulation(r[:, 0], r[:, 1], lc)
        if fld.get('lcfs') is not None:
            inside = Path(A(fld['lcfs'])).contains_points(r[:, :2], radius=1.E-6)
            tri.set_mask(~np.any(inside[lc], axis=1))
        return ax.tricontour(tri, A(fld[key]), levels=levels, **kw)

    def frame(ax, fld, lcfs_color='k'):
        l = A(fld['lim']); ax.plot(l[:, 0], l[:, 1], color='0.5', lw=1.0)
        if fld.get('lcfs') is not None:
            c = A(fld['lcfs']); ax.plot(c[:, 0], c[:, 1], color=lcfs_color, lw=1.6)
        dr = 0.05*(l[:, 0].max()-l[:, 0].min())
        ax.set_xlim(l[:, 0].min()-dr, l[:, 0].max()+dr); ax.set_ylim(l[:, 1].min()-dr, l[:, 1].max()+dr)
        ax.set_aspect('equal'); ax.grid(False); ax.set_xlabel('R [m]')

    # ---- cut-cell map vs traced-q map
    def fig_map(d, name, title):
        m = d['phi_n']['map']
        x = A(m['x']); pm = A(m['phi_map']); jm = A(m['jac_map'])
        ps = A(m['psi_ref']); pr = A(m['phi_ref']); jr = A(m['jac_ref']); ok = np.isfinite(jr)
        fig, ax = plt.subplots(2, 2, figsize=(10, 7))
        a = ax[0, 0]
        a.plot(ps, pr, 'k-', lw=3, alpha=0.3, label='traced-q map (q traced on ~1700 surfaces)')
        a.plot(x, pm, 'C0--', label='cut-cell map (get_torflux_map)')
        a.set_xlabel(r'$\hat\psi$'); a.set_ylabel(r'$\hat\Phi$'); a.legend(fontsize=8); a.set_title(r'$\hat\Phi(\hat\psi)$')
        a = ax[0, 1]
        a.plot(x, pm-np.interp(x, ps, pr), 'C0')
        a.set_xlabel(r'$\hat\psi$'); a.set_ylabel(r'$\hat\Phi_{cut}-\hat\Phi_{traced}$'); a.set_title('Map difference')
        a.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        a = ax[1, 0]
        a.plot(ps[ok], jr[ok], 'k-', lw=3, alpha=0.3, label='traced-q map')
        a.plot(x, jm, 'C0--', label='cut-cell map')
        a.set_xlabel(r'$\hat\psi$'); a.set_ylabel(r'$d\hat\Phi/d\hat\psi$'); a.legend(fontsize=8)
        a.set_title(r'Jacobian $d\hat\Phi/d\hat\psi = q/Q$, $Q=\int_0^1 q\,d\hat\psi$')
        a = ax[1, 1]
        msk = (x > 0.0) & (x < 0.995)
        a.semilogy(x[msk], np.abs(jm[msk]/np.interp(x[msk], ps[ok], jr[ok])-1.0), 'C0', label='Jacobian, cut-cell vs traced (relative)')
        ib = np.abs(A(m['psi_back'])-x); ib[ib == 0] = 1.E-17
        a.semilogy(x, ib, 'C3', label=r'cut-cell inverse: $|\hat\psi-\hat\psi(\hat\Phi(\hat\psi))|$')
        a.set_xlabel(r'$\hat\psi$'); a.set_ylabel('difference'); a.legend(fontsize=8); a.set_title('Jacobian and inverse-map')
        fig.suptitle(title+'\n'+r'$\Phi_{{edge}}$: cut-cell {0:.6g} Wb, traced q {1:.6g} Wb, $\int B_t\,dA$ {2:.6g} Wb'.format(
            m['phi_edge_map'], m['phi_edge_ref'], m['tflux']), fontsize=10)
        fig.tight_layout(); save(fig, name)
        return {'map_err': float(np.max(np.abs(pm-np.interp(x, ps, pr)))), 'inv_err': float(np.max(ib)),
                'Q_rel': m['Q_map']/m['Q_ref']-1.0}

    def fig_roundtrip(d, name, title):
        fig = plt.figure(figsize=(13, 7.5))
        gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.2, 1.2])
        a = fig.add_subplot(gs[:, 0])
        lev = np.linspace(0.1, 0.9, 9)
        for c, ls in (('psi_n', '-'), ('phi_n', '--'), ('phi_n_relabel', ':')):
            contours(a, d[c]['field'], lev, colors=CC[c], linestyles=ls, linewidths=1.0)
            a.plot([], [], color=CC[c], ls=ls, label=LBL[c])
        frame(a, d['psi_n']['field'])
        a.plot([], [], 'k', lw=1.6, label='LCFS'); a.set_ylabel('Z [m]'); a.legend(fontsize=7, loc='lower left')
        a.set_title(r'$\hat\psi$ = 0.1..0.9 contours')
        s0 = d['psi_n']['snap']
        a = fig.add_subplot(gs[0, 1])
        for c, ls in (('psi_n', '-'), ('phi_n', '--'), ('phi_n_relabel', ':')):
            s = d[c]['snap']; a.plot(A(s['psi_n']), A(s['q']), color=CC[c], ls=ls, label=LBL[c])
        a.set_xlabel(r'$\hat\psi$'); a.set_ylabel('q'); a.legend(fontsize=8); a.set_title('Safety factor')
        a = fig.add_subplot(gs[1, 1])
        for c in ('phi_n', 'phi_n_relabel'):
            s = d[c]['snap']; a.semilogy(A(s['psi_n']), np.abs(A(s['q'])/A(s0['q'])-1.0), color=CC[c], label=LBL[c])
        a.set_xlabel(r'$\hat\psi$'); a.set_ylabel(r'$|\Delta q|/q$'); a.legend(fontsize=8); a.set_title('q difference from reference')
        for k, (lab, f) in enumerate(((r"$FF'$", lambda s: A(s['F'])*A(s['Fp'])), (r"$P'$", lambda s: A(s['Pp'])))):
            a = fig.add_subplot(gs[k, 2])
            for c, ls in (('psi_n', '-'), ('phi_n', '--'), ('phi_n_relabel', ':')):
                s = d[c]['snap']; a.plot(A(s['prof_psi']), f(s), color=CC[c], ls=ls)
            a.set_xlabel(r'$\hat\psi$'); a.set_ylabel(lab); a.set_title(lab+r'$(\hat\psi)$ seen by the solver')
        rows = ['{0}: its {1}, {2:.2f} s'.format(c, d[c]['its'], d[c]['t']) for c in ('psi_n', 'phi_n', 'phi_n_relabel')]
        fig.suptitle(title+'\n'+'   |   '.join(rows), fontsize=10)
        fig.tight_layout(); save(fig, name)

    def fig_coords(d, name):
        fld = d['phi_n']['field']; m = d['phi_n']['map']
        x = A(m['x']); pm = A(m['phi_map']); pf = A(fld['psi_n'])
        fld = dict(fld, phi_n=np.where(pf <= 1.0, np.interp(np.clip(pf, 0, 1), x, pm), pf))
        fig, ax = plt.subplots(1, 3, figsize=(12, 5.2), gridspec_kw={'width_ratios': [1, 1, 1.4]})
        lev = np.linspace(0.1, 0.9, 9)
        for a, key, t in ((ax[0], 'psi_n', r'Equal $\hat\psi$ = 0.1..0.9'), (ax[1], 'phi_n', r'Equal $\hat\Phi$ = 0.1..0.9')):
            contours(a, fld, lev, key=key, colors='C0', linewidths=1.0)
            frame(a, fld); a.set_title(t)
        ax[0].set_ylabel('Z [m]'); ax[0].plot([], [], 'k', lw=1.6, label='LCFS'); ax[0].legend(fontsize=8, loc='lower left')
        s = d['phi_n']['snap']
        a = ax[2]
        a.plot(x, pm, 'C0', label=r'$\hat\Phi(\hat\psi)$ (cut-cell map)')
        a.set_xlabel(r'$\hat\psi$'); a.set_ylabel(r'$\hat\Phi$')
        b = a.twinx(); b.plot(A(s['psi_n']), A(s['q']), 'C1'); b.set_ylabel('q', color='C1'); b.grid(False)
        a.legend(loc='upper left', fontsize=8); a.set_title(r'$d\hat\Phi/d\hat\psi = q/Q$')
        fig.tight_layout(); save(fig, name)

    def fig_nodes(d, name):
        ns = sorted(int(k) for k in d)
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
        for c, m_ in (('phi_n', 'o'), ('phi_n_relabel', 's')):
            ax[0].loglog(ns, [d[str(n)][c]['psi_err'] for n in ns], m_+'-', color=CC[c], label=LBL[c])
            ax[1].loglog(ns, [d[str(n)][c]['q_err_mid'] for n in ns], m_+'-', color=CC[c], label=LBL[c]+r', $0.01\leq\hat\psi\leq0.95$')
            ax[1].loglog(ns, [d[str(n)][c]['q_err'] for n in ns], m_+':', color=CC[c], mfc='none', alpha=0.6,
                         label=LBL[c]+r', $0.001\leq\hat\psi\leq0.999$')
        e0 = d[str(ns[0])]['phi_n']['psi_err']
        for a in ax:
            a.loglog(ns, [e0*(ns[0]/n)**2 for n in ns], ':', color='0.5', label=r'$N^{-2}$')
            a.set_xlabel('profile nodes N'); a.legend(fontsize=8)
        ax[0].set_ylabel(r'$\|\psi-\psi_{ref}\|/\|\psi_{ref}\|$'); ax[0].set_title(r'$\psi$ difference')
        ax[1].set_ylabel(r'max $|\Delta q|/q$'); ax[1].set_title('q difference')
        fig.suptitle("ITER: FF' and p' are 'linterp' profiles, piecewise linear in their own coordinate. The reference is linear\n"
                     r"between $\hat\psi$ nodes; the re-solves are linear between the same surfaces' $\hat\Phi$ nodes, which the solver"
                     "\nmaps to $\\hat\\psi$ each iteration. A line in $\\hat\\Phi$ is a curve in $\\hat\\psi$, so the solves differ by O(h$^2$).",
                     fontsize=9)
        fig.tight_layout(); save(fig, name)

    def fig_boot(d, name):
        fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
        a = ax[0]
        a.plot(A(d['x']), A(d['Te'])/1.E3, 'C1', label=r'on $\hat\psi$ nodes')
        a.plot(A(d['x_phi']), A(d['Te'])/1.E3, 'C1--', label=r'on $\hat\Phi$ nodes')
        a.set_xlabel(r'$\hat\psi$ or $\hat\Phi$'); a.set_ylabel(r'$T_e$ [keV]'); a.legend(fontsize=8, loc='upper right')
        a.set_title('Same $T_e$ profile, two coordinate representations', fontsize=10)
        a.text(0.03, 0.03, "1. solve_bootstrap, coord='psi_n', nodes x\n2. trace q of that equilibrium $\\rightarrow\\hat\\Phi(\\hat\\psi)$\n"
               "3. $\\hat\\Phi$ nodes = $\\hat\\Phi$(x), same values\n4. solve_bootstrap, coord='phi_n'",
               transform=a.transAxes, fontsize=7.5, va='bottom', bbox=dict(boxstyle='round', fc='white', ec='0.6'))
        a = ax[1]
        for k, c, ls in (('psi', 'k', '-'), ('phi', 'C0', '--')):
            r = d[k]
            a.plot(A(r['psi_n']), A(r['j_bs_final'])/1.E6, color=c, ls=ls, label=r"$j_{BS}$, coord='%s_n'" % k)
            a.plot(A(r['psi_n']), A(r['total_j_phi'])/1.E6, color=c, ls=ls, alpha=0.5, lw=0.8)
        a.set_xlabel(r'$\hat\psi$ (from each run)'); a.set_ylabel(r'MA/m$^2$'); a.legend(fontsize=8)
        a.set_title(r'Bootstrap (bold) and total $j_\phi$ (thin)')
        a = ax[2]
        s0 = d['snap_psi']; s1 = d['snap_phi']
        a.plot(A(s0['psi_n']), A(s0['q']), 'k', label='psi_n'); a.plot(A(s1['psi_n']), A(s1['q']), 'C0--', label='phi_n')
        a.set_xlabel(r'$\hat\psi$'); a.set_ylabel('q'); a.legend(fontsize=8)
        a.set_title('q: rel ψ diff {0:.1e}, max Δq/q {1:.1e}'.format(d['diff']['psi_err'], d['diff']['q_err']))
        fig.suptitle("ITER solve_bootstrap: all profiles on $\\hat\\Phi$ vs on $\\hat\\psi$; wall time {0:.1f} s vs {1:.1f} s".format(
            d['phi']['t'], d['psi']['t']), fontsize=10)
        fig.tight_layout(); save(fig, name)

    # ---- q backends
    bs = {}
    for o in orders:
        for t in threads:
            try:
                bs[(o, t)] = load_json(out_dir, 'bench_o{0}_t{1}.json'.format(o, t))
            except FileNotFoundError:
                pass

    def fig_qacc(name):
        cols = [o for o in orders if (o, max(threads)) in bs]
        fig, ax = plt.subplots(2, len(cols), figsize=(5.5*len(cols), 7.5), squeeze=False)
        for j, o in enumerate(cols):
            a_ = bs[(o, max(threads))]['acc']; x = A(a_['psi_n'])
            series = (('cut-cell', 'k', A(a_['rel_cut'])), ('tracer, default tol', 'C0', A(a_['rel_tr']['1e-08'])),
                      (r'tracer, high res (tol$\times$0.01)', 'C1', A(a_['rel_tr']['1e-10'])))
            for i, xx, xl in ((0, x, r'$\hat\psi$'), (1, 1.0-x, r'$1-\hat\psi$ (distance from LCFS)')):
                a = ax[i, j]
                for l, c, e in series:
                    a.plot(xx, np.abs(e)+1.E-16, '.-', ms=2, lw=0.7, color=c, label=l)
                a.set_yscale('log'); a.set_ylim(1.E-13, 1.E-2); a.set_xlabel(xl); a.set_ylabel(r'$|q/q_{ref}-1|$')
                if i == 1:
                    a.set_xscale('log')
                a.set_title('order {0}: '.format(o)+('relative q error' if i == 0 else 'near the separatrix'))
            ax[0, j].legend(fontsize=8, loc='upper left')
        fig.suptitle(r'q accuracy on ITER (diverted), 235 surfaces. $q_{ref}$: field-line tracer at tol $10^{-12}$ ($10^{-14}$ near LCFS)'
                     '\non the same FE field; errors compare q/F, so F cancels', fontsize=10)
        fig.tight_layout(); save(fig, name)

    def fig_pareto(name):
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
        for (o, t), b in sorted(bs.items()):
            a_ = b['acc']; tols = sorted(a_['rel_tr'], key=float, reverse=True)[:-1]
            lab = 'order {0}, {1} thr'.format(o, t)
            for i, stat in enumerate((np.max, np.median)):
                ax[i].loglog([a_['t_tr'][k] for k in tols], [stat(np.abs(A(a_['rel_tr'][k]))) for k in tols], 'o-',
                             **mk(col_tr[t], o), label='tracer tol sweep, '+lab)
                ax[i].loglog([a_['t_cut']], [stat(np.abs(A(a_['rel_cut'])))], '*', ms=15, **mk(col_cut[t], o), label='cut-cell, '+lab)
                if (o, t) == (min(orders), max(threads)):
                    for k, txt in (('1e-05', 'tol 1e-5'), ('1e-08', 'default'), ('1e-10', 'high res')):
                        ax[i].annotate(txt, (a_['t_tr'][k], stat(np.abs(A(a_['rel_tr'][k])))), textcoords='offset points',
                                       xytext=(5, 5), fontsize=8)
        for i, t in enumerate(('max', 'median')):
            ax[i].set_xlabel('wall time per call, 235 surfaces [s]'); ax[i].set_ylabel('{0} relative q error'.format(t))
            ax[i].set_title('Cost vs accuracy ({0} over surfaces)'.format(t))
        h, l = ax[0].get_legend_handles_labels()
        fig.legend(h, l, loc='lower center', ncol=4, fontsize=8, bbox_to_anchor=(0.5, -0.12))
        fig.suptitle('Colour: method and thread count; hollow: order {0}, solid: order {1}'.format(min(orders), max(orders)), fontsize=9)
        fig.tight_layout(); save(fig, name)

    def fig_timing(name):
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
        for i, (key, title) in enumerate((('tr_def', 'default tracer'), ('tr_high', 'high-res tracer'))):
            for (o, t), b in sorted(bs.items()):
                tm = b['timing']; lab = 'order {0}, {1} thr'.format(o, t)
                ax[i].loglog(tm['nr'], tm['cut'], '*-', ms=9, **mk(col_cut[t], o), label='cut-cell, '+lab)
                ax[i].loglog(tm['nr'], tm[key], 'o-', **mk(col_tr[t], o), label='tracer, '+lab)
            ax[i].set_xlabel('number of surfaces'); ax[i].set_title('Cut-cell vs '+title)
        ax[0].set_ylabel('wall time per call [s]'); ax[0].legend(fontsize=7, ncol=2)
        fig.suptitle(r'Cost of one q/F evaluation (ITER, surfaces uniform in $\hat\psi\in[0.01,0.99]$). '
                     'Hollow: order {0}, solid: order {1}'.format(min(orders), max(orders)), fontsize=9)
        fig.tight_layout(); save(fig, name)

    def fig_solve(name):
        keys = ('psi_n', 'cutcell', 'tracer_default', 'tracer_high')
        cols = {'psi_n': col_psi, 'cutcell': col_cut, 'tracer_default': col_tr, 'tracer_high': col_tr}
        lbl = ['$\\hat\\psi$ profiles\n(no map)', "'phi_n'\ncut-cell map", "'phi_n'\ntracer map", "'phi_n'\ntracer map\n(high res)"]
        cfg = sorted(bs); w = 0.8/len(cfg)
        fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
        for j, (o, t) in enumerate(cfg):
            s = bs[(o, t)]['solve']
            for i, k in enumerate(keys):
                c = cols[k][t]; xb = i+(j-(len(cfg)-1)/2)*w
                ax[0].bar(xb, s[k]['t'], w, color=c if fill(o) else 'white', edgecolor=c, lw=1.5)
                ax[0].annotate(str(s[k]['its']), (xb, s[k]['t']), ha='center', va='bottom', fontsize=7)
        ax[0].set_xticks(np.arange(4)); ax[0].set_xticklabels(lbl, fontsize=8); ax[0].set_yscale('log')
        ax[0].set_ylabel('wall time per solve [s]')
        ax[0].set_title('Full ITER solve. Bars per group: {0}\n(numbers: nonlinear iterations; hollow: order {1}, solid: order {2})'.format(
            ', '.join('o{0}/{1}thr'.format(o, t) for o, t in cfg), min(orders), max(orders)), fontsize=9)
        rows = []
        for (o, t) in cfg:
            s = bs[(o, t)]['solve']
            rows.append(['order {0}, {1} thr'.format(o, t), '{0:.1e}'.format(s['cut_vs_trhigh']['psi_err']),
                         '{0:.1e}'.format(s['trdef_vs_trhigh']['psi_err']), '{0:.2f}x'.format(s['cutcell']['t']/s['psi_n']['t']),
                         '{0:.1f}x'.format(s['tracer_default']['t']/s['cutcell']['t'])])
        ax[1].axis('off')
        tb = ax[1].table(cellText=rows, colLabels=['case', 'rel |Δψ|\ncut-cell vs\ntracer high', 'rel |Δψ|\ntracer default\nvs tracer high',
                                                    "cut-cell 'phi_n'\n/ psi_n time", 'tracer default\n/ cut-cell time'],
                         loc='center', cellLoc='center')
        tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(1.0, 2.4)
        ax[1].set_title('Converged equilibria from each map backend')
        fig.tight_layout(); save(fig, name)

    summary = {}
    if os.path.exists(os.path.join(out_dir, 'data', 'map_iter.json')):
        it = load_json(out_dir, 'map_iter.json')
        summary['map_iter'] = fig_map(it, 'map_check_iter.png', 'ITER (diverted): cut-cell map vs map from traced q')
        fig_roundtrip(it, 'roundtrip_iter.png', r"ITER (diverted): $\hat\psi$ profiles re-expressed on $\hat\Phi$ (40 nodes) and re-solved")
        fig_coords(it, 'coordinates_iter.png')
    if os.path.exists(os.path.join(out_dir, 'data', 'map_nodes.json')):
        fig_nodes(load_json(out_dir, 'map_nodes.json'), 'node_convergence_iter.png')
    if os.path.exists(os.path.join(out_dir, 'data', 'map_ltx.json')):
        lt = load_json(out_dir, 'map_ltx.json')
        summary['map_ltx'] = fig_map(lt, 'map_check_ltx.png', 'LTX (limited): cut-cell map vs map from traced q')
        fig_roundtrip(lt, 'roundtrip_ltx.png', r"LTX (limited): $\hat\psi$ profiles re-expressed on $\hat\Phi$ (50 nodes) and re-solved")
    if os.path.exists(os.path.join(out_dir, 'data', 'map_boot.json')):
        fig_boot(load_json(out_dir, 'map_boot.json'), 'bootstrap_iter.png')
    if bs:
        fig_qacc('q_accuracy.png'); fig_pareto('q_cost_vs_accuracy.png')
        fig_timing('q_timing_vs_surfaces.png'); fig_solve('solve_timing.png')
    print(json.dumps(summary, indent=1))


# ============================================================ driver
if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', default=os.path.join(os.getcwd(), 'torflux_analysis'), help='output directory')
    p.add_argument('--steps', nargs='+', default=['bench', 'iter', 'nodes', 'ltx', 'boot', 'plots'])
    p.add_argument('--orders', nargs='+', type=int, default=[2, 3], help='FE orders for the q backend comparison')
    p.add_argument('--threads', nargs='+', type=int, default=[1, 4], help='thread counts (1 and/or 4)')
    p.add_argument('--worker', nargs='+', help=argparse.SUPPRESS)
    args = p.parse_args()
    out_dir = os.path.abspath(args.out)
    if args.worker:
        case = args.worker[0]
        if case == 'bench':
            case_bench(out_dir, int(args.worker[1]), int(args.worker[2]))
        else:
            {'iter': case_iter, 'nodes': case_nodes, 'ltx': case_ltx, 'boot': case_boot}[case](out_dir)
        sys.exit(0)
    for sub in ('data', 'figures'):
        os.makedirs(os.path.join(out_dir, sub), exist_ok=True)
    jobs = []
    for step in args.steps:
        if step == 'bench':
            jobs += [['bench', str(o), str(t)] for o in args.orders for t in args.threads]
        elif step != 'plots':
            jobs.append([step])
    for job in jobs:
        print('==== ' + ' '.join(job), flush=True)
        t0 = time.perf_counter()
        res = subprocess.run([sys.executable, os.path.abspath(__file__), '--out', out_dir, '--worker'] + job)
        print('==== {0}: exit {1}, {2:.0f} s'.format(' '.join(job), res.returncode, time.perf_counter()-t0), flush=True)
    if 'plots' in args.steps:
        make_plots(out_dir, args.orders, args.threads)
