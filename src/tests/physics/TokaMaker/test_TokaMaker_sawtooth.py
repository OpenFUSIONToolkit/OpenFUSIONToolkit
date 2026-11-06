#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: sawtooth reset (saw_q_s) of the q profile, in 1-D and in the ITER internal bootstrap solve
import pytest
import numpy as np
from tokamaker_testing import mp_run, validate_dict, setup_ITER_bootstrap_internal
from OpenFUSIONToolkit._interface import oftpy_dump_cov


def run_ITER_bootstrap_saw_internal(mesh_resolution, fe_order, mp_q):
    '''Sawtooth reset (saw_q_s) in the internal bootstrap solve'''
    setup = setup_ITER_bootstrap_internal(mesh_resolution, fe_order)
    if setup is None:
        mp_q.put(None)
        return
    mygs, psi_sample, n_sample, Ip_target = (setup[k] for k in ('mygs', 'psi_sample', 'n_sample', 'Ip_target'))
    common = dict(
        ffp_prof={'type': 'jphi-split-bootstrap', 'x': psi_sample, 'y': setup['inductive_jphi']},
        te_prof={'type': 'linterp', 'x': psi_sample, 'y': setup['Te'] / 1e3},
        ne_prof={'type': 'linterp', 'x': psi_sample, 'y': setup['ne']},
        ti_prof={'type': 'linterp', 'x': psi_sample, 'y': setup['Ti'] / 1e3},
        ni_prof={'type': 'linterp', 'x': psi_sample, 'y': setup['ni']},
        Zeff=setup['Zeff_val'], Ip_target=Ip_target, scale_jbs=1.0,
    )

    def close(a, b, rtol):
        return np.allclose(a, b, rtol=rtol, atol=rtol*np.max(np.abs(b)))

    def sum_ok(p):
        return close(p['jphi_total'], p['jphi_ind'] + p['jphi_bs'] + p['jphi_pd'] + p['jphi_fixed'] + p['jphi_saw'], 1e-10)

    try:
        base = mygs.solve_bootstrap(**common)
        _, q_base, _, _, _, _ = mygs.get_q(npsi=40)
        # An input jphi_saw_prof with the reset off acts as a fixed current
        bump = 1.0e5 * np.exp(-((psi_sample - 0.2) / 0.1)**2)
        p_fix = mygs.solve_bootstrap(jphi_fixed_prof={'x': psi_sample, 'y': bump}, **common)
        mygs.boot_reset()  # jphi_fixed_prof would otherwise be kept
        p_saw = mygs.solve_bootstrap(jphi_saw_prof={'x': psi_sample, 'y': bump}, **common)
        # Both converge j_BS to djbs_tol (1e-4) and may stop an iteration apart
        if not close(p_saw['jphi_total'], p_fix['jphi_total'], 2e-4):
            raise AssertionError("jphi_saw_prof (saw off) total differs from the same jphi_fixed_prof")
        if not close(p_saw['jphi_saw'], np.interp(p_saw['psi_n'], psi_sample, bump), 1e-6) or np.any(p_saw['jphi_fixed'] != 0.0):
            raise AssertionError("j_saw / jphi_fixed with jphi_saw_prof (saw off) wrong")
        if not (sum_ok(p_fix) and sum_ok(p_saw)):
            raise AssertionError("jphi_total != j_ind + j_bs + jphi_fixed + jphi_saw")
        # jphi_saw survives an equilibrium save/load: re-solving the reloaded state keeps j_saw
        # (zeroed before the reload, since replace_eq loads into a copy of the live equilibrium)
        mygs._tMaker_equil.save_TokaMaker('ITER_saw_roundtrip.h5')
        mygs.set_kinetic_profiles(jphi_saw_prof={'x': np.array([0.0, 1.0]), 'y': np.zeros(2)})
        mygs.replace_eq(source_file='ITER_saw_roundtrip.h5')
        mygs.solve()
        p_rel = mygs.get_boot_profs()
        if not (close(p_rel['jphi_saw'], p_saw['jphi_saw'], 1e-6) and close(p_rel['jphi_total'], p_saw['jphi_total'], 1e-3)):
            raise AssertionError("jphi_saw lost in the equilibrium save/load")
        # q_s below min q: no reset, same result as saw off
        mygs.boot_reset()  # jphi_saw_prof would otherwise be kept
        p_low = mygs.solve_bootstrap(saw_q_s=0.5*np.min(q_base), **common)
        if np.any(p_low['jphi_saw'] != 0.0) or p_low['saw_n_dips'] != 0:
            raise AssertionError("saw_q_s below min q changed j_saw")
        for key in ('jphi_total', 'jphi_ind', 'jphi_bs'):
            if not close(p_low[key], base[key], 2e-4):
                raise AssertionError(f"saw_q_s below min q: '{key}' differs from saw off")
        # q_s above the unconstrained q0 with the local rule (the reversed-shear base puts the axis inside the
        # reset region): q0 -> q_s, monotone inside rho_out, Ip on target
        q_s, dq = float(q_base[0]) + 0.2, 0.03
        p_on = mygs.solve_bootstrap(saw_q_s=q_s, saw_dq=dq, diagnose_bs=False, **common)
        psi_q, q_on, _, _, _, _ = mygs.get_q(npsi=40)
        print(f"saw: q0 base {q_base[0]:.4f} -> {q_on[0]:.4f} (q_s {q_s:.4f}), rho_m {p_on['saw_rho_m']:.3f}, "
              f"rho_out {p_on['saw_rho_out']:.3f}, n_dips {p_on['saw_n_dips']}, max|j_saw| {np.max(np.abs(p_on['jphi_saw'])):.3e}")
        if abs(q_on[0] - q_s) > 0.03:
            raise AssertionError(f"saw: q0 {q_on[0]:.4f} not at q_s {q_s:.4f}")
        inner = q_on < q_s + dq
        inner[np.argmax(~inner):] = False
        if np.any(np.diff(q_on[inner]) < -2e-3):
            raise AssertionError("saw: q not monotone inside the mixing radius")
        ip = mygs.get_stats()['Ip']
        if abs(ip/Ip_target - 1.0) > 1e-3:
            raise AssertionError(f"saw: Ip {ip:.6e} off target {Ip_target:.6e}")
        if not sum_ok(p_on) or p_on['saw_n_dips'] != 1 or not (0.0 < p_on['saw_rho_m'] < 1.0) \
                or abs(p_on['saw_rho_out'] - p_on['saw_rho_m']) > 1e-12:
            raise AssertionError("saw: profile sum / rho_m / rho_out / n_dips wrong")
        # j_saw vanishes outside rho_out (rho_tor_norm from int q dpsi)
        def rho_tor(psi_q, q, psi):
            phi = np.concatenate(([0.0], np.cumsum(0.5*(q[1:] + q[:-1])*np.diff(psi_q))))
            return np.sqrt(np.interp(psi, psi_q, phi/phi[-1]))
        rho_of_psi = rho_tor(psi_q, q_on, p_on['psi_n'])
        outside = rho_of_psi > p_on['saw_rho_out'] + 0.05
        if np.max(np.abs(p_on['jphi_saw'][outside])) > 1e-3*np.max(np.abs(p_on['jphi_saw'])):
            raise AssertionError("saw: j_saw nonzero outside rho_out")
        # A second solve from the converged reset state must reset again (an early iterate with no dip must
        # not freeze the reset off)
        mygs.solve_bootstrap(saw_q_s=q_s, saw_dq=dq, **common)
        _, q_on2, _, _, _, _ = mygs.get_q(npsi=40)
        if abs(q_on2[0] - q_s) > 0.03:
            raise AssertionError(f"saw: repeated solve q0 {q_on2[0]:.4f} not at q_s {q_s:.4f}")
        # Node-dense axis (psi_N = rho^2, as on a FUSE core_profiles grid): inner cumulative-area bins hold no
        # quadrature points; the reset must stay finite and still reach q_s
        x2 = np.linspace(0.0, 1.0, n_sample)**2
        dense = {k: (dict(v, x=x2, y=np.interp(x2, v['x'], v['y'])) if isinstance(v, dict) else v)
                 for k, v in common.items()}
        p_d = mygs.solve_bootstrap(saw_q_s=q_s, saw_dq=dq, **dense)
        _, q_d, _, _, _, _ = mygs.get_q(npsi=40)
        print(f"saw dense axis: q0 {q_d[0]:.4f} (q_s {q_s:.4f}), rho_m {p_d['saw_rho_m']:.3f}")
        if not np.all(np.isfinite(p_d['jphi_saw'])) or abs(q_d[0] - q_s) > 0.03:
            raise AssertionError("saw: dense-axis grid reset not finite or q0 not at q_s")
        # Reversed shear with q_s between the off-axis min q and q0: the local rule raises only the dip (q0 kept,
        # j_saw zero inside the region's inner end) while fuse resets from the axis (q0 -> q_s)
        q_s = 0.5*(float(q_base[0]) + float(np.min(q_base)))
        out = {}
        for rule in (2, 1):
            p_r = mygs.solve_bootstrap(saw_q_s=q_s, saw_dq=dq, saw_rule=rule, **common)
            psi_r, q_r, _, _, _, _ = mygs.get_q(npsi=100)
            out[rule] = (p_r, psi_r, q_r)
            print(f"saw reversed shear, rule {rule}: q0 {q_base[0]:.4f} -> {q_r[0]:.4f}, min q {np.min(q_base):.4f} -> "
                  f"{np.min(q_r):.4f} (q_s {q_s:.4f}), rho_m {p_r['saw_rho_m']:.3f}, rho_out {p_r['saw_rho_out']:.3f}, "
                  f"n_dips {p_r['saw_n_dips']}, Ip/target-1 {mygs.get_stats()['Ip']/Ip_target - 1.0:.1e}")
            if not sum_ok(p_r) or p_r['saw_n_dips'] != 1 or abs(np.min(q_r) - q_s) > 0.01:
                raise AssertionError(f"saw reversed shear, rule {rule}: sum / n_dips / min q wrong")
        p_l, psi_l, q_l = out[2]
        # local: q0 is not lowered (it may rise slightly as alpha and j_BS respond to the raised dip)
        if q_l[0] < max(q_base[0] - 0.01, q_s + dq) or abs(out[1][2][0] - q_s) > 0.01:
            raise AssertionError("saw reversed shear: local must not lower q0, fuse must lower it to q_s")
        rho_l = rho_tor(psi_l, q_l, p_l['psi_n'])
        rho_in = rho_tor(psi_l, q_l, psi_l[np.argmax(q_l < q_s + dq)])
        away = (rho_l < rho_in - 0.05) | (rho_l > p_l['saw_rho_out'] + 0.05)
        if not (0.1 < rho_in < p_l['saw_rho_out']) or np.max(np.abs(p_l['jphi_saw'][away])) > 1e-3*np.max(np.abs(p_l['jphi_saw'])):
            raise AssertionError("saw reversed shear: local j_saw not confined to [rho_in, rho_out]")
        # A loaded bootstrap state is re-solvable as saved: a bare solve after loading reproduces the call with
        # every keyword input and option passed explicitly
        eC = 1.602176634e-19
        p_kin = eC*(setup['ne']*setup['Te'] + setup['ni']*setup['Ti'])
        kw = dict(jphi_fixed_prof={'x': psi_sample, 'y': 0.5*bump}, jphi_saw_prof={'x': psi_sample, 'y': 0.5*bump},
                  p_fixed_prof={'x': psi_sample, 'y': p_kin[0]*(0.1*np.exp(-((psi_sample - 0.3)/0.15)**2) + 1.0e-3)},
                  saw_q_s=q_s, saw_dq=dq, saw_rule=2, taper_edge_jphi=True)
        mygs.solve_bootstrap(**kw, **common)
        mygs._tMaker_equil.save_TokaMaker('ITER_boot_state.h5')
        mygs.boot_reset()
        mygs.replace_eq(source_file='ITER_boot_state.h5')
        p_kept = mygs.solve_bootstrap(**common)
        mygs.boot_reset()
        mygs.replace_eq(source_file='ITER_boot_state.h5')
        p_expl = mygs.solve_bootstrap(**kw, **common)
        for key in ('jphi_total', 'jphi_ind', 'jphi_bs', 'jphi_fixed', 'jphi_saw', 'jphi_pd'):
            if not close(p_kept[key], p_expl[key], 1e-10):
                raise AssertionError(f"loaded state: '{key}' differs from the explicit call")
        mygs.boot_reset()
        if mygs._tMaker_equil._boot_inputs is not None or mygs._tMaker_equil._boot_ops['saw_q_s'] != 0.0 \
                or mygs._tMaker_equil._boot_ops['taper_edge_jphi']:
            raise AssertionError("boot_reset() did not restore the defaults")
        p_def = mygs.solve_bootstrap(**common)
        if np.any(p_def['jphi_fixed'] != 0.0) or np.any(p_def['jphi_saw'] != 0.0):
            raise AssertionError("boot_reset(): jphi_fixed / j_saw not back to zero")
    except Exception as e:
        print(e)
        mp_q.put(None)
        return
    mp_q.put([{'saw_checks': 1.0}])
    oftpy_dump_cov()


@pytest.mark.slow
@pytest.mark.parametrize("order", (2,))
def test_ITER_bootstrap_saw_internal(order):
    '''Sawtooth q-profile reset inside the internal bootstrap solve on ITER: a sawtooth current profile with the reset
    off acts as a fixed current and survives save/load, no reset when q_s is below min q, q0 raised to q_s with the
    compensating current confined to the mixing region for both reset rules and a node-dense axis grid, and a saved
    state re-solves as saved.'''
    results = mp_run(run_ITER_bootstrap_saw_internal, (1.0, order), timeout=600)
    assert validate_dict(results, {'saw_checks': 1.0})


def run_saw_reset_1d(mp_q):
    import traceback
    try:
        _check_saw_reset_1d()
        mp_q.put('ok')
    except Exception:
        mp_q.put(traceback.format_exc())


def test_saw_reset_1d():
    '''The 1-D sawtooth reset on an analytic cylinder (I = rho^2/q, A = rho^2): mixing radius, monotone and
    current-neutral q reset for on-axis and off-axis dips, merged and split reset regions, both reset rules, and
    rejection of non-finite or inconsistent inputs.'''
    result = mp_run(run_saw_reset_1d, (), timeout=60)
    assert result == 'ok', result


def _check_saw_reset_1d():
    from OpenFUSIONToolkit.TokaMaker.bootstrap import _saw_reset_1d
    rho = np.linspace(0.0, 1.0, 201)
    A, c1 = rho**2, np.ones_like(rho)
    q_s, dq = 1.025, 0.03
    g = lambda c, w: np.exp(-((rho - c)/w)**2)
    reset = lambda q, rule: _saw_reset_1d(rho, q, A/q, A, c1, q_s, dq, rule=rule)
    # Single on-axis dip, rule 1 (fuse)
    q = 0.9 + 2.5*rho**2
    q_new, dj, sc = reset(q, 1)
    assert sc['n_dips'] == 1 and sc['w'] == 1.0
    assert np.isclose(sc['rho_s'], np.sqrt((q_s - 0.9)/2.5), atol=2e-3)
    assert np.isclose(sc['rho_m'], np.sqrt((q_s + dq - 0.9)/2.5), atol=2e-3)
    assert sc['rho_out'] == sc['rho_m'] and np.all(sc['wr'] == np.where(rho < sc['rho_m'], 1.0, 0.0))
    assert np.isclose(q_new[0], q_s) and np.all(np.diff(q_new) >= -1e-12)
    jm = np.searchsorted(rho, sc['rho_m'])
    assert np.all(q_new[jm:] == q[jm:]) and np.all(dj[jm+1:] == 0.0)
    # Ip-neutral up to the O(h) gradient error at the kink of dI at rho_m (the solver's alpha absorbs it)
    assert abs(np.trapezoid(dj, A)) < 5e-3*np.trapezoid(np.abs(dj), A)
    # Rule 2 (local) on a monotone on-axis dip is rule 1
    q_l, dj_l, sc_l = reset(q, 2)
    assert np.allclose(q_l, q_new, rtol=1e-14, atol=0.0) and np.allclose(dj_l, dj, rtol=1e-12, atol=1e-12*np.abs(dj).max())
    assert sc_l['n_dips'] == 1 and sc_l['rho_m'] == sc['rho_m'] and np.isclose(sc_l['rho_out'], sc['rho_out'])
    # q above q_s everywhere: no reset
    for rule in (1, 2):
        q_new2, dj2, sc2 = reset(q + 0.2, rule)
        assert sc2['n_dips'] == 0 and np.all(dj2 == 0.0) and np.all(q_new2 == q + 0.2) and sc2['rho_out'] == 0.0
    # Non-finite q (a failed surface): no reset, no abort, and the failure is returned
    q_nan = q.copy()
    q_nan[100] = np.nan
    q_new2, dj2, sc2 = reset(q_nan, 2)
    assert np.all(dj2 == 0.0) and np.array_equal(q_new2, q_nan, equal_nan=True) and sc2['ierr'] != 0
    assert sc['ierr'] == 0
    # Inputs of different lengths are rejected
    try:
        _saw_reset_1d(rho[:-1], q, A/q, A, c1, q_s, dq)
    except ValueError:
        pass
    else:
        raise AssertionError('_saw_reset_1d accepted inputs of different lengths')
    # Off-axis dip only (q0 > q_s + dq): local leaves q0 and everything outside [rho_in, rho_out] alone,
    # dI = 0 at both ends (Ip-neutral); fuse lowers q0 to q_s. rho_m is the same axis-anchored radius.
    q4 = 1.10 + 0.6*rho**2 - 0.17*g(0.35, 0.08)
    out = {r: reset(q4, r) for r in (1, 2)}
    q_l, dj_l, sc_l = out[2]
    i_in = np.flatnonzero(q4 < q_s + dq)[0]
    rho_in = np.interp(q_s + dq, q4[i_in-1:i_in+1][::-1], rho[i_in-1:i_in+1][::-1])
    inside = (rho > rho_in) & (rho < sc_l['rho_out'])
    assert sc_l['n_dips'] == 1 and 0.2 < rho_in < 0.35 < sc_l['rho_out'] < 0.5
    assert np.all(q_l[~inside] == q4[~inside]) and np.all(sc_l['wr'][~inside] == 0.0) and np.all(sc_l['wr'][inside] == 1.0)
    assert np.isclose(q_l.min(), q_s) and np.all(q_l[inside] >= q_s - 1e-12)
    dI = (A/q4)*(q4/q_l - 1.0)
    assert np.all(dI[~inside] == 0.0)
    near = np.convolve(inside, np.ones(3), mode='same') > 0   # gradient stencil reaches one node out
    assert np.all(dj_l[~near] == 0.0)
    assert abs(np.trapezoid(dj_l, A)) < 5e-3*np.trapezoid(np.abs(dj_l), A)
    assert np.isclose(out[1][0][0], q_s) and out[1][2]['rho_m'] == sc_l['rho_m'] and sc_l['rho_m'] > sc_l['rho_out'] - 1e-12
    # Two dips with the hump between them below q_s + dq: one region, q >= q_s across it, the hump
    # partly kept (merged reset blended toward the split one by the hump height)
    q2 = 1.035 + 0.12*rho**2 - 0.12*g(0.0, 0.08) - 0.05*g(0.4, 0.06)
    q_l, dj_l, sc_l = reset(q2, 2)
    assert sc_l['n_dips'] == 1 and sc_l['rho_out'] > 0.45 and np.isclose(q_l[0], q_s)
    assert np.all(q_l[rho < sc_l['rho_out']] >= q_s - 1e-12)
    assert np.all(q_l[rho >= sc_l['rho_out']] == q2[rho >= sc_l['rho_out']])
    # Hump height through q_s + dq: the merged region splits in two. Unblended, q_new jumps by dq there;
    # blended, only the parabola-vs-node hump-top mismatch remains (~1 % of dq), and q_new >= q_s
    qh = lambda a: reset(0.995 + a*g(0.3, 0.06) + 0.6*rho**6, 2)
    lo, hi = 0.035, 0.075
    assert qh(lo)[2]['n_dips'] == 1 and qh(hi)[2]['n_dips'] == 2
    for _ in range(50):
        lo, hi = ((lo + hi)/2, hi) if qh((lo + hi)/2)[2]['n_dips'] == 1 else (lo, (lo + hi)/2)
    assert np.abs(qh(hi)[0] - qh(lo)[0]).max() < 0.05*dq
    assert all(np.all(r[0][rho < r[2]['rho_out']] >= q_s - 1e-12) for r in map(qh, np.linspace(0.035, 0.075, 41)))
    # Hump above q_s + dq: two regions, the hump between them untouched
    q3 = 1.06 + 0.6*rho**2 - 0.12*g(0.0, 0.1) - 0.12*g(0.35, 0.06)
    out = {r: reset(q3, r) for r in (1, 2)}
    q_l, dj_l, sc_l = out[2]
    hump = (rho > 0.15) & (rho < 0.3) & (q3 >= q_s + dq)
    assert sc_l['n_dips'] == 2 and np.any(hump) and np.all(q_l[hump] == q3[hump]) and np.all(sc_l['wr'][hump] == 0.0)
    assert np.isclose(sc_l['rho_out'], out[1][2]['rho_m']) and sc_l['rho_m'] == out[1][2]['rho_m']
    assert np.all(q_l[rho < sc_l['rho_out']] >= q_s - 1e-12) and np.isclose(q_l[0], q_s)
    # Region reaching the axis, off-axis minimum: q0 in (q_s, q_s + dq) is kept; q0 < q_s is raised to a flat q_s
    for q0, keep in ((1.04, True), (1.0, False)):
        q5 = q0 + 0.6*rho**2 - 0.1*g(0.2, 0.08)
        q_l, dj_l, sc_l = reset(q5, 2)
        assert sc_l['n_dips'] == 1 and np.isclose(q_l.min(), q_s) and np.all(q_l[rho < sc_l['rho_out']] >= q_s - 1e-12)
        if keep:
            assert np.isclose(q_l[0], q5[0], rtol=1e-14) and np.all(q_l[rho < 0.15] <= q5[0] + 1e-12)
        else:
            assert np.allclose(q_l[rho < 0.15], q_s, rtol=1e-14)
        assert np.all(q_l[rho >= sc_l['rho_out']] == q5[rho >= sc_l['rho_out']])
    # As the outer dip crosses q_s (hump above q_s + dq), rule 1 jumps from the inner to the outer mixing
    # radius; rule 2 adds the outer region through its weight ramp, continuously
    bg = 1.06 + 0.6*rho**2 - 0.12*g(0.0, 0.1)
    norms = {}
    for rule in (1, 2):
        nj = []
        for d in np.linspace(0.104, 0.113, 46):   # outer dip minimum crosses q_s near d = 0.1085
            qd = bg - d*g(0.35, 0.06)
            nj.append(np.sqrt(np.trapezoid(reset(qd, rule)[1]**2, A)))
        norms[rule] = np.abs(np.diff(nj)).max() / np.max(nj)
    assert norms[1] > 0.5 and norms[2] < 0.05
