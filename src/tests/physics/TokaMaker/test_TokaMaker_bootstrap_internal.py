#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: self-consistent bootstrap current from the internal solve_bootstrap() on ITER
import pytest
import numpy as np
from tokamaker_testing import mp_run, validate_dict, setup_ITER_bootstrap_internal
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.util import mu0, eC


def _check_boot_state_restored(mygs, restore, label):
    '''Check that restore() (a replace_eq call) restores the bootstrap options shadow and profiles'''
    ops_ref, profs_ref = dict(mygs._tMaker_equil._boot_ops), mygs.get_boot_profs()
    mygs._tMaker_equil._boot_ops['scale_jbs'] = -999.0  # the restore must overwrite the shadow
    restore()
    ops, profs = mygs._tMaker_equil._boot_ops, mygs.get_boot_profs()
    if ops is None or profs is None or profs_ref is None:
        raise AssertionError(f"_boot_ops or boot_profs missing around replace_eq({label}=...)")
    for key, ref in ops_ref.items():
        if not np.isclose(ops[key], ref, rtol=0.0, atol=1e-10):
            raise AssertionError(f"_boot_ops['{key}'] = {ops[key]} != {ref} after replace_eq({label}=...)")
    for key, ref in profs_ref.items():
        if key not in profs:
            raise AssertionError(f"boot_profs['{key}'] missing after replace_eq({label}=...)")
        if not np.allclose(profs[key], ref, rtol=1e-3):
            rel = np.max(np.abs(profs[key] - ref) / (np.abs(ref) + 1e-30))
            raise AssertionError(f"boot_profs['{key}'] differs after replace_eq({label}=...) (max rel diff {rel:.3e})")


def run_ITER_bootstrap_case_internal(mesh_resolution, fe_order, mp_q):
    setup = setup_ITER_bootstrap_internal(mesh_resolution, fe_order)
    if setup is None:
        mp_q.put(None)
        return
    mygs, myOFT, psi_sample, n_sample, Ip_target, Zeff_val, ne, Te, ni, Ti, inductive_jphi = (
        setup[k] for k in ('mygs', 'myOFT', 'psi_sample', 'n_sample', 'Ip_target', 'Zeff_val',
                           'ne', 'Te', 'ni', 'Ti', 'inductive_jphi'))

    try:
        mygs.solve_bootstrap(
            ffp_prof={'type': 'jphi-split-bootstrap', 'x': psi_sample, 'y': inductive_jphi},
            te_prof={'type': 'linterp', 'x': psi_sample, 'y': Te / 1e3},
            ne_prof={'type': 'linterp', 'x': psi_sample, 'y': ne},
            ti_prof={'type': 'linterp', 'x': psi_sample, 'y': Ti / 1e3},
            ni_prof={'type': 'linterp', 'x': psi_sample, 'y': ni},
            Zeff=Zeff_val,
            Ip_target=Ip_target,
            scale_jbs=1.0,
            diagnose_bs=True,
        )
    except Exception:
        mp_q.put(None)
        return

    eq_info = mygs.get_stats(li_normalization='ITER')
    psi_i, F_i, Fp_i, P_i, Pp_i = mygs.get_profiles(npsi=n_sample, psi_pad=1e-3)
    _, q_i, ravgs_i, _, _, _     = mygs.get_q(npsi=n_sample, psi_pad=1e-3)
    jtor_i = F_i * Fp_i * ravgs_i['<1/R>'] / mu0 + Pp_i * ravgs_i['<R>']

    eq_info['jphi_axis'] = float(jtor_i[0])
    eq_info['jphi_max']  = float(np.max(np.abs(jtor_i)))
    eq_info['q_axis']    = float(q_i[0])
    sample_idx = np.round(np.linspace(0, len(jtor_i) - 1, 10)).astype(int)
    eq_info['jphi_prof'] = [float(jtor_i[i]) for i in sample_idx]

    # Conversion of the bootstrap <J.B> to j_phi (doc_tokamaker_current_conventions eqs. A7, A8), checked against
    # flux-surface geometry traced from the converged equilibrium, away from the axis and edge
    try:
        bp = mygs.get_boot_profs()
        inner = (bp['psi_n'] > 0.05) & (bp['psi_n'] < 0.95)
        psi_c = bp['psi_n'][inner]
        _, F_c, Fp_c, _, Pp_c = mygs.get_profiles(psi=psi_c)
        _, _, ravgs_c, _, _, _ = mygs.get_q(psi=psi_c)
        _, _, _, modb_c = mygs.sauter_fc(psi=psi_c)
        R_c, invR_c, B2_c = ravgs_c['<R>'], ravgs_c['<1/R>'], modb_c[1]
        jdotB = bp['jdotB_bs_raw'][inner]
        field_aligned = jdotB * F_c * invR_c / B2_c
        p_term = Pp_c * (R_c - F_c**2 * invR_c / B2_c)
        scale = np.max(np.abs(bp['jphi_bs_raw'][inner]))
        j_bs_raw = bp['jphi_bs_raw'][inner] + bp['jphi_pd_raw'][inner]
        err = np.max(np.abs(j_bs_raw - (field_aligned + p_term))) / scale
        # The flipped pressure term must be clearly worse, so the check is sensitive to its sign
        err_flip = np.max(np.abs(j_bs_raw - (field_aligned - p_term))) / scale
        print(f"jphi_bs_raw vs A7: rel err {err:.3e} (pressure term flipped: {err_flip:.3e}), "
              f"max|P'G|/max|jphi_bs_raw| = {np.max(np.abs(p_term))/scale:.3e}")
        if not (err < 1.0e-2 and err_flip > 5.0 * err):
            raise AssertionError(f"jphi_bs_raw does not match eq. A7 (rel err {err:.3e}, flipped {err_flip:.3e})")
        # The stored pressure-driven part is the A7 P'G term
        err_pd = np.max(np.abs(bp['jphi_pd_raw'][inner] - p_term)) / scale
        print(f"jphi_pd_raw vs P'G: rel err {err_pd:.3e}")
        if err_pd > 1.0e-2:
            raise AssertionError(f"jphi_pd_raw does not match P'(<R> - F^2<1/R>/<B^2>) (rel err {err_pd:.3e})")
        # A8: the equilibrium's own <J.B> = F P' + F'<B^2>/mu0 equals the sum of its components'
        jdotB_eq = F_c * Pp_c + Fp_c * B2_c / mu0
        jdotB_sum = jdotB + (bp['jphi_ind'][inner] + bp['jphi_fixed'][inner]) * B2_c / (F_c * invR_c)
        err_par = np.max(np.abs(jdotB_eq - jdotB_sum)) / np.max(np.abs(jdotB_eq))
        print(f"<J.B> balance (A8): rel err {err_par:.3e}")
        if err_par > 2.0e-2:
            raise AssertionError(f"equilibrium <J.B> != sum of components (rel err {err_par:.3e})")
        # A9d: I_p integrated from TokaMaker's flux-surface j_phi equals the finite element I_p
        psi_f = np.linspace(1.0e-4, 1.0 - 1.0e-4, 1001)
        _, F_f, Fp_f, _, Pp_f = mygs.get_profiles(psi=psi_f)
        _, _, rv_f, _, _, _ = mygs.get_q(psi=psi_f)
        J_f = F_f * Fp_f * rv_f['<1/R>'] / mu0 + Pp_f * rv_f['<R>']
        psi_phys = mygs.psi_bounds[0] + psi_f * (mygs.psi_bounds[1] - mygs.psi_bounds[0])
        w = rv_f['dV/dPsi'] / (2.0 * np.pi)
        ip_exact = abs(np.trapezoid(w * (J_f * rv_f['<1/R^2>'] / rv_f['<1/R>']
                   + Pp_f * (1.0 - rv_f['<R>'] * rv_f['<1/R^2>'] / rv_f['<1/R>'])), psi_phys))
        ip_qtmp = abs(np.trapezoid(w * J_f / rv_f['<R>'], psi_phys))
        print(f"I_p: FEM {eq_info['Ip']:.6e}, A9d {ip_exact:.6e} ({ip_exact/eq_info['Ip']-1:+.3e}), "
              f"old <R><1/R> measure {ip_qtmp:.6e} ({ip_qtmp/eq_info['Ip']-1:+.3e})")
        if abs(ip_exact / eq_info['Ip'] - 1.0) > 2.0e-3:
            raise AssertionError(f"A9d I_p {ip_exact:.6e} != FEM I_p {eq_info['Ip']:.6e}")
    except Exception:
        import traceback
        traceback.print_exc()
        mp_q.put(None)
        return

    # Bootstrap options and profiles survive save/load and copy_eq, and replace_eq updates the Python-side copy
    save_file = 'ITER_boot_ops_test.h5'
    try:
        _, F_s, Fp_s, _, _ = mygs.get_profiles(npsi=50)
        mygs._tMaker_equil.save_TokaMaker(save_file)
        _check_boot_state_restored(mygs, lambda: mygs.replace_eq(source_file=save_file), 'source_file')
        # F*F' round-trips exactly because the saved profile knots are restored, not recomputed
        _, F_l, Fp_l, _, _ = mygs.get_profiles(npsi=50)
        ffp_err = np.max(np.abs(F_l*Fp_l - F_s*Fp_s)) / np.max(np.abs(F_s*Fp_s))
        if ffp_err > 1e-10:
            raise AssertionError(f"F*F' changed through save/load (max rel diff {ffp_err:.3e})")
        mygs_copied = mygs.copy_eq()
        _check_boot_state_restored(mygs, lambda: mygs.replace_eq(source_eq=mygs_copied), 'source_eq')
    except Exception:
        import traceback
        traceback.print_exc()
        mp_q.put(None)
        return

    # A Zeff profile must be used, not just its scalar value
    try:
        zeff_common_kwargs = dict(
            ffp_prof={'type': 'jphi-split-bootstrap', 'x': psi_sample, 'y': inductive_jphi},
            te_prof={'type': 'linterp', 'x': psi_sample, 'y': Te / 1e3},
            ne_prof={'type': 'linterp', 'x': psi_sample, 'y': ne},
            ti_prof={'type': 'linterp', 'x': psi_sample, 'y': Ti / 1e3},
            ni_prof={'type': 'linterp', 'x': psi_sample, 'y': ni},
            Ip_target=Ip_target,
            scale_jbs=1.0,
            diagnose_bs=False,
        )
        profs_scalar = mygs.solve_bootstrap(Zeff=Zeff_val, **zeff_common_kwargs)
        profs_linear = mygs.solve_bootstrap(
            Zeff={'x': psi_sample, 'y': np.linspace(1.0, 2.5, n_sample)},
            **zeff_common_kwargs,
        )
        j_scalar = profs_scalar['jphi_bs_raw']
        j_linear = profs_linear['jphi_bs_raw']
        magnitude = 0.5 * (np.abs(j_scalar) + np.abs(j_linear))
        rel_diff = np.where(magnitude > 0, (j_linear - j_scalar) / magnitude, 0.0)
        print("\nZeff scalar vs linear jphi_bs_raw relative difference (j_linear-j_scalar)/|mean|:")
        print(rel_diff)
        if np.allclose(j_scalar, j_linear):
            raise AssertionError(
                "Bootstrap profiles with scalar Zeff and linearly-increasing Zeff profile "
                "are identical; expected them to differ."
            )
    except Exception as e:
        print(e)
        mp_q.put(None)
        return

    # A fixed current profile is added unscaled and displaces inductive current; it persists until boot_reset()
    try:
        jfix_y = 2.0e5 * np.exp(-((psi_sample - 0.5) / 0.1)**2)
        profs_fixed = mygs.solve_bootstrap(
            Zeff=Zeff_val,
            jphi_fixed_prof={'type': 'linterp', 'x': psi_sample, 'y': jfix_y},
            **zeff_common_kwargs,
        )
        if not np.allclose(profs_fixed['jphi_fixed'], np.interp(profs_fixed['psi_n'], psi_sample, jfix_y),
                           rtol=1e-6, atol=1e-6*jfix_y.max()):
            raise AssertionError("boot_profs['jphi_fixed'] does not match the input jphi_fixed_prof")
        j_sum = profs_fixed['jphi_ind'] + profs_fixed['jphi_bs'] + profs_fixed['jphi_pd'] + profs_fixed['jphi_fixed']
        if not np.allclose(profs_fixed['jphi_total'], j_sum, rtol=1e-10, atol=1e-6):
            raise AssertionError("jphi_total != jphi_ind + jphi_bs + jphi_pd + jphi_fixed")
        if not np.sum(profs_fixed['jphi_ind']) < np.sum(profs_scalar['jphi_ind']):
            raise AssertionError("jphi_fixed did not reduce the inductive current share")
        profs_kept = mygs.solve_bootstrap(Zeff=Zeff_val, **zeff_common_kwargs)
        if not np.array_equal(profs_kept['jphi_fixed'], profs_fixed['jphi_fixed']):
            raise AssertionError("jphi_fixed not kept by solve_bootstrap without jphi_fixed_prof")
        try:  # The kept profile is on psi_N, so it cannot be reused with toroidal flux coordinates
            mygs.solve_bootstrap(Zeff=Zeff_val, coord='phi_n', **zeff_common_kwargs)
        except ValueError as err:
            if 'kept coord' not in str(err):
                raise
        else:
            raise AssertionError("kept jphi_fixed_prof reused on a different coord")
        mygs.boot_reset()
        profs_reset = mygs.solve_bootstrap(Zeff=Zeff_val, **zeff_common_kwargs)
        if np.any(profs_reset['jphi_fixed'] != 0.0):
            raise AssertionError("jphi_fixed persisted after boot_reset()")
    except Exception as e:
        print(e)
        mp_q.put(None)
        return

    # A fixed pressure [Pa] adds to the kinetic pressure, matching pres_prof = kinetic + p_fixed; the two options
    # are mutually exclusive and negative values are rejected
    try:
        p_kin = eC * (ne * Te + ni * Ti)
        # Offset keeps pres_prof safely above the kinetic pressure (roundoff) at the edge
        pf = p_kin[0] * (0.2 * np.exp(-((psi_sample - 0.3) / 0.15)**2) + 1.0e-3)
        profs_pfix = mygs.solve_bootstrap(Zeff=Zeff_val, p_fixed_prof={'x': psi_sample, 'y': pf}, **zeff_common_kwargs)
        _, _, _, P_pfix, _ = mygs.get_profiles(npsi=n_sample)
        P_ax_pfix = np.max(mygs.get_profiles(psi=np.array([0.0, 1.0]))[3])
        profs_pres = mygs.solve_bootstrap(Zeff=Zeff_val, pres_prof={'x': psi_sample, 'y': p_kin + pf}, **zeff_common_kwargs)
        _, _, _, P_pres, _ = mygs.get_profiles(npsi=n_sample)
        # Both solves converge jphi_bs to djbs_tol (1e-4) and may stop an iteration apart
        for key in profs_pres:
            if not np.allclose(profs_pfix[key], profs_pres[key], rtol=2e-4, atol=2e-4*np.max(np.abs(profs_pres[key]))):
                raise AssertionError(f"p_fixed_prof vs pres_prof=kinetic+p_fixed: boot_profs['{key}'] differ")
        if not np.allclose(P_pfix, P_pres, rtol=1e-6, atol=1e-6*np.max(P_pres)):
            raise AssertionError("p_fixed_prof vs pres_prof=kinetic+p_fixed: pressure profiles differ")
        # TokaMaker holds P = 0 at the LCFS: the axis value is the kinetic + p_fixed drop
        p_tot = p_kin + pf
        if not np.isclose(P_ax_pfix, p_tot[0] - p_tot[-1], rtol=1e-4):
            raise AssertionError(f"axis pressure {P_ax_pfix:.6e} != kinetic + p_fixed drop {p_tot[0] - p_tot[-1]:.6e}")
        mygs.boot_reset()
        mygs.solve_bootstrap(Zeff=Zeff_val, **zeff_common_kwargs)
        if np.isclose(np.max(mygs.get_profiles(psi=np.array([0.0, 1.0]))[3]), P_ax_pfix, rtol=1e-4):
            raise AssertionError("p_fixed_prof did not change the axis pressure")
        for bad_kwargs, label in (({'pres_prof': {'x': psi_sample, 'y': p_kin + pf}}, "with pres_prof"),
                                  ({}, "negative")):
            try:
                y_bad = -pf if label == "negative" else pf
                mygs.solve_bootstrap(Zeff=Zeff_val, p_fixed_prof={'x': psi_sample, 'y': y_bad}, **bad_kwargs, **zeff_common_kwargs)
            except ValueError:
                continue
            raise AssertionError(f"p_fixed_prof {label} did not raise ValueError")
    except Exception as e:
        print(e)
        mp_q.put(None)
        return

    # Isolating the edge bootstrap current acts on the field-aligned part only: the pressure-driven part is
    # unchanged and the field-aligned part is flat in the core
    try:
        profs_iso = mygs.solve_bootstrap(Zeff=Zeff_val, isolate_edge_jbs=True, **zeff_common_kwargs)
        if not np.allclose(profs_iso['jphi_pd'], profs_iso['jphi_pd_raw'], rtol=1e-12, atol=1e-6):
            raise AssertionError("isolate_edge_jbs changed jphi_pd")
        core = profs_iso['jphi_bs'][profs_iso['psi_n'] < 0.3]
        if np.ptp(core) > 1e-10 * np.max(np.abs(profs_iso['jphi_bs'])):
            raise AssertionError(f"isolated jphi_bs not flat in the core (spread {np.ptp(core):.3e} A/m^2)")
        j_sum = profs_iso['jphi_ind'] + profs_iso['jphi_bs'] + profs_iso['jphi_pd'] + profs_iso['jphi_fixed']
        if not np.allclose(profs_iso['jphi_total'], j_sum, rtol=1e-10, atol=1e-6):
            raise AssertionError("isolated: jphi_total != jphi_ind + jphi_bs + jphi_pd + jphi_fixed")
    except Exception as e:
        print(e)
        mp_q.put(None)
        return

    mp_q.put([eq_info])
    oftpy_dump_cov()


ITER_bootstrap_internal_eq_dict = {
    'Ip':        13000003.24961143,
    'kappa':     1.8710874971500062,
    'R_geo':     6.2229574417901645,
    'a_geo':     1.9785979198829686,
    'q_0':       1.3748401369017682,
    'q_95':      3.504083741872661,
    # P_ax = p(axis) - p(LCFS): pax targets the pressure drop carried by TokaMaker's P'
    'P_ax':      7.232028E+05,
    'beta_pol':  8.320959E+01,
    'beta_tor':  2.435303E+00,
    'jphi_axis': 1010416.1199979713,
    'jphi_max':  1173298.5135544762,
    'q_axis':    1.410565702638605,
    'jphi_prof': [1010416.1199979713, 1155117.92412601,   1162817.6369466858,
                  1072123.3558973635,  898492.2764938722,  683025.6545245483,
                   449301.00087938644, 249118.67207980127, 173771.68344307062,
                   136306.8890837856],
}


@pytest.mark.slow
@pytest.mark.parametrize("order", (2,))
def test_ITER_bootstrap_internal(order):
    '''Self-consistent bootstrap current solve with TokaMaker's internal solve_bootstrap() on ITER H-mode profiles:
    the equilibrium and current profile match reference values, the current decomposition is consistent with the
    flux-surface geometry, and the bootstrap state survives save/load and copy. Also covers Zeff profiles, fixed
    current and pressure profiles, and edge bootstrap isolation.

    References: A. Redl et al., Phys. Plasmas 28, 022502 (2021); doc_tokamaker_current_conventions.
    '''
    results = mp_run(run_ITER_bootstrap_case_internal, (1.0, order), timeout=300)
    assert validate_dict(results, ITER_bootstrap_internal_eq_dict)
