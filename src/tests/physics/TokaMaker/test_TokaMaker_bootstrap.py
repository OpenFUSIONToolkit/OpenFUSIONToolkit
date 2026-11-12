#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: bootstrap current from the Python solve_with_bootstrap() workflow and Redl model on ITER
import pytest
import numpy as np
from tokamaker_testing import mp_run, validate_dict, solve_ITER_reference, ITER_Hmode_profiles
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.util import eC
from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun


def run_ITER_bootstrap_case(mesh_resolution, fe_order, mp_q):
    from OpenFUSIONToolkit.TokaMaker.bootstrap import solve_with_bootstrap
    mygs, _ = solve_ITER_reference(mesh_resolution, fe_order)
    if mygs is None:
        mp_q.put(None)
        return
    Ip_target = 15.6E6
    n_sample = 257
    psi_sample = np.linspace(0.0, 1.0, n_sample)
    psi_pad = 1.E-3
    jphi_prof = create_power_flux_fun(len(psi_sample), 2.25, 2.5)
    inductive_jphi = jphi_prof['y']
    ne, Te = ITER_Hmode_profiles(n_sample)
    ni = ne.copy()       # Assuming quasineutrality
    Ti = Te.copy()       # Assuming isothermal
    Zeff = np.full(n_sample, 1.7)
    try:
        bs_results = solve_with_bootstrap(
            mygs,
            ne, Te, ni, Ti, Zeff,
            Ip_target,
            inductive_jphi,
            scale_jbs=1.0,
            isolate_edge_jbs=False,
            psi_pad=psi_pad,
            iterations=2,
            diagnostic_plots=False,
            use_python_solve=True,
        )
    except Exception as e:
        print("Bootstrap solve failed: {0}".format(e))
        mp_q.put(None)
        return
    eq_info = mygs.get_stats(li_normalization='ITER')
    # Bootstrap j_phi including its pressure-driven part (doc_tokamaker_current_conventions eq. A7)
    jphi_bs_raw = bs_results['jphi_bs_raw'] + bs_results['jphi_pd_raw']
    jphi_total = bs_results['jphi_total']
    jphi_ind = bs_results['jphi_ind']
    eq_info['jphi_bs_max'] = float(np.max(np.abs(jphi_bs_raw)))
    eq_info['jphi_bs_axis'] = float(jphi_bs_raw[0])
    eq_info['jphi_axis'] = float(jphi_total[0])
    eq_info['jphi_max'] = float(np.max(np.abs(jphi_total)))
    eq_info['jphi_ind_axis'] = float(jphi_ind[0])
    # Bootstrap fraction estimated in psi_N rather than by area, sufficient for regression
    bs_frac = np.trapezoid(jphi_bs_raw, psi_sample) / np.trapezoid(jphi_total, psi_sample) \
              if np.trapezoid(jphi_total, psi_sample) != 0 else 0.0
    eq_info['bs_fraction'] = float(bs_frac)
    mp_q.put([eq_info])
    oftpy_dump_cov()


ITER_bootstrap_eq_dict = {
    'Ip': 15599997.261988742,
    'kappa': 1.8746806271900014,
    'R_geo': 6.222524490498655,
    'a_geo': 1.9807072056151687,
    'q_0': 1.036565329524374,
    'q_95': 2.863029898626678,
    # P_ax = p(axis) - p(LCFS): pax targets the pressure drop carried by TokaMaker's P'
    'P_ax': 7.232048E+05,
    'jphi_bs_max': 219163.55709184994,
    'jphi_bs_axis': 4696.213397541225,
    'jphi_axis': 1397130.088093367,
    'jphi_max': 1507460.349786304,
    # returned jphi_ind at the axis: alpha * seed (the Python path's jphi_ind excludes scale_j0)
    'jphi_ind_axis': 1.392583E+06,
    'bs_fraction': 0.17910707149956598,
}

@pytest.mark.parametrize("order", (2,))
def test_ITER_bootstrap(order):
    '''Bootstrap current iteration with the Python solve_with_bootstrap() workflow on the ITER equilibrium with
    H-mode profiles and the Redl model; equilibrium quantities, current density profile values, and the bootstrap
    fraction are compared to reference values.

    Reference: A. Redl et al., Phys. Plasmas 28, 022502 (2021).
    '''
    results = mp_run(run_ITER_bootstrap_case, (1.0, order), timeout=300)
    assert validate_dict(results, ITER_bootstrap_eq_dict)


def run_Redl_jBS_case(mesh_resolution, fe_order, mp_q):
    from OpenFUSIONToolkit.TokaMaker.bootstrap import redl_bootstrap, calculate_ln_lambda
    mygs, _ = solve_ITER_reference(mesh_resolution, fe_order)
    if mygs is None:
        mp_q.put(None)
        return
    # Geometry, gradients, Coulomb logarithms, and collisionalities are computed here as in solve_with_bootstrap(),
    # so redl_bootstrap() is tested on its own with the same inputs as test_ITER_bootstrap
    n_psi = 257
    psi_N = np.linspace(0.0, 1.0, n_psi)
    psi_pad = 1.E-3

    ne, Te = ITER_Hmode_profiles(n_psi)
    ni = ne.copy()
    Ti = Te.copy()
    Zeff = np.full(n_psi, 1.7)
    # Flux-surface quantities are ill-defined exactly on the axis and LCFS
    psi_eval = np.clip(psi_N, psi_pad, 1.0 - psi_pad)
    _, f, _, _, pp_eq = mygs.get_profiles(psi=psi_eval)
    _, fc, r_avgs, b_avgs, eps = mygs.sauter_fc(psi=psi_eval, return_eps=True)
    ft = 1 - fc
    _, qvals, ravgs_q, _, _, _ = mygs.get_q(psi=psi_eval)
    R_avg = ravgs_q['<R>']
    fsa = mygs.get_fsa(psi=psi_eval)
    R_geo = 0.5*(fsa['R_min'] + fsa['R_max'])
    psi_range = mygs.psi_bounds[1] - mygs.psi_bounds[0]
    dn_e_dpsi = np.gradient(ne, psi_N, edge_order=2) / psi_range
    dT_e_dpsi = np.gradient(Te, psi_N, edge_order=2) / psi_range
    dn_i_dpsi = np.gradient(ni, psi_N, edge_order=2) / psi_range
    dT_i_dpsi = np.gradient(Ti, psi_N, edge_order=2) / psi_range
    ln_le, ln_lii = calculate_ln_lambda(
        Te, Ti, ne, ni, Zeff,
        electron_lnLambda_model='NRL',
        ion_lnLambda_model='Zavg',
    )
    # Normalized collisionalities (Sauter et al. 1999, Eqs. 18b and 18c)
    Zdom = 1.0
    Zavg = ne / ni
    Zion = (Zdom**2 * Zavg * Zeff)**0.25
    nu_i_star = (4.90e-18 * np.abs(qvals) * R_geo * ni
                 * Zion**4 * ln_lii / (Ti**2 * eps**1.5))
    nu_e_star = (6.921e-18 * np.abs(qvals) * R_geo * ne
                 * Zeff * ln_le / (Te**2 * eps**1.5))
    try:
        j_BS_neo, coeffs = redl_bootstrap(
            psi_N=psi_N, Te=Te, Ti=Ti, ne=ne, ni=ni,
            pe=eC*(ne*Te), pi=eC*(ni*Ti),
            Zeff=Zeff, R=R_geo, q=qvals, eps=eps, fT=ft, I_psi=f,
            dT_e_dpsi=dT_e_dpsi, dT_i_dpsi=dT_i_dpsi,
            dn_e_dpsi=dn_e_dpsi, dn_i_dpsi=dn_i_dpsi,
            ln_lambda_e=ln_le, ln_lambda_ii=ln_lii,
            nu_e_star_override=nu_e_star,
            nu_i_star_override=nu_i_star,
            use_legacy_L34=False,
            use_sign_q=False,
            formula_form='jboot1',
        )
    except Exception as e:
        print("redl_bootstrap failed: {0}".format(e))
        mp_q.put(None)
        return
    # Convert <J.B> to j_phi [A/m^2] (doc_tokamaker_current_conventions eq. A7)
    inv_R, B2 = ravgs_q['<1/R>'], b_avgs[1]
    jphi_bs = j_BS_neo * f * inv_R / B2 + pp_eq * (R_avg - f**2 * inv_R / B2)
    jphi_bs = np.nan_to_num(jphi_bs, nan=0.0)
    results = {}
    results['jphi_bs_max'] = float(np.max(np.abs(jphi_bs)))
    results['jphi_bs_axis'] = float(jphi_bs[0])
    results['jphi_bs_edge'] = float(jphi_bs[-1])
    results['L31_axis'] = float(coeffs['L31'][0])
    results['L32_axis'] = float(coeffs['L32'][0])
    results['alpha_axis'] = float(coeffs['alpha'][0])
    results['nu_e_star_axis'] = float(coeffs['nu_e_star'][0])
    results['nu_i_star_axis'] = float(coeffs['nu_i_star'][0])
    mp_q.put([results])
    oftpy_dump_cov()


Redl_jBS_eq_dict = {
    'jphi_bs_max': 195234.5633663803,
    'jphi_bs_axis': 5182.40772662609,
    'jphi_bs_edge': 88568.31833389413,
    'L31_axis': 0.10849481301194958,
    'L32_axis': -0.014581267675566556,
    'alpha_axis': -0.6150406171285213,
    'nu_e_star_axis': 0.3420647094964808,
    'nu_i_star_axis': 0.2975729435421919,
}


@pytest.mark.parametrize("order", (2,))
def test_Redl_jBS(order):
    '''Redl bootstrap current model evaluated directly on the ITER equilibrium with H-mode profiles: current
    density, transport coefficients, and collisionalities on axis are compared to reference values.

    References: O. Sauter et al., Phys. Plasmas 6, 2834 (1999); A. Redl et al., Phys. Plasmas 28, 022502 (2021).
    '''
    results = mp_run(run_Redl_jBS_case, (1.0, order), timeout=300)
    assert validate_dict(results, Redl_jBS_eq_dict)


@pytest.mark.coverage
@pytest.mark.parametrize("use_python_solve", (True, False))
def test_bootstrap_x_validation(use_python_solve):
    '''solve_with_bootstrap() rejects invalid profile grids (wrong length, unsorted, non-finite, not spanning [0,1],
    too few points, or padding too large) and handles the deprecated psi_N argument, before any solve'''
    from OpenFUSIONToolkit.TokaMaker.bootstrap import solve_with_bootstrap
    n = 65
    ne = np.full(n, 1.0e20)
    Te = np.full(n, 2.0e3)
    Zeff = np.full(n, 1.7)
    def call(x=None, psi_pad=1.E-3, psi_N=None, m=n):
        return solve_with_bootstrap(None, ne[:m], Te[:m], ne[:m].copy(), Te[:m].copy(), Zeff[:m],
                                    1.0e6, inductive_jphi=np.ones(m), psi_pad=psi_pad, x=x, psi_N=psi_N,
                                    use_python_solve=use_python_solve)
    good = np.linspace(0.0, 1.0, n)
    dup = good.copy()
    dup[32] = dup[31]
    unsorted_grid = good.copy()
    unsorted_grid[[10, 11]] = good[[11, 10]]
    nan_grid = good.copy()
    nan_grid[5] = np.nan
    with pytest.raises(ValueError, match="same length"):
        call(np.linspace(0.0, 1.0, n-1))
    with pytest.raises(ValueError, match="strictly increasing"):
        call(dup)
    with pytest.raises(ValueError, match="strictly increasing"):
        call(unsorted_grid)
    with pytest.raises(ValueError, match="non-finite"):
        call(nan_grid)
    with pytest.raises(ValueError, match="span"):
        call(np.linspace(-0.1, 1.0, n))
    with pytest.raises(ValueError, match="span"):
        call(np.linspace(0.3, 1.0, n))
    with pytest.raises(ValueError, match="span"):
        call(np.linspace(0.0, 0.5, n))
    with pytest.raises(ValueError, match="at least 5 points"):
        call(None, m=4)
    if use_python_solve:  # psi_pad is only used by the Python solver
        with pytest.raises(ValueError, match="exceeds the first/last"):
            call(good, psi_pad=0.5)
    # The deprecated psi_N alias warns, then is validated like x
    with pytest.warns(DeprecationWarning, match="psi_N"):
        with pytest.raises(ValueError, match="same length"):
            call(psi_N=good[:-1])
    with pytest.raises(ValueError, match="x only"):
        call(good, psi_N=good)


@pytest.mark.coverage
@pytest.mark.parametrize("use_python_solve", (True, False))
def test_bootstrap_deprecated_names(use_python_solve):
    '''Deprecated solve_with_bootstrap() argument and result names warn but still work, and passing both old and
    new argument names is an error'''
    from OpenFUSIONToolkit.TokaMaker.bootstrap import solve_with_bootstrap, _BootstrapResults
    n = 65
    ne = np.full(n, 1.0e20)
    Te = np.full(n, 2.0e3)
    bad_grid = np.linspace(0.0, 1.0, n)
    bad_grid[[10, 11]] = bad_grid[[11, 10]]
    with pytest.warns(DeprecationWarning, match="scale_jBS"):
        with pytest.raises(ValueError, match="strictly increasing"):
            solve_with_bootstrap(None, ne, Te, ne.copy(), Te.copy(), 1.7, 1.0e6, inductive_jphi=np.ones(n),
                                 scale_jBS=1.0, psi_N=bad_grid, use_python_solve=use_python_solve)
    with pytest.raises(TypeError, match="scale_jBS"):
        solve_with_bootstrap(None, ne, Te, ne.copy(), Te.copy(), 1.7, 1.0e6, inductive_jphi=np.ones(n),
                             scale_jbs=2.0, scale_jBS=2.0, use_python_solve=use_python_solve)
    results = _BootstrapResults({'jphi_total': 1})
    with pytest.warns(DeprecationWarning, match="total_j_phi"):
        assert results['total_j_phi'] == 1
    with pytest.warns(DeprecationWarning, match="total_j_phi"):
        assert 'total_j_phi' in results
