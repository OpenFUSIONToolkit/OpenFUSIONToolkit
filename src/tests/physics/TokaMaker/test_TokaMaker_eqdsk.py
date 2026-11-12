#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: g-file and p-file readers/writers and COCOS conversion in TokaMaker.eqdsk
import os
import pytest
import numpy as np
from tokamaker_testing import validate_dict, load_fixture
test_dir = os.path.abspath(os.path.dirname(__file__))

# Expected values come from the OMFIT readers, from which TokaMaker.eqdsk was adapted, run once offline on the
# same files (see generate_eqdsk_expected.py). Raw file contents match to round-off; derived flux-surface
# quantities use slightly different algorithms (eg. triangularity definition), so eqdsk_gfile_tol.json relaxes
# their tolerances.
eqdsk_gfile_expected_dict = load_fixture('eqdsk_gfile_expected.json')
eqdsk_gfile_tol_dict      = load_fixture('eqdsk_gfile_tol.json')


def run_eqdsk_gfile_case():
    '''Raw and derived quantities from ITER_test.eqdsk, keyed as in eqdsk_gfile_expected.json'''
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    eq = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    raw = eq._raw
    results = {}

    for key in ['NW','NH','RDIM','ZDIM','RCENTR','RLEFT','ZMID',
                'RMAXIS','ZMAXIS','SIMAG','SIBRY','BCENTR','CURRENT']:
        results[key] = float(raw[key])

    for key in ['FPOL','PRES','FFPRIM','PPRIME','QPSI']:
        arr = np.asarray(raw[key])
        results[f'{key}_0']   = float(arr[0])
        results[f'{key}_mid'] = float(arr[len(arr)//2])
        results[f'{key}_end'] = float(arr[-1])
        results[f'{key}_sum'] = float(np.sum(arr))
        results[f'{key}_len'] = len(arr)

    psi = np.asarray(raw['PSIRZ'])
    results['PSIRZ_min'] = float(np.min(psi))
    results['PSIRZ_max'] = float(np.max(psi))
    results['PSIRZ_sum'] = float(np.sum(psi))

    for key in ['RBBBS','ZBBBS','RLIM','ZLIM']:
        arr = np.asarray(raw[key])
        results[f'{key}_min'] = float(np.min(arr))
        results[f'{key}_max'] = float(np.max(arr))
        results[f'{key}_len'] = len(arr)

    geo = eq.geometry
    avg = eq.averages
    mid = eq.midplane
    li  = eq.li
    betas = eq.betas
    results['Ip_enclosed_edge'] = float(avg['ip'][-1])
    results['q_axis']           = float(eq.q_profile[0])
    results['kappa_edge']       = float(geo['kappa'][-1])
    results['delta_edge']       = float(geo['delta'][-1])
    results['a_edge']           = float(geo['a'][-1])
    results['R_geo_edge']       = float(geo['R'][-1])
    results['vol_edge']         = float(geo['vol'][-1])
    results['cxArea_edge']      = float(geo['cxArea'][-1])
    results['li_1']             = float(li['li(1)'])
    results['li_3']             = float(li['li(3)'])
    results['R_mid_edge']       = float(mid['R'][-1])
    results['Btot_mid_edge']    = float(mid['Btot'][-1])
    results['beta_t']           = float(betas.get('beta_t', 0.0))
    results['beta_p']           = float(betas.get('beta_p', 0.0))
    results['beta_n']           = float(betas.get('beta_n', 0.0))
    results['cocos']            = int(eq.cocos)

    # Flux-surface averages are on psi_N_levels, which differs from the g-file grid if nlevels is set
    psi_N_lvl  = eq.psi_N_levels
    j_tor      = eq.j_tor_averaged         # <Jt/R>/<1/R>  standard convention
    j_tor_dir  = eq.j_tor_averaged_direct  # -(p'<R> + FF'<1/R>/mu0)*(2pi)^exp_Bp
    q          = eq.q_profile
    for pct in (10, 25, 50, 75, 90):
        psin_target = pct / 100.0
        results[f'j_tor_psiN_{pct:02d}']        = float(np.interp(psin_target, psi_N_lvl, j_tor))
        results[f'j_tor_direct_psiN_{pct:02d}'] = float(np.interp(psin_target, psi_N_lvl, j_tor_dir))
        results[f'q_psiN_{pct:02d}']            = float(np.interp(psin_target, psi_N_lvl, q))
    return results


def test_eqdsk_gfile():
    '''gEQDSK reader: raw file contents and derived flux-surface quantities (q, shape, internal inductance, beta,
    averaged toroidal current) match the OMFIT reader'''
    results = run_eqdsk_gfile_case()
    assert validate_dict([results], eqdsk_gfile_expected_dict,
                         tol_dict=eqdsk_gfile_tol_dict)


eqdsk_pfile_expected_dict = load_fixture('eqdsk_pfile_expected.json')


def run_eqdsk_pfile_case():
    '''Profile samples from D3Dlike_Hmode_test.peqdsk, keyed as in eqdsk_pfile_expected.json'''
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_pfile
    pf = read_pfile(os.path.join(test_dir, 'D3Dlike_Hmode_test.peqdsk'))
    results = {}

    for key in pf.keys:
        if key == 'N Z A':
            nza = pf[key]
            results['N Z A_N_0'] = float(nza['N'][0])
            results['N Z A_Z_0'] = float(nza['Z'][0])
            results['N Z A_A_0'] = float(nza['A'][0])
            results['N Z A_len'] = len(nza['N'])
            continue
        d = np.asarray(pf._get_data(key))
        psin = np.asarray(pf.psinorm_for(key))
        results[f'{key}_len'] = len(d)
        results[f'{key}_0']   = float(d[0])
        results[f'{key}_mid'] = float(d[len(d)//2])
        results[f'{key}_end'] = float(d[-1])
        results[f'{key}_min'] = float(np.min(d))
        results[f'{key}_max'] = float(np.max(d))
        results[f'{key}_sum'] = float(np.sum(d))
        results[f'{key}_psinorm_end'] = float(psin[-1])
    return results


def test_eqdsk_pfile():
    '''Osborne p-file reader: every profile and its psi_N grid match the OMFIT reader'''
    results = run_eqdsk_pfile_case()
    assert validate_dict([results], eqdsk_pfile_expected_dict)


# COCOS indices 1-8 and 11-18 (O. Sauter and S. Yu. Medvedev, Comput. Phys. Commun. 184, 293 (2013))
_COCOS_VALID = tuple(list(range(1, 9)) + list(range(11, 19)))


def _cocos_params_reference(cocos_index):
    '''Sign conventions for a COCOS index (Sauter and Medvedev), written independently of eqdsk.py'''
    if cocos_index < 1 or cocos_index > 18 or cocos_index in (9, 10):
        raise ValueError(f'Invalid COCOS index: {cocos_index}')
    exp_Bp = 0 if cocos_index < 10 else 1
    base = cocos_index if cocos_index < 10 else cocos_index - 10
    return dict(
        sigma_Bp    = +1 if base in (1, 2, 5, 6) else -1,
        sigma_RpZ   = +1 if base in (1, 2, 7, 8) else -1,
        sigma_rhotp = +1 if base in (1, 3, 5, 7) else -1,
        exp_Bp      = exp_Bp,
    )


def _expected_cocos_factors(cocos_in, cocos_out):
    '''Factor applied to each raw g-file field when converting from `cocos_in` to `cocos_out` (Sauter and
    Medvedev, Eqs. 14 and 23)'''
    cc_in  = _cocos_params_reference(cocos_in)
    cc_out = _cocos_params_reference(cocos_out)
    sBp    = cc_out['sigma_Bp']    * cc_in['sigma_Bp']
    sRpZ   = cc_out['sigma_RpZ']   * cc_in['sigma_RpZ']
    srhotp = cc_out['sigma_rhotp'] * cc_in['sigma_rhotp']
    exp_eff = cc_out['exp_Bp'] - cc_in['exp_Bp']
    twopi_exp = (2.0 * np.pi) ** exp_eff
    psi_fac  = sRpZ * sBp * twopi_exp
    dpsi_fac = sRpZ * sBp / twopi_exp
    bt_fac   = sRpZ
    ip_fac   = sRpZ
    q_fac    = srhotp
    return {
        'SIMAG':   psi_fac,
        'SIBRY':   psi_fac,
        'PSIRZ':   psi_fac,
        'PPRIME':  dpsi_fac,
        'FFPRIM':  dpsi_fac,
        'FPOL':    bt_fac,
        'BCENTR':  bt_fac,
        'CURRENT': ip_fac,
        'QPSI':    q_fac,
    }


@pytest.mark.parametrize('cocos_out', _COCOS_VALID)
def test_eqdsk_cocos_roundtrip(cocos_out):
    '''Converting the COCOS 1 test file to each COCOS convention and back reproduces the original data'''
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    eq_ref = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    eq_rt  = (eq_ref.cocosify(cocos_out, copy=True)
                    .cocosify(eq_ref.cocos, copy=True))
    assert eq_rt.cocos == eq_ref.cocos
    for key in ('SIMAG', 'SIBRY', 'BCENTR', 'CURRENT',
                'FPOL', 'PRES', 'PPRIME', 'FFPRIM', 'QPSI', 'PSIRZ'):
        assert np.allclose(eq_rt._raw[key], eq_ref._raw[key],
                           rtol=1e-12, atol=1e-12), f'{key} mismatch at COCOS {cocos_out}'


@pytest.mark.parametrize('cocos_out', _COCOS_VALID)
def test_eqdsk_cocos_signs(cocos_out):
    '''Converting to each COCOS convention scales every flux, field, current, and q entry by the sign and 2 pi
    factors predicted from the COCOS definitions'''
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    eq_ref = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    eq_n   = eq_ref.cocosify(cocos_out, copy=True)
    assert eq_n.cocos == cocos_out
    factors = _expected_cocos_factors(eq_ref.cocos, cocos_out)
    for key, factor in factors.items():
        expected = np.asarray(eq_ref._raw[key]) * factor
        assert np.allclose(eq_n._raw[key], expected,
                           rtol=1e-12, atol=1e-12), (
            f'{key} expected scale {factor:+g} at COCOS {cocos_out}')


@pytest.mark.parametrize('bad', [0, 9, 10, 19, -1])
def test_eqdsk_cocos_invalid(bad):
    '''Invalid COCOS indices raise ValueError both when reading and when converting'''
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    with pytest.raises(ValueError):
        read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'), cocos=bad)
    eq = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    with pytest.raises(ValueError):
        eq.cocosify(bad, copy=True)


def test_eqdsk_flip_Bt_Ip_roundtrip():
    '''Reversing the toroidal field and plasma current twice restores the original data'''
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    eq_ref = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    eq_rt  = eq_ref.flip_Bt_Ip(copy=True).flip_Bt_Ip(copy=True)
    for key in ('BCENTR', 'FPOL', 'CURRENT', 'SIMAG', 'SIBRY',
                'PSIRZ', 'PPRIME', 'FFPRIM'):
        assert np.allclose(eq_rt._raw[key], eq_ref._raw[key],
                           rtol=1e-12, atol=1e-12), f'{key} not invariant under flip^2'


def test_eqdsk_save_load_roundtrip(tmp_path):
    '''Writing a gEQDSK file and reading it back reproduces the data to the precision of the fixed-format
    text'''
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    eq = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    out_path = tmp_path / 'ITER_roundtrip.eqdsk'
    eq.save(str(out_path))
    eq_rt = read_geqdsk(str(out_path))
    for key in ('NW', 'NH', 'RDIM', 'ZDIM', 'RCENTR', 'RLEFT', 'ZMID',
                'RMAXIS', 'ZMAXIS', 'SIMAG', 'SIBRY', 'BCENTR', 'CURRENT'):
        assert np.isclose(eq._raw[key], eq_rt._raw[key], rtol=1e-6, atol=0), key
    for key in ('FPOL', 'PRES', 'PPRIME', 'FFPRIM', 'QPSI', 'PSIRZ',
                'RBBBS', 'ZBBBS', 'RLIM', 'ZLIM'):
        assert np.allclose(eq._raw[key], eq_rt._raw[key],
                           rtol=1e-6, atol=1e-10), key


def test_eqdsk_cocos_bytes_roundtrip():
    '''Serializing a gEQDSK equilibrium to bytes and back reproduces the data'''
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk, GEQDSKEquilibrium
    eq = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    raw_bytes = eq.to_bytes()
    eq_rt = GEQDSKEquilibrium.from_bytes(raw_bytes, cocos=eq.cocos)
    for key in ('FPOL', 'PRES', 'PSIRZ', 'QPSI'):
        assert np.allclose(eq._raw[key], eq_rt._raw[key],
                           rtol=1e-6, atol=1e-10), key


def test_pfile_bytes_roundtrip():
    '''Serializing a p-file to bytes and back reproduces the profiles'''
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_pfile, PFile
    pf = read_pfile(os.path.join(test_dir, 'D3Dlike_Hmode_test.peqdsk'))
    pf_rt = PFile.from_bytes(pf.to_bytes())
    for key in ('ne', 'te', 'ni', 'ti'):
        d_ref = np.asarray(pf._get_data(key))
        d_rt  = np.asarray(pf_rt._get_data(key))
        assert np.allclose(d_ref, d_rt, rtol=1e-6, atol=1e-10), key
