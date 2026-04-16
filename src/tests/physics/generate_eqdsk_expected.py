#!/usr/bin/env python3
"""Regenerate the hardcoded expected-value dicts for the TokaMaker eqdsk tests.

Reads `ITER_test.eqdsk` and `D3Dlike_Hmode_test.peqdsk` with the OMFIT
`OMFITgeqdsk` / `OMFITpFile` readers (which the TokaMaker `eqdsk` module
was adapted from) and prints two Python dicts on stdout.  Paste the output
into `test_TokaMaker.py` as `eqdsk_gfile_expected_dict` and
`eqdsk_pfile_expected_dict` if the test input files change.

Raw scalar/profile quantities come from OMFIT's raw parser and should match
the TokaMaker reader to floating-point precision.  Derived flux-surface
quantities (q, geometry, li, betas, j_tor samples) come from OMFIT's
`fluxSurfaces` tracer and the TokaMaker reader is expected to be "close
enough" (within ~0.5-10% depending on the quantity -- see per-key
tolerances in the test).

Requires `omfit_classes` from an environment with scipy<1.14 and numpy<2.0
(the standard `omfit_env` conda env works).  Compat shims below allow the
raw g-file/p-file parsers to run under newer scipy/numpy but the
fluxSurfaces tracer requires the legacy versions.
"""

import os
import sys
import importlib.util
import warnings
import pprint

warnings.filterwarnings('ignore')

# --- Compat shims (no-ops if already on legacy scipy/numpy) ---
import numpy as np
if not hasattr(np, 'NaN'):
    np.NaN = np.nan
if not hasattr(np, 'RankWarning'):
    try:
        np.RankWarning = np.exceptions.RankWarning
    except AttributeError:
        pass

import scipy.integrate
if not hasattr(scipy.integrate, 'cumtrapz'):
    scipy.integrate.cumtrapz = scipy.integrate.cumulative_trapezoid

HERE = os.path.abspath(os.path.dirname(__file__))
GFILE = os.path.join(HERE, 'ITER_test.eqdsk')
PFILE = os.path.join(HERE, 'D3Dlike_Hmode_test.peqdsk')


def gen_gfile_expected():
    from omfit_classes.omfit_eqdsk import OMFITgeqdsk
    g = OMFITgeqdsk(GFILE)
    g.load()
    fs = g['fluxSurfaces']
    avg = fs['avg']
    geo = fs['geo']
    info = fs['info']
    mid = fs['midplane']
    li_info = info['internal_inductance']

    out = {}
    # --- Raw scalars (from OMFIT raw parser) ---
    for key in ['NW','NH','RDIM','ZDIM','RCENTR','RLEFT','ZMID',
                'RMAXIS','ZMAXIS','SIMAG','SIBRY','BCENTR','CURRENT']:
        out[key] = float(g[key])

    # --- 1-D profiles (raw) ---
    for key in ['FPOL','PRES','FFPRIM','PPRIME','QPSI']:
        arr = np.asarray(g[key])
        out[f'{key}_0']   = float(arr[0])
        out[f'{key}_mid'] = float(arr[len(arr)//2])
        out[f'{key}_end'] = float(arr[-1])
        out[f'{key}_sum'] = float(np.sum(arr))
        out[f'{key}_len'] = len(arr)

    # --- 2-D PSIRZ ---
    psi = np.asarray(g['PSIRZ'])
    out['PSIRZ_min'] = float(np.min(psi))
    out['PSIRZ_max'] = float(np.max(psi))
    out['PSIRZ_sum'] = float(np.sum(psi))

    # --- Boundary and limiter ---
    for key in ['RBBBS','ZBBBS','RLIM','ZLIM']:
        arr = np.asarray(g[key])
        out[f'{key}_len'] = len(arr)
        if len(arr) > 0:
            out[f'{key}_min'] = float(np.min(arr))
            out[f'{key}_max'] = float(np.max(arr))

    # --- Derived FSA quantities (from OMFIT fluxSurfaces) ---
    # Note: OMFIT's geo dict uses key 'kap' for elongation (we map to 'kappa_edge')
    out['Ip_enclosed_edge'] = float(avg['ip'][-1])
    out['q_axis']           = float(avg['q'][0])
    out['q_edge']           = float(avg['q'][-1])
    out['kappa_edge']       = float(geo['kap'][-1])
    out['delta_edge']       = float(geo['delta'][-1])
    out['a_edge']           = float(geo['a'][-1])
    out['R_geo_edge']       = float(geo['R'][-1])
    out['vol_edge']         = float(geo['vol'][-1])
    out['cxArea_edge']      = float(geo['cxArea'][-1])
    out['li_1']             = float(li_info['li_(1)_EFIT'])
    out['li_3']             = float(li_info['li_(3)_IMAS'])
    out['R_mid_edge']       = float(mid['R'][-1])
    bp_edge  = float(mid['Bp'][-1])
    bt_edge  = float(mid['Bt'][-1])
    out['Btot_mid_edge']    = float(np.sqrt(bp_edge**2 + bt_edge**2))
    out['beta_t']           = float(avg['beta_t'][-1])
    out['beta_p']           = float(avg['beta_p'][-1])
    out['beta_n']           = float(avg['beta_n'][-1])
    out['cocos']            = 1

    # --- j_tor and q profile samples at interior psi_N positions ---
    # We record both j_tor variants: standard <Jt/R>/<1/R> and the "direct"
    # GS form <Jt>_GS = -sigma_Bp * (p'*<R> + FF'*<1/R>/mu0) * (2pi)^exp_Bp
    # (the latter is the TokaMaker default).
    import scipy.constants as _sc
    om_psin      = np.asarray(geo['psin'])
    om_Jt_over_R = np.asarray(avg['Jt/R'])
    om_invR      = np.asarray(avg['1/R'])
    om_R         = np.asarray(avg['R'])
    om_pprime    = np.asarray(avg['PPRIME'])
    om_ffprim    = np.asarray(avg['FFPRIM'])
    om_q         = np.asarray(avg['q'])
    om_jt_std    = om_Jt_over_R / om_invR
    # COCOS 1: sigma_Bp = +1, exp_Bp = 0 -> leading coefficient -1
    om_jt_direct = -(om_pprime * om_R + om_ffprim * om_invR / _sc.mu_0)
    for pct in (10, 25, 50, 75, 90):
        psin_target = pct / 100.0
        out[f'j_tor_psiN_{pct:02d}']        = float(np.interp(psin_target, om_psin, om_jt_std))
        out[f'j_tor_direct_psiN_{pct:02d}'] = float(np.interp(psin_target, om_psin, om_jt_direct))
        out[f'q_psiN_{pct:02d}']            = float(np.interp(psin_target, om_psin, om_q))
    return out


def gen_pfile_expected():
    from omfit_classes.omfit_osborne import OMFITpFile
    pf = OMFITpFile(PFILE)
    pf.load()

    out = {}
    for key in list(pf.keys()):
        v = pf[key]
        if key == 'N Z A':
            out['N Z A_N_0'] = float(v['N'][0])
            out['N Z A_Z_0'] = float(v['Z'][0])
            out['N Z A_A_0'] = float(v['A'][0])
            out['N Z A_len'] = len(v['N'])
            continue
        if hasattr(v, 'keys') and 'data' in list(v.keys()):
            d = np.asarray(v['data'])
            psin = np.asarray(v['psinorm'])
            out[f'{key}_len'] = len(d)
            out[f'{key}_0']   = float(d[0])
            out[f'{key}_mid'] = float(d[len(d)//2])
            out[f'{key}_end'] = float(d[-1])
            out[f'{key}_min'] = float(np.min(d))
            out[f'{key}_max'] = float(np.max(d))
            out[f'{key}_sum'] = float(np.sum(d))
            out[f'{key}_psinorm_end'] = float(psin[-1])
    return out


if __name__ == '__main__':
    print('# ----- eqdsk_gfile_expected_dict (paste into test_TokaMaker.py) -----')
    print('eqdsk_gfile_expected_dict = \\')
    pprint.pprint(gen_gfile_expected(), width=100)
    print()
    print('# ----- eqdsk_pfile_expected_dict (paste into test_TokaMaker.py) -----')
    print('eqdsk_pfile_expected_dict = \\')
    pprint.pprint(gen_pfile_expected(), width=100)
