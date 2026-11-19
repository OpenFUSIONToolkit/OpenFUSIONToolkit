#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Regression tests for the TokaMaker VARYPED equilibrium scan.'''
import os

import numpy as np
import pytest

from tokamaker_testing import mp_run
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker, eqdsk, meshing, varyped


# Coil current targets of the D3Dlike VARYPED example (F coils in A-turns)
D3DLIKE_TARGETS = {
    'ECOILA': -0.977888676757812/61.0, 'ECOILB': -0.962711173828125/61.0,
    'F1A': 0.115971984375, 'F2A': 0.05980789453125, 'F3A': -0.03076001171875,
    'F4A': -0.0401077421875, 'F5A': 0.000723009033203125, 'F6A': -0.1178582568125,
    'F7A': 0.0341264296875, 'F8A': -0.05660116015625, 'F9A': 0.236625375,
    'F1B': 0.128368578125, 'F2B': 0.0763701328125, 'F3B': -0.02234583203125,
    'F4B': -0.13314096875, 'F5B': 0.000399709045410156, 'F6B': -0.15356990625,
    'F7B': 0.0415082109375, 'F8B': -0.05138975390625, 'F9B': 0.252380265625,
}


def run_D3Dlike_scan_point(out_root, reverse, output_cocos, mp_q):
    '''Scan one D3Dlike VARYPED point and compare the written files.'''
    g_eq = eqdsk.read_geqdsk('TkMkr_D3Dlike_Hmode_BSamp=1.0.geqdsk')
    if reverse:
        g_eq = g_eq.flip_Bt_Ip(copy=True)
    expected_g = g_eq.cocosify(output_cocos, copy=True)
    p_eq = eqdsk.read_pfile('D3Dlike_Hmode_test.peqdsk')
    psi_ptot = np.asarray(p_eq.psinorm_for('ptot'), dtype=float)
    p_eq = p_eq.remap(np.append(psi_ptot, [1.02, 1.05]))

    myOFT = OFT_env(nthreads=2)
    mygs = TokaMaker(myOFT)
    mygs.settings.maxits = 200
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = meshing.load_gs_mesh('DIIID_mesh.h5')
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.setup(order=2, F0=abs(g_eq.R_center*g_eq.B_center))
    mygs.set_coil_vsc({'F9A': 1.0, 'F9B': -1.0})
    isoflux_pts, weights = varyped.gfile_boundary_constraints(g_eq)
    mygs.set_isoflux_constraints(isoflux_pts, weights)
    reg_terms = []
    for name, target in D3DLIKE_TARGETS.items():
        target *= 1.E6
        if name.startswith('F'):
            target /= coil_dict[name]['nturns']
            reg_terms.append(mygs.coil_reg_term({name: 1.0}, target=target, weight=coil_dict[name]['nturns']))
        else:
            reg_terms.append(mygs.coil_reg_term({name: 1.0}, target=target, weight=1.0))
    reg_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E-2))
    mygs.set_coil_reg(reg_terms=reg_terms)

    results = varyped.equilibrium_scan(mygs, p_eq, g_eq, [1.5], path_to_output=out_root,
                                       file_output=True, result_output=False,
                                       output_cocos=output_cocos)
    point = results[1]
    if point['status'] == 'failed':
        mp_q.put({'status': point['status']})
        return
    out_dir = results[0]['out_dir']
    g_out = eqdsk.read_geqdsk(os.path.join(out_dir, point['gfile_name']), cocos=results[0]['cocos'])
    p_out = eqdsk.read_pfile(os.path.join(out_dir, point['pfile_name']))
    g_pres = np.asarray(g_out.pres, dtype=float)
    g_psi = np.asarray(g_out.psi_N, dtype=float)
    p_ptot = np.asarray(p_out['ptot']['data'], dtype=float)*1.E3
    p_psi = np.asarray(p_out.psinorm_for('ptot'), dtype=float)
    core = g_psi <= 0.9
    mp_q.put({
        'status': point['status'],
        'psep': point['psep'],
        'g_pres_sep': g_pres[-1],
        'p_ptot_sep': p_ptot[-1],
        'p_psi_max': p_psi[-1],
        'core_dev': np.max(np.abs(np.interp(g_psi[core], p_psi, p_ptot)/g_pres[core] - 1.0)),
        'expected_ip': expected_g.Ip,
        'expected_bt': expected_g.B_center,
        'ip_target': point['Ip_target'],
        'ip_out': g_out.Ip,
        'bt_out': g_out.B_center,
        'jtor': g_out.j_tor_averaged_direct[20],
        'q_file': g_out.qpsi[20],
        'q_trace': g_out.q_profile[20],
    })


@pytest.mark.parametrize('reverse,output_cocos', [(False, 2), (True, 7)])
def test_D3Dlike_scan_point(tmp_path, reverse, output_cocos):
    '''Written pressures agree and g-file field and current directions are consistent.'''
    result = mp_run(run_D3Dlike_scan_point, (str(tmp_path), reverse, output_cocos), timeout=600)
    assert result is not None
    assert result['status'] in ('converged', 'partial')
    # p-file points outside the separatrix are dropped
    assert result['p_psi_max'] == 1.0
    # Untruncated g-file: its boundary pressure is the separatrix pressure
    assert abs(result['g_pres_sep']/result['psep'] - 1.0) < 1.E-3
    # Output p-file ptot is the scaled input ptot, as solved
    assert abs(result['p_ptot_sep']/result['psep'] - 1.0) < 1.E-3
    assert result['core_dev'] < 1.E-2
    assert np.isclose(result['ip_target'], abs(result['expected_ip']))
    assert np.sign(result['ip_out']) == np.sign(result['expected_ip']) == np.sign(result['jtor'])
    assert np.sign(result['bt_out']) == np.sign(result['expected_bt'])
    assert np.sign(result['q_file']) == np.sign(result['q_trace']) != 0
