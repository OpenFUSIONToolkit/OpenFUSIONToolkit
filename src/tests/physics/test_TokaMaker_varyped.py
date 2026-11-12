#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Regression tests for the TokaMaker VARYPED equilibrium scan.'''
import os
import shutil
import sys
import time
import multiprocessing

import numpy as np

test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..', 'python')))
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


def mp_run(target, args, timeout=30):
    '''! Run a test case in a separate process (one OFT_env per process)

    @param target Function to run, which puts its result on the queue
    @param args Arguments to `target`, before the queue
    @param timeout Time limit [s]
    @result Result put on the queue by `target`, or `None` on timeout or error
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
    else: # Reached timeout
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


def run_D3Dlike_scan_point(mp_q):
    '''! Scan one VARYPED point of the D3Dlike example and compare the written files

    The input p-file is extended past the separatrix (psi_N = 1.02, 1.05), as
    p-files with scrape-off-layer points are, to check that the scan drops them.

    @param mp_q Queue for the result
    '''
    g_eq = eqdsk.read_geqdsk('TkMkr_D3Dlike_Hmode_BSamp=1.0.geqdsk')
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

    out_root = os.path.join(test_dir, 'varyped_test_output')
    shutil.rmtree(out_root, ignore_errors=True)
    os.makedirs(out_root)
    try:
        results = varyped.equilibrium_scan(mygs, p_eq, g_eq, [1.5], path_to_output=out_root,
                                           file_output=True, result_output=False,
                                           pin_edge_pprime=True)
    except Exception as err:
        print(err)
        mp_q.put(None)
        return
    point = results[1]
    if point['status'] == 'failed':
        mp_q.put({'status': point['status']})
        return
    out_dir = results[0]['out_dir']
    g_out = eqdsk.read_geqdsk(os.path.join(out_dir, point['gfile_name']))
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
    })


def test_D3Dlike_scan_point():
    '''! Written g-file and p-file pressures agree with each other and the separatrix pressure'''
    result = mp_run(run_D3Dlike_scan_point, (), timeout=600)
    assert result is not None
    assert result['status'] in ('converged', 'partial')
    # p-file points outside the separatrix are dropped
    assert result['p_psi_max'] == 1.0
    # Untruncated g-file: its boundary pressure is the separatrix pressure
    assert abs(result['g_pres_sep']/result['psep'] - 1.0) < 1.E-3
    # Output p-file ptot is the scaled input ptot, as solved
    assert abs(result['p_ptot_sep']/result['psep'] - 1.0) < 1.E-3
    assert result['core_dev'] < 1.E-2
