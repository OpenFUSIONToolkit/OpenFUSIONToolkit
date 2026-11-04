#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: LTX free-boundary equilibrium, wall modes, and linear stability
import os
import json
import pytest
import numpy as np
from tokamaker_testing import mp_run, validate_dict
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain, save_gs_mesh, load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import create_isoflux, create_power_flux_fun


def run_LTX_case(fe_order,test_type,mp_q):
    def create_mesh():
        with open('LTX_geom.json','r') as fid:
            LTX_geom = json.load(fid)
        plasma_dx = 0.02
        coil_dx = 0.02
        vv_dx = 0.015
        vac_dx = 0.05
        gs_mesh = gs_Domain()
        #
        gs_mesh.define_region('air',vac_dx,'boundary')
        gs_mesh.define_region('plasma',plasma_dx,'plasma')
        gs_mesh.define_region('shellU',vv_dx,'conductor',eta=4.E-7,noncontinuous=True)
        gs_mesh.define_region('shellL',vv_dx,'conductor',eta=4.E-7,noncontinuous=True)
        for i, vv_segment in enumerate(LTX_geom['vv']):
            gs_mesh.define_region('vv{0}'.format(i),vv_dx,'conductor',eta=vv_segment[1])
        for key, coil in LTX_geom['coils'].items():
            if key.startswith('OH'):
                gs_mesh.define_region(key,coil_dx,'coil',nTurns=coil['nturns'],coil_set='OH')
            else:
                gs_mesh.define_region(key,coil_dx,'coil',nTurns=coil['nturns'])
        #
        gs_mesh.add_polygon(LTX_geom['limiter'],'plasma',parent_name='air')
        gs_mesh.add_polygon(LTX_geom['shell'],'shellU',parent_name='air')
        shell_lower = np.array(LTX_geom['shell'].copy()); shell_lower[:,1] *= -1.0
        gs_mesh.add_polygon(shell_lower,'shellL',parent_name='air')
        for i, vv_segment in enumerate(LTX_geom['vv']):
            gs_mesh.add_polygon(vv_segment[0],'vv{0}'.format(i),parent_name='air')
        for key, coil in LTX_geom['coils'].items():
            gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name='air')
        #
        mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
        coil_dict = gs_mesh.get_coils()
        cond_dict = gs_mesh.get_conductors()
        save_gs_mesh(mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict,'LTX_mesh.h5')
    if not os.path.exists('LTX_mesh.h5'):
        try:
            create_mesh()
        except Exception as e:
            print(e)
            mp_q.put(None)
            return
    # Run EQ
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('LTX_mesh.h5')
    mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    mygs.setup(order=fe_order,F0=0.10752)
    #
    if test_type == 'eig':
        eig_vals, _ = mygs.compute_wall_modes(10)
        mp_q.put([{'Tau_w': eig_vals[:5]}])
        oftpy_dump_cov()
        return
    #
    mygs.set_coil_vsc({'INTERNALU': 1.0, 'INTERNALL': -1.0})
    #
    Ip_target = 8.0E4
    mygs.set_targets(Ip=Ip_target,Ip_ratio=2.0)
    isoflux_pts = create_isoflux(20,0.40,0.0,0.22,1.5,0.1)
    mygs.set_isoflux(isoflux_pts)
    # Set regularization weights
    disable_list = ('YELLOW',)
    regularization_terms = []
    for name in mygs.coil_sets:
        if name[:-1] in disable_list:
            regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E4))
            continue
        if name == 'OH': # OH coil has no mirror
            regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E-1))
            continue
        elif name[-1] == 'L':
            continue
        regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E-1))
        regularization_terms.append(mygs.coil_reg_term({name: 1.0, name[:-1]+'L': -1.0},target=0.0,weight=1.E2))
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E-4))
    mygs.set_coil_reg(reg_terms=regularization_terms)
    #
    ffp_prof = create_power_flux_fun(50,1.5,2.0)
    pp_prof = create_power_flux_fun(50,4.0,1.0)
    mygs.set_profiles(ffp_prof=ffp_prof,pp_prof=pp_prof)
    #
    mygs.init_psi(0.42,0.0,0.15,1.5,0.6)
    mygs.settings.pm=True
    mygs.update_settings()
    mygs.solve()
    if test_type == 'stab':
        eig_vals, _ = mygs.compute_linear_stability(1.E3,10,False)
        mp_q.put([{'gamma': eig_vals[:5]}])
        oftpy_dump_cov()
        return
    #
    psi_last = mygs.get_psi(False)
    mygs.set_psi_dt(psi_last,5.E-3)
    Ip_target = 9.0E4
    mygs.set_targets(Ip=Ip_target,Ip_ratio=2.0)
    mygs.solve()
    mygs.save_eqdsk('test.eqdsk')
    #
    mp_q.put([mygs.get_stats()])
    oftpy_dump_cov()

# Test runners for LTX test cases
@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_LTX_eig(order):
    exp_dict = {
        'Tau_w': [5.152566E-03, 3.953030E-03, 2.536384E-03, 2.172948E-03, 1.853882E-03]
    }
    results = mp_run(run_LTX_case,(order,'eig'))
    assert validate_dict(results,exp_dict)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_LTX_stability(order):
    exp_dict = {
        'gamma': [234.1051, -214.4196, -282.0877, -388.7592, -388.7592]
    }
    results = mp_run(run_LTX_case,(order,'stab'))
    assert validate_dict(results,exp_dict)

LTX_eq_dict = {
    'Ip': 90002.51679781199,
    'Ip_centroid': [ 4.05458767e-01, None],
    'kappa': 1.525595596236063,
    'kappaU': 1.5256161060199729,
    'kappaL': 1.5255750864521527,
    # 'delta': 0.1274386709874723,
    # 'deltaU': 0.13292909306765727,
    # 'deltaL': 0.1219482489072871,
    'R_geo': 0.39198831687969443,
    'a_geo': 0.2378387910900877,
    'vol': 0.6507554668762836,
    'q_0': 1.3280540982807334,
    'q_95': 5.901997513881755,
    'P_ax': 1720.958666106632,
    'W_MHD': 563.2852958452944,
    'beta_pol': 41.4052788464515,
    'dflux': 0.0009601908294886685,
    'tflux': 0.08551496495989133,
    'l_i': 1.0271521431711803,
    'beta_tor': 1.9276444168027145,
    'beta_n': 1.3972402635015146
}

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_LTX_eq(order):
    results = mp_run(run_LTX_case,(order,''))
    assert validate_dict(results,LTX_eq_dict)
