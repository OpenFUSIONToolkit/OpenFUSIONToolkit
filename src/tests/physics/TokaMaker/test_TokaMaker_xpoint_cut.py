#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: q near the LCFS with a second X-point close to the active one (DIII-D)
import pytest
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from tokamaker_testing import mp_run
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh


def run_DIIID_xpoint_pair_case(fe_order,mp_q):
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('DIIID_mesh.h5')
    mygs.setup_model(mesh_pts, mesh_lc, reg=mesh_reg, cond_dict=cond_dict, coil_dict=coil_dict, order=fe_order, F0=2.96)
    mygs.set_coil_vsc({'F9A': 1.0, 'F9B': -1.0})
    mygs.set_targets(Ip=1.19E6,pax=6.0E4)
    # Two lower X-points 9 cm apart, which skews the core axis of the active one by ~45 degrees from O - X
    x_points = np.array([[1.183,-1.227],[1.11,-1.17]])
    lcfs = np.array([[2.287,0.067],[2.275,0.205],[2.241,0.343],[2.187,0.487],[2.101,0.636],[1.972,0.78],
                     [1.793,0.908],[1.568,0.968],[1.346,0.904],[1.195,0.725],[1.128,0.512],[1.094,0.323],
                     [1.077,0.151],[1.076,-0.017],[1.083,-0.193],[1.09,-0.405],[1.117,-0.686],[1.542,-0.97],
                     [1.792,-0.768],[1.959,-0.607],[2.077,-0.465],[2.164,-0.331],[2.227,-0.201],[2.271,-0.069]])
    mygs.set_isoflux_constraints(np.vstack((lcfs,x_points)))
    mygs.set_saddle_constraints(x_points)
    reg_terms = [mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E-2) for name in mygs.coil_sets]
    reg_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E-2))
    mygs.set_coil_reg(reg_terms=reg_terms)
    psi_prof = np.array([0.0,0.1,0.3,0.5,0.7,0.85,0.925,0.95,0.975,1.0])
    ffp_prof = np.array([1.0,0.571,0.355,0.233,0.155,0.112,0.071,-0.002,-0.205,0.0])
    pp_prof = np.array([1.0,1.868,1.758,1.235,0.881,0.75,0.91,1.642,3.084,0.0])
    mygs.set_profiles(ffp_prof={'type': 'linterp', 'x': psi_prof, 'y': ffp_prof},
                      pp_prof={'type': 'linterp', 'x': psi_prof, 'y': pp_prof})
    mygs.init_psi(1.7,0.0,0.6,1.7,0.4)
    mygs.solve()
    # q near the LCFS from get_q and from a field line followed once around poloidally from the outboard midplane
    psi_q = np.array([0.99,0.9999])
    _,qvals,_,_,_,_ = mygs.get_q(psi=psi_q)
    B_field = mygs.get_field_eval('B')
    psi_field = mygs.get_field_eval('PSI')
    R0, Z0 = mygs.o_point
    def field_line(s,y): # (R, Z, toroidal angle, poloidal angle) vs poloidal arc length
        B = B_field.eval(y[:2])
        Bp = np.hypot(B[0],B[2])
        dR, dZ = B[0]/Bp, B[2]/Bp
        return [dR, dZ, B[1]/(y[0]*Bp), ((y[0]-R0)*dZ-(y[1]-Z0)*dR)/((y[0]-R0)**2+(y[1]-Z0)**2)]
    one_turn = lambda s,y: abs(y[3])-2.0*np.pi
    one_turn.terminal = True
    q_trace = []
    for psi in psi_q:
        psi_surf = mygs.psi_bounds[1]+psi*(mygs.psi_bounds[0]-mygs.psi_bounds[1])
        R_start = brentq(lambda R: psi_field.eval(np.array([R,Z0]))[0]-psi_surf,R0,2.36,xtol=1.E-12) # 2.36 m: inside the outboard limiter
        sol = solve_ivp(field_line,[0.0,50.0],[R_start,Z0,0.0,0.0],events=one_turn,rtol=1.E-10,atol=1.E-12)
        q_trace.append(abs(sol.y_events[0][0][2])/(2.0*np.pi))
    mp_q.put([{'diverted': mygs.diverted, 'nx_points': len(mygs.get_xpoints()[0]), 'q': qvals, 'q_trace': np.array(q_trace)}])
    oftpy_dump_cov()


@pytest.mark.coverage
def test_DIIID_xpoint_pair():
    results = mp_run(run_DIIID_xpoint_pair_case,(2,))
    assert results is not None
    results = results[0]
    assert results['diverted'] and (results['nx_points'] > 1)
    # 1% as q diverges at the LCFS (main: +6.5% at psi_N = 0.9999)
    assert np.allclose(results['q'],results['q_trace'],rtol=1.E-2)
