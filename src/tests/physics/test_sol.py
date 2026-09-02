import pytest

from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.reconstruction import reconstruction
from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk, read_mhdin, read_kfile

import numpy as np

@pytest.mark.coverage
def test_current_consistent():
    myOFT = OFT_env(nthreads=2)
    mygs = TokaMaker(myOFT)

    mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('./src/examples/TokaMaker/DIIID/DIIID_mesh.h5')
    eqdsk = read_eqdsk('./src/examples/TokaMaker/DIIID/g192185.02440')
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    mygs.setup(order = 2, F0 = eqdsk['rcentr']*eqdsk['bcentr'])

    mygs.set_coil_vsc({'F9A': 1.0, 'F9B': -1.0})
    mygs.set_targets(Ip = eqdsk['ip'], pax=eqdsk['pres'][0])

    isoflux_pts = eqdsk['rzout'].copy()
    mygs.set_isoflux_constraints(isoflux_pts)

    target_currents = {
        'ECOILA': -16030.961914062493,
        'ECOILB': -15782.150390625,
        'F1A': 1999.5169719827586,
        'F2A': 1031.1705953663793,
        'F3A': -530.3450296336207,
        'F4A': -691.5127963362069,
        'F5A': 12.465672986260776,
        'F6A': -2142.877414772727,
        'F7A': 620.4805397727273,
        'F8A': -975.8820716594828,
        'F9A': 4302.279545454546,
        'F1B': 2213.2513469827586,
        'F2B': 1316.7264278017242,
        'F3B': -385.2729660560345,
        'F4B': -2295.533943965517,
        'F5B': 6.891535265692344,
        'F6B': -2792.180113636364,
        'F7B': 754.6947443181818,
        'F8B': -886.030239762931,
        'F9B': 4588.732102272727,
    }

    # Set regularization weights
    regularization_terms = []
    for name, target_current in target_currents.items():
        # Set specific target currents from input equilibrium and different weights depending on the coil set
        if name.startswith('ECOIL'):
            regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=target_current,weight=61.0))
        elif name.startswith('F'):
            if name.startswith('F5'):
                regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=target_current,weight=1.E2))
            else:
                regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=target_current,weight=1.E0))
    # Set zero target current and small weight on virtual VSC to allow up-down adjustment
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E-2))

    # Pass regularization terms to TokaMaker
    mygs.set_coil_reg(reg_terms=regularization_terms)

    # === TEST SOL ===
    ffprim = eqdsk['ffprim']
    pprime = eqdsk['pprime']
    ffprim -= ffprim[-1]
    ffprim = ffprim / ffprim[0]
    pprime -= pprime[-1]
    pprime = pprime / pprime[0]

    x = np.linspace(0.0, 1.0, len(ffprim))
    x_sol = np.linspace(1.01, 1.3, 25)

    x = np.linspace(0.0, 1.0, len(ffprim))
    x_sol = np.linspace(1.01, 1.3, 25)
    x = np.append(x, x_sol)
    pprime = np.append(pprime, len(x_sol) * [0.0])
    ffprim = np.append(ffprim, len(x_sol) * [0.0]) 

    def supergaussian(psi_n, centr, std=0.05):
        return np.exp(-(psi_n - centr) ** 4 / (2 * std**4))

    ffp_bump = 0.1 * supergaussian(x, 1.15, std=0.12)
    ffprim += ffp_bump
    pp_bump = 0.1 * supergaussian(x, 1.15, std=0.12)
    pprime += pp_bump
    mygs.set_profiles(ffp_prof={'type': 'linterp', 'y': ffprim, 'x': x},pp_prof={'type': 'linterp', 'y': pprime, 'x': x}, f_SOL=True)

    R0 = 1.7
    Z0 = 0.0
    a = 0.5
    kappa = 1.3
    delta = -0.4

    mygs.init_psi(R0, Z0, a, kappa, delta)
    _ = mygs.solve()

    eq_stats = mygs.get_stats()
    stats_Ip = eq_stats['Ip']

    j_dens = mygs.calc_jtor_plasma()
    jtot = 0
    for cell in range(mygs.nc):
        if mygs.reg[cell] not in [1, 3]:
            continue
        idx1, idx2, idx3 = mygs.lc[cell]
        v1 = mygs.r[idx1]
        v2 = mygs.r[idx2]
        v3 = mygs.r[idx3]

        area = 0.5 * np.linalg.norm(np.cross(v2-v1, v3-v1))
        jtot += (j_dens[idx1] + j_dens[idx2] + j_dens[idx3]) * area / 3.0

    err = np.abs((jtot - stats_Ip) / stats_Ip)
    assert err < 0.001