'''Why solve with profiles on normalised toroidal flux (Phi_N)?

During experimental reconstruction, quantities such as pressure are often reconstructed on the
psi_N surfaces of a given equilibrium (e.g. EFIT02) via tomographic inversion.

These psi_N surfaces are sensitive to the shape of the current profile, which is difficult to
reconstruct. If we re-solve the equilibrium with an improved current-profile model, the psi_N
surfaces can shift relative to the lab frame. Quantities mapped on psi_N will then deviate from
their original experimental reconstructions.

Normalised toroidal flux (Phi_N) is much less sensitive to changes in the current profile than
psi_N. Tokamaks employ a large toroidal field supplied by external coils, so the
toroidal flux inside a surface is set mainly by the area it encloses. The poloidal flux inside a
surface, by contrast, is produced by the plasma current inside it. The two are linked by
dPhi = q dpsi, so a new current profile changes q and the psi_N -> Phi_N map, while Phi_N stays
tied to the geometry. If we pin flux-surface quantities to Phi_N during the re-solve, they deviate
less in real space than quantities pinned to psi_N. They still move somewhat, because the shapes
of the inner flux surfaces also depend on the current profile.

This script generates pressure_anomaly_rz.png to illustrate this point. A reference equilibrium
with an EFIT02-like current profile is solved with both current and pressure profiles mapped on
psi_N. It is then re-solved with the same pressure profile in three ways (panels of
pressure_anomaly_rz.png, each showing the change in pressure relative to the reference):
  (a) the same EFIT02-like current profile, both profiles mapped on Phi_N. This control should
      reproduce the reference; its residual is the numerical floor (interpolation, Phi_N map).
  (b) a more physics-informed, sawtoothing H-mode current profile (see profiles_compare.png),
      both profiles mapped on psi_N: significant deviation.
  (c) the same H-mode current profile, both profiles mapped on Phi_N: less deviation.

Not part of the automated tests. Usage (run from this directory, ~1 min):
    python tokamaker_torflux_motivation.py [--out DIR]
'''
import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(HERE, '..', '..', 'python')))
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk

# ============================================================ inputs
GEQDSK_FILE = os.path.join(HERE, 'D3Dlike_Hmode_baseline.geqdsk')  # LCFS, F0 and pressure
MESH_FILE = os.path.join(HERE, '..', '..', 'examples', 'TokaMaker', 'DIIID', 'DIIID_mesh.h5')
IP_TARGET = 1200000.0  # [A]
EDGE_RAMP = 0.05       # FF' and p' ramped linearly to zero over the outer 5% of psi_N

# Normalised j_phi = FF'<1/R>/mu0 + p'<R> on PSI_NODES (psi_N, 0 at the axis)
PSI_NODES = np.linspace(0.0, 1.0, 51)
# EFIT02-like: DIII-D 174956 EFIT02 gEQDSK (t = 1.98 s)
JPHI_EFIT = np.array([
    1.0000, 0.9854, 0.9566, 0.9284, 0.9006, 0.8734, 0.8466, 0.8204, 0.7947, 0.7695, 0.7448, 0.7206,
    0.6969, 0.6737, 0.6510, 0.6289, 0.6072, 0.5861, 0.5654, 0.5453, 0.5256, 0.5065, 0.4878, 0.4697,
    0.4520, 0.4349, 0.4183, 0.4022, 0.3865, 0.3714, 0.3567, 0.3426, 0.3290, 0.3158, 0.3032, 0.2910,
    0.2794, 0.2683, 0.2576, 0.2474, 0.2377, 0.2285, 0.2198, 0.2115, 0.2037, 0.1963, 0.1891, 0.1823,
    0.1754, 0.1680, 0.1635])
# H-mode: D3D-like baseline gEQDSK, flattened inside psi_N = 0.2 (sawteeth) with a smoothstep to 0.3
JPHI_HMODE = np.array([
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.9970,
    0.9796, 0.9435, 0.8954, 0.8535, 0.8236, 0.7937, 0.7640, 0.7345, 0.7054, 0.6767, 0.6487, 0.6214,
    0.5950, 0.5695, 0.5450, 0.5217, 0.4996, 0.4787, 0.4591, 0.4408, 0.4237, 0.4079, 0.3932, 0.3795,
    0.3666, 0.3541, 0.3424, 0.3320, 0.3236, 0.3174, 0.3142, 0.3157, 0.3237, 0.3432, 0.3975, 0.5475,
    0.6625, 0.3594, 0.1410])

# Coil current targets [A-turns] for regularisation, as in examples/TokaMaker/DIIID
COIL_TARGETS = {
    'ECOILA': -16030.961914062493, 'ECOILB': -15782.150390625,
    'F1A': 1999.5169719827586, 'F2A': 1031.1705953663793, 'F3A': -530.3450296336207,
    'F4A': -691.5127963362069, 'F5A': 12.465672986260776, 'F6A': -2142.877414772727,
    'F7A': 620.4805397727273, 'F8A': -975.8820716594828, 'F9A': 4302.279545454546,
    'F1B': 2213.2513469827586, 'F2B': 1316.7264278017242, 'F3B': -385.2729660560345,
    'F4B': -2295.533943965517, 'F5B': 6.891535265692344, 'F6B': -2792.180113636364,
    'F7B': 754.6947443181818, 'F8B': -886.030239762931, 'F9B': 4588.732102272727,
}


# ============================================================ TokaMaker
def coil_weight(name):
    if name.startswith('ECOIL'):
        return 61.0
    return 1.E2 if name.startswith('F5') else 1.E0


def solve(mygs, eqdsk, ffp_prof, pp_prof, name, out_dir):
    '''Free-boundary solve with the gEQDSK LCFS, Ip and pax targets; returns pressure at mesh nodes.'''
    mygs.reset()
    mygs.settings.maxits = 200
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(MESH_FILE)
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.setup(order=2)
    mygs.set_coil_vsc({'F9A': 1.0, 'F9B': -1.0})
    reg_terms = [mygs.coil_reg_term({k: 1.0}, target=v, weight=coil_weight(k)) for k, v in COIL_TARGETS.items()]
    reg_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E-2))
    mygs.set_coil_reg(reg_terms=reg_terms)

    mygs.set_targets(Ip=IP_TARGET, pax=eqdsk['pres'][0])
    mygs.set_isoflux_constraints(eqdsk['rzout'].copy())
    mygs.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof, foffset=eqdsk['rcentr']*eqdsk['bcentr'])
    mygs.init_psi(1.7, 0.0, 0.5, 1.3, 0.4)  # R0, Z0, a, kappa, delta
    mygs.solve()

    fig, ax = plt.subplots(1,1)
    mygs.plot_machine(fig,ax,coil_colormap='seismic',coil_symmap=True,coil_scale=1.E-6,coil_clabel=r'$I_C$ [MA]')
    mygs.plot_psi(fig,ax,vacuum_nlevels=4)
    ax.plot(eqdsk['rzout'][:,0], eqdsk['rzout'][:,1], 'r--')
    ax.set_title(name)
    fig.savefig(os.path.join(out_dir, f'machine_{name}.png'), dpi=150)
    plt.close(fig)
    return mygs.get_field_eval('P').eval(mygs.r[:,:2])[:,0]


def phi_N_map(mygs):
    '''Phi_N(psi_N) = int_0^psi_N q / int_0^1 q of the current solve, as (psi_N, Phi_N) arrays.'''
    psi, q, _, _, _, _ = mygs.get_q(npsi=400, psi_pad=0.001)
    psi = np.r_[0.0, psi, 1.0]
    q = np.abs(np.r_[q[0], q, q[-1]])
    phi = np.r_[0.0, np.cumsum(0.5*(q[1:]+q[:-1])*np.diff(psi))]
    return psi, phi/phi[-1]


def psi_profiles(jphi, pp):
    return ({'type': 'jphi-linterp', 'x': PSI_NODES, 'y': jphi},
            {'type': 'linterp', 'x': PSI_NODES, 'y': pp})


def phi_profiles(jphi, pp, phi_nodes):
    '''jphi as usual, on the Phi_N image of PSI_NODES; pp is dp/dPhi_N.'''
    return ({'type': 'jphi-linterp', 'x': phi_nodes, 'y': jphi, 'coord': 'phi_n_relabel'},
            {'type': 'linterp', 'x': phi_nodes, 'y': pp, 'coord': 'phi_n'})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(HERE, 'torflux_motivation_out'))
    out_dir = ap.parse_args().out
    os.makedirs(out_dir, exist_ok=True)

    # ============================================================ profiles
    eqdsk = read_eqdsk(GEQDSK_FILE)
    psi_g = np.linspace(0.0, 1.0, eqdsk['nr'])  # gEQDSK psi_N grid
    ramp = np.clip((1.0-PSI_NODES)/EDGE_RAMP, 0.0, 1.0)
    jphi_efit = ramp*JPHI_EFIT
    jphi_hmode = ramp*JPHI_HMODE
    pp_psi = ramp*np.interp(PSI_NODES, psi_g, np.gradient(eqdsk['pres'], psi_g))  # dp/dpsi_N (scale set by pax)

    # ============================================================ solves
    mygs = TokaMaker(OFT_env(nthreads=2))
    p_rz = {'ref': solve(mygs, eqdsk, *psi_profiles(jphi_efit, pp_psi), 'ref', out_dir)}

    # Pin the reference profiles to Phi_N using the reference equilibrium's map
    psi_m, phi_m = phi_N_map(mygs)
    phi_nodes = np.interp(PSI_NODES, psi_m, phi_m)
    phi_g = np.interp(psi_g, psi_m, phi_m)
    pp_phi = ramp*np.interp(phi_nodes, phi_g, np.gradient(eqdsk['pres'], phi_g))  # dp/dPhi_N

    fig, ax = plt.subplots(1,1)
    ax.plot(phi_nodes, jphi_efit, label=r'$j_\phi$, EFIT02-like')
    ax.plot(phi_nodes, jphi_hmode, '--', label=r'$j_\phi$, H-mode, sawtoothing')
    ax.plot(phi_g, eqdsk['pres']/eqdsk['pres'][0], 'k:', label=r'$p/p_{axis}$')
    ax.set_xlabel(r'$\Phi_N$ (reference equilibrium)')
    ax.set_ylabel('normalised profile')
    ax.legend(fontsize=8)
    ax.grid(True)
    fig.savefig(os.path.join(out_dir, 'profiles_compare.png'), dpi=150)
    plt.close(fig)

    cases = [
        ('a', r'(a) EFIT02-like $j_\phi$, on $\Phi_N$', phi_profiles(jphi_efit, pp_phi, phi_nodes)),
        ('b', r'(b) H-mode $j_\phi$, on $\psi_N$', psi_profiles(jphi_hmode, pp_psi)),
        ('c', r'(c) H-mode $j_\phi$, on $\Phi_N$', phi_profiles(jphi_hmode, pp_phi, phi_nodes)),
    ]
    for k, _, (ffp_prof, pp_prof) in cases:
        p_rz[k] = solve(mygs, eqdsk, ffp_prof, pp_prof, k, out_dir)

    # ============================================================ pressure anomaly
    p_ref = p_rz['ref']
    dp = {k: (p_rz[k]-p_ref)/p_ref.max() for k, _, _ in cases}
    vmax = max(np.abs(v).max() for v in dp.values())

    # Area-weighted means over cells inside either plasma (cell value = mean of its nodes)
    x_c, y_c = mygs.r[mygs.lc,0], mygs.r[mygs.lc,1]
    area = 0.5*np.abs((x_c[:,1]-x_c[:,0])*(y_c[:,2]-y_c[:,0]) - (x_c[:,2]-x_c[:,0])*(y_c[:,1]-y_c[:,0]))
    stats = {}
    for k, _, _ in cases:
        m = np.any((p_ref[mygs.lc] > 0) | (p_rz[k][mygs.lc] > 0), axis=1)
        w = area[m]/area[m].sum()
        stats[k] = (np.abs(dp[k]).max(), np.sum(w*np.abs(dp[k][mygs.lc[m]]).mean(axis=1)), np.sum(w*dp[k][mygs.lc[m]].mean(axis=1)))
        print(f'({k}): max|dp| = {stats[k][0]:.3e}, mean|dp| = {stats[k][1]:.3e}, mean dp = {stats[k][2]:+.3e}  (normalised to max p_ref)')

    p_lev = np.linspace(0.1, 0.9, 9)*p_ref.max()
    rz_pl = mygs.r[np.unique(mygs.lc[mygs.reg == 1]), :2]  # plasma-region nodes, set the view
    tri = (mygs.r[:,0], mygs.r[:,1], mygs.lc)
    fig, axs = plt.subplots(1, 3, figsize=(13, 7), sharex=True, sharey=True)
    for ax, (k, label, _) in zip(axs, cases):
        cf = ax.tripcolor(*tri, dp[k], shading='gouraud', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        ax.tricontour(*tri, p_ref, levels=p_lev, colors='k', linewidths=0.8)
        ax.tricontour(*tri, p_rz[k], levels=p_lev, colors='g', linewidths=0.8, linestyles='--')
        ax.set_title(f'{label}\nmax|dp| = {stats[k][0]:.2e}\nmean|dp| = {stats[k][1]:.2e}\nmean dp = {stats[k][2]:+.2e}', fontsize=10)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('R [m]')
    axs[0].set_xlim(rz_pl[:,0].min()-0.05, rz_pl[:,0].max()+0.05)
    axs[0].set_ylim(rz_pl[:,1].min()-0.05, rz_pl[:,1].max()+0.05)
    axs[0].set_ylabel('Z [m]')
    axs[0].plot([], [], 'k-', label='p (reference)')
    axs[0].plot([], [], 'g--', label='p (re-solved)')
    axs[0].legend(loc='lower left', fontsize=8)
    fig.colorbar(cf, ax=axs, shrink=0.8, label=r'$(p - p_{ref})/\max(p_{ref})$')
    fig.savefig(os.path.join(out_dir, 'pressure_anomaly_rz.png'), dpi=150, bbox_inches='tight')
    np.savez(os.path.join(out_dir, 'pressure_rz.npz'), r=mygs.r[:,:2], lc=mygs.lc, **p_rz)


if __name__ == '__main__':
    main()
