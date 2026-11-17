"""Figures for the breakdown / start-up example (machine, traces, breakdown maps, equilibria)."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.path import Path

from startup_model import eval_fields, cell_geometry, MU0


class MeshGeometry:
    """Triangulations, limiter and coil/segment cells of a MachineModel, derived once."""

    def __init__(self, model):
        gs = model.gs
        self.R, self.Z = gs.r[:, 0].copy(), gs.r[:, 1].copy()
        self.lc = np.asarray(gs.lc)
        self.lim = max(gs.lim_contours, key=len)[:, :2]
        cen = np.column_stack([self.R[self.lc].mean(axis=1), self.Z[self.lc].mean(axis=1)])
        self.tri = mtri.Triangulation(self.R, self.Z, self.lc)
        self.tri_in = mtri.Triangulation(self.R, self.Z, self.lc, mask=~Path(self.lim).contains_points(cen))
        P, LC, REG = model.mesh
        self.P, self.LC = P[:, :2], LC
        area, cc = cell_geometry(P, LC)
        self.coil_cells = [np.where(REG == model.coil_dict[n]['reg_id'])[0] for n in model.coil_names]
        self.coil_centre = np.array([cc[c].mean(axis=0) for c in self.coil_cells])
        self.seg_cells = [s['cells'] for s in model.segments]
        self.seg_area = np.array([area[c].sum() for c in self.seg_cells])
        self.seg_region = np.array([s['region'] for s in model.segments])
        self.seg_of_cell = np.concatenate([[j] * len(c) for j, c in enumerate(self.seg_cells)]).astype(int)
        self.coil_of_cell = np.concatenate([[i] * len(c) for i, c in enumerate(self.coil_cells)]).astype(int)
        self.seg_polys = self.P[LC[np.concatenate(self.seg_cells)]]
        self.coil_polys = self.P[LC[np.concatenate(self.coil_cells)]]
        allp = np.vstack([self.seg_polys.reshape(-1, 2), self.coil_polys.reshape(-1, 2)])
        lo, hi = allp.min(axis=0), allp.max(axis=0)
        pad = 0.05 * (hi - lo).max()
        self.extent = (max(lo[0] - pad, 0.0), hi[0] + pad, lo[1] - pad, hi[1] + pad)


def _structure(ax, g, I_seg=None, Jlim=None, used=None):
    """Draw the structure (coloured by current density if I_seg is given) and the coils."""
    if I_seg is None:
        regions = list(dict.fromkeys(g.seg_region))
        base = np.array([plt.get_cmap('tab10')(regions.index(r))[:3] for r in g.seg_region])
        shade = np.where(np.arange(len(g.seg_cells)) % 2, 0.75, 1.0)[:, None]
        fc = np.c_[(base * shade)[g.seg_of_cell], np.ones(len(g.seg_of_cell))]
        ax.add_collection(PolyCollection(g.seg_polys, facecolors=fc, edgecolors=fc, linewidths=0.3))
        h = None
    else:
        J = I_seg / g.seg_area / 1e6
        h = PolyCollection(g.seg_polys, array=J[g.seg_of_cell], cmap='PuOr_r', edgecolors='face', linewidths=0.2)
        h.set_clim(-Jlim, Jlim)
        ax.add_collection(h)
    used = np.ones(len(g.coil_cells), bool) if used is None else used
    cc = np.array(['0.35' if used[i] else '0.8' for i in g.coil_of_cell])
    ax.add_collection(PolyCollection(g.coil_polys, facecolors=cc, edgecolors=cc, linewidths=0.3))
    ax.plot(g.lim[:, 0], g.lim[:, 1], 'k', lw=1.2)
    ax.set_aspect('equal'); ax.set_xlim(g.extent[:2]); ax.set_ylim(g.extent[2:])
    return h


def plot_machine(model, g, used, region, channel, shapes):
    """Coils (driven / open), structure segments by region, limiter, breakdown region, early
    channel and the ramp's target shapes."""
    fig, ax = plt.subplots(figsize=(9, 8))
    _structure(ax, g, used=used)
    for i, n in enumerate(model.coil_names):
        R, Z = g.coil_centre[i]
        ax.text(R, Z, n, fontsize=7, ha='center', va='center', color='w' if used[i] else 'k')
    R0, Z0, r = region
    ax.add_patch(plt.Circle((R0, Z0), r, fill=False, ls='--', ec='C3', lw=1.5, label='breakdown region'))
    ax.add_patch(plt.Circle(channel[:2], channel[2], fill=False, ec='C1', lw=1.5, label='early plasma channel'))
    for i, p in enumerate(shapes):
        ax.plot(np.r_[p[:, 0], p[0, 0]], np.r_[p[:, 1], p[0, 1]], color='C4', lw=0.7, alpha=0.7,
                label=None if i else 'ramp target shapes')
    regions = list(dict.fromkeys(g.seg_region))
    for rg in regions:
        ax.plot([], [], 's', color=plt.get_cmap('tab10')(regions.index(rg)),
                label=f'{rg} ({int((g.seg_region == rg).sum())} segments)')
    ax.plot([], [], 's', color='0.35', label='coil, driven'); ax.plot([], [], 's', color='0.8', label='coil, open')
    ax.plot([], [], 'k', label='limiter')
    ax.set(xlabel='R [m]', ylabel='Z [m]', title='ARC V3A start-up model')
    ax.legend(fontsize=8, loc='center left', bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    return fig


def _shade(ax, window):
    ax.axvspan(*window, color='gold', alpha=0.35, lw=0)


def plot_stage1(s1, res, t_lim=(-2.05, 0.5)):
    """Stage-1 traces: coil currents, supply voltages, eddy currents, breakdown-region field, E_phi
    and Lloyd quantity, plasma current and flux."""
    c, m, t = s1.c, s1.m, res['t']
    bd = s1.scen['breakdown']
    used = np.abs(c.T).sum(axis=1) > 0
    lim = m.coil_lim[:, 1]
    fig, ax = plt.subplots(2, 3, figsize=(16, 8), sharex=True)
    ax = ax.ravel()
    cols = plt.get_cmap('tab20')(np.linspace(0, 1, 20))
    for i in np.flatnonzero(used):
        ax[0].plot(t, res['I_coil'][:, i] / lim[i], color=cols[i % 20], lw=1.2, label=m.coil_names[i])
    ax[0].set(ylabel='$I / I_{max}$', title='coil currents', ylim=(-1.08, 1.08))
    ax[0].legend(fontsize=6, ncol=2)
    tt = np.repeat(t, 2)[1:-1]
    for j in range(c.ncirc):
        i = np.flatnonzero(c.T[:, j])[0]
        ax[1].plot(tt, np.repeat(res['V'][:, j], 2) / c.V_max[j], color=cols[i % 20], lw=1)
    ax[1].set(ylabel='$V / V_{max}$', title='supply voltages (constant on each step)', ylim=(-1.08, 1.08))
    ax[2].plot(t, res['I_struct'].sum(axis=1) / 1e6, 'C5', label='total structure current')
    ax[2].plot(t, s1.s * res['Ip'] / 1e6, 'C3', label='plasma current $I_p$')
    ax[2].set(ylabel='MA', title='eddy and plasma currents'); ax[2].legend(fontsize=8)
    ax[3].semilogy(t, res['B_pts'].max(axis=1) * 1e3, 'C0', label='max over region')
    ax[3].semilogy(t, res['B_pts'].mean(axis=1) * 1e3, 'C0--', label='mean over region')
    ax[3].set(ylabel='mT', title='poloidal field $|B_p|$ in the breakdown region'); ax[3].legend(fontsize=8)
    ax[4].plot(t[1:], res['E_pts'][1:].min(axis=1), 'C2', label='min $sE_\\phi$ over region [V/m]')
    ax[4].axhline(bd['E_min'], color='C2', ls=':')
    a2 = ax[4].twinx()
    a2.semilogy(t[1:], np.maximum(res['lloyd_pts'][1:].min(axis=1), 1e-1), 'k', lw=1)
    a2.axhline(bd['margin'] * bd['E_star'], color='k', ls=':')
    a2.set_ylabel('min $E_\\phi B_\\phi / B_p$ [V/m] (black)')
    ax[4].set(ylim=(-0.2, 3), title='breakdown criteria (dotted: required)'); ax[4].legend(fontsize=8, loc='upper left')
    ax[5].plot(t, res['flux_ref_Wb'], 'C4')
    ax[5].set(ylabel='Wb', title='flux available at the region centre $2\\pi s\\,\\psi$')
    for a in ax:
        _shade(a, bd['window']); a.set_xlim(*t_lim); a.grid(alpha=0.3)
    for a in ax[3:]:
        a.set_xlabel('t [s]')
    fig.tight_layout()
    return fig


def nodal_flux(model, circuit, x, plasma=None, Ip=0.0):
    """Nodal flux of a circuit state: coils + structure (+ the channel if plasma is given)."""
    psi = model.psi_u.T @ (circuit.U @ x)
    if plasma is not None and Ip:
        psi = psi + Ip * plasma.psi_unit(model.gs.r[:, :2])
    return psi


def field_maps(model, circuit, x, xdot, s, R_lim, Z_lim, n=61):
    """|B_p|, E_phi and the Lloyd quantity s E_phi B_phi / |B_p| on an R-Z grid (channel included)."""
    Rg, Zg = np.meshgrid(np.linspace(*R_lim, n), np.linspace(*Z_lim, n))
    pts = np.column_stack([Rg.ravel(), Zg.ravel()])
    f = eval_fields(model.gs, np.vstack([model.psi_u.T @ (circuit.U @ x), model.psi_u.T @ (circuit.U @ xdot)]), pts)
    if circuit.plasma is not None:
        gp = circuit.plasma.field_at(pts)
        for key, gg in zip(('psi', 'Br', 'Bz'), gp):
            f[key] = f[key] + np.outer([x[circuit.ip], xdot[circuit.ip]], gg)
    Bp = np.hypot(f['Br'][0], f['Bz'][0]).reshape(Rg.shape)
    E = (-f['psi'][1] / pts[:, 0]).reshape(Rg.shape)
    return {'R': Rg, 'Z': Zg, 'Bp': Bp, 'E': E, 'psi': f['psi'][0].reshape(Rg.shape),
            'lloyd': s * E * np.abs(model.Bphi(Rg)) / np.maximum(Bp, 1e-12)}


def plot_breakdown_map(ax, mp, g, region, threshold, title=''):
    """|B_p| (colour, log scale), flux contours (grey) and the contour where the Lloyd criterion
    is exactly met (black): inside it the gas breaks down."""
    h = ax.pcolormesh(mp['R'], mp['Z'], mp['Bp'] * 1e3, norm=LogNorm(0.1, 100), cmap='viridis', shading='gouraud')
    ax.contour(mp['R'], mp['Z'], mp['psi'], 15, colors='w', linewidths=0.5, alpha=0.6)
    ax.contour(mp['R'], mp['Z'], mp['lloyd'], [threshold], colors='k', linewidths=2)
    R0, Z0, r = region
    ax.add_patch(plt.Circle((R0, Z0), r, fill=False, ls='--', ec='C3', lw=1.5))
    ax.plot(g.lim[:, 0], g.lim[:, 1], 'k', lw=1)
    ax.set_xlim(mp['R'].min(), mp['R'].max()); ax.set_ylim(mp['Z'].min(), mp['Z'].max())
    ax.set_aspect('equal'); ax.set(xlabel='R [m]', ylabel='Z [m]', title=title)
    return h


def plot_snapshot(ax, g, psi, I_seg, Jlim, used, lcfs=None, iso=None, title=''):
    """Whole machine at one time: structure current density, coils, flux contours, and (ramp
    slices) the last closed flux surface and the target boundary points."""
    h = _structure(ax, g, I_seg=I_seg, Jlim=Jlim, used=used)
    inside = Path(g.lim).contains_points(np.column_stack([g.R, g.Z]))
    lev = np.linspace(psi[inside].min(), psi[inside].max(), 21)
    ax.tricontour(g.tri, psi, levels=lev, colors='0.4', linewidths=0.5)
    if lcfs is not None:
        ax.tricontour(g.tri_in, psi, levels=[lcfs], colors='C3', linewidths=1.8)
    if iso is not None:
        ax.plot(iso[:, 0], iso[:, 1], 'C1+', ms=6, mew=1.5)
    ax.set_title(title, fontsize=10)
    return h


def plot_pulse(s1, ramp, res):
    """The whole start-up: stage 1 (to the hand-off, dotted line) then the GS ramp."""
    c, m = s1.c, s1.m
    t1, t2 = res['t'], ramp.t
    used = np.abs(c.T).sum(axis=1) > 0
    lim = m.coil_lim[:, 1]
    I2 = (ramp.c.T @ ramp.x[:, :ramp.c.ncirc].T).T
    cols = plt.get_cmap('tab20')(np.linspace(0, 1, 20))
    fig, ax = plt.subplots(2, 2, figsize=(15, 8), sharex=True)
    ax = ax.ravel()
    ax[0].plot(t1, s1.s * res['Ip'] / 1e6, 'C3', label='stage 1: channel current')
    ax[0].plot(t2, ramp.Ip / 1e6, 'C3o-', ms=3, label='stage 2: prescribed $I_p$')
    ax[0].plot(t2[1:], [i['Ip'] / 1e6 for i in ramp.tm_info[1:]], 'kx', label='TokaMaker equilibria')
    ax[0].set(ylabel='MA', title='plasma current'); ax[0].legend(fontsize=8)
    for i in np.flatnonzero(used):
        ax[1].plot(np.r_[t1, t2[1:]], np.r_[res['I_coil'][:, i], I2[1:, i]] / lim[i], color=cols[i % 20],
                   lw=1.2, label=m.coil_names[i])
    ax[1].set(ylabel='$I / I_{max}$', title='coil currents', ylim=(-1.08, 1.08)); ax[1].legend(fontsize=6, ncol=2)
    tt = np.repeat(np.r_[t1, t2[1:]], 2)[1:-1]
    VV = np.vstack([res['V'], ramp.V])
    for j in range(c.ncirc):
        i = np.flatnonzero(c.T[:, j])[0]
        ax[2].plot(tt, np.repeat(VV[:, j], 2) / c.V_max[j], color=cols[i % 20], lw=1)
    ax[2].set(ylabel='$V / V_{max}$', title='supply voltages', ylim=(-1.08, 1.08))
    ax[3].plot(t1, res['flux_ref_Wb'], 'C4', label='stage 1: flux at the region centre')
    ax[3].plot(t2, 2 * np.pi * ramp.Fb, 'C4o-', ms=3, label="stage 2: plasma boundary flux (Ohm's law)")
    ax[3].set(ylabel='Wb', title='flux available to the pulse'); ax[3].legend(fontsize=8)
    for a in ax:
        _shade(a, s1.scen['breakdown']['window']); a.axvline(t2[0], color='k', ls=':', lw=1)
        a.set_xlim(-0.6, t2[-1] + 0.02); a.grid(alpha=0.3)
    for a in ax[2:]:
        a.set_xlabel('t [s]')
    fig.tight_layout()
    return fig


def lcfs_distance(g, psi, psi_lcfs, pts, center):
    """Distance [m] from each point in pts to the contour psi = psi_lcfs inside the limiter (the last
    closed flux surface); of several pieces, the one around center (the magnetic axis) is used."""
    fig, ax = plt.subplots()
    polys = ax.tricontour(g.tri_in, psi, levels=[psi_lcfs]).get_paths()[0].to_polygons(closed_only=False)
    plt.close(fig)
    inside = [q for q in polys if len(q) > 2 and Path(q).contains_point(center)]
    q = max(inside or polys, key=len)
    a, b = q[:-1], q[1:]                                 # segments of the contour polyline
    ab = b - a
    tt = np.clip(np.einsum('ijk,jk->ij', pts[:, None, :] - a[None], ab) / np.maximum((ab ** 2).sum(axis=1), 1e-30), 0, 1)
    return np.linalg.norm(pts[:, None, :] - (a[None] + tt[..., None] * ab[None]), axis=2).min(axis=1)
