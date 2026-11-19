#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
"""Figures for the breakdown / start-up example: the machine, stage-1 and whole-pulse traces,
breakdown maps, whole-machine snapshots and an animated GIF of the start-up. Also provides
lcfs_distance, the distance from target points to TokaMaker's last closed flux surface."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.ticker as mticker
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from matplotlib.figure import Figure
from matplotlib.path import Path

from startup_model import eval_fields, cell_geometry


class MeshGeometry:
    """Triangulations, limiter and coil/segment cells of a MachineModel, derived once.

    On TokaMaker's mesh: R, Z, lc (nodes, triangles), lim (longest limiter contour [n, 2]),
    cell_centre, node_in_lim, tri (whole mesh) and tri_in (triangles inside the limiter).
    On the model mesh: coil_cells / seg_cells (cell indices per coil / segment), coil_centre,
    seg_area, seg_region, region_colour (tab10 colour per region), coil_polys / seg_polys (their
    triangles [n, 3, 2]), coil_of_cell / seg_of_cell (owner of each of those triangles) and
    extent (R_min, R_max, Z_min, Z_max) framing every coil and segment."""

    def __init__(self, model):
        gs = model.gs
        self.R, self.Z = gs.r[:, 0].copy(), gs.r[:, 1].copy()
        self.lc = np.asarray(gs.lc)
        self.lim = max(gs.lim_contours, key=len)[:, :2]
        lim_path = Path(self.lim)
        self.cell_centre = np.column_stack([self.R[self.lc].mean(axis=1), self.Z[self.lc].mean(axis=1)])
        self.node_in_lim = lim_path.contains_points(np.column_stack([self.R, self.Z]))
        self.tri = mtri.Triangulation(self.R, self.Z, self.lc)
        self.tri_in = mtri.Triangulation(self.R, self.Z, self.lc, mask=~lim_path.contains_points(self.cell_centre))

        P, LC, REG = model.mesh
        area, centroid = cell_geometry(P, LC)
        self.coil_cells = [np.where(REG == model.coil_dict[n]['reg_id'])[0] for n in model.coil_names]
        self.coil_centre = np.array([centroid[c].mean(axis=0) for c in self.coil_cells])
        self.seg_cells = [s['cells'] for s in model.segments]
        self.seg_area = np.array([area[c].sum() for c in self.seg_cells])
        self.seg_region = np.array([s['region'] for s in model.segments])
        self.region_colour = {r: plt.get_cmap('tab10')(i) for i, r in enumerate(dict.fromkeys(self.seg_region))}
        self.seg_of_cell = np.repeat(np.arange(len(self.seg_cells)), [len(c) for c in self.seg_cells])
        self.coil_of_cell = np.repeat(np.arange(len(self.coil_cells)), [len(c) for c in self.coil_cells])
        self.seg_polys = P[:, :2][LC[np.concatenate(self.seg_cells)]]
        self.coil_polys = P[:, :2][LC[np.concatenate(self.coil_cells)]]

        all_pts = np.vstack([self.seg_polys.reshape(-1, 2), self.coil_polys.reshape(-1, 2)])
        lo, hi = all_pts.min(axis=0), all_pts.max(axis=0)
        pad = 0.05 * (hi - lo).max()
        self.extent = (max(lo[0] - pad, 0.0), hi[0] + pad, lo[1] - pad, hi[1] + pad)


def plot_machine(model, g, used, region, channel, shapes):
    """Coils (driven dark, open light), structure segments by region (alternate segments darker),
    limiter, breakdown region (R0, Z0, radius), early channel (Rc, Zc, a) and the ramp's target
    shapes (list of [n, 2] boundaries, drawn closed)."""
    fig, ax = plt.subplots(figsize=(9, 8))
    base = np.array([g.region_colour[r][:3] for r in g.seg_region])
    shade = np.where(np.arange(len(g.seg_cells)) % 2, 0.75, 1.0)[:, None]
    fc = np.c_[(base * shade)[g.seg_of_cell], np.ones(len(g.seg_of_cell))]
    ax.add_collection(PolyCollection(g.seg_polys, facecolors=fc, edgecolors=fc, linewidths=0.3))
    coil_fc = np.array(['0.35' if used[i] else '0.8' for i in g.coil_of_cell])
    ax.add_collection(PolyCollection(g.coil_polys, facecolors=coil_fc, edgecolors=coil_fc, linewidths=0.3))
    ax.plot(g.lim[:, 0], g.lim[:, 1], 'k', lw=1.2)
    for i, name in enumerate(model.coil_names):
        ax.text(*g.coil_centre[i], name, fontsize=7, ha='center', va='center', color='w' if used[i] else 'k')
    ax.add_patch(plt.Circle(region[:2], region[2], fill=False, ls='--', ec='C3', lw=1.5, label='breakdown region'))
    ax.add_patch(plt.Circle(channel[:2], channel[2], fill=False, ec='C1', lw=1.5, label='early plasma channel'))
    for i, p in enumerate(shapes):
        ax.plot(np.r_[p[:, 0], p[0, 0]], np.r_[p[:, 1], p[0, 1]], color='C4', lw=0.7, alpha=0.7,
                label=None if i else 'ramp target shapes')
    for rg, col in g.region_colour.items():                 # legend-only entries
        ax.plot([], [], 's', color=col, label=f'{rg} ({int((g.seg_region == rg).sum())} segments)')
    ax.plot([], [], 's', color='0.35', label='coil, driven')
    ax.plot([], [], 's', color='0.8', label='coil, open')
    ax.plot([], [], 'k', label='limiter')
    ax.set(xlim=g.extent[:2], ylim=g.extent[2:], aspect='equal', xlabel='R [m]', ylabel='Z [m]',
           title='ARC V3A start-up model')
    ax.legend(fontsize=8, loc='center left', bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    return fig


def _coil_and_supply_panels(ax_I, ax_V, s1, t, I_coil, V, V_title):
    """Driven-coil currents (ax_I) and supply voltages (ax_V, staircases: V[k] holds on
    [t[k], t[k+1]]), each normalised to its limit. Both coils of an up/down pair (names differing
    only in a trailing U/L) share a colour, the lower one dashed; a supply takes the style of the
    first coil it drives."""
    m, circ = s1.m, s1.c
    stem = [n[:-1] if n[-1:] in ('U', 'L') else n for n in m.coil_names]
    stems = list(dict.fromkeys(stem))
    cmap = plt.get_cmap('tab10' if len(stems) <= 10 else 'tab20')
    style = [{'color': cmap(stems.index(s) % cmap.N), 'ls': '--' if n.endswith('L') else '-'}
             for n, s in zip(m.coil_names, stem)]
    for i in np.flatnonzero(np.abs(circ.T).sum(axis=1) > 0):
        ax_I.plot(t, I_coil[:, i] / m.coil_lim[i, 1], lw=1.2, label=m.coil_names[i], **style[i])
    ax_I.set(ylabel='$I / I_{max}$', title='coil currents', ylim=(-1.08, 1.08))
    ax_I.legend(fontsize=6, ncol=2)
    for j in range(circ.ncirc):
        ax_V.plot(np.repeat(t, 2)[1:-1], np.repeat(V[:, j], 2) / circ.V_max[j], lw=1,
                  **style[np.flatnonzero(circ.T[:, j])[0]])
    ax_V.set(ylabel='$V / V_{max}$', title=V_title, ylim=(-1.08, 1.08))


def plot_stage1(s1, res, t_lim=None):
    """Stage-1 traces: coil currents, supply voltages, eddy currents, breakdown-region field, E_phi
    and Lloyd quantity, plasma current and flux. t_lim: time axis (default: the whole stage)."""
    t, bd = res['t'], s1.scen['breakdown']
    t_lim = (t[0] - 0.05, t[-1]) if t_lim is None else t_lim
    fig, ax = plt.subplots(2, 3, figsize=(16, 8), sharex=True)
    ax = ax.ravel()
    _coil_and_supply_panels(ax[0], ax[1], s1, t, res['I_coil'], res['V'], 'supply voltages (constant on each step)')
    ax[2].plot(t, res['I_struct'].sum(axis=1) / 1e6, 'C5', label='total structure current')
    ax[2].plot(t, s1.s * res['Ip'] / 1e6, 'C3', label='plasma current $I_p$')
    ax[2].set(ylabel='MA', title='eddy and plasma currents')
    ax[2].legend(fontsize=8)
    ax[3].semilogy(t, res['B_pts_ext'].max(axis=1) * 1e3, 'C0', label='max over region')
    ax[3].semilogy(t, res['B_pts_ext'].mean(axis=1) * 1e3, 'C0--', label='mean over region')
    ax[3].set(ylabel='mT', title='applied $|B_p|$ (coils + structure) in the breakdown region')
    ax[3].legend(fontsize=8)
    ax[4].plot(t[1:], res['E_pts'][1:].min(axis=1), 'C2', label='min $sE_\\phi$ over region [V/m]')
    ax[4].axhline(bd['E_min'], color='C2', ls=':')
    ax_lloyd = ax[4].twinx()
    ax_lloyd.semilogy(t[1:], np.maximum(res['lloyd_pts'][1:].min(axis=1), 1e-1), 'k', lw=1)
    ax_lloyd.axhline(bd['margin'] * bd['E_star'], color='k', ls=':')
    ax_lloyd.set_ylabel('min $E_\\phi B_\\phi / B_p$ [V/m] (black)')
    ax[4].set(ylim=(-0.2, 3), title='breakdown criteria (dotted: required)')
    ax[4].legend(fontsize=8, loc='upper left')
    ax[5].plot(t, res['flux_bound_Wb'], 'C4')
    ax[5].set(ylabel='Wb', title='flux at the plasma boundary (inboard contact point) $2\\pi\\psi$')
    for i, a in enumerate(ax):
        a.axvspan(*bd['window'], color='gold', alpha=0.35, lw=0)
        a.set_xlim(*t_lim)
        a.grid(alpha=0.3)
        if i >= 3:
            a.set_xlabel('t [s]')
    fig.tight_layout()
    return fig


def plot_pulse(s1, ramp, res, t_lim=None):
    """The whole start-up: stage 1 (to the hand-off, dotted line) then the GS ramp. t_lim: time
    axis (default: from 0.6 s before the breakdown window to the end of the ramp)."""
    window = s1.scen['breakdown']['window']
    t_lim = (window[0] - 0.6, ramp.t[-1] + 0.02) if t_lim is None else t_lim
    t1, t2 = res['t'], ramp.t
    fig, ax = plt.subplots(2, 2, figsize=(15, 8), sharex=True)
    ax = ax.ravel()
    ax[0].plot(t1, s1.s * res['Ip'] / 1e6, 'C3', label='stage 1: channel current')
    ax[0].plot(t2, ramp.Ip / 1e6, 'C3o-', ms=3, label='stage 2: prescribed $I_p$')
    ax[0].plot(t2[1:], [info['Ip'] / 1e6 for info in ramp.tm_info[1:]], 'kx', label='TokaMaker equilibria')
    ax[0].set(ylabel='MA', title='plasma current')
    ax[0].legend(fontsize=8)
    I_coil_ramp = (ramp.c.T @ ramp.x[:, :ramp.c.ncirc].T).T
    _coil_and_supply_panels(ax[1], ax[2], s1, np.r_[t1, t2[1:]],      # the ramp's first slice is the hand-off
                            np.vstack([res['I_coil'], I_coil_ramp[1:]]), np.vstack([res['V'], ramp.V]),
                            'supply voltages')
    ax[3].plot(t1, res['flux_bound_Wb'], 'C4', label='stage 1: channel boundary flux')
    ax[3].plot(t2, 2 * np.pi * ramp.Fb, 'C4o-', ms=3, label="stage 2: plasma boundary flux (Ohm's law)")
    ax[3].set(ylabel='Wb', title='flux available to the pulse')
    ax[3].legend(fontsize=8)
    for i, a in enumerate(ax):
        a.axvspan(*window, color='gold', alpha=0.35, lw=0)
        a.axvline(t2[0], color='k', ls=':', lw=1)
        a.set_xlim(*t_lim)
        a.grid(alpha=0.3)
        if i >= 2:
            a.set_xlabel('t [s]')
    fig.tight_layout()
    return fig


def nodal_flux(model, circuit, x, plasma=None, Ip=0.0):
    """Nodal flux of a state x: coils + structure (+ the channel if plasma is given). circuit: anything
    with the units <- state map U, e.g. a CircuitModel with the channel added, or a StartupProblem."""
    psi = model.psi_u.T @ (circuit.U @ x)
    if plasma is not None and Ip:
        psi = psi + Ip * plasma.psi_unit(model.gs.r[:, :2])
    return psi


def field_maps(model, circuit, x, xdot, s, R_lim, Z_lim, n=61):
    """|B_p|, E_phi and the Lloyd quantity s E_phi B_phi / |B_p| on an R-Z grid (channel included).
    circuit: anything with U, plasma and ip, e.g. a CircuitModel with the channel added, or a
    StartupProblem (needed when I_p is prescribed and the channel is not a circuit state).

    Returns a dict of [n, n] arrays: R, Z, Bp, E, psi, Br, Bz and lloyd."""
    Rg, Zg = np.meshgrid(np.linspace(*R_lim, n), np.linspace(*Z_lim, n))
    pts = np.column_stack([Rg.ravel(), Zg.ravel()])
    # row 0: the state x, row 1: its time derivative xdot
    fields = eval_fields(model.gs, np.vstack([model.psi_u.T @ (circuit.U @ x),
                                              model.psi_u.T @ (circuit.U @ xdot)]), pts)
    if circuit.plasma is not None:
        for key, val in zip(('psi', 'Br', 'Bz'), circuit.plasma.field_at(pts)):
            fields[key] = fields[key] + np.outer([x[circuit.ip], xdot[circuit.ip]], val)
    Bp = np.hypot(fields['Br'][0], fields['Bz'][0]).reshape(Rg.shape)
    E = (-fields['psi'][1] / pts[:, 0]).reshape(Rg.shape)
    return {'R': Rg, 'Z': Zg, 'Bp': Bp, 'E': E, 'psi': fields['psi'][0].reshape(Rg.shape),
            'Br': fields['Br'][0].reshape(Rg.shape), 'Bz': fields['Bz'][0].reshape(Rg.shape),
            'lloyd': s * E * np.abs(model.Bphi(Rg)) / np.maximum(Bp, 1e-12)}


def plot_breakdown_map(ax, mp, g, region, threshold, title=''):
    """|B_p| (colour, log scale), flux contours (white) and the contour where the Lloyd criterion
    is exactly met (black): inside it the gas breaks down. mp: a dict from field_maps."""
    h = ax.pcolormesh(mp['R'], mp['Z'], mp['Bp'] * 1e3, norm=LogNorm(0.1, 100), cmap='viridis', shading='gouraud')
    ax.contour(mp['R'], mp['Z'], mp['psi'], 15, colors='w', linewidths=0.5, alpha=0.6)
    ax.contour(mp['R'], mp['Z'], mp['lloyd'], [threshold], colors='k', linewidths=2)
    ax.add_patch(plt.Circle(region[:2], region[2], fill=False, ls='--', ec='C3', lw=1.5))
    ax.plot(g.lim[:, 0], g.lim[:, 1], 'k', lw=1)
    ax.set(xlim=(mp['R'].min(), mp['R'].max()), ylim=(mp['Z'].min(), mp['Z'].max()), aspect='equal',
           xlabel='R [m]', ylabel='Z [m]', title=title)
    return h


def plot_snapshot(ax, g, psi, I_seg, Jlim, used, lcfs=None, iso=None, title=''):
    """Whole machine at one time: structure current density [MA/m^2] on +-Jlim (the returned
    collection, for a colour bar), coils (open ones light grey), flux contours and, on ramp slices,
    the last closed flux surface psi = lcfs and the target boundary points iso."""
    h = PolyCollection(g.seg_polys, array=(I_seg / g.seg_area / 1e6)[g.seg_of_cell], cmap='PuOr_r',
                       edgecolors='face', linewidths=0.2)
    h.set_clim(-Jlim, Jlim)
    ax.add_collection(h)
    coil_fc = np.array(['0.35' if used[i] else '0.8' for i in g.coil_of_cell])
    ax.add_collection(PolyCollection(g.coil_polys, facecolors=coil_fc, edgecolors=coil_fc, linewidths=0.3))
    ax.plot(g.lim[:, 0], g.lim[:, 1], 'k', lw=1.2)
    ax.set(xlim=g.extent[:2], ylim=g.extent[2:], aspect='equal')
    psi_in = psi[g.node_in_lim]
    ax.tricontour(g.tri, psi, levels=np.linspace(psi_in.min(), psi_in.max(), 21), colors='0.4', linewidths=0.5)
    if lcfs is not None:                                 # the closed piece around the plasma only
        poly = _lcfs_poly(g, psi, lcfs, iso.mean(axis=0) if iso is not None else g.lim.mean(axis=0))
        if poly is not None:
            ax.plot(poly[:, 0], poly[:, 1], 'C3', lw=1.8)
    if iso is not None:
        ax.plot(iso[:, 0], iso[:, 1], 'C1+', ms=6, mew=1.5)
    ax.set_title(title, fontsize=10)
    return h


def _lcfs_pieces(g, psi, psi_lcfs):
    """Every piece of the contour psi = psi_lcfs inside the limiter, as [n, 2] polylines."""
    path = Figure().add_subplot().tricontour(g.tri_in, psi, levels=[psi_lcfs]).get_paths()[0]
    verts, codes = path.vertices, path.codes
    if codes is None:
        return [verts] if len(verts) > 1 else []
    # split at MOVETO by hand: Path.to_polygons simplifies paths of >= 128 vertices (pixel threshold)
    starts = np.flatnonzero(codes == Path.MOVETO).tolist() + [len(verts)]
    pieces = []
    for i0, i1 in zip(starts[:-1], starts[1:]):
        piece = verts[i0:i1].copy()
        if codes[i1 - 1] == Path.CLOSEPOLY:           # closed piece: its last vertex is a placeholder
            piece[-1] = piece[0]
        if len(piece) > 1:
            pieces.append(piece)
    return pieces


def _lcfs_poly(g, psi, psi_lcfs, center):
    """For drawing: the piece of the last closed flux surface around center (the magnetic axis),
    or the longest piece if none encloses it (a limited boundary cut by the wall)."""
    pieces = _lcfs_pieces(g, psi, psi_lcfs)
    enclosing = [q for q in pieces if len(q) > 2 and Path(q).contains_point(center)]
    return max(enclosing or pieces, key=len) if pieces else None


def lcfs_distance(g, psi, psi_lcfs, pts, center=None):
    """Distance [m] from each point in pts to the last closed flux surface psi = psi_lcfs: the nearest
    point over every contour piece inside the limiter (a limited boundary can be cut into several
    pieces by the wall, none of which need enclose the axis). Returns inf for every point if the
    contour is empty. center is not used; it is kept for existing callers."""
    pieces = _lcfs_pieces(g, psi, psi_lcfs)
    if not pieces:
        return np.full(len(pts), np.inf)
    a = np.vstack([q[:-1] for q in pieces])              # every line segment a -> b of every piece
    ab = np.vstack([q[1:] for q in pieces]) - a
    frac = np.clip(np.einsum('ijk,jk->ij', pts[:, None, :] - a[None], ab)
                   / np.maximum((ab ** 2).sum(axis=1), 1e-30), 0, 1)
    return np.linalg.norm(pts[:, None, :] - (a[None] + frac[..., None] * ab[None]), axis=2).min(axis=1)


# ============================================================================ start-up animation
# plain journal style, applied only inside make_startup_gif
_RC = {'font.family': 'serif', 'font.serif': ['STIXGeneral', 'DejaVu Serif'], 'mathtext.fontset': 'stix',
       'font.size': 8.5, 'axes.labelsize': 8.5, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'font.weight': 'normal',
       'axes.labelweight': 'normal', 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in',
       'xtick.top': True, 'ytick.right': True, 'xtick.minor.visible': True, 'ytick.minor.visible': True,
       'xtick.major.size': 3.0, 'ytick.major.size': 3.0, 'xtick.minor.size': 1.6, 'ytick.minor.size': 1.6,
       'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
       'lines.linewidth': 1.1, 'lines.markeredgewidth': 0.6, 'axes.grid': False}


def startup_frames(model, g, s1, res, ramp, zoom, n_grid=61):
    """Everything the animation draws, one entry per time slice: stage 1 (premagnetisation to the
    hand-off) then the TokaMaker ramp. Stage-1 slices carry |B_p| and the Lloyd quantity on an
    R-Z grid over the zoom window (R_min, R_max, Z_min, Z_max): the applied field (coils +
    structure) up to the end of the breakdown window, then the stray field, i.e. the applied field
    minus the vertical field that the channel's radial force balance needs. Ramp slices carry the
    last closed flux surface. Returns a dict for make_startup_gif."""
    rg, pl, bd = s1.scen['region'], s1.scen['plasma'], s1.scen['breakdown']
    frames = []
    for k, t in enumerate(s1.t):
        x_ext, xdot_ext = res['x'][k].copy(), res['xdot'][k].copy()
        x_ext[s1.ip] = xdot_ext[s1.ip] = 0.0             # applied field only: coils + structure
        mp = field_maps(model, s1, x_ext, xdot_ext, s1.s, zoom[:2], zoom[2:], n=n_grid)
        phase = 'channel' if t > s1.t_p0 + 1e-9 else ('breakdown' if t > bd['window'][0] + 1e-9 else 'premag')
        frames.append({'t': t, 'Ip': res['Ip'][k], 'I_coil': res['I_coil'][k], 'I_seg': res['I_struct'][k],
                       'psi': nodal_flux(model, s1, res['x'][k], s1.plasma, res['Ip'][k]), 'phase': phase,
                       'R': mp['R'], 'Z': mp['Z'],     # channel: stray field, applied minus the required B_v
                       'Bp': np.hypot(mp['Br'], mp['Bz'] + s1.kf * s1.s * res['Ip'][k]) if phase == 'channel'
                       else mp['Bp'],
                       'lloyd': mp['lloyd'] if 0 < k and phase != 'channel' else None})
    nc = ramp.c.ncirc
    for k in range(1, len(ramp.t)):                       # slice 0 is the hand-off, already shown
        info = ramp.tm_info[k]
        frames.append({'t': ramp.t[k], 'Ip': ramp.Ip[k], 'I_coil': ramp.c.T @ ramp.x[k][:nc], 'I_seg': ramp.x[k][nc:],
                       'psi': ramp.psi_tot[k], 'phase': 'ramp', 'psi_bounds': info['psi_bounds'], 'iso': ramp.iso[k],
                       'lcfs': _lcfs_poly(g, ramp.psi_tot[k], info['psi_bounds'][0], info['o_point'])})
    traces = {'t': np.r_[s1.t, ramp.t[1:]], 'Ip': np.r_[res['Ip'], ramp.Ip[1:]],
              't1': s1.t, 'Bmax': res['B_pts_ext'].max(axis=1), 'flux1': res['flux_bound_Wb'],
              't2': ramp.t, 'flux2': 2 * np.pi * ramp.Fb, 'stray': s1.channel_stray()}
    return {'frames': frames, 'traces': traces, 'used': np.abs(s1.c.T).sum(axis=1) > 0,
            'coil_lim': np.abs(model.coil_lim[:, 1]), 'region': (rg['R'], rg['Z'], rg['radius']),
            'channel': (pl['Rc'], pl['Zc'], pl['a']), 'window': bd['window'],
            'threshold': bd['margin'] * bd['E_star'], 'zoom': zoom}


def make_startup_gif(d, g, path, dpi=100, colors=96, ms_frame=160, ms_hold=1500):
    """Animated GIF of the start-up from startup_frames: (a) machine, (b) plasma region,
    (c)-(e) plasma current, largest |B_p| in the breakdown region, and flux at the plasma boundary
    (the channel edge in stage 1, TokaMaker's boundary on the ramp). colors: GIF palette size;
    ms_frame / ms_hold: duration of each frame / of the first and last frames [ms]. Returns path."""
    from PIL import Image
    BLUE, RED, LIGHT = '#0077BB', '#CC3311', '0.85'
    frames, traces, zoom = d['frames'], d['traces'], d['zoom']
    # fixed colour scales and trace limits over the whole animation
    J = np.array([f['I_seg'] for f in frames]) / g.seg_area / 1e3      # structure current density [kA/m^2]
    Jlim = np.percentile(np.abs(J), 99)
    Jnorm, Inorm, Bnorm = TwoSlopeNorm(0, -Jlim, Jlim), Normalize(-1, 1), LogNorm(0.1, 100)
    Ip_kA, Bmax_mT = traces['Ip'] / 1e3, 1e3 * traces['Bmax']
    flux = np.r_[traces['flux1'], traces['flux2']]
    ylims = [(-0.03 * Ip_kA.max(), 1.05 * Ip_kA.max()),
             (10 ** np.floor(np.log10(max(Bmax_mT.min(), 1e-3))), 10 ** np.ceil(np.log10(max(Bmax_mT.max(), 1e-2)))),
             (flux.min() - 0.05 * np.ptp(flux), flux.max() + 0.05 * np.ptp(flux))]
    images = []
    with plt.rc_context(_RC):
        for frame, J_f in zip(frames, J):
            fig = Figure(figsize=(7.0, 4.4), dpi=dpi)
            canvas = FigureCanvasAgg(fig)
            gsp = fig.add_gridspec(3, 3, width_ratios=[1.0, 1.05, 1.15], left=0.06, right=0.985, bottom=0.31,
                                   top=0.9, wspace=0.42, hspace=0.12)
            ax_m, ax_z = fig.add_subplot(gsp[:, 0]), fig.add_subplot(gsp[:, 1])
            ax_tr = [fig.add_subplot(gsp[i, 2]) for i in range(3)]

            # (a) whole machine: structure current density, coil current / limit, poloidal flux
            ax_m.add_collection(PolyCollection(g.seg_polys, array=J_f[g.seg_of_cell], cmap='RdBu_r', norm=Jnorm,
                                               edgecolors='face', linewidths=0.2))
            rel = frame['I_coil'] / d['coil_lim']
            fc = [plt.get_cmap('RdBu_r')(Inorm(rel[i])) if d['used'][i] else LIGHT for i in g.coil_of_cell]
            ax_m.add_collection(PolyCollection(g.coil_polys, facecolors=fc, edgecolors=fc, linewidths=0.2))
            psi = frame['psi']
            psi_in = psi[g.node_in_lim]
            ax_m.tricontour(g.tri, psi, levels=np.linspace(psi_in.min() - 0.5 * np.ptp(psi_in),
                                                           psi_in.max() + 0.5 * np.ptp(psi_in), 30),
                            colors='0.55', linewidths=0.3)
            ax_m.plot(g.lim[:, 0], g.lim[:, 1], 'k', lw=0.6)
            ax_m.add_patch(plt.Rectangle((zoom[0], zoom[2]), zoom[1] - zoom[0], zoom[3] - zoom[2], fill=False,
                                         ec='k', lw=0.5, ls=':'))
            ax_m.set(xlim=g.extent[:2], ylim=g.extent[2:], xlabel='$R$ [m]', ylabel='$Z$ [m]', aspect='equal')

            # (b) zoom window: |B_p| and the breakdown contour in stage 1, flux surfaces on the ramp
            if frame['phase'] != 'ramp':
                grid_pts = np.column_stack([frame['R'].ravel(), frame['Z'].ravel()])
                outside = ~Path(g.lim).contains_points(grid_pts).reshape(frame['R'].shape)
                bp = np.ma.masked_where(outside, np.clip(frame['Bp'] * 1e3, 0.1001, 99.9))
                ax_z.contourf(frame['R'], frame['Z'], bp, levels=np.logspace(-1, 2, 25), norm=Bnorm, cmap='Greys')
                if frame['lloyd'] is not None:
                    lloyd = np.ma.masked_where(outside, frame['lloyd'])
                    if lloyd.min() < d['threshold'] < lloyd.max():
                        ax_z.contour(frame['R'], frame['Z'], lloyd, levels=[d['threshold']], colors=RED,
                                     linewidths=1.0)
                ax_z.add_patch(plt.Circle(d['region'][:2], d['region'][2], fill=False, ls='--', ec=BLUE, lw=0.9))
                if frame['phase'] == 'channel':
                    ax_z.add_patch(plt.Circle(d['channel'][:2], d['channel'][2], fill=False, ec=RED, lw=1.2))
            else:
                ax_z.tricontour(g.tri_in, psi, levels=np.linspace(psi_in.min(), psi_in.max(), 25), colors='0.75',
                                linewidths=0.3)
                poly = frame['lcfs']
                if poly is not None:
                    # normalised flux surfaces, drawn only on the triangles inside the LCFS
                    psi_bound, psi_axis = frame['psi_bounds']       # TokaMaker order: [LCFS, axis]
                    tri_plasma = mtri.Triangulation(g.R, g.Z, g.lc, mask=~Path(poly).contains_points(g.cell_centre))
                    ax_z.tricontour(tri_plasma, (psi - psi_axis) / (psi_bound - psi_axis),
                                    levels=np.linspace(0.1, 0.9, 9), colors='0.4', linewidths=0.5)
                    ax_z.plot(poly[:, 0], poly[:, 1], color=RED, lw=1.2)
                ax_z.plot(frame['iso'][:, 0], frame['iso'][:, 1], 'o', mfc='none', mec=BLUE, ms=2.5)
            ax_z.plot(g.lim[:, 0], g.lim[:, 1], 'k', lw=0.8)
            ax_z.set(xlim=zoom[:2], ylim=zoom[2:], xlabel='$R$ [m]', aspect='equal')
            phase_label = {'premag': 'premagnetisation and swing', 'breakdown': 'breakdown window',
                           'channel': 'early plasma channel', 'ramp': 'free-boundary ramp'}[frame['phase']]
            fig.text(0.5, 0.965, f'$t$ = {frame["t"]:+.3f} s     $I_p$ = {frame["Ip"] / 1e3:5.0f} kA     {phase_label}',
                     ha='center', va='center', fontsize=9)

            # (c)-(e) traces with a cursor at the current time
            ax_tr[0].plot(traces['t'], traces['Ip'] / 1e3, color='k')
            ax_tr[0].set(ylabel='$I_p$ [kA]')
            stray = traces['stray']
            after = stray['t'] >= d['window'][1] - 1e-9
            ax_tr[1].semilogy(traces['t1'], 1e3 * traces['Bmax'], color='k', label='applied, breakdown region')
            ax_tr[1].semilogy(stray['t'][after], 1e3 * stray['stray'][after], color=RED, label='stray, channel')
            ax_tr[1].set(ylabel='max $|B_p|$ [mT]')
            ax_tr[1].legend(loc='upper right', fontsize=6.5, frameon=False, handlelength=1.2, borderaxespad=0.2)
            ax_tr[1].yaxis.set_major_formatter(mticker.LogFormatterMathtext())
            ax_tr[2].plot(traces['t1'], traces['flux1'], color='k')
            ax_tr[2].plot(traces['t2'], traces['flux2'], color=RED)
            ax_tr[2].set(ylabel='flux [Wb]', xlabel='time [s]')
            for i, (axt, yl) in enumerate(zip(ax_tr, ylims)):
                axt.axvspan(*d['window'], color=LIGHT, lw=0)
                axt.axvline(frame['t'], color=BLUE, lw=1.0)
                axt.set_xlim(traces['t'][0] - 0.05, traces['t'][-1] + 0.05)
                axt.set_ylim(*yl)
                if i < 2:
                    axt.tick_params(labelbottom=False)

            # panel labels; (b) sits above its map, which fills the panel
            ax_m.text(0.03, 0.985, '(a)', transform=ax_m.transAxes, ha='left', va='top')
            ax_z.text(0.0, 1.01, '(b)', transform=ax_z.transAxes, ha='left', va='bottom')
            for axt, lab in zip(ax_tr, 'cde'):
                axt.text(0.02, 0.93, f'({lab})', transform=axt.transAxes, ha='left', va='top')

            # colour bars: (grid cell, norm, colour map, label, bottom of the bar in figure units)
            B_label = ('stray $|B_p|$ (applied minus required $B_v$) [mT]' if frame['phase'] in ('channel', 'ramp')
                       else 'applied $|B_p|$ [mT]')
            for cell, norm, cmap, label, y0 in ((gsp[:, 0], Inorm, 'RdBu_r', 'coil current / limit', 0.175),
                                                (gsp[:, 0], Jnorm, 'RdBu_r', 'structure current density [kA/m$^2$]', 0.065),
                                                (gsp[:, 1], Bnorm, 'Greys', B_label, 0.175)):
                p = cell.get_position(fig)
                cb = fig.colorbar(plt.cm.ScalarMappable(norm, cmap), orientation='horizontal',
                                  cax=fig.add_axes([p.x0 + 0.02, y0, p.width - 0.04, 0.02]))
                cb.set_label(label, fontsize=7.5, labelpad=1)
                cb.ax.tick_params(labelsize=7, which='both', top=False, pad=1)
                if norm is Bnorm:                        # log scale
                    cb.set_ticks([0.1, 1, 10, 100])
                    cb.set_ticklabels(['0.1', '1', '10', '100'])
                else:
                    cb.ax.xaxis.set_major_locator(mticker.MaxNLocator(5))
            canvas.draw()
            images.append(Image.fromarray(np.asarray(canvas.buffer_rgba())[..., :3])
                          .convert('P', palette=Image.ADAPTIVE, colors=colors))
    durations = [ms_frame] * len(images)
    durations[0] = durations[-1] = ms_hold
    images[0].save(path, save_all=True, append_images=images[1:], duration=durations, loop=0, optimize=True)
    return path
