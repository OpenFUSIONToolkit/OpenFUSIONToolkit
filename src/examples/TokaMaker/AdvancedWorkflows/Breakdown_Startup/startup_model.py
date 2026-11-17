"""Machine, circuit and early-plasma models for breakdown and start-up design with TokaMaker.

Contents
    segment_conductors  split every passive conducting region of a TokaMaker mesh into segments,
                        each of which TokaMaker then treats as a 1-turn coil with a 1/R current
                        distribution (the eddy-current model of the structure)
    MachineModel        inductance matrix, resistances, unit flux responses and coil limits,
                        all computed by TokaMaker
    CircuitModel        M dx/dt + R x = S V for the driven coils + structure (+ plasma channel),
                        with an exact zero-order-hold discretisation
    ChannelPlasma       the early plasma as a uniform-current ring: one more circuit element
    FieldTracer, townsend_window, paschen_Emin
                        connection length and Townsend avalanche check of a breakdown

Conventions: coil currents are per-turn amps (as TokaMaker's set_coil_currents takes them);
psi is TokaMaker's poloidal flux per radian [Wb/rad], so the loop voltage is -2 pi dpsi/dt and
E_phi = -(dpsi/dt) / R; B_R = -(dpsi/dZ)/R, B_Z = (dpsi/dR)/R.
"""
import numpy as np
import scipy.sparse as sp
from scipy.linalg import expm
from matplotlib.path import Path

MU0 = 4e-7 * np.pi
LN_LAMBDA = 15.0        # Coulomb logarithm, used for the Spitzer resistivity in both stages


# ============================================================================ field evaluation
def eval_fields(gs, psi_vecs, pts, second=False, h=5.0e-3):
    """psi, B_R, B_Z of each nodal flux vector in psi_vecs at the points pts, using TokaMaker's own
    finite-element interpolation. Returns a dict of arrays [n_vec, n_pts]; with second=True also
    dBz_dR (central difference of the interpolated gradient, step h)."""
    psi_vecs = np.atleast_2d(psi_vecs)
    pts = np.atleast_2d(np.asarray(pts, float))
    nv, npt = psi_vecs.shape[0], pts.shape[0]
    out = {k: np.zeros((nv, npt)) for k in ('psi', 'Br', 'Bz')}
    if second:
        out['dBz_dR'] = np.zeros((nv, npt))

    def B_at(fg, x):
        g = fg.eval(x)[:2]
        return -g[1] / x[0], g[0] / x[0]

    for i in range(nv):
        psi_i = np.ascontiguousarray(psi_vecs[i])
        fe, fg = gs.get_field_eval('eval', values=psi_i), gs.get_field_eval('grad', values=psi_i)
        for j in range(npt):
            x = pts[j]
            out['psi'][i, j] = fe.eval(x)[0]
            out['Br'][i, j], out['Bz'][i, j] = B_at(fg, x)
            if second:
                out['dBz_dR'][i, j] = (B_at(fg, x + [h, 0.0])[1] - B_at(fg, x - [h, 0.0])[1]) / (2 * h)
    return out


def disk_points(Rc, Zc, radius, n_rings=3, n_per_ring=8):
    """Centre plus concentric rings of points covering a disk (the breakdown region)."""
    pts = [[Rc, Zc]]
    for k in range(1, n_rings + 1):
        r = radius * k / n_rings
        th = 2 * np.pi * (np.arange(n_per_ring * k) + 0.5 * (k % 2)) / (n_per_ring * k)
        pts += [[Rc + r * np.cos(t), Zc + r * np.sin(t)] for t in th]
    return np.array(pts)


# ============================================================================ structure segments
def cell_geometry(P, LC):
    """Area and centroid of every triangle of the mesh."""
    v1, v2 = P[LC[:, 1]] - P[LC[:, 0]], P[LC[:, 2]] - P[LC[:, 0]]
    return 0.5 * np.abs(v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]), P[LC].mean(axis=1)


def _connected_pieces(groups, LC, min_cells):
    """Split each group of cells into node-connected pieces (a patch must be ONE closed toroidal
    loop, not two disjoint ones), then merge pieces smaller than min_cells into a piece they touch."""
    from scipy.sparse.csgraph import connected_components
    allc = np.concatenate(groups)
    nn = int(LC.max()) + 1
    inc = sp.csr_matrix((np.ones(LC[allc].size), (np.repeat(np.arange(len(allc)), LC.shape[1]),
                                                  LC[allc].ravel())), shape=(len(allc), nn))
    adj = (inc @ inc.T).tocsr()                          # cells sharing at least one node
    pos = {c: i for i, c in enumerate(allc)}
    pieces = []
    for g in groups:
        li = np.array([pos[c] for c in g])
        ncomp, lab = connected_components(adj[li][:, li], directed=False)
        pieces += [g[lab == j] for j in range(ncomp)]
    big = [c for c in pieces if len(c) >= min_cells]
    for c in (c for c in pieces if len(c) < min_cells):
        touch = set(allc[adj[np.array([pos[x] for x in c])].indices])
        cand = [j for j, b in enumerate(big) if touch.intersection(b)]
        if cand:
            big[cand[0]] = np.concatenate([big[cand[0]], c])
        else:
            big.append(c)
    return big


def segment_conductors(P, LC, REG, cond_dict, n_segments, center, min_cells=2, seed=12345):
    """Split every conducting region (eta > 0) of a TokaMaker mesh into segments.

    Each region is cut into n_segments[name] equal-size patches by k-means clustering of the cell
    centroids (so a thick region, e.g. a blanket, is also split through its thickness); patches
    are ordered by poloidal angle about `center`. Every segment is later declared to TokaMaker as
    a 1-turn coil carrying a 1/R current distribution (uniform E_phi), and its resistance follows
    from the region's resistivity: R = 2 pi eta / sum(A_cell / R_cell).

    Returns (reg, segments, seg_coils, vac_dict):
        reg       new region id per mesh cell (each segment gets its own id)
        segments  list of dicts: name, region, reg_id, cells, eta, R_ohm
        seg_coils TokaMaker coil_dict entries for the segments
        vac_dict  the non-conducting regions (TokaMaker cond_dict without the segmented ones)
    """
    from scipy.cluster.vq import kmeans2
    area, cen = cell_geometry(P, LC)
    reg = REG.copy()
    next_id = int(REG.max()) + 1
    cond = {k: v for k, v in cond_dict.items() if v.get('eta', -1.0) > 0}
    vac_dict = {k: v for k, v in cond_dict.items() if not (v.get('eta', -1.0) > 0)}
    segments, seg_coils = [], {}
    for name, info in cond.items():
        idx = np.where(REG == info['reg_id'])[0]
        if len(idx) == 0:
            continue
        nk = min(int(n_segments[name]), len(idx))
        _, lab = kmeans2(cen[idx], nk, minit='++', seed=seed, iter=50)
        order = np.argsort([np.arctan2(*(cen[idx[lab == j]].mean(axis=0) - center)[::-1])
                            if np.any(lab == j) else 9.0 for j in range(nk)])
        groups = [idx[lab == j] for j in order if np.any(lab == j)]
        for j, cells in enumerate(_connected_pieces(groups, LC, min_cells)):
            rid = info['reg_id'] if j == 0 else next_id
            next_id += (j > 0)
            reg[cells] = rid
            R_ohm = 2 * np.pi * info['eta'] / np.sum(area[cells] / cen[cells, 0])     # J ~ 1/R
            seg = {'name': f'{name}#{j:03d}', 'region': name, 'reg_id': int(rid), 'cells': cells,
                   'eta': float(info['eta']), 'R_ohm': float(R_ohm)}
            segments.append(seg)
            seg_coils[seg['name']] = {'reg_id': int(rid), 'nturns': 1}
            if 'allow_xpoints' in info:
                seg_coils[seg['name']]['allow_xpoints'] = info['allow_xpoints']
    return reg, segments, seg_coils, vac_dict


# ============================================================================ machine model
class MachineModel:
    """Coils + segmented structure on a TokaMaker instance that was set up with the segment coils.

    Attributes
        coil_names, seg_names, unit_names   units = coils followed by segments
        L        [nu, nu]  mutual inductances (TokaMaker's Lcoils, per-turn basis)
        R        [nu]      resistances (superconducting coils: 0; segments: from eta)
        psi_u    [nu, np]  nodal flux of a unit current in each unit (vacuum)
        coil_lim [nc, 2]   per-turn current limits from a current-density limit J_max
    """

    def __init__(self, gs, mesh, coil_dict, segments, F0, J_max=40.0e6, log=print):
        P, LC, REG = mesh
        self.gs, self.mesh, self.coil_dict, self.segments, self.F0 = gs, mesh, coil_dict, segments, F0
        self.coil_names = list(coil_dict.keys())
        self.seg_names = [s['name'] for s in segments]
        self.unit_names = self.coil_names + self.seg_names
        self.nturns = np.array([coil_dict[c]['nturns'] for c in self.coil_names], float)
        # an equilibrium object must exist before current distributions / field evaluations
        gs.set_coil_currents({self.coil_names[0]: 1.0})
        gs.solve(vacuum=True)
        rn = gs.r[:, 0]
        dist = np.where(rn > 1e-8, 1.0 / np.maximum(rn, 1e-8), 0.0)
        for s in self.seg_names:                         # uniform E_phi inside each segment
            gs.set_coil_current_dist(s, dist, normalize=True)
        ids = [gs.coil_sets[n]['id'] for n in self.unit_names]
        L = np.array(gs.Lcoils)[np.ix_(ids, ids)]
        self.L = 0.5 * (L + L.T)
        log(f'computing the vacuum flux of a unit current in each of {len(self.unit_names)} units ...')
        self.psi_u = np.array([self.unit_psi({n: 1.0}) for n in self.unit_names])
        area, _ = cell_geometry(P, LC)
        self.coil_area = np.array([area[REG == coil_dict[c]['reg_id']].sum() for c in self.coil_names])
        Imax = J_max * self.coil_area / self.nturns
        self.coil_lim = np.column_stack([-Imax, Imax])
        self.R = np.concatenate([np.zeros(self.nc), [s['R_ohm'] for s in segments]])

    def unit_psi(self, currents):
        """Static vacuum flux [Wb/rad] of a dict of unit currents (all other units zero)."""
        gs = self.gs
        gs.set_coil_currents(currents)
        gs.set_psi(np.zeros(gs.np))
        return gs.vac_solve().get_psi(normalized=False)

    @property
    def nc(self):
        return len(self.coil_names)

    @property
    def ns(self):
        return len(self.seg_names)

    def Bphi(self, R):
        return self.F0 / np.asarray(R)

    def structure_time_constants(self):
        """L/R eigenvalues of the structure alone [s], slowest first."""
        nc = self.nc
        return np.sort(np.linalg.eigvals(np.linalg.solve(np.diag(self.R[nc:]), self.L[nc:, nc:])).real)[::-1]


# ============================================================================ circuit model
class CircuitModel:
    """M dx/dt + R x = S V.

    State x = [I_circ (one current per driven coil, per-turn A); I_seg (structure segments);
    I_p (plasma channel, once add_plasma is called)]. Coils listed in `exclude` (substring match)
    are left open: no current, not in the state. V_max: supply voltage limit [V]. Between time slices the supply voltages V are held
    constant, and the dynamics are integrated exactly with a matrix exponential (zoh), so there is
    no time-step error and superconducting (R = 0) circuits need no special treatment."""

    R_OPEN = 1.0e3          # plasma 'resistance' before the channel exists (effectively open)

    def __init__(self, machine, V_max, exclude=()):
        m = machine
        nc, ns = m.nc, m.ns
        keep = [j for j, n in enumerate(m.coil_names) if not any(e in n for e in exclude)]
        self.circ_names = [m.coil_names[j] for j in keep]
        self.ncirc = len(keep)
        T = np.zeros((nc, self.ncirc))
        T[keep, np.arange(self.ncirc)] = 1.0
        self.T = T                                       # coil <- circuit incidence
        L = m.L
        Lcc, Lcs, Lss = L[:nc, :nc], L[:nc, nc:], L[nc:, nc:]
        self.M = np.block([[T.T @ Lcc @ T, T.T @ Lcs], [Lcs.T @ T, Lss]])
        self.Rm = np.diag(np.concatenate([T.T @ m.R[:nc], m.R[nc:]]))
        self.S = np.vstack([np.eye(self.ncirc), np.zeros((ns, self.ncirc))])
        self.U = np.block([[T, np.zeros((nc, ns))], [np.zeros((ns, self.ncirc)), np.eye(ns)]])   # units <- state
        self.nx = self.ncirc + ns
        self.V_max = np.full(self.ncirc, float(V_max))  # supply voltage limit per circuit [V]
        self.plasma = None
        self.Minv = np.linalg.inv(self.M)
        self._zoh = {}

    def add_plasma(self, plasma):
        """Append the plasma channel current I_p as one more state. The channel has no supply: it
        is driven only by the changing flux of the coils and structure, L_p dIp/dt + R_p Ip +
        sum_u M_pu dI_u/dt = 0. Its resistance is set interval by interval in zoh(dt, Rp)."""
        n0 = self.nx
        Mpx = plasma.M_u @ self.U
        M = np.zeros((n0 + 1, n0 + 1)); M[:n0, :n0] = self.M
        M[n0, :n0] = M[:n0, n0] = Mpx; M[n0, n0] = plasma.L_p
        R = np.zeros((n0 + 1, n0 + 1)); R[:n0, :n0] = self.Rm; R[n0, n0] = self.R_OPEN
        self.M, self.Rm = M, R
        self.S = np.vstack([self.S, np.zeros((1, self.ncirc))])
        self.U = np.hstack([self.U, np.zeros((self.U.shape[0], 1))])
        self.nx, self.ip, self.plasma = n0 + 1, n0, plasma
        self.Minv = np.linalg.inv(self.M)
        self._zoh = {}

    def output_matrix(self, G_units, G_plasma=None):
        """Per-unit response matrix [n_out, nu] (and per-ampere plasma response) -> [n_out, nx]."""
        C = G_units @ self.U
        if G_plasma is not None and self.plasma is not None:
            C[:, self.ip] = G_plasma
        return C

    def zoh(self, dt, Rp=None):
        """Exact step x_k = A x_{k-1} + B V_k for voltages held constant over an interval dt
        (Rp: plasma resistance on the interval). exp([[Ac, Bc], [0, 0]] dt) = [[A, B], [0, I]]."""
        key = (round(float(dt), 12), None if Rp is None else float(f'{Rp:.6e}'))
        if key not in self._zoh:
            n, m = self.nx, self.ncirc
            Rm = self.Rm
            if Rp is not None and self.plasma is not None:
                Rm = Rm.copy(); Rm[self.ip, self.ip] = Rp
            Z = np.zeros((n + m, n + m))
            Z[:n, :n], Z[:n, n:] = -self.Minv @ Rm * dt, self.Minv @ self.S * dt
            E = expm(Z)
            self._zoh[key] = (E[:n, :n], E[:n, n:])
        return self._zoh[key]


# ============================================================================ early plasma channel
def _disk_quadrature(Rc, Zc, a, n_r, n_t):
    pts, w = [], []
    edges = np.linspace(0, a, n_r + 1)
    for i in range(n_r):
        r = 0.5 * (edges[i] + edges[i + 1])
        for j in range(n_t):
            th = 2 * np.pi * (j + 0.5 * (i % 2)) / n_t
            pts.append([Rc + r * np.cos(th), Zc + r * np.sin(th)])
            w.append(np.pi * (edges[i + 1] ** 2 - edges[i] ** 2) / n_t)
    w = np.array(w)
    return np.array(pts), w / w.sum()


class ChannelPlasma:
    """The early plasma as a toroidal ring of uniform current: minor radius a, centred at (Rc, Zc).

        L_p = mu0 Rc (ln(8 Rc/a) - 2 + l_i/2)       self-inductance of a ring
        M_u = 2 pi <psi_u>_channel                  mutual inductance with every coil and segment
        R_p = eta(Te) 2 pi Rc / (pi a^2)            Spitzer resistance,
        eta = 5.2e-5 Zeff lnL / Te^1.5  [Ohm m, Te in eV]  parallel Spitzer (NRL Plasma Formulary)
    """

    def __init__(self, machine, Rc, Zc, a, li=0.5, Zeff=2.0, lnL=LN_LAMBDA):
        self.Rc, self.Zc, self.a, self.li, self.Zeff, self.lnL = Rc, Zc, a, li, Zeff, lnL
        self.qp, self.qw = _disk_quadrature(Rc, Zc, a, n_r=8, n_t=24)
        psi = eval_fields(machine.gs, machine.psi_u, self.qp)['psi']
        self.M_u = 2 * np.pi * psi @ self.qw
        self.L_p = MU0 * Rc * (np.log(8 * Rc / a) - 2 + li / 2)

    def resistance(self, Te_eV):
        eta = 5.2e-5 * self.Zeff * self.lnL / np.maximum(Te_eV, 0.1) ** 1.5
        return eta * 2 * np.pi * self.Rc / (np.pi * self.a ** 2)

    def field_at(self, pts, h=1e-3):
        """psi, B_R, B_Z per ampere of channel current at pts (TokaMaker sign convention)."""
        pts = np.atleast_2d(pts)
        dR = (self.psi_unit(pts + [h, 0]) - self.psi_unit(pts - [h, 0])) / (2 * h)
        dZ = (self.psi_unit(pts + [0, h]) - self.psi_unit(pts - [0, h])) / (2 * h)
        return self.psi_unit(pts), -dZ / pts[:, 0], dR / pts[:, 0]

    def psi_unit(self, x):
        """Flux per ampere at points x: a sum of filaments (TokaMaker's eval_green), where each
        filament stands for a small uniform-current cell, so points inside the channel see a
        finite, smooth field instead of the filament's log singularity."""
        from OpenFUSIONToolkit.TokaMaker.util import eval_green
        x = np.atleast_2d(np.asarray(x, float))
        delta = self.a * np.sqrt(self.qw)                # radius of each filament's cell
        q0 = self.qp[0]
        if not hasattr(self, '_cR'):                     # psi ~ -cR ln(d) near a filament
            d1, d2 = 1e-4, 1e-3
            g = -eval_green(np.array([[q0[0] + d1, q0[1]], [q0[0] + d2, q0[1]]]), q0)
            self._cR = -(g[1] - g[0]) / np.log(d2 / d1)
        tot = np.zeros(len(x))
        for q, w, dl in zip(self.qp, self.qw, delta):
            d = np.hypot(x[:, 0] - q[0], x[:, 1] - q[1])
            f = -eval_green(x, q)
            inn = d < dl
            if inn.any():                                # uniform-wire core inside the cell
                di = np.maximum(d[inn], 1e-12)
                f[inn] += -(self._cR * q[0] / q0[0]) * (np.log(dl) - np.log(di) + 0.5 * (di ** 2 / dl ** 2 - 1))
            tot += w * f
        return tot


def Te_burnthrough(t, t_bd, Te0=5.0, Te1=100.0, tau=0.1):
    """Prescribed electron temperature [eV] of the channel: Te0 at the end of breakdown, rising
    to Te1 with time constant tau (a placeholder for a burn-through model)."""
    t = np.asarray(t, float)
    return np.where(t < t_bd, Te0, Te0 + (Te1 - Te0) * (1 - np.exp(-(t - t_bd) / tau)))


# ============================================================================ breakdown physics
TORR = 133.322          # Pa
# first Townsend coefficient of H2, alpha = A p exp(-B p / E); close to the values tabulated by
# Raizer, Gas Discharge Physics (1991), A = 5 /(cm Torr), B = 130 V/(cm Torr)
TOWNSEND_A = 510.0      # 1/(m Torr)
TOWNSEND_B = 1.25e4     # V/(m Torr)


class FieldTracer:
    """Follows magnetic field lines (poloidal field from a nodal flux + toroidal field F0/R) from a
    point in both directions until they leave the wall polygon: L_c = L_forward + L_backward."""

    def __init__(self, gs, psi_nodal, F0, wall, Bp_min=1e-6):
        self.fg = gs.get_field_eval('grad', values=np.ascontiguousarray(psi_nodal))
        self.F0, self.wall, self.Bp_min = abs(F0), Path(wall[:, :2]), Bp_min

    def _rhs(self, x, sgn):
        g = self.fg.eval(x)[:2]
        br, bz, bt = -g[1] / x[0], g[0] / x[0], self.F0 / x[0]
        bp = max(np.hypot(br, bz), self.Bp_min)
        return np.array([sgn * br / bp, sgn * bz / bp]), np.sqrt(bp ** 2 + bt ** 2) / bp

    def _rk4(self, x, sgn, dl):
        k1, w1 = self._rhs(x, sgn)
        k2, w2 = self._rhs(x + 0.5 * dl * k1, sgn)
        k3, w3 = self._rhs(x + 0.5 * dl * k2, sgn)
        k4, w4 = self._rhs(x + dl * k3, sgn)
        return x + dl / 6 * (k1 + 2 * k2 + 2 * k3 + k4), dl / 6 * (w1 + 2 * w2 + 2 * w3 + w4)

    def trace(self, x0, sgn, dl=0.01, L_max=1e5, n_max=2000):
        """Length [m] along B to the wall (poloidal step dl; L_max if the line never leaves)."""
        x, L = np.array(x0, float), 0.0
        for _ in range(n_max):
            xn, dL = self._rk4(x, sgn, dl)
            if not self.wall.contains_point(xn):         # bisect the wall-crossing step
                lo, hi = 0.0, dl
                for _ in range(20):
                    mid = 0.5 * (lo + hi)
                    lo, hi = (mid, hi) if self.wall.contains_point(self._rk4(x, sgn, mid)[0]) else (lo, mid)
                return min(L + self._rk4(x, sgn, lo)[1], L_max)
            L += dL
            if L > L_max:
                return L_max
            x = xn
        return L_max

    def connection_length(self, x0, **kw):
        return self.trace(x0, +1, **kw) + self.trace(x0, -1, **kw)


def townsend_window(E, L, A=TOWNSEND_A, B=TOWNSEND_B):
    """Range of prefill pressures [Pa] for which an electron avalanche runs away (alpha L >= 1)."""
    p = np.logspace(-7, 3, 10000)
    pt = p / TORR
    ok = A * pt * L * np.exp(-B * pt / max(E, 1e-9)) >= 1.0
    return (float(p[ok].min()), float(p[ok].max())) if ok.any() else None


def paschen_Emin(p, L, A=TOWNSEND_A, B=TOWNSEND_B):
    """Minimum E_phi [V/m] for an avalanche at pressure p [Pa] over a length L [m]."""
    pt = np.asarray(p) / TORR
    x = A * pt * L
    return np.where(x > 1, B * pt / np.log(np.maximum(x, 1 + 1e-12)), np.inf)
