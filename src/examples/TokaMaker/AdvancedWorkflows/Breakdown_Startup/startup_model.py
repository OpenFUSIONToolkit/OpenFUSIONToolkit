#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
"""Machine, circuit and early-plasma models for breakdown and start-up design with TokaMaker.

    segment_conductors  split the passive conductors into segments, each a 1-turn coil with a 1/R
                        current distribution (the eddy-current model of the structure)
    MachineModel        inductances, resistances, unit flux maps and coil limits from TokaMaker
    CircuitModel        M dx/dt + R x = S V for the supplies, structure (and plasma channel)
    feed_forward        structure (and channel) currents driven by given coil currents (validation)
    ChannelPlasma       the early plasma as a uniform-current ring: one more circuit element
    FieldTracer, townsend_window, paschen_Emin   connection length and Townsend avalanche check

SI units (m, s, A, V) except electron temperatures [eV] and the Townsend constants [Torr based].
Coil currents are per-turn [A]; psi is TokaMaker's poloidal flux per radian [Wb/rad], so the loop
voltage is -2 pi dpsi/dt, E_phi = -(dpsi/dt) / R, B_R = -(dpsi/dZ) / R and B_Z = (dpsi/dR) / R.
"""
import numpy as np
import scipy.sparse as sp
from scipy.cluster.vq import kmeans2
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components
from matplotlib.path import Path

MU0 = 4e-7 * np.pi      # vacuum permeability [H/m]
LN_LAMBDA = 15.0        # Coulomb logarithm, used for the Spitzer resistivity in both stages


# ============================================================================ geometry and fields
def eval_fields(gs, psi_vecs, pts, second=False, h=5.0e-3):
    """psi [Wb/rad], B_R and B_Z [T] of nodal flux vectors psi_vecs [n_vec, n_nodes] at points
    pts [n_pts, 2] (R, Z) [m], by TokaMaker's finite-element interpolation. If second, also
    dBz_dR [T/m] by a central difference of the gradient with step h [m].
    Returns a dict of arrays [n_vec, n_pts]."""
    psi_vecs = np.atleast_2d(psi_vecs)
    pts = np.atleast_2d(np.asarray(pts, float))
    R, dR = pts[:, 0], np.array([h, 0.0])
    out = {k: np.zeros((len(psi_vecs), len(pts))) for k in ('psi', 'Br', 'Bz', 'dBz_dR')[:4 if second else 3]}
    for i, psi_i in enumerate(psi_vecs):
        psi_i = np.ascontiguousarray(psi_i)
        out['psi'][i] = gs.get_field_eval('eval', values=psi_i).eval(pts)[:, 0]
        f_grad = gs.get_field_eval('grad', values=psi_i)
        grad = f_grad.eval(pts)                          # columns: dpsi/dR, dpsi/dZ
        out['Br'][i], out['Bz'][i] = -grad[:, 1] / R, grad[:, 0] / R
        if second:
            out['dBz_dR'][i] = (f_grad.eval(pts + dR)[:, 0] / (R + h)
                                - f_grad.eval(pts - dR)[:, 0] / (R - h)) / (2 * h)
    return out


def disk_points(Rc, Zc, radius, n_rings=3, n_per_ring=8):
    """Points [n_pts, 2] covering a disk of radius [m] centred at (Rc, Zc) [m]: the centre plus
    n_rings rings, ring k holding n_per_ring * k points (odd rings rotated by half a spacing)."""
    pts = [[Rc, Zc]]
    for k in range(1, n_rings + 1):
        r, n_k = radius * k / n_rings, n_per_ring * k
        theta = 2 * np.pi * (np.arange(n_k) + 0.5 * (k % 2)) / n_k
        pts += [[Rc + r * np.cos(t), Zc + r * np.sin(t)] for t in theta]
    return np.array(pts)


def _disk_quadrature(Rc, Zc, a, n_r, n_t):
    """Midpoint quadrature of a disk of radius a [m] centred at (Rc, Zc) [m]: n_r annuli of equal
    width, n_t points each (odd annuli rotated by half a spacing).
    Returns (points [n_r * n_t, 2], weights [n_r * n_t]) with area-fraction weights summing to 1."""
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


def cell_geometry(P, LC):
    """(area [m^2] [n_cells], centroid [m] [n_cells, 2]) of the mesh triangles (nodes P, cells LC)."""
    v1, v2 = P[LC[:, 1]] - P[LC[:, 0]], P[LC[:, 2]] - P[LC[:, 0]]
    return 0.5 * np.abs(v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]), P[LC].mean(axis=1)


# ============================================================================ structure segments
def segment_conductors(P, LC, REG, cond_dict, n_segments, center, min_cells=2, seed=12345):
    """Split every conducting region (eta > 0) of a TokaMaker mesh into segments.

    Each region is cut into n_segments[name] patches by k-means of the cell centroids (so thick
    regions are also split through their thickness), ordered by poloidal angle about center
    (R, Z) [m]. Patches that are not node-connected are split (a segment must be ONE toroidal
    loop) and pieces with < min_cells cells are merged into the first large piece they touch.
    Each segment is a 1-turn coil with J ~ 1/R, so R_seg = 2 pi eta / sum(A_cell / R_cell) [Ohm].

    Args: P, LC, REG mesh nodes, cells and region id per cell; cond_dict TokaMaker conductor dict
    {name: {'reg_id', 'eta' [Ohm m], ...}}; n_segments {region name: count}; seed of k-means.

    Returns (reg, segments, seg_coils, vac_dict):
        reg        new region id per cell (the first segment of a region keeps the region's id)
        segments   list of dicts: name, region, reg_id, cells, eta [Ohm m], R_ohm [Ohm]
        seg_coils  TokaMaker coil_dict entries of the segments
        vac_dict   the non-conducting entries of cond_dict
    """
    area, centroid = cell_geometry(P, LC)
    reg = REG.copy()
    next_id = int(REG.max()) + 1
    vac_dict = {k: v for k, v in cond_dict.items() if not (v.get('eta', -1.0) > 0)}
    segments, seg_coils = [], {}
    for name, info in cond_dict.items():
        region_cells = [] if name in vac_dict else np.where(REG == info['reg_id'])[0]
        if len(region_cells) == 0:
            continue
        n_k = min(int(n_segments[name]), len(region_cells))
        _, label = kmeans2(centroid[region_cells], n_k, minit='++', seed=seed, iter=50)
        groups = [region_cells[label == j] for j in range(n_k) if np.any(label == j)]
        dR, dZ = (np.array([centroid[c].mean(axis=0) for c in groups]) - center).T
        groups = [groups[j] for j in np.argsort(np.arctan2(dZ, dR))]
        # connected pieces of each patch, as rows of all_cells; adj = cells sharing a node
        all_cells = np.concatenate(groups)
        n_all = len(all_cells)
        incidence = sp.csr_matrix((np.ones(LC[all_cells].size),
                                   (np.repeat(np.arange(n_all), LC.shape[1]), LC[all_cells].ravel())),
                                  shape=(n_all, int(LC.max()) + 1))
        adj = (incidence @ incidence.T).tocsr()
        pieces, i0 = [], 0
        for cells in groups:
            rows = np.arange(i0, i0 + len(cells))
            i0 += len(cells)
            n_comp, lab = connected_components(adj[rows][:, rows], directed=False)
            pieces += [rows[lab == j] for j in range(n_comp)]
        large = [p for p in pieces if len(p) >= min_cells]
        for p in pieces:
            if len(p) < min_cells:
                touching = set(adj[p].indices)
                j = next((j for j, big in enumerate(large) if touching.intersection(big)), None)
                if j is None:
                    large.append(p)
                else:
                    large[j] = np.concatenate([large[j], p])
        for j, rows in enumerate(large):
            cells = all_cells[rows]
            rid = info['reg_id'] if j == 0 else next_id + j - 1
            reg[cells] = rid
            R_ohm = 2 * np.pi * info['eta'] / np.sum(area[cells] / centroid[cells, 0])
            seg_name = f'{name}#{j:03d}'
            segments.append({'name': seg_name, 'region': name, 'reg_id': int(rid), 'cells': cells,
                             'eta': float(info['eta']), 'R_ohm': float(R_ohm)})
            seg_coils[seg_name] = {'reg_id': int(rid), 'nturns': 1,
                                   **{k: info[k] for k in ('allow_xpoints',) if k in info}}
        next_id += len(large) - 1
    return reg, segments, seg_coils, vac_dict


# ============================================================================ machine model
class MachineModel:
    """Driven coils + segmented structure of a TokaMaker instance set up with the segment coils.

    Args: gs TokaMaker instance; mesh (P, LC, REG); coil_dict the driven coils (no segments), as
    from load_gs_mesh: entries sharing a 'coil_set' are one coil (in series, one current per turn);
    segments from segment_conductors; F0 = R B_phi [T m]; J_max coil current-density limit
    [A/m^2]; coil_lim optional {coil: (lo, hi)} per-turn current limits [A] overriding J_max;
    log progress-message callable.

    Attributes: coil_names (coil sets), seg_names, unit_names (units = coils then segments); nc, ns;
    nturns [nc] (net turns of each set); L [nu, nu] mutual inductances [H] (per-turn basis); R [nu]
    resistances [Ohm] (coils 0: superconducting); psi_u [nu, n_nodes] vacuum flux [Wb/rad] of 1 A per
    turn in each unit; coil_cells (cells of each coil) and coil_turns (turns carried by each of those
    cells, each element's turns spread over its area); coil_area [nc] [m^2]; coil_lim [nc, 2] per-turn
    current limits (lo, hi) [A].
    """

    def __init__(self, gs, mesh, coil_dict, segments, F0, J_max=40.0e6, coil_lim=None, log=print):
        P, LC, REG = mesh
        self.gs, self.mesh, self.coil_dict, self.segments, self.F0 = gs, mesh, coil_dict, segments, F0
        sets = {}
        for name, c in coil_dict.items():
            sets.setdefault(c.get('coil_set', name), []).append(c)
        self.coil_names = list(sets)
        self.seg_names = [s['name'] for s in segments]
        self.unit_names = self.coil_names + self.seg_names
        self.nc, self.ns = len(self.coil_names), len(self.seg_names)
        area = cell_geometry(P, LC)[0]
        self.nturns = np.zeros(self.nc)
        self.coil_cells, self.coil_turns = [], []
        for i, n in enumerate(self.coil_names):
            cells, turns = [], []
            for c in sets[n]:                            # every element (mesh region) of the coil set
                x = np.where(REG == c['reg_id'])[0]
                cells.append(x)
                turns.append(c['nturns'] * area[x] / area[x].sum())   # the element's turns spread over its area
                self.nturns[i] += c['nturns']
            self.coil_cells.append(np.concatenate(cells))
            self.coil_turns.append(np.concatenate(turns))
        gs.set_coil_currents({self.coil_names[0]: 1.0})  # an equilibrium must exist before the
        gs.solve(vacuum=True)                            # current distributions / field evaluations
        r_nodes = gs.r[:, 0]
        inv_r = np.where(r_nodes > 1e-8, 1.0 / np.maximum(r_nodes, 1e-8), 0.0)
        for s in self.seg_names:                         # J ~ 1/R: uniform E_phi in each segment
            gs.set_coil_current_dist(s, inv_r, normalize=True)
        ids = [gs.coil_sets[n]['id'] for n in self.unit_names]
        L = np.array(gs.Lcoils)[np.ix_(ids, ids)]
        self.L = 0.5 * (L + L.T)                         # exact symmetry
        log(f'computing the vacuum flux of a unit current in each of {len(self.unit_names)} units ...')
        self.psi_u = np.array([self.unit_psi({n: 1.0}) for n in self.unit_names])
        self.coil_area = np.array([area[c].sum() for c in self.coil_cells])
        I_max = J_max * self.coil_area / self.nturns
        self.coil_lim = np.column_stack([-I_max, I_max])
        for name, lim in (coil_lim or {}).items():
            self.coil_lim[self.coil_names.index(name)] = lim
        self.R = np.concatenate([np.zeros(self.nc), [s['R_ohm'] for s in segments]])

    def unit_psi(self, currents):
        """Vacuum flux [Wb/rad] at the nodes of {unit: current [A per turn]} (other units zero)."""
        self.gs.set_coil_currents(currents)
        self.gs.set_psi(np.zeros(self.gs.np))
        return self.gs.vac_solve().get_psi(normalized=False)

    def Bphi(self, R):
        """Toroidal field [T] at major radius R [m]."""
        return self.F0 / np.asarray(R)

    def structure_time_constants(self):
        """L/R eigenvalues [s] of the structure alone, slowest first."""
        nc = self.nc
        return np.sort(np.linalg.eigvals(np.linalg.solve(np.diag(self.R[nc:]), self.L[nc:, nc:])).real)[::-1]


# ============================================================================ circuit model
class CircuitModel:
    """Lumped circuit M dx/dt + R x = S V of the supplies, structure and plasma channel.

    State x = [I_circ (per circuit, per-turn A); I_seg (segments, A); I_p (A, after add_plasma)].
    Circuits: default one per coil, coils matching a substring in `exclude` left open; or
    circuits = {name: {coil: +1 (series) or -1 (anti-series), ...}}; coils in no circuit are open.
    V_max, V_min (default -V_max) [V] and extra series resistance R_circ [Ohm] are scalars,
    per-circuit arrays or {circuit: value} dicts (V_max must cover every circuit; circuits left
    out of V_min get -V_max, out of R_circ get 0).
    scheme 'zoh' (default): exact for voltages constant over each step (no time-step error, R = 0
    circuits need no special treatment); 'be': backward Euler (M + R dt) x_k = M x_{k-1} + S V_k dt,
    for comparison (over-damps eddy currents for steps of order the wall time).

    Attributes: M [H], Rm [Ohm] (diagonal), S supply input matrix, T [nc, ncirc] coil <- circuit,
    U [nu, nx] unit currents <- state, V_lo, V_hi [V] per circuit, V_max = max |V| per circuit,
    nx state size, ip index of I_p.
    """

    R_OPEN = 1.0e3          # plasma 'resistance' [Ohm] before the channel exists (effectively open)

    def __init__(self, machine, V_max, exclude=(), circuits=None, V_min=None, R_circ=0.0, scheme='zoh'):
        if scheme not in ('zoh', 'be'):
            raise ValueError(f"scheme must be 'zoh' or 'be', not {scheme!r}")
        m, nc, ns = machine, machine.nc, machine.ns
        if circuits is None:                             # one circuit per coil, except the excluded ones
            circuits = {}
            for n in m.coil_names:
                if not any(e in n for e in exclude):
                    circuits[n] = {n: 1.0}
        names = self.circ_names = list(circuits)
        self.ncirc = len(names)
        T = self.T = np.zeros((nc, self.ncirc))
        for j, name in enumerate(names):
            for coil, sign in circuits[name].items():
                T[m.coil_names.index(coil), j] = sign
        if isinstance(V_max, dict):
            missing = [n for n in names if n not in V_max]
            if missing:
                raise ValueError(f'no supply limit given for circuits {missing}')
            V_max = [V_max[n] for n in names]
        self.V_hi = np.full(self.ncirc, V_max, float)
        if isinstance(V_min, dict):
            V_min = [V_min.get(n, -v) for n, v in zip(names, self.V_hi)]
        self.V_lo = -self.V_hi if V_min is None else np.full(self.ncirc, V_min, float)
        self.V_max = np.maximum(np.abs(self.V_lo), np.abs(self.V_hi))
        if isinstance(R_circ, dict):
            R_circ = [R_circ.get(n, 0.0) for n in names]
        Lcc, Lcs, Lss = m.L[:nc, :nc], m.L[:nc, nc:], m.L[nc:, nc:]
        self.M = np.block([[T.T @ Lcc @ T, T.T @ Lcs], [Lcs.T @ T, Lss]])
        # circuit resistance: its coils' (whatever their winding sense) + R_circ
        R_c = np.abs(T).T @ m.R[:nc] + np.full(self.ncirc, R_circ, float)
        self.Rm = np.diag(np.concatenate([R_c, m.R[nc:]]))
        self.S = np.vstack([np.eye(self.ncirc), np.zeros((ns, self.ncirc))])
        self.U = np.block([[T, np.zeros((nc, ns))], [np.zeros((ns, self.ncirc)), np.eye(ns)]])
        self.nx, self.scheme, self.plasma = self.ncirc + ns, scheme, None
        self.Minv = np.linalg.inv(self.M)
        self._zoh = {}                                   # cache of step matrices

    def add_plasma(self, plasma):
        """Append the plasma channel current I_p as a state with no supply:
        L_p dI_p/dt + R_p I_p + sum_u M_pu dI_u/dt = 0, R_p set per step by zoh(dt, Rp) (default R_OPEN)."""
        M_px = plasma.M_u @ self.U
        self.M = np.block([[self.M, M_px[:, None]], [M_px, plasma.L_p]])
        self.Rm = np.diag(np.append(np.diag(self.Rm), self.R_OPEN))
        self.S, self.U = np.pad(self.S, ((0, 1), (0, 0))), np.pad(self.U, ((0, 0), (0, 1)))
        self.ip, self.plasma = self.nx, plasma
        self.nx += 1
        self.Minv = np.linalg.inv(self.M)
        self._zoh = {}

    def output_matrix(self, G_units, G_plasma=None):
        """Per-state output matrix C [n_out, nx] from a per-unit response G_units [n_out, nu] (and
        the per-ampere plasma response G_plasma [n_out])."""
        C = G_units @ self.U
        if G_plasma is not None and self.plasma is not None:
            C[:, self.ip] = G_plasma
        return C

    def zoh(self, dt, Rp=None, forcing=None):
        """Step x_k = A x_{k-1} + B V_k (+ Bw w_k) over dt [s], V (and w) constant on the step.

        Rp: plasma resistance [Ohm] on the step (None: R_OPEN). forcing: optional [nx, q] matrix
        through which q extra 'voltages' w enter (e.g. a prescribed plasma's flux linkage).
        Returns (A, B), or (A, B, Bw) with forcing; 'zoh': expm([[-M^-1 R, M^-1 [S, F]], [0, 0]] dt)
        = [[A, [B, Bw]], [0, I]]. Cached on (dt, Rp, shape and bytes of forcing).
        """
        key = (round(float(dt), 12), None if Rp is None else float(f'{Rp:.6e}'),
               None if forcing is None else (forcing.shape, hash(forcing.tobytes())))
        if key not in self._zoh:
            n, m = self.nx, self.ncirc
            Rm = self.Rm.copy()
            if Rp is not None and self.plasma is not None:
                Rm[self.ip, self.ip] = Rp
            B_in = self.S if forcing is None else np.hstack([self.S, forcing])
            if self.scheme == 'be':
                Mi = np.linalg.inv(self.M + Rm * dt)
                A, B = Mi @ self.M, Mi @ B_in * dt
            else:
                Z = np.zeros((n + B_in.shape[1],) * 2)
                Z[:n, :n], Z[:n, n:] = -self.Minv @ Rm * dt, self.Minv @ B_in * dt
                E = expm(Z)
                A, B = E[:n, :n], E[:n, n:]
            self._zoh[key] = (A, B[:, :m], B[:, m:])[:2 if forcing is None else 3]
        return self._zoh[key]


def feed_forward(machine, t, I_coil, plasma=None, Ip=None, Rp=None):
    """Feed-forward circuit (validation against a discharge): given coil currents I_coil [nt, nc]
    (per turn, linear in time between the times t [s]) drive the structure,
    M_ss dI_s/dt + R_s I_s = -M_sc dI_c/dt, stepped exactly (a CircuitModel with no supply
    circuits; the coil current rates enter through its forcing). The structure starts at rest.

    Plasma channel (a ChannelPlasma), optional: with Ip [nt] [A] its current is prescribed (replay:
    its changing flux linkage also drives the structure); with Rp [nt-1] [Ohm], its resistance on each
    step, it is a circuit state, L_p dI_p/dt + R_p I_p + sum_u M_pu dI_u/dt = 0 (prediction).

    Returns X [nt, nc + ns + 1]: coil (as given), segment and channel currents [A]."""
    nc = machine.nc
    # the unknowns are the segment currents only: a circuit model with no supply circuits
    circ = CircuitModel(machine, 1.0, circuits={})
    M_xc = machine.L[nc:, :nc]                            # state <- coil mutual inductances [H]
    W = np.diff(I_coil, axis=0) / np.diff(t)[:, None]     # forcing: coil current rates [A/s], constant per step
    if Rp is not None:
        # prediction: the channel is one more state, also driven by the coil current rates
        circ.add_plasma(plasma)
        M_xc = np.vstack([M_xc, plasma.M_u[:nc]])
    elif Ip is not None:
        # replay: the channel current is known, so its rate is one more forcing of the segments
        M_xc = np.hstack([M_xc, plasma.M_u[nc:, None]])
        W = np.hstack([W, (np.diff(Ip) / np.diff(t))[:, None]])
    x = np.zeros((len(t), circ.nx))                       # structure at rest at t[0]
    for k in range(len(t) - 1):
        # exact step of M dx/dt + R x = -M_xc w with the forcing w constant on the step
        A, _, Bw = circ.zoh(t[k + 1] - t[k], None if Rp is None else Rp[k], forcing=-M_xc)
        x[k + 1] = A @ x[k] + Bw @ W[k]
    if Rp is not None:                                    # predicted: the last state
        I_p = x[:, -1]
    elif Ip is not None:                                  # replayed
        I_p = Ip
    else:                                                 # no channel
        I_p = np.zeros(len(t))
    return np.column_stack([I_coil, x[:, :machine.ns], I_p])


# ============================================================================ early plasma channel
class ChannelPlasma:
    """The early plasma as a toroidal ring of uniform current, minor radius a [m], centred at
    (Rc, Zc) [m]; li internal inductance, Zeff effective charge, lnL Coulomb logarithm.

        L_p = mu0 Rc (ln(8 Rc/a) - 2 + l_i/2)   self-inductance [H]
        M_u = 2 pi <psi_u>_channel              mutual inductance with every coil and segment [H]
        R_p = eta 2 pi Rc / (pi a^2)            resistance [Ohm], with the parallel Spitzer
        eta = 5.2e-5 Zeff lnL / Te^1.5          resistivity [Ohm m, Te in eV] (NRL Formulary)
    """

    def __init__(self, machine, Rc, Zc, a, li=0.5, Zeff=2.0, lnL=LN_LAMBDA):
        self.Rc, self.Zc, self.a, self.li, self.Zeff, self.lnL = Rc, Zc, a, li, Zeff, lnL
        self.qp, self.qw = _disk_quadrature(Rc, Zc, a, n_r=8, n_t=24)
        self.M_u = 2 * np.pi * eval_fields(machine.gs, machine.psi_u, self.qp)['psi'] @ self.qw
        self.L_p = MU0 * Rc * (np.log(8 * Rc / a) - 2 + li / 2)

    def resistance(self, Te_eV):
        """Channel resistance [Ohm] at electron temperature Te_eV [eV] (floored at 0.1 eV)."""
        eta = 5.2e-5 * self.Zeff * self.lnL / np.maximum(Te_eV, 0.1) ** 1.5
        return eta * 2 * np.pi * self.Rc / (np.pi * self.a ** 2)

    def field_at(self, pts, h=1e-3):
        """psi [Wb/rad], B_R, B_Z [T] per ampere of channel current at pts [n_pts, 2]; gradient by
        central difference with step h [m]."""
        pts = np.atleast_2d(pts)
        dpsi_dR = (self.psi_unit(pts + [h, 0]) - self.psi_unit(pts - [h, 0])) / (2 * h)
        dpsi_dZ = (self.psi_unit(pts + [0, h]) - self.psi_unit(pts - [0, h])) / (2 * h)
        return self.psi_unit(pts), -dpsi_dZ / pts[:, 0], dpsi_dR / pts[:, 0]

    def psi_unit(self, x):
        """Flux [Wb/rad] per ampere of channel current at points x [n, 2]: a sum of filaments
        (eval_green), one per quadrature cell, whose log singularity inside the cell is replaced by
        the flux of a straight uniform wire of the cell's area (finite, smooth inside the channel)."""
        from OpenFUSIONToolkit.TokaMaker.util import eval_green
        x = np.atleast_2d(np.asarray(x, float))
        q0 = self.qp[0]
        if not hasattr(self, '_cR'):
            # near a filament psi ~ -cR ln(d): fit cR at the first filament (cR scales with its R)
            d1, d2 = 1e-4, 1e-3
            g = -eval_green(np.array([[q0[0] + d1, q0[1]], [q0[0] + d2, q0[1]]]), q0)
            self._cR = -(g[1] - g[0]) / np.log(d2 / d1)
        total = np.zeros(len(x))
        for q, w, rc in zip(self.qp, self.qw, self.a * np.sqrt(self.qw)):   # rc: radius of the cell's area
            d = np.hypot(x[:, 0] - q[0], x[:, 1] - q[1])
            f = -eval_green(x, q)
            inside = d < rc
            if inside.any():                             # uniform-wire core
                di = np.maximum(d[inside], 1e-12)
                f[inside] += -(self._cR * q[0] / q0[0]) * (np.log(rc) - np.log(di)
                                                           + 0.5 * (di ** 2 / rc ** 2 - 1))
            total += w * f
        return total


def Te_burnthrough(t, t_bd, Te0=5.0, Te1=100.0, tau=0.1):
    """Prescribed channel Te [eV] at times t [s]: Te0 until the end of breakdown t_bd [s], then
    rising towards Te1 with time constant tau [s] (placeholder for a burn-through model)."""
    t = np.asarray(t, float)
    return np.where(t < t_bd, Te0, Te0 + (Te1 - Te0) * (1 - np.exp(-(t - t_bd) / tau)))


# ============================================================================ breakdown physics
TORR = 133.322          # 1 Torr in Pa
# First Townsend coefficient of H2, alpha = A p exp(-B p / E) with p in Torr; close to Raizer, Gas
# Discharge Physics (1991): A = 5 /(cm Torr), B = 130 V/(cm Torr).
TOWNSEND_A = 510.0      # 1/(m Torr)
TOWNSEND_B = 1.25e4     # V/(m Torr)


class FieldTracer:
    """Traces field lines of nodal flux psi_nodal [Wb/rad] plus B_phi = |F0|/R (F0 [T m]) from a
    point until they leave the wall polygon [n, >=2] (R, Z) [m]. Bp_min [T] floors B_pol so the
    step stays finite near nulls."""

    def __init__(self, gs, psi_nodal, F0, wall, Bp_min=1e-6):
        self.fg = gs.get_field_eval('grad', values=np.ascontiguousarray(psi_nodal))
        self.F0, self.wall, self.Bp_min = abs(F0), Path(wall[:, :2]), Bp_min

    def _rk4(self, x, sgn, dl):
        """One RK4 step of poloidal length dl [m] along sgn * B_pol: (new position, length along B [m])."""
        k, w = [], []
        for c in (0.0, 0.5, 0.5, 1.0):
            y = x + c * dl * k[-1] if k else x
            g = self.fg.eval(y)[:2]
            br, bz, bt = -g[1] / y[0], g[0] / y[0], self.F0 / y[0]
            bp = max(np.hypot(br, bz), self.Bp_min)
            k.append(np.array([sgn * br / bp, sgn * bz / bp]))
            w.append(np.sqrt(bp ** 2 + bt ** 2) / bp)    # dL / dl_pol = |B| / B_pol
        return x + dl / 6 * (k[0] + 2 * k[1] + 2 * k[2] + k[3]), dl / 6 * (w[0] + 2 * w[1] + 2 * w[2] + w[3])

    def trace(self, x0, sgn, dl=0.01, L_max=1e5, n_max=2000):
        """Length [m] along B from x0 to the wall in direction sgn (+1/-1), poloidal step dl [m];
        L_max if longer or not out after n_max steps."""
        x, L = np.array(x0, float), 0.0
        for _ in range(n_max):
            x_new, dL = self._rk4(x, sgn, dl)
            if not self.wall.contains_point(x_new):      # bisect the wall-crossing step
                lo, hi = 0.0, dl
                for _ in range(20):
                    mid = 0.5 * (lo + hi)
                    if self.wall.contains_point(self._rk4(x, sgn, mid)[0]):
                        lo = mid
                    else:
                        hi = mid
                return min(L + self._rk4(x, sgn, lo)[1], L_max)
            L += dL
            if L > L_max:
                return L_max
            x = x_new
        return L_max

    def connection_length(self, x0, **kw):
        """Connection length [m] through x0: traces in both directions."""
        return self.trace(x0, +1, **kw) + self.trace(x0, -1, **kw)


def townsend_window(E, L, A=TOWNSEND_A, B=TOWNSEND_B):
    """(p_min, p_max) [Pa] for which an avalanche runs away (alpha L >= 1) at field E [V/m] over
    connection length L [m], scanned on a log grid 1e-7..1e3 Pa; None if there is none."""
    p = np.logspace(-7, 3, 10000)
    p_torr = p / TORR
    ok = A * p_torr * L * np.exp(-B * p_torr / max(E, 1e-9)) >= 1.0
    return (float(p[ok].min()), float(p[ok].max())) if ok.any() else None


def paschen_Emin(p, L, A=TOWNSEND_A, B=TOWNSEND_B):
    """Minimum E_phi [V/m] for an avalanche (alpha L = 1) at pressure p [Pa] over length L [m];
    inf where A p L <= 1."""
    p_torr = np.asarray(p) / TORR
    x = A * p_torr * L
    return np.where(x > 1, B * p_torr / np.log(np.maximum(x, 1 + 1e-12)), np.inf)
