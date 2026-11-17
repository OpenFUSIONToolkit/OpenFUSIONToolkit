#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
"""Optimisation of a tokamak start-up: breakdown and early plasma (stage 1), then a ramp of
TokaMaker free-boundary equilibria (stage 2), optimised together.

Contents
    ConicProblem    a small sparse conic-program builder on the Clarabel solver
    time_grid       the stage-1 time grid
    StartupProblem  stage 1 as ONE convex program: premagnetisation currents and every supply
                    voltage, with the breakdown criteria, the coil/supply limits and the early
                    plasma channel's force balance, stability and current target as constraints
    GSRamp          stage 2: TokaMaker equilibria coupled to the same circuit model by Picard
                    iteration; each convex program also re-optimises stage 1 (joint optimisation)
"""
import numpy as np
import scipy.sparse as sp
from scipy.linalg import expm

from startup_model import (eval_fields, disk_points, cell_geometry, ChannelPlasma, Te_burnthrough,
                           LN_LAMBDA, MU0)


# ============================================================================ conic programs
def solved(status):
    return status in ('Solved', 'AlmostSolved')


class Var:
    def __init__(self, name, start, n, scale):
        self.name, self.start, self.n = name, start, n
        self.scale = np.broadcast_to(np.asarray(scale, float), (n,)).copy()

    @property
    def sl(self):
        return slice(self.start, self.start + self.n)


class ConicProblem:
    """    minimise   1/2 sum_i w_i^2 ||A_i x - b_i||^2 + q^T x
           subject to linear equalities, linear inequalities and second-order cones.

    Variables are declared in named blocks with a typical size (scale); constraints are written in
    physical units as dicts {Var: matrix} and the builder rescales rows and columns before calling
    Clarabel. Least-squares terms become extra residual variables with a diagonal quadratic cost."""

    def __init__(self):
        self.vars, self.nvar = {}, 0
        self._eq, self._ineq, self._soc, self._quad, self._lin = [], [], [], [], []

    def clone(self, objective=True):
        """Copy with the same variables and terms, to which more can be added (objective=False:
        without the linear objective terms)."""
        P = ConicProblem()
        P.vars, P.nvar = dict(self.vars), self.nvar
        P._eq, P._ineq, P._soc, P._quad = list(self._eq), list(self._ineq), list(self._soc), list(self._quad)
        P._lin = list(self._lin) if objective else []
        return P

    def var(self, name, n, scale=1.0):
        v = Var(name, self.nvar, int(n), scale)
        self.vars[name] = v
        self.nvar += v.n
        return v

    @staticmethod
    def _m(terms):
        m = None
        for v, M in terms.items():
            M = sp.csr_matrix(M)
            m = M.shape[0] if m is None else m
            assert M.shape == (m, v.n), (v.name, M.shape, (m, v.n))
        return m

    def eq(self, terms, rhs, tag=''):
        """sum_v M_v v = rhs"""
        self._eq.append((terms, np.broadcast_to(np.asarray(rhs, float), (self._m(terms),)).copy(), tag))

    def le(self, terms, rhs, tag=''):
        """sum_v M_v v <= rhs"""
        self._ineq.append((terms, np.broadcast_to(np.asarray(rhs, float), (self._m(terms),)).copy(), tag))

    def soc(self, t_terms, t0, v_terms, v0, tag=''):
        """|| sum M_v v + v0 ||_2 <= sum a_v v + t0"""
        mv = self._m(v_terms)
        self._soc.append((t_terms, float(t0), v_terms, np.broadcast_to(np.asarray(v0, float), (mv,)).copy(), tag))

    def lsq(self, terms, target, weight, tag=''):
        """adds (weight^2 / 2) || sum M_v v - target ||^2 to the objective"""
        m = self._m(terms)
        self._quad.append((terms, np.broadcast_to(np.asarray(target, float), (m,)).copy(), float(weight), tag))

    def linear(self, var, coeff):
        """adds coeff . var to the objective"""
        self._lin.append((var, np.broadcast_to(np.asarray(coeff, float), (var.n,)).copy()))

    def _row_block(self, terms, nrows):
        out = sp.csr_matrix((nrows, self.nvar))
        for v, M in terms.items():
            M = sp.csr_matrix(M) @ sp.diags(v.scale)                 # column scaling
            out = out + sp.hstack([sp.csr_matrix((nrows, v.start)), M,
                                   sp.csr_matrix((nrows, self.nvar - v.start - v.n))], format='csr')
        out = out.tocsr()                                # drop numerical dust below 1e-10 x row max
        if out.nnz:
            rmax = np.asarray(abs(out).max(axis=1).todense()).ravel()
            rows = np.repeat(np.arange(out.shape[0]), np.diff(out.indptr))
            out.data[np.abs(out.data) < 1e-10 * rmax[rows]] = 0.0
            out.eliminate_zeros()
        return out

    def solve(self):
        import clarabel
        nres = sum(q[1].size for q in self._quad)
        ntot = self.nvar + nres
        pad = lambda M: sp.hstack([M, sp.csr_matrix((M.shape[0], nres))], format='csr')
        rows_A, rows_b, cones = [], [], []

        def add_rows(A, b, kind, tag):                   # drop empty rows, then row-normalise
            rm = np.asarray(abs(A).max(axis=1).todense()).ravel()
            keep = rm != 0
            bad = (np.abs(b[~keep]) > 0) if kind == 'eq' else (b[~keep] < 0)
            if np.any(bad):
                raise ValueError(f'constraint "{tag}" has empty rows that are infeasible')
            rows_A.append(sp.diags(1 / rm[keep]) @ A[keep]); rows_b.append(b[keep] / rm[keep])
            return int(keep.sum())

        P_diag = np.zeros(ntot)
        r0, nz = self.nvar, 0
        for terms, target, w, tag in self._quad:         # residuals r = w (A x - target)
            m = target.size
            R = sp.csr_matrix((np.full(m, -1.0), (np.arange(m), r0 + np.arange(m) - self.nvar)), shape=(m, nres))
            nz += add_rows(sp.hstack([self._row_block(terms, m) * w, R], format='csr'), w * target, 'eq', tag)
            P_diag[r0:r0 + m] = 1.0
            r0 += m
        for terms, rhs, tag in self._eq:
            nz += add_rows(pad(self._row_block(terms, rhs.size)), rhs, 'eq', tag)
        if nz:
            cones.append(clarabel.ZeroConeT(nz))
        nn = 0
        for terms, rhs, tag in self._ineq:
            nn += add_rows(pad(self._row_block(terms, rhs.size)), rhs, 'le', tag)
        if nn:
            cones.append(clarabel.NonnegativeConeT(nn))
        for t_terms, t0, v_terms, v0, tag in self._soc:
            A = sp.vstack([-pad(self._row_block(t_terms, 1)), pad(self._row_block(v_terms, v0.size))], format='csr')
            b = np.concatenate([[t0], -v0])
            sc = max(abs(A).max(), 1e-300)
            rows_A.append(A / sc); rows_b.append(b / sc)
            cones.append(clarabel.SecondOrderConeT(A.shape[0]))
        A = sp.vstack(rows_A, format='csc')
        b = np.concatenate(rows_b)
        q = np.zeros(ntot)
        for v, c in self._lin:
            q[v.sl] += c * v.scale
        Pm = sp.triu(sp.diags(P_diag, format='csc'), format='csc')
        st = clarabel.DefaultSettings()
        st.verbose = False
        sol = clarabel.DefaultSolver(Pm, q, A, b, cones, st).solve()
        if 'NumericalError' in str(sol.status) or 'InsufficientProgress' in str(sol.status):
            # retry with a much wider equilibration range
            for k, val in {'equilibrate_max_iter': 50, 'equilibrate_min_scaling': 1e-8,
                           'equilibrate_max_scaling': 1e8, 'max_iter': 400}.items():
                setattr(st, k, val)
            sol = clarabel.DefaultSolver(Pm, q, A, b, cones, st).solve()
        self.status, self._x = str(sol.status), np.array(sol.x)
        return self.status

    def value(self, var):
        return self._x[var.sl] * var.scale

    def linear_terms(self):
        """The linear objective as {Var: coefficient row [1, n]}, in physical units."""
        lin = {}
        for v, q in self._lin:
            lin[v] = lin.get(v, 0) + q[None, :]
        return lin


# ============================================================================ stage 1
def time_grid(tc):
    """t_start --dt_pre--> t_fine_start --dt_fine--> t_fine_end --dt_ramp--> t_end"""
    def seg(a, b, d):
        if b - a <= 1e-12:
            return np.zeros(0)
        return np.linspace(a, b, max(1, int(np.ceil((b - a) / d - 1e-9))) + 1)[1:]
    t0, tf, tfe, t1 = tc['t_start'], tc['t_fine_start'], tc['t_fine_end'], tc['t_end']
    return np.concatenate([[t0], seg(t0, tf, tc['dt_pre']), seg(tf, tfe, tc['dt_fine']), seg(tfe, t1, tc['dt_ramp'])])


class StartupProblem:
    """Stage 1: premagnetisation -> breakdown -> early plasma channel up to the hand-off.

    Decision variables: the premagnetisation current I0 of every driven circuit (structure at rest
    at t_start) and the supply voltages V_1..V_N, constant on each interval (t_{k-1}, t_k].
    Because the circuit is linear and stepped exactly, every state is an affine map of them,
        x_k = Gi_k I0 + Gv_k V,
    so every coil current, eddy current, field and electric field is linear in (I0, V) and the
    whole design is one convex program. E_phi is the interval average -(psi_k - psi_{k-1})/(R dt).
    """

    def __init__(self, machine, circuit, scen):
        self.m, self.c, self.scen = machine, circuit, scen
        self.s = 1.0                                     # sign of I_p and E_phi: positive I_p only
        self.t = time_grid(scen['time'])
        rg = scen['region']
        self.pts = disk_points(rg['R'], rg['Z'], rg['radius'], rg['n_rings'], rg['n_per_ring'])
        self.ref = np.array([[rg['R'], rg['Z']]])       # flux objective point: the region centre
        pl = scen['plasma']
        self.plasma = ChannelPlasma(machine, pl['Rc'], pl['Zc'], pl['a'], li=pl['li'], Zeff=pl['Zeff'])
        circuit.add_plasma(self.plasma)
        np_ = len(self.pts)
        G = eval_fields(machine.gs, machine.psi_u, np.vstack([self.pts, self.ref]))
        gp = self.plasma.field_at(self.pts)              # the channel's own field at the region
        self.Cpsi = circuit.output_matrix(G['psi'][:, :np_].T, gp[0])
        self.CBr = circuit.output_matrix(G['Br'][:, :np_].T, gp[1])
        self.CBz = circuit.output_matrix(G['Bz'][:, :np_].T, gp[2])
        self.Cref = circuit.output_matrix(G['psi'][:, np_:].T)          # external flux only
        Gc = eval_fields(machine.gs, machine.psi_u, np.array([[pl['Rc'], pl['Zc']]]), second=True)
        self.Cbz_c = circuit.output_matrix(Gc['Bz'].T)   # external fields at the channel centre
        self.Cbr_c = circuit.output_matrix(Gc['Br'].T)
        self.Cdbz_c = circuit.output_matrix(Gc['dBz_dR'].T)
        self.Cpsi_c = circuit.output_matrix(Gc['psi'].T)
        self.Rp = self.pts[:, 0]
        self.Bphi = np.abs(machine.Bphi(self.Rp))

    def build(self):
        c, m, scen, t = self.c, self.m, self.scen, self.t
        N, nx, nv = len(t) - 1, c.nx, c.ncirc
        Ilim = np.abs(m.coil_lim).max(axis=1)
        Iscale = np.array([Ilim[np.abs(c.U[:m.nc, j]) > 0].max() for j in range(nv)])
        P = ConicProblem()
        I0 = P.var('I0', nv, Iscale)
        V = P.var('V', N * nv, np.tile(c.V_max, N))
        self.P, self.I0v, self.V = P, I0, V

        # --- plasma resistance per interval: open circuit before t_p0, Spitzer(Te(t)) after
        pl = scen['plasma']
        self.t_p0 = scen['breakdown']['window'][1]       # the channel exists after the window
        tm = 0.5 * (t[1:] + t[:-1])
        self.Te = Te_burnthrough(tm, self.t_p0, pl['Te0'], pl['Te1'], pl['tau_Te'])
        self.Rplas = np.where(t[:-1] >= self.t_p0 - 1e-12, self.plasma.resistance(self.Te), c.R_OPEN)

        # --- propagate the affine state maps x_k = Gi_k I0 + Gv_k V
        Gi = np.zeros((nx, nv)); Gi[:nv] = np.eye(nv)
        Gv = np.zeros((nx, N * nv))
        self.Gi, self.Gv = [Gi.copy()], [Gv.copy()]
        for k in range(1, N + 1):
            A, B = c.zoh(t[k] - t[k - 1], self.Rplas[k - 1])
            Gi, Gv = A @ Gi, A @ Gv
            Gv[:, (k - 1) * nv:k * nv] += B
            self.Gi.append(Gi.copy()); self.Gv.append(Gv.copy())

        # --- coil current limits at every slice, supply voltage limits on every interval
        Tu = c.U[:m.nc][np.abs(c.T).sum(axis=1) > 0]
        used = np.abs(c.T).sum(axis=1) > 0
        lo, hi = m.coil_lim[used, 0], m.coil_lim[used, 1]
        Ii = np.vstack([Tu @ self.Gi[k] for k in range(N + 1)])
        Iv = np.vstack([Tu @ self.Gv[k] for k in range(N + 1)])
        P.le({I0: Ii, V: Iv}, np.tile(hi, N + 1), tag='I_max')
        P.le({I0: -Ii, V: -Iv}, -np.tile(lo, N + 1), tag='I_min')
        P.le({V: sp.eye(N * nv)}, np.tile(c.V_max, N), tag='V_max')
        P.le({V: -sp.eye(N * nv)}, np.tile(c.V_max, N), tag='V_min')

        # --- breakdown window: E_phi floor and the Lloyd criterion at every region point
        bd = scen['breakdown']
        w0, w1 = bd['window']
        self.K = [k for k in range(1, N + 1) if t[k - 1] >= w0 - 1e-12 and t[k] <= w1 + 1e-12]
        np_ = len(self.pts)
        for k in self.K:
            Ck = -self.s * self.Cpsi / self.Rp[:, None] / (t[k] - t[k - 1])     # s E_phi on interval k
            Ei, Ev = Ck @ (self.Gi[k] - self.Gi[k - 1]), Ck @ (self.Gv[k] - self.Gv[k - 1])
            P.le({I0: -Ei, V: -Ev}, -bd['E_min'] * np.ones(np_), tag='E_floor')
            # ||B_p|| <= B_phi s E_phi / (m E*), against B at BOTH ends of the interval
            for j in (k - 1, k):
                Bri, Brv = self.CBr @ self.Gi[j], self.CBr @ self.Gv[j]
                Bzi, Bzv = self.CBz @ self.Gi[j], self.CBz @ self.Gv[j]
                for p in range(np_):
                    a = self.Bphi[p] / (bd['margin'] * bd['E_star'])
                    P.soc({I0: a * Ei[p:p + 1], V: a * Ev[p:p + 1]}, 0.0,
                          {I0: np.vstack([Bri[p], Bzi[p]]), V: np.vstack([Brv[p], Bzv[p]])}, [0, 0], tag='lloyd')

        # --- early plasma channel
        eip = np.zeros((1, nx)); eip[0, c.ip] = 1.0
        Kp = [k for k in range(1, N + 1) if t[k] > self.t_p0 + 1e-12]
        for k in Kp:                                     # I_p never falls once it exists
            P.le({I0: -self.s * eip @ (self.Gi[k] - self.Gi[k - 1]),
                  V: -self.s * eip @ (self.Gv[k] - self.Gv[k - 1])}, [0.0], tag='Ip_mono')
        P.eq({I0: self.s * eip @ self.Gi[N], V: self.s * eip @ self.Gv[N]}, [pl['Ip_target']], tag='Ip_target')
        Rc, a = pl['Rc'], pl['a']
        # Shafranov: B_v,req = -kf I_p, kf = mu0/(4 pi Rc) (ln(8 Rc/a) + beta_p + l_i/2 - 3/2)
        kf = MU0 / (4 * np.pi * Rc) * (np.log(8 * Rc / a) + pl.get('beta_p', 0.0) + pl['li'] / 2 - 1.5)
        fb = pl['force_balance']
        nlo, nhi = pl['decay_index']
        for k in Kp:
            Gi_, Gv_ = self.Gi[k], self.Gv[k]
            ipi, ipv = eip @ Gi_, eip @ Gv_
            if t[k] <= pl['t_burn_end'] + 1e-12:        # external E_phi at the channel during burn-through
                ce = -self.s * self.Cpsi_c / (Rc * (t[k] - t[k - 1]))
                P.le({I0: -ce @ (Gi_ - self.Gi[k - 1]), V: -ce @ (Gv_ - self.Gv[k - 1])}, [-pl['E_min_burn']], tag='E_burn')
            # radial force balance |B_z,ext + kf I_p| <= frac kf |I_p| + B_tol, and |B_R,ext| likewise
            frac, Btol = fb['frac'], fb['B_tol']
            ei, ev = self.Cbz_c @ Gi_ + kf * ipi, self.Cbz_c @ Gv_ + kf * ipv
            P.le({I0: ei - frac * kf * self.s * ipi, V: ev - frac * kf * self.s * ipv}, [Btol], tag='Bv_hi')
            P.le({I0: -ei - frac * kf * self.s * ipi, V: -ev - frac * kf * self.s * ipv}, [Btol], tag='Bv_lo')
            bi, bv = self.Cbr_c @ Gi_, self.Cbr_c @ Gv_
            P.le({I0: bi - frac * kf * self.s * ipi, V: bv - frac * kf * self.s * ipv}, [Btol], tag='Br_hi')
            P.le({I0: -bi - frac * kf * self.s * ipi, V: -bv - frac * kf * self.s * ipv}, [Btol], tag='Br_lo')
            # decay index n = -(Rc/B_z) dB_z/dR in [nlo, nhi], multiplied through by B_z (sign -s)
            qi, qv = -Rc * self.Cdbz_c @ Gi_, -Rc * self.Cdbz_c @ Gv_
            zi, zv = self.Cbz_c @ Gi_, self.Cbz_c @ Gv_
            lo_, hi_ = (nhi, nlo) if self.s > 0 else (nlo, nhi)
            P.le({I0: lo_ * zi - qi, V: lo_ * zv - qv}, [0.0], tag='n_lo')
            P.le({I0: qi - hi_ * zi, V: qv - hi_ * zv}, [0.0], tag='n_hi')

        # --- weak regularisation: smooth voltages (dV) and small voltages (V), relative to V_max
        reg = scen['regularization']
        Dv = sp.kron(sp.diags([-np.ones(N - 1), np.ones(N - 1)], [0, 1], shape=(N - 1, N)), sp.diags(1.0 / c.V_max))
        P.lsq({V: Dv}, 0.0, reg['dV'], tag='reg_dV')
        P.lsq({V: sp.kron(sp.eye(N), sp.diags(1.0 / c.V_max))}, 0.0, reg['V'], tag='reg_V')

        # --- objective: maximise the flux s psi(t_N) left at the region centre
        for v_, M_ in self.handoff_terms().items():
            P.linear(v_, -self.s * (self.Cref @ M_).ravel())
        return P

    def solve(self):
        st = self.P.solve()
        if not solved(st):
            raise RuntimeError(f'stage-1 convex program not solved: {st}')
        self.load(self.P)
        return st

    def handoff_terms(self):
        """The hand-off (last-slice) state as {variable: matrix}: x_N = Gi_N I0 + Gv_N V."""
        return {self.I0v: self.Gi[-1], self.V: self.Gv[-1]}

    def flux_handoff(self):
        """Flux left for the pulse at the hand-off, 2 pi s psi_ref(t_N) [Wb]."""
        return float(2 * np.pi * self.s * (self.Cref @ self.x[-1])[0])

    def load(self, P):
        """Read I0 and V from a solved problem containing this stage's variables (its own, or the
        joint stage-1 + ramp problem) and rebuild the state trajectory."""
        N, nv = len(self.t) - 1, self.c.ncirc
        i0, vv = P.value(self.I0v), P.value(self.V)
        self.v = vv.reshape(N, nv)
        self.x = np.array([self.Gi[k] @ i0 + self.Gv[k] @ vv for k in range(N + 1)])

    def results(self):
        """Trajectories: currents, voltages, region fields, Lloyd quantity, channel quantities."""
        c, m, t, x = self.c, self.m, self.t, self.x
        xdot = np.zeros_like(x)
        xdot[1:] = np.diff(x, axis=0) / np.diff(t)[:, None]           # interval averages
        E = -(self.Cpsi @ xdot.T) / self.Rp[:, None]                   # [np, N+1]
        B = np.hypot(self.CBr @ x.T, self.CBz @ x.T)
        return {'t': t, 'x': x, 'xdot': xdot, 'V': self.v, 'I_coil': (c.T @ x[:, :c.ncirc].T).T,
                'I_struct': x[:, c.ncirc:c.ncirc + m.ns], 'Ip': x[:, c.ip], 'Te': self.Te,
                'flux_ref_Wb': 2 * np.pi * self.s * (self.Cref @ x.T).ravel(),
                'E_pts': (self.s * E).T, 'B_pts': B.T,
                'lloyd_pts': (self.s * E * self.Bphi[:, None] / np.maximum(B, 1e-12)).T,
                'Bz_ext_c': (self.Cbz_c @ x.T).ravel(), 'Br_ext_c': (self.Cbr_c @ x.T).ravel(),
                'n_c': -self.scen['plasma']['Rc'] * (self.Cdbz_c @ x.T).ravel()
                       / np.where(np.abs(self.Cbz_c @ x.T).ravel() > 0, (self.Cbz_c @ x.T).ravel(), np.nan),
                'K': self.K}


# ============================================================================ stage 2
def interp_keys(keys, t):
    """[[t, value, ...], ...] keyframes -> values at the times t (linear interpolation)."""
    k = np.asarray(keys, float)
    return np.array([np.interp(t, k[:, 0], k[:, j]) for j in range(1, k.shape[1])]).T


class UnitSampler:
    """Flux linkage of a nodal flux with every coil and segment: Phi_u = 2 pi n_u <psi>_u, averaged
    over the unit's cells with its current distribution (uniform for coils, 1/R for segments)."""

    def __init__(self, machine):
        P, LC, REG = machine.mesh
        area, cen = cell_geometry(P, LC)
        pts, rows, cols, vals, npt = [], [], [], [], 0
        groups = [(np.where(REG == machine.coil_dict[n]['reg_id'])[0], machine.nturns[i], False)
                  for i, n in enumerate(machine.coil_names)]
        groups += [(s['cells'], 1.0, True) for s in machine.segments]
        for i, (cells, nt, inv_r) in enumerate(groups):
            w = area[cells] / cen[cells, 0] if inv_r else area[cells].copy()
            pts.append(cen[cells]); rows += [i] * len(cells)
            cols += list(range(npt, npt + len(cells))); vals += list(nt * w / w.sum()); npt += len(cells)
        self.pts = np.vstack(pts)
        self.S = sp.csr_matrix((vals, (rows, cols)), shape=(len(groups), npt))

    def linkage(self, gs, psi_nodal):
        return 2 * np.pi * (self.S @ eval_fields(gs, psi_nodal[None, :], self.pts)['psi'][0])


class GSRamp:
    """Stage 2: a ramp of TokaMaker free-boundary equilibria from the stage-1 hand-off, with
    I_p(t) and the target boundary shape prescribed per time slice, optimised jointly with stage 1.

    Picard iteration:
      1. TokaMaker solves every slice for its target shape, with the coil and structure currents
         pinned to the current trajectory -> the plasma's own flux psi_p per slice.
      2. One convex program, plasma flux held fixed: all ramp voltages AND all stage-1 variables.
         The plasma enters the circuit as a known source (its flux linkage with every unit
         changes in time, inducing currents) and the isoflux shape targets as least-squares
         terms; a 0-D Ohm's law sets the boundary flux. Stage 1 keeps all its constraints, and
         its flux objective becomes a floor (within flux_slack_Wb of its stand-alone optimum).
      3. Repeat until the currents, TokaMaker's response and the plasma flux stop changing.
    """

    def __init__(self, machine, circuit, scen, stage1, log=print):
        """circuit: a plasma-free CircuitModel with the same circuits as stage 1.
        stage1: the SOLVED StartupProblem the ramp starts from (and re-optimises)."""
        self.m, self.c, self.scen, self.log = machine, circuit, scen, log
        rp = scen['ramp']
        tc = rp['time']
        n = int(round((tc['t_end'] - tc['t_start']) / tc['dt']))
        self.t = np.linspace(tc['t_start'], tc['t_end'], n + 1)
        self.Ip = interp_keys(rp['Ip'], self.t).ravel()
        self.shape = interp_keys(rp['shape'], self.t)            # R0, Z0, a, kappa, delta
        self.s = 1.0                                             # positive I_p only, as in stage 1
        if abs(stage1.t[-1] - self.t[0]) > 1e-9 or abs(stage1.s * stage1.x[-1, stage1.c.ip] - self.Ip[0]) > 0.02 * self.Ip[0]:
            raise ValueError('stage 1 must end where the ramp starts (same time and I_p)')
        self.sampler = UnitSampler(machine)
        from OpenFUSIONToolkit.TokaMaker.util import create_isoflux
        self.iso = [create_isoflux(rp['n_iso'], *sh) for sh in self.shape]
        # flux at each slice's target points per unit current: fixed, so computed once
        self.G_iso = [circuit.output_matrix(eval_fields(machine.gs, machine.psi_u, iso)['psi'].T) for iso in self.iso]
        self.Bu = circuit.U.T                    # plasma forcing enters as unit 'voltages'
        self._zoh = {}
        # joint optimisation with stage 1
        self.stage1 = stage1
        F_star = stage1.flux_handoff()
        lin = stage1.P.linear_terms()
        L_star = float(sum((q @ stage1.P.value(v)).sum() for v, q in lin.items()))
        self.joint = {'P0': stage1.P.clone(objective=False), 'flux_star_Wb': F_star, 'lin': lin,
                      'L_cap': L_star + rp['joint']['flux_slack_Wb'] / (2 * np.pi),
                      'X0': {v: M[:circuit.nx] for v, M in stage1.handoff_terms().items()}}
        self.x0 = stage1.x[-1][:circuit.nx].copy()

    def zoh(self, dt):
        """Exact step with voltages V and unit 'voltages' w (plasma forcing) held constant."""
        key = round(float(dt), 12)
        if key not in self._zoh:
            c = self.c
            n, m, q = c.nx, c.ncirc, self.Bu.shape[1]
            Z = np.zeros((n + m + q, n + m + q))
            Z[:n, :n] = -c.Minv @ c.Rm * dt
            Z[:n, n:n + m] = c.Minv @ c.S * dt
            Z[:n, n + m:] = c.Minv @ self.Bu * dt
            E = expm(Z)
            self._zoh[key] = (E[:n, :n], E[:n, n:n + m], E[:n, n + m:])
        return self._zoh[key]

    # ----------------------------------------------------------- TokaMaker update
    def tm_update(self, I_units, w_pin, first=False):
        """Solve every slice with TokaMaker, coil + structure currents pinned to I_units by a
        regularisation of weight w_pin; store the plasma's own flux psi_p = psi - psi_external."""
        from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
        gs, m, rp = self.m.gs, self.m, self.scen['ramp']
        gs.set_profiles(ffp_prof=create_power_flux_fun(40, *rp['ffp']), pp_prof=create_power_flux_fun(40, *rp['pp']))
        gs.set_coil_bounds({n: [-1e9, 1e9] for n in m.unit_names})
        gs.set_saddle_constraints(None)
        gs.settings.maxits = 300
        gs.update_settings()
        nfail, prev_I, dpsi = 0, None, []
        for k in range(len(self.t)):
            self.stale[k] = True
            Ik = I_units[k] if (not first or prev_I is None) else prev_I   # first pass: march in time
            units = dict(zip(m.unit_names, Ik))
            gs.set_targets(Ip=float(abs(self.Ip[k])), Ip_ratio=rp['Ip_ratio'])
            gs.set_isoflux_constraints(self.iso[k], weights=np.full(len(self.iso[k]), rp['iso_weight']))
            gs.set_psi_dt(np.zeros(gs.np), 1.0e9)        # static equilibria (no eddy-current term)
            if self.psi_tot[k] is not None:
                gs.set_psi(self.psi_tot[k])
            elif first and k > 0 and self.psi_tot[k - 1] is not None:
                gs.set_psi(self.psi_tot[k - 1])
            else:
                gs.init_psi(*self.shape[k])
            ok = False
            for attempt, wfac in enumerate([1.0, 0.1, 0.01]):     # loosen the pin if TokaMaker fails
                gs.set_coil_reg(reg_terms=[gs.coil_reg_term({n: 1.0}, target=float(v), weight=w_pin * wfac)
                                           for n, v in units.items()])
                gs.set_coil_currents(units)
                if attempt > 0:
                    gs.init_psi(*self.shape[k])
                try:
                    gs.solve()
                    st = gs.get_stats(li_normalization='iter')
                    ok = True
                    break
                except Exception as e:
                    err = e
            if ok and attempt > 0:
                self.log(f'    TokaMaker slice {k} (t={self.t[k]:.2f} s) solved only with the coil pin loosened {1 / wfac:.0f}x')
            if not ok:
                nfail += 1
                self.log(f'    TokaMaker slice {k} (t={self.t[k]:.2f} s) failed: {err}')
                continue
            cc, _ = gs.get_coil_currents()
            Itm = np.array([cc[n] for n in m.unit_names])
            prev_I = Itm
            psi = gs.get_psi(normalized=False)
            self.psi_tot[k] = psi
            new = psi - m.psi_u.T @ Itm
            if self.psi_p[k] is not None:
                dpsi.append(np.abs(new - self.psi_p[k]).max() / max(np.abs(new).max(), 1e-12))
            alpha = rp['relax']                          # under-relax the plasma flux
            self.psi_p[k] = new if (first or self.psi_p[k] is None) else (1 - alpha) * self.psi_p[k] + alpha * new
            self.stale[k] = False
            self.tm_info[k] = {'psi_bounds': [float(x) for x in gs.psi_bounds], 'Ip': float(st['Ip']),
                               'li3': float(st['l_i']), 'o_point': [float(v) for v in gs.o_point],
                               'L_i': float(MU0 * gs.o_point[0] * st['l_i'] / 2),
                               'ext_dev_rel': self.rel_change(Itm - Ik, Ik), 'pin_factor': wfac}
        self.dpsi_rel = float(max(dpsi)) if dpsi and not nfail else np.inf
        return nfail

    # ----------------------------------------------------------- joint convex program
    def solve_qp(self):
        c, m, J, s1 = self.c, self.m, self.joint, self.stage1
        t, N, nv = self.t, len(self.t) - 1, c.ncirc
        rp = self.scen['ramp']
        Phi = np.array([self.sampler.linkage(m.gs, p) for p in self.psi_p])        # [N+1, nu]
        isoA = [eval_fields(m.gs, self.psi_p[k][None, :], self.iso[k])['psi'][0] for k in range(N + 1)]
        G_iso = self.G_iso
        P = J['P0'].clone()                              # stage 1's variables and constraints
        V = P.var('V_ramp', N * nv, np.tile(c.V_max, N))
        Fb = P.var('Fb', N + 1, 1.0)                     # boundary flux per slice [Wb/rad]
        # x_k = Ph_k x_0 + h_k + Gv_k V, with x_0 = X0 . (stage-1 variables)
        Gv, Ph, h = np.zeros((c.nx, N * nv)), np.eye(c.nx), np.zeros(c.nx)
        self.Gv, self.Ph, self.h = [Gv.copy()], [Ph.copy()], [h.copy()]
        for k in range(1, N + 1):
            A, B, Bw = self.zoh(t[k] - t[k - 1])
            w = -(Phi[k] - Phi[k - 1]) / (t[k] - t[k - 1])         # plasma flux change -> forcing
            Gv = A @ Gv
            Gv[:, (k - 1) * nv:k * nv] += B
            Ph, h = A @ Ph, A @ h + Bw @ w
            self.Gv.append(Gv.copy()); self.Ph.append(Ph.copy()); self.h.append(h.copy())

        def aff(C, k):
            """(terms, constant) of the output C x_k."""
            d = {V: C @ self.Gv[k]}
            CP = C @ self.Ph[k]
            for v, M in J['X0'].items():
                d[v] = CP @ M
            return d, C @ self.h[k]
        neg = lambda d: {v: -M for v, M in d.items()}
        used = np.abs(c.T).sum(axis=1) > 0
        Tu = c.U[:m.nc][used]
        lo, hi = m.coil_lim[used, 0], m.coil_lim[used, 1]
        Isc = np.abs(m.coil_lim[used]).max(axis=1)
        for k in range(1, N + 1):
            d, c0 = aff(Tu, k)
            P.le(d, hi - c0, tag='I_max')
            P.le(neg(d), -(lo - c0), tag='I_min')
        P.le({V: np.eye(N * nv)}, np.tile(c.V_max, N), tag='V_max')
        P.le({V: -np.eye(N * nv)}, np.tile(c.V_max, N), tag='V_min')
        # isoflux: psi_ext + psi_p at the target boundary points = Fb_k, weighted so that a flux error
        # equal to a boundary displacement of shape_tol costs 1: dpsi ~ R0 B_p,edge displacement
        for k in range(1, N + 1):
            R0, Bp = self._edge_field(k)
            e = np.zeros((len(self.iso[k]), N + 1)); e[:, k] = -1.0
            d, c0 = aff(G_iso[k], k)
            d[Fb] = e
            P.lsq(d, -(c0 + isoA[k]), 1.0 / (R0 * Bp * rp['shape_tol']), tag='isoflux')
        # smooth voltages: through the ramp, across the hand-off, and back over the whole stage-1
        # channel phase (every interval after the breakdown window)
        D = lambda n: sp.diags([-np.ones(n - 1), np.ones(n - 1)], [0, 1], shape=(n - 1, n))
        Vs = sp.diags(1.0 / c.V_max)
        P.lsq({V: sp.kron(D(N), Vs)}, 0.0, rp['dV'], tag='reg_dV')
        V1, N1 = s1.V, s1.V.n // nv
        E1 = sp.hstack([sp.eye(nv), sp.csr_matrix((nv, (N - 1) * nv))]).tocsr()
        EL = sp.hstack([sp.csr_matrix((nv, (N1 - 1) * nv)), sp.eye(nv)]).tocsr()
        P.lsq({V: Vs @ E1, V1: -Vs @ EL}, 0.0, rp['dV'], tag='smooth_handoff')
        # row j of D is V_{j+1} - V_j (0-based intervals): keep the steps whose earlier interval
        # (t_j, t_{j+1}] already carries the channel
        steps = [j for j in range(N1 - 1) if s1.t[j] >= s1.t_p0 - 1e-12]
        rows = [j * nv + i for j in steps for i in range(nv)]      # each step is a block of nv rows
        P.lsq({V1: sp.kron(D(N1), Vs).tocsr()[rows]}, 0.0, rp['dV'], tag='smooth_channel')
        # stage 1's flux objective becomes a floor: hand-off flux >= stand-alone optimum - slack
        P.le(J['lin'], [J['L_cap']], tag='objective_floor')
        P.lsq({V: sp.kron(sp.eye(N), Vs)}, 0.0, rp['V_reg'], tag='reg_V')
        # proximal pull toward the previous Picard iterate (removes drift between equivalent optima)
        if self.x_prev is not None:
            for k in range(1, N + 1):
                d, c0 = aff(Tu / Isc[:, None], k)
                P.lsq(d, (Tu @ self.x_prev[k]) / Isc - c0, rp['prox_weight'], tag='prox')
        # 0-D Ohm's law of the plasma at its boundary, in energy form V_b I_p = R_p I_p^2 + d/dt(L_i I_p^2/2):
        #   2 pi (Fb_k - Fb_{k-1}) = -[R_p I_p dt + (L_i,k I_k^2 - L_i,k-1 I_k-1^2) / (I_k + I_k-1)]
        # with Fb_0 the boundary flux of the hand-off state
        e0 = np.zeros((1, N + 1)); e0[0, 0] = 1.0
        d0 = {v: -np.mean(G_iso[0] @ M, axis=0, keepdims=True) for v, M in J['X0'].items()}
        d0[Fb] = e0
        P.eq(d0, [float(np.mean(isoA[0]))], tag='Fb0')
        self.Rp_gs, self.Li_gs = self.plasma_RL()
        for k in range(1, N + 1):
            I1, I0_ = abs(self.Ip[k]), abs(self.Ip[k - 1])
            ind = (self.Li_gs[k] * I1 ** 2 - self.Li_gs[k - 1] * I0_ ** 2) / max(I1 + I0_, 1e-9)
            res = 0.5 * (self.Rp_gs[k] + self.Rp_gs[k - 1]) * 0.5 * (I1 + I0_) * (t[k] - t[k - 1])
            e = np.zeros((1, N + 1)); e[0, k] = 1.0; e[0, k - 1] = -1.0
            P.eq({Fb: e}, [-self.s * (res + ind) / (2 * np.pi)], tag='ohm_0d')
        st = P.solve()
        self.P, self.status = P, st
        if solved(st):
            s1.load(P)                                   # stage 1 re-optimised: new hand-off state
            self.x0 = s1.x[-1][:c.nx].copy()
            vv = P.value(V)
            self.V, self.Fb = vv.reshape(N, nv), P.value(Fb)
            self.x = np.array([self.Ph[k] @ self.x0 + self.h[k] + self.Gv[k] @ vv for k in range(N + 1)])
            iso_res = [float(np.abs(G_iso[k] @ self.x[k] + isoA[k] - self.Fb[k]).max()) for k in range(1, N + 1)]
            self.shape_err_m = [r / np.prod(self._edge_field(k)) for k, r in enumerate(iso_res, start=1)]
        return st

    def _edge_field(self, k):
        """(R0, B_p at the plasma edge) of slice k's target shape."""
        R0, _, a, kap, _ = self.shape[k]
        return R0, MU0 * max(abs(self.Ip[k]), 1e3) / (2 * np.pi * a * np.sqrt((1 + kap ** 2) / 2))

    def plasma_RL(self):
        """Spitzer resistance from the T_e keyframes and internal inductance L_i = mu0 R l_i(3) / 2
        from each slice's TokaMaker equilibrium."""
        rp = self.scen['ramp']
        Te = interp_keys(rp['Te'], self.t).ravel()
        R0, a, kap = self.shape[:, 0], self.shape[:, 2], self.shape[:, 3]
        eta = 5.2e-5 * rp['Zeff'] * LN_LAMBDA / np.maximum(Te, 1.0) ** 1.5
        return eta * 2 * np.pi * R0 / (np.pi * a ** 2 * kap), np.array([i['L_i'] for i in self.tm_info])

    def rel_change(self, dI, I_ref):
        """Largest relative change: coils relative to their limits, structure relative to the
        largest segment current."""
        nc = self.m.nc
        dI, I_ref = np.atleast_2d(dI), np.atleast_2d(I_ref)
        lim = np.abs(self.m.coil_lim).max(axis=1)
        return max(float((np.abs(dI[:, :nc]) / lim).max()),
                   float(np.abs(dI[:, nc:]).max() / max(np.abs(I_ref[:, nc:]).max(), 1.0)))

    def unit_currents(self, x):
        return (self.c.U @ x.T).T

    def voltage_steps(self):
        """Largest |V_j - V_j-1| / V_max per circuit from the end of the breakdown window to the end
        of the ramp (stage-1 channel phase, hand-off, ramp)."""
        s1 = self.stage1
        k0 = int(np.searchsorted(s1.t, s1.t_p0 - 1e-12))
        return np.abs(np.diff(np.vstack([s1.v[k0:], self.V]), axis=0)).max(axis=0) / self.c.V_max

    # ----------------------------------------------------------- Picard driver
    def run(self, max_iter=15, tol=2e-3):
        """Converged when every TokaMaker slice solves, the currents change by < tol, TokaMaker's
        external currents stay within 5 tol of the convex program's, and the raw plasma flux
        changes by < 5 tol."""
        rp, n = self.scen['ramp'], len(self.t)
        self.x_prev = None
        self.psi_p, self.psi_tot, self.tm_info = [None] * n, [None] * n, [None] * n
        self.stale, self.dpsi_rel = [True] * n, np.inf
        I_prev = self.unit_currents(np.repeat(self.x0[None, :], n, axis=0))
        self.history = []
        gs = self.m.gs
        maxits0 = gs.settings.maxits
        try:
            self.tm_update(I_prev, w_pin=rp['seed_weight'], first=True)
            for it in range(max_iter):
                if any(p is None for p in self.psi_p):
                    self.log('  a TokaMaker slice is missing; stopping'); return False
                st = self.solve_qp()
                if not solved(st):
                    self.log(f'  Picard {it}: convex program {st}'); return False
                I_new = self.unit_currents(self.x)
                self.x_prev = self.x.copy()
                dI = self.rel_change(I_new - I_prev, I_new)
                nf = self.tm_update(I_new, w_pin=rp['pin_weight'])
                devs = [i['ext_dev_rel'] for k, i in enumerate(self.tm_info) if i and k > 0 and not self.stale[k]]
                dev = max(devs) if devs else np.inf
                h = {'iter': it, 'status': st, 'dI_rel': dI, 'tm_dev': dev, 'tm_fail': nf, 'dpsi_rel': self.dpsi_rel,
                     'shape_err_max_cm': 100 * max(self.shape_err_m), 'flux_handoff_Wb': self.stage1.flux_handoff(),
                     'dV_max': float(self.voltage_steps().max())}
                self.history.append(h)
                self.log(f"  Picard {it}: {st:8s} current change {dI:.1e}, TokaMaker deviation {dev:.1e}, "
                         f"plasma flux change {self.dpsi_rel:.1e}, max shape error {h['shape_err_max_cm']:.2f} cm")
                I_prev = I_new
                if nf == 0 and dI < tol and dev < 5 * tol and self.dpsi_rel < 5 * tol:
                    return True
            return False
        finally:
            # undo the ramp's settings: constraints, targets, coil bounds, time step, iteration limit
            # (the coil regularisation cannot be cleared through the API; the next user sets their own)
            gs.set_isoflux_constraints(None)
            gs.set_targets()
            gs.set_coil_bounds(None)
            gs.set_psi_dt(np.zeros(gs.np), 1.0e9)
            gs.settings.maxits = maxits0
            gs.update_settings()

    def handoff_jump(self):
        """Change of the plasma's flux linkage with the coils and structure at the hand-off, where
        stage 1's uniform ring is replaced by the first TokaMaker equilibrium (same I_p). The ramp
        does not model this jump. If it happened at once, flux conservation (M dx = -U^T dPhi)
        would kick the currents by dx; the kicks are returned relative to the coil limits and to
        the largest structure current at the hand-off."""
        s1, nc = self.stage1, self.m.nc
        Phi_ring = s1.plasma.M_u * s1.x[-1, s1.c.ip]
        Phi_gs = self.sampler.linkage(self.m.gs, self.psi_p[0])
        dPhi = Phi_gs - Phi_ring
        dI = self.unit_currents((-self.c.Minv @ (self.c.U.T @ dPhi))[None, :])[0]
        I0 = self.unit_currents(self.x0[None, :])[0]
        return {'dPhi_rel': float(np.abs(dPhi).max() / np.abs(Phi_gs).max()),
                'coil_kick_rel': float((np.abs(dI[:nc]) / np.abs(self.m.coil_lim).max(axis=1)).max()),
                'struct_kick_rel': float(np.abs(dI[nc:]).max() / np.abs(I0[nc:]).max())}

    def summary(self):
        return {'flux_handoff_standalone_Wb': self.joint['flux_star_Wb'],
                'flux_handoff_Wb': self.stage1.flux_handoff(),
                'flux_boundary_handoff_Wb': float(2 * np.pi * self.Fb[0]),
                'flux_boundary_end_Wb': float(2 * np.pi * self.Fb[-1]),
                'ramp_flux_cost_Wb': float(2 * np.pi * (self.Fb[0] - self.Fb[-1])),
                'shape_err_max_cm': float(100 * max(self.shape_err_m)),
                'dV_max_after_breakdown': float(self.voltage_steps().max()),
                'picard_iterations': len(self.history)}
