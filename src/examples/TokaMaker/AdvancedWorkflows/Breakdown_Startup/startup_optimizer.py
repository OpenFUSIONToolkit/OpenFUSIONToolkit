#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
"""Start-up optimisation: breakdown and early plasma (stage 1), then a ramp of TokaMaker
free-boundary equilibria (stage 2), optimised together.

    ConicProblem    sparse conic-program builder on the Clarabel solver
    time_grid       stage-1 time grid
    StartupProblem  stage 1 as one convex program in the premagnetisation currents and supply voltages
    GSRamp          stage 2: Picard iteration between TokaMaker and a joint (stage 1 + ramp) convex program

The scenario keys each stage reads (defaults, units) are listed in the StartupProblem and GSRamp
docstrings. Units: SI (A, V, T, m, s, Ohm); psi in Wb/rad (TokaMaker's convention, flux in Wb is
2 pi psi); electron temperatures in eV.
"""
import time

import numpy as np
import scipy.sparse as sp

from startup_model import (eval_fields, disk_points, cell_geometry, ChannelPlasma, Te_burnthrough,
                           _disk_quadrature, LN_LAMBDA, MU0)


# ============================================================================ conic programs
# Clarabel settings tried in turn after NumericalError / InsufficientProgress / MaxIterations / an
# unacceptable AlmostSolved: wider equilibration; stronger static regularisation with more iterative
# refinement; looser tolerances.
SOLVER_LADDER = (
    {'equilibrate_max_iter': 50, 'equilibrate_min_scaling': 1e-8, 'equilibrate_max_scaling': 1e8,
     'max_iter': 400},
    {'static_regularization_constant': 1e-6, 'iterative_refinement_max_iter': 30,
     'iterative_refinement_reltol': 1e-14, 'max_iter': 400},
    {'tol_feas': 1e-6, 'tol_gap_abs': 1e-6, 'tol_gap_rel': 1e-6, 'max_iter': 600,
     'equilibrate_max_iter': 50, 'equilibrate_min_scaling': 1e-8, 'equilibrate_max_scaling': 1e8},
)


class Var:
    """A block of n decision variables at columns start..start+n; scale: typical size of each
    (scalar or length n), used to scale the solver's columns."""

    def __init__(self, name, start, n, scale):
        self.name, self.start, self.n = name, start, n
        self.scale, self.sl = np.full(n, scale, float), slice(start, start + n)


class ConicProblem:
    """    minimise   1/2 sum_i w_i^2 ||A_i x - b_i||^2 + q^T x
           subject to linear equalities, linear inequalities and second-order cones.

    Constraints are dicts {Var: matrix} in physical units; rows and columns are rescaled before the
    Clarabel solve. Least-squares terms become residual variables with a diagonal quadratic cost.
    After a solve, self.violations holds the largest violation of each constraint tag (physical units)."""

    def __init__(self, settings=None):
        """settings: extra Clarabel settings for the first attempt, e.g. {'direct_solve_method': 'faer'}."""
        self.vars, self.nvar, self.settings = {}, 0, dict(settings or {})
        self._eq, self._ineq, self._soc, self._quad, self._lin = [], [], [], [], []

    def clone(self, objective=True):
        """Copy with the same variables and terms, to add more to (objective=False: no linear objective)."""
        P = ConicProblem(self.settings)
        P.vars, P.nvar = dict(self.vars), self.nvar
        P._eq, P._ineq, P._soc, P._quad = list(self._eq), list(self._ineq), list(self._soc), list(self._quad)
        P._lin = list(self._lin) if objective else []
        return P

    def var(self, name, n, scale=1.0):
        """Declare a block of n variables of typical size scale; returns its Var."""
        v = self.vars[name] = Var(name, self.nvar, int(n), scale)
        self.nvar += v.n
        return v

    @staticmethod
    def _rhs(terms, x):
        """x as a float vector with one entry per row of the terms {Var: matrix} (block shapes checked)."""
        m = sp.csr_matrix(next(iter(terms.values()))).shape[0]
        for v, M in terms.items():
            assert sp.csr_matrix(M).shape == (m, v.n), (v.name, sp.csr_matrix(M).shape, (m, v.n))
        return np.full(m, x, float)

    def eq(self, terms, rhs, tag=''):
        """sum_v M_v v = rhs"""
        self._eq.append((terms, self._rhs(terms, rhs), tag))

    def le(self, terms, rhs, tag=''):
        """sum_v M_v v <= rhs"""
        self._ineq.append((terms, self._rhs(terms, rhs), tag))

    def soc(self, t_terms, t0, v_terms, v0, tag=''):
        """|| sum_v M_v v + v0 ||_2 <= sum_v a_v v + t0   (t_terms {Var: [1, n]}, v_terms {Var: [m, n]})"""
        self._soc.append((t_terms, float(t0), v_terms, self._rhs(v_terms, v0), tag))

    def lsq(self, terms, target, weight, tag=''):
        """Add (weight^2 / 2) || sum_v M_v v - target ||^2 to the objective."""
        self._quad.append((terms, self._rhs(terms, target), float(weight), tag))

    def linear(self, var, coeff):
        """Add coeff . var to the objective."""
        self._lin.append((var, np.full(var.n, coeff, float)))

    def _row_block(self, terms, nrows, ncol):
        """{Var: matrix} -> [nrows, ncol] CSR in the solver's scaled columns, dust below 1e-10 x row max dropped."""
        out = sp.csr_matrix((nrows, ncol))
        for v, M in terms.items():
            out = out + sp.hstack([sp.csr_matrix((nrows, v.start)), sp.csr_matrix(M) @ sp.diags(v.scale),
                                   sp.csr_matrix((nrows, ncol - v.start - v.n))], format='csr')
        out = out.tocsr()
        if out.nnz:
            rmax = np.asarray(abs(out).max(axis=1).todense()).ravel()
            rows = np.repeat(np.arange(out.shape[0]), np.diff(out.indptr))
            out.data[np.abs(out.data) < 1e-10 * rmax[rows]] = 0.0
            out.eliminate_zeros()
        return out

    def solve(self):
        """Assemble and solve with Clarabel, walking down SOLVER_LADDER on numerical trouble. Returns the
        final status (also self.status); self.solve_log lists the attempts."""
        import clarabel
        nres = sum(q[1].size for q in self._quad)        # least-squares residual variables
        ntot = self.nvar + nres
        rows_A, rows_b, cones = [], [], []

        def add_rows(A, b, is_eq, tag):                  # drop empty rows (if feasible), row-normalise
            rm = np.asarray(abs(A).max(axis=1).todense()).ravel()
            keep = rm != 0
            if np.any(np.abs(b[~keep]) > 0 if is_eq else b[~keep] < 0):
                raise ValueError(f'constraint "{tag}" has empty rows that are infeasible')
            rows_A.append(sp.diags(1 / rm[keep]) @ A[keep])
            rows_b.append(b[keep] / rm[keep])
            return int(keep.sum())

        # zero cone: residual definitions r = w (A x - target), then the equalities
        n_zero, r0 = 0, 0
        for terms, target, w, tag in self._quad:
            m = target.size
            R = sp.csr_matrix((np.full(m, -1.0), (np.arange(m), r0 + np.arange(m))), shape=(m, nres))
            n_zero += add_rows(sp.hstack([self._row_block(terms, m, self.nvar) * w, R], format='csr'),
                               w * target, True, tag)
            r0 += m
        for terms, rhs, tag in self._eq:
            n_zero += add_rows(self._row_block(terms, rhs.size, ntot), rhs, True, tag)
        if n_zero:
            cones.append(clarabel.ZeroConeT(n_zero))
        n_nonneg = sum(add_rows(self._row_block(terms, rhs.size, ntot), rhs, False, tag)
                       for terms, rhs, tag in self._ineq)
        if n_nonneg:
            cones.append(clarabel.NonnegativeConeT(n_nonneg))
        for t_terms, t0, v_terms, v0, tag in self._soc:  # each cone scaled by its largest coefficient
            A = sp.vstack([-self._row_block(t_terms, 1, ntot), self._row_block(v_terms, v0.size, ntot)], format='csr')
            sc = max(abs(A).max(), 1e-300)
            rows_A.append(A / sc)
            rows_b.append(np.concatenate([[t0], -v0]) / sc)
            cones.append(clarabel.SecondOrderConeT(A.shape[0]))
        A, b = sp.vstack(rows_A, format='csc'), np.concatenate(rows_b)
        q = np.zeros(ntot)
        for v, c in self._lin:
            q[v.sl] += c * v.scale
        Pm = sp.triu(sp.diags(np.r_[np.zeros(self.nvar), np.ones(nres)], format='csc'), format='csc')
        self.solve_log = []
        for step, extra in enumerate(({},) + SOLVER_LADDER):
            st = clarabel.DefaultSettings()
            st.verbose = False
            for k, val in {**self.settings, **extra}.items():
                setattr(st, k, val)
            sol = clarabel.DefaultSolver(Pm, q, A, b, cones, st).solve()
            self.status, self._x = str(sol.status), np.array(sol.x)
            r = A @ self._x - b                          # residuals of the normalised rows
            self._eq_res = float(np.abs(r[:n_zero]).max(initial=0.0))
            self._ineq_res = float(np.maximum(r[n_zero:n_zero + n_nonneg], 0.0).max(initial=0.0))
            self.violations = self.evaluate({v: self.value(v) for v in self.vars.values()})
            note = (f'equality residual {self._eq_res:.1e}',) if self.status == 'AlmostSolved' else ()
            self.solve_log.append((step, self.status) + note)
            # stop when accepted, or on a definite answer (e.g. infeasible)
            if self.accepted() or not any(f in self.status for f in
                                          ('NumericalError', 'InsufficientProgress', 'MaxIterations', 'AlmostSolved')):
                break
        return self.status

    def accepted(self, tol=1e-6):
        """Solved, or AlmostSolved with every equality and inequality met to tol in the normalised rows."""
        return self.status == 'Solved' or (self.status == 'AlmostSolved' and self._eq_res <= tol
                                           and self._ineq_res <= tol)

    def evaluate(self, values):
        """Largest violation of every constraint tag (physical units) for variable values {Var: array}:
        |lhs - rhs| for equalities, the amount by which inequalities and cones fail (0 when met)."""
        def lhs(terms):
            return sum(np.asarray(sp.csr_matrix(M) @ values[v]).ravel() for v, M in terms.items())
        out = {}
        for tag, val in ([(tag, np.abs(lhs(T) - r).max()) for T, r, tag in self._eq]
                         + [(tag, np.maximum(lhs(T) - r, 0.0).max()) for T, r, tag in self._ineq]
                         + [(tag, max(np.linalg.norm(lhs(Vt) + v0) - (lhs(T)[0] + t0), 0.0))
                            for T, t0, Vt, v0, tag in self._soc]):
            out[tag] = max(out.get(tag, 0.0), float(val))
        return out

    def value(self, var):
        """Solution of a variable block, in physical units."""
        return self._x[var.sl] * var.scale

    def linear_terms(self):
        """The linear objective as {Var: coefficient row [1, n]}, in physical units."""
        lin = {}
        for v, q in self._lin:
            lin[v] = lin.get(v, 0) + q[None, :]
        return lin


# ============================================================================ stage 1
def time_grid(tc):
    """Stage-1 slice times [s]: t_start --dt_pre--> t_fine_start --dt_fine--> t_fine_end --dt_ramp--> t_end,
    each segment split into equal steps no longer than its dt."""
    ts, out = (tc['t_start'], tc['t_fine_start'], tc['t_fine_end'], tc['t_end']), [[float(tc['t_start'])]]
    for a, b, d in zip(ts[:-1], ts[1:], (tc['dt_pre'], tc['dt_fine'], tc['dt_ramp'])):
        if b - a > 1e-12:
            out.append(np.linspace(a, b, max(1, int(np.ceil((b - a) / d - 1e-9))) + 1)[1:])
    return np.concatenate(out)


class StartupProblem:
    """Stage 1: premagnetisation -> breakdown -> early plasma channel up to the hand-off.

    Variables: premagnetisation currents I0 of the driven circuits (structure at rest at t_start) and
    supply voltages V_1..V_N, constant on each interval (t_{k-1}, t_k]. The circuit is linear and stepped
    exactly, so every state is affine in them, x_k = Gi_k I0 + Gv_k V + h_k (h_k != 0 only with a
    prescribed I_p), and the design is one convex program. E_phi = -(psi_k - psi_{k-1}) / (R dt).

    Constraints: coil current, supply voltage and slew limits; in the breakdown window an E_phi floor and
    the Lloyd criterion over the breakdown region; for the channel, radial force balance and decay index
    at its centre, an E_phi floor during burn-through, the stray-field limit, monotone I_p and the
    hand-off current. Objective: maximise the flux at the plasma boundary (the channel's inboard contact
    point) at the hand-off.

    Scenario keys (required unless a default is given):
      time            t_start, t_fine_start, t_fine_end, t_end [s]; dt_pre, dt_fine, dt_ramp [s]: largest
                      step in each segment (time_grid)
      region          breakdown disk R, Z, radius [m]; n_rings, n_per_ring: its points (disk_points)
      breakdown       window [t0, t1] [s] (the channel exists after t1); E_min [V/m]: E_phi floor;
                      E_star [V/m], margin [-]: Lloyd criterion E_phi B_phi / B_p >= margin E_star
      plasma          Rc, Zc, a [m]: channel, limited on the inboard wall at Rc - a; li, Zeff [-];
                      beta_p [-] (default 0.0) in the Shafranov B_v; Te0, Te1 [eV], tau_Te [s]: burn-through
                      T_e(t); Ip_target [A]: I_p at t_end (unused with Ip_prescribed);
                      Ip_prescribed (default None): [[t, I_p], ...] keyframes [s, A];
                      force_balance {'frac' [-], 'B_tol' [T]}: |B_z,ext - B_v,req|, |B_R,ext| <= frac kf I_p + B_tol;
                      decay_index [n_lo, n_hi] [-]; E_min_burn [V/m] up to t_burn_end [s]: E_phi floor at Rc;
                      stray (default None: no limit) {'self_frac' [-], 't_end' [s]}
      regularization  dV, V [-]: weights of the voltage smoothing / size penalties (on V / V_max)
      slew            (default None: no limit) supply slew limit [V_max per s]; also used by GSRamp
    """

    def __init__(self, machine, circuit, scen, settings=None):
        """circuit: a plasma-free CircuitModel. The channel is added to it as one more state (current set by
        the circuit equations, Spitzer resistance); with Ip_prescribed it is not a circuit state: its
        current is a known last state component and its changing flux linkage a known forcing.
        settings: extra Clarabel settings (inherited by the ramp's joint programs)."""
        self.m, self.c, self.scen, self.settings = machine, circuit, scen, settings
        self.s = 1.0                                     # sign of I_p and E_phi: positive I_p only
        self.t = time_grid(scen['time'])
        rg, pl = scen['region'], scen['plasma']
        self.pts = disk_points(rg['R'], rg['Z'], rg['radius'], rg['n_rings'], rg['n_per_ring'])
        self.plasma = ChannelPlasma(machine, pl['Rc'], pl['Zc'], pl['a'], li=pl['li'], Zeff=pl['Zeff'])
        self.prescribed = pl.get('Ip_prescribed') is not None
        # state: circuit states, then I_p; U maps state -> unit currents (zero I_p column)
        self.ip, self.nx = circuit.nx, circuit.nx + 1
        self.U = np.hstack([circuit.U, np.zeros((circuit.U.shape[0], 1))])
        if self.prescribed:   # M_xx dx/dt + R x = S V - M_xp dI_p/dt: forcing _F times w = dI_p/dt
            self.Ip_presc = interp_keys(pl['Ip_prescribed'], self.t).ravel()
            self._F = -(self.plasma.M_u @ circuit.U)[:, None]
        else:
            circuit.add_plasma(self.plasma)

        def out(G_units, G_plasma=None):                 # per-unit (+ per-ampere channel) response -> [n_out, nx]
            C = G_units @ self.U
            if G_plasma is not None:
                C[:, self.ip] = G_plasma
            return C
        gs, pu = machine.gs, machine.psi_u
        G, gp = eval_fields(gs, pu, self.pts), self.plasma.field_at(self.pts)
        self.Cpsi, self.CBr, self.CBz = out(G['psi'].T, gp[0]), out(G['Br'].T, gp[1]), out(G['Bz'].T, gp[2])
        # boundary flux: total flux at the channel's contact point on the inboard wall (as in the ramp)
        contact = np.array([[pl['Rc'] - pl['a'], pl['Zc']]])
        self.Cbound = out(eval_fields(gs, pu, contact)['psi'].T, self.plasma.psi_unit(contact))
        Gc = eval_fields(gs, pu, np.array([[pl['Rc'], pl['Zc']]]), second=True)   # applied fields at the centre
        self.Cbz_c, self.Cbr_c, self.Cdbz_c, self.Cpsi_c = (out(Gc[k].T) for k in ('Bz', 'Br', 'dBz_dR', 'psi'))
        # stray-field points: rings across the channel plus 16 on its edge; applied B_R, B_Z, total psi
        th = 2 * np.pi * np.arange(16) / 16
        p_s = np.vstack([_disk_quadrature(pl['Rc'], pl['Zc'], pl['a'], 3, 8)[0],
                         np.column_stack([pl['Rc'] + pl['a'] * np.cos(th), pl['Zc'] + pl['a'] * np.sin(th)])])
        Gs = eval_fields(gs, pu, p_s)
        self.Cbr_s, self.Cbz_s = out(Gs['Br'].T), out(Gs['Bz'].T)
        self.Cpsi_s, self.R_s = out(Gs['psi'].T, self.plasma.psi_unit(p_s)), p_s[:, 0]
        self.Rp = self.pts[:, 0]
        self.Bphi = np.abs(machine.Bphi(self.Rp))

    def aff(self, C, k):
        """Output C x_k as ({Var: matrix}, constant), x_k = Gi_k I0 + Gv_k V + h_k."""
        return {self.I0v: C @ self.Gi[k], self.V: C @ self.Gv[k]}, C @ self.h[k]

    def daff(self, C, k):
        """C (x_k - x_{k-1}) as ({Var: matrix}, constant)."""
        (d1, c1), (d0, c0) = self.aff(C, k), self.aff(C, k - 1)
        return {v: d1[v] - d0[v] for v in d1}, c1 - c0

    def build(self):
        """Assemble the stage-1 convex program (self.P, also returned)."""
        c, m, scen, t, s = self.c, self.m, self.scen, self.t, self.s
        pl, bd = scen['plasma'], scen['breakdown']
        N, nx, nv = len(t) - 1, self.nx, c.ncirc
        Ilim = np.abs(m.coil_lim).max(axis=1)
        P = self.P = ConicProblem(self.settings)
        I0 = self.I0v = P.var('I0', nv, [Ilim[np.abs(c.U[:m.nc, j]) > 0].max() for j in range(nv)])
        V = self.V = P.var('V', N * nv, np.tile(c.V_max, N))
        aff, daff = self.aff, self.daff

        # plasma resistance per interval: open circuit until the channel exists (t_p0), then Spitzer(T_e)
        self.t_p0 = bd['window'][1]
        self.Te = Te_burnthrough(0.5 * (t[1:] + t[:-1]), self.t_p0, pl['Te0'], pl['Te1'], pl['tau_Te'])
        self.Rplas = np.where(t[1:] > self.t_p0 + 1e-12,       # same steps as Kp below
                               self.plasma.resistance(self.Te), c.R_OPEN)

        # affine state maps; the circuit steps its c.nx states (all of them, or all but a prescribed I_p)
        n0 = c.nx
        Gi, Gv, h = np.eye(nx, nv), np.zeros((nx, N * nv)), np.zeros(nx)
        if self.prescribed:
            h[self.ip] = self.Ip_presc[0]
        self.Gi, self.Gv, self.h = [Gi.copy()], [Gv.copy()], [h.copy()]
        for k in range(1, N + 1):
            dt = t[k] - t[k - 1]
            if self.prescribed:
                A, B, Bw = c.zoh(dt, forcing=self._F)
                h[:n0] = A @ h[:n0] + Bw[:, 0] * (self.Ip_presc[k] - self.Ip_presc[k - 1]) / dt
                h[self.ip] = self.Ip_presc[k]
            else:
                A, B = c.zoh(dt, self.Rplas[k - 1])
            Gi[:n0], Gv[:n0] = A @ Gi[:n0], A @ Gv[:n0]
            Gv[:n0, (k - 1) * nv:k * nv] += B
            self.Gi.append(Gi.copy())
            self.Gv.append(Gv.copy())
            self.h.append(h.copy())

        # coil current limits at every slice, supply voltage limits on every interval
        used = np.abs(c.T).sum(axis=1) > 0               # coils driven by some circuit
        Tu, lo, hi = self.U[:m.nc][used], m.coil_lim[used, 0], m.coil_lim[used, 1]
        Ii, Iv = np.vstack([Tu @ G for G in self.Gi]), np.vstack([Tu @ G for G in self.Gv])
        Ic = np.concatenate([Tu @ hk for hk in self.h])
        P.le({I0: Ii, V: Iv}, np.tile(hi, N + 1) - Ic, tag='I_max')
        P.le({I0: -Ii, V: -Iv}, -np.tile(lo, N + 1) + Ic, tag='I_min')
        P.le({V: sp.eye(N * nv)}, np.tile(c.V_hi, N), tag='V_max')
        P.le({V: -sp.eye(N * nv)}, -np.tile(c.V_lo, N), tag='V_min')

        # breakdown window: E_phi >= E_min at every region point
        w0, w1 = bd['window']
        self.K = [k for k in range(1, N + 1) if t[k - 1] >= w0 - 1e-12 and t[k] <= w1 + 1e-12]
        for k in self.K:
            d, c0 = daff(s * self.Cpsi / self.Rp[:, None] / (t[k] - t[k - 1]), k)
            P.le(d, -bd['E_min'] + c0, tag='E_floor')
        # Lloyd-type cones at both ends of each step k, with E_phi of the step:
        #   ||(B_R, B_Z + kb I_p)|| <= B_phi s E_phi / (margin E_star) + ks s I_p
        # 'lloyd': breakdown region, total field, kb = ks = 0. 'stray': points across the channel up to
        # stray t_end, applied field minus the Shafranov B_v,req = -kf I_p (kb = kf), plus a fraction of the
        # plasma's own edge field (ks = self_frac mu0 / (2 pi a))
        eip = np.eye(nx)[[self.ip]]                      # [1, nx]: selects I_p
        Kp = [k for k in range(1, N + 1) if t[k] > self.t_p0 + 1e-12]       # steps with a channel
        Rc, a = pl['Rc'], pl['a']
        self.kf = kf = MU0 / (4 * np.pi * Rc) * (np.log(8 * Rc / a) + pl.get('beta_p', 0.0) + pl['li'] / 2 - 1.5)
        sy = pl.get('stray')
        cones = [('lloyd', self.K, self.Cpsi, self.CBr, self.CBz, self.Rp, 0.0, 0.0)]
        if sy:
            cones.append(('stray', [k for k in Kp if t[k] <= sy['t_end'] + 1e-12], self.Cpsi_s, self.Cbr_s,
                          self.Cbz_s, self.R_s, kf, sy['self_frac'] * MU0 / (2 * np.pi * a)))
        for tag, Ks, Cpsi, CBr, CBz, R, kb, ks in cones:
            lloyd = np.abs(m.Bphi(R)) / (bd['margin'] * bd['E_star'])
            for k in Ks:
                dE, cE = daff(-s * Cpsi / R[:, None] / (t[k] - t[k - 1]), k)        # s E_phi on step k
                for j in (k - 1, k):
                    (dI, cI), (dBr, cBr), (dBz, cBz) = aff(ks * s * eip, j), aff(CBr, j), aff(CBz + kb * eip, j)
                    for p in range(len(R)):
                        P.soc({v: lloyd[p] * dE[v][p:p + 1] + dI[v] for v in dE}, lloyd[p] * cE[p] + cI[0],
                              {v: np.vstack([dBr[v][p], dBz[v][p]]) for v in dBr}, [cBr[p], cBz[p]], tag=tag)

        # early channel
        if not self.prescribed:                          # I_p never falls, and is Ip_target at t_end
            for k in Kp:
                d, c0 = daff(-s * eip, k)
                P.le(d, -c0, tag='Ip_mono')
            d, c0 = aff(s * eip, N)
            P.eq(d, [pl['Ip_target'] - c0[0]], tag='Ip_target')
        # force balance |B_z,ext + kf I_p|, |B_R,ext| <= frac kf |I_p| + B_tol; decay index
        # n = -(Rc/B_z) dB_z/dR in [n_lo, n_hi], multiplied through by B_z (sign -s); tags name the bound
        frac, Btol = pl['force_balance']['frac'], pl['force_balance']['B_tol']
        nlo, nhi = pl['decay_index']
        lo_, hi_ = (nhi, nlo) if s > 0 else (nlo, nhi)
        band, Cv, q = frac * kf * s * eip, self.Cbz_c + kf * eip, -Rc * self.Cdbz_c
        rows = [(Cv - band, Btol, 'Bv_hi'), (-Cv - band, Btol, 'Bv_lo'), (self.Cbr_c - band, Btol, 'Br_hi'),
                (-self.Cbr_c - band, Btol, 'Br_lo'), (lo_ * self.Cbz_c - q, 0.0, 'n_hi' if s > 0 else 'n_lo'),
                (q - hi_ * self.Cbz_c, 0.0, 'n_lo' if s > 0 else 'n_hi')]
        for k in Kp:
            if t[k] <= pl['t_burn_end'] + 1e-12:         # E_phi floor at the centre during burn-through
                d, c0 = daff(s * self.Cpsi_c / (Rc * (t[k] - t[k - 1])), k)
                P.le(d, [-pl['E_min_burn'] - c0[0]], tag='E_burn')
            for C, b, tag in rows:
                d, c0 = aff(C, k)
                P.le(d, [b - c0[0]], tag=tag)

        # weak regularisation (smooth and small V / V_max); slew limit |V_k - V_{k-1}| <= slew V_max dt_k
        reg, Vs = scen['regularization'], sp.diags(1.0 / c.V_max)
        D = sp.diags([-1.0, 1.0], [0, 1], shape=(N - 1, N))            # row j: x_{j+1} - x_j
        P.lsq({V: sp.kron(D, Vs)}, 0.0, reg['dV'], tag='reg_dV')
        P.lsq({V: sp.kron(sp.eye(N), Vs)}, 0.0, reg['V'], tag='reg_V')
        if scen.get('slew'):
            Ds = sp.kron(D, sp.eye(nv))
            cap = np.repeat(scen['slew'] * np.diff(t)[1:], nv) * np.tile(c.V_max, N - 1)
            P.le({V: Ds}, cap, tag='slew_up')
            P.le({V: -Ds}, cap, tag='slew_dn')

        # objective: maximise s psi at the plasma boundary at the hand-off
        for v, M in self.handoff_terms().items():
            P.linear(v, -s * (self.Cbound @ M).ravel())
        return P

    def solve(self):
        """Solve the program from build() and load the solution; raises if it is not accepted."""
        st = self.P.solve()
        if not self.P.accepted():
            raise RuntimeError(f'stage-1 convex program not solved: {st} (solver attempts {self.P.solve_log})')
        self.load(self.P)
        return st

    def handoff_terms(self):
        """Variable part of the hand-off state x_N as {Var: matrix} (x_N = Gi_N I0 + Gv_N V + h_N)."""
        return {self.I0v: self.Gi[-1], self.V: self.Gv[-1]}

    def channel_stray(self):
        """Stray field across the early channel at every slice: the applied poloidal field minus the B_v
        force balance needs, largest over the stray-field points, and the plasma's own edge field
        mu0 I_p / (2 pi a); while the stray field is larger the channel has no closed flux surfaces.
        Returns {'t' [s], 'stray' [T], 'self' [T], 'ratio' [-]}: ratio is the largest stray field / its
        limit as the stray-field constraint applies it (NaN where it does not apply)."""
        pl, bd, t, x = self.scen['plasma'], self.scen['breakdown'], self.t, self.x
        sy = pl.get('stray') or {'self_frac': 1.0, 't_end': np.inf}
        Ip = x[:, self.ip]                               # (the stray maps have no I_p column)
        stray = np.hypot(self.Cbr_s @ x.T, self.Cbz_s @ x.T + self.kf * Ip[None, :])     # [n_pts, N+1]
        lloyd = np.abs(self.m.Bphi(self.R_s)) / (bd['margin'] * bd['E_star'])
        ks = sy['self_frac'] * MU0 * self.s / (2 * np.pi * pl['a'])
        ratio = np.full(len(t), np.nan)                  # each slice: larger ratio of the steps it ends
        for k in range(1, len(t)):
            if self.t_p0 + 1e-12 < t[k] <= sy['t_end'] + 1e-12:
                E = -self.s * (self.Cpsi_s @ (x[k] - x[k - 1])) / self.R_s / (t[k] - t[k - 1])
                for j in (k - 1, k):
                    ratio[j] = np.fmax(ratio[j], (stray[:, j] / np.maximum(lloyd * E + ks * Ip[j], 1e-12)).max())
        return {'t': t, 'stray': stray.max(axis=0), 'self': MU0 * np.abs(Ip) / (2 * np.pi * pl['a']), 'ratio': ratio}

    def flux_handoff(self):
        """Flux left at the hand-off: 2 pi s psi at the plasma boundary [Wb]."""
        return float(2 * np.pi * self.s * (self.Cbound @ self.x[-1])[0])

    def load(self, P):
        """Read I0 and V from a solved program with this stage's variables (its own or the joint one);
        rebuild the states self.x [N+1, nx] and voltages self.v [N, ncirc]."""
        i0, vv = P.value(self.I0v), P.value(self.V)
        self.v = vv.reshape(len(self.t) - 1, self.c.ncirc)
        self.x = np.array([Gi @ i0 + Gv @ vv + h for Gi, Gv, h in zip(self.Gi, self.Gv, self.h)])

    def results(self):
        """Trajectories: currents, voltages, region fields, Lloyd quantity, channel quantities."""
        c, m, t, x = self.c, self.m, self.t, self.x
        xdot = np.zeros_like(x)
        xdot[1:] = np.diff(x, axis=0) / np.diff(t)[:, None]          # interval averages
        E = -(self.Cpsi @ xdot.T) / self.Rp[:, None]                   # [n_pts, N+1]
        B = np.hypot(self.CBr @ x.T, self.CBz @ x.T)
        x_ext = x.copy()
        x_ext[:, self.ip] = 0.0                                        # coils + structure only
        Bz_c = (self.Cbz_c @ x.T).ravel()
        n_c = -self.scen['plasma']['Rc'] * (self.Cdbz_c @ x.T).ravel() / np.where(np.abs(Bz_c) > 0, Bz_c, np.nan)
        return {'t': t, 'x': x, 'xdot': xdot, 'V': self.v, 'I_coil': (c.T @ x[:, :c.ncirc].T).T,
                'I_struct': x[:, c.ncirc:c.ncirc + m.ns], 'Ip': x[:, self.ip], 'Te': self.Te,
                'flux_bound_Wb': 2 * np.pi * self.s * (self.Cbound @ x.T).ravel(),
                'E_pts': (self.s * E).T, 'B_pts': B.T, 'B_pts_ext': np.hypot(self.CBr @ x_ext.T, self.CBz @ x_ext.T).T,
                'lloyd_pts': (self.s * E * self.Bphi[:, None] / np.maximum(B, 1e-12)).T,
                'Bz_ext_c': Bz_c, 'Br_ext_c': (self.Cbr_c @ x.T).ravel(), 'n_c': n_c, 'K': self.K}


# ============================================================================ stage 2
def interp_keys(keys, t):
    """[[t, value, ...], ...] keyframes -> values at the times t (linear interpolation), [len(t), n_values]."""
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
                  for i, n in enumerate(machine.coil_names)] + [(s['cells'], 1.0, True) for s in machine.segments]
        for i, (cells, nt, inv_r) in enumerate(groups):
            w = area[cells] / cen[cells, 0] if inv_r else area[cells].copy()
            pts.append(cen[cells])
            rows += [i] * len(cells)
            cols += list(range(npt, npt + len(cells)))
            vals += list(nt * w / w.sum())
            npt += len(cells)
        self.pts = np.vstack(pts)
        self.S = sp.csr_matrix((vals, (rows, cols)), shape=(len(groups), npt))

    def linkage(self, gs, psi_nodal):
        """Flux linkage [Wb] of the nodal flux psi_nodal [Wb/rad] with every unit (coils, then segments)."""
        return 2 * np.pi * (self.S @ eval_fields(gs, psi_nodal[None, :], self.pts)['psi'][0])


class GSRamp:
    """Stage 2: a ramp of TokaMaker free-boundary equilibria from the stage-1 hand-off, I_p(t) and the
    boundary shape prescribed per slice, optimised jointly with stage 1 by Picard iteration:
      1. TokaMaker solves every slice for its target shape, coil and structure currents pinned to the
         current trajectory -> the plasma's own flux psi_p per slice.
      2. One convex program in all ramp voltages AND all stage-1 variables, psi_p held fixed: the plasma
         is a known source in the circuit (its changing flux linkage with every unit), the isoflux targets
         are least-squares terms, a 0-D Ohm's law sets the boundary flux, and wall points away from the
         target's contact point stay outside the boundary. Stage 1 keeps its constraints (the slew limit
         also covers the ramp and hand-off); its flux objective becomes a floor (flux_slack_Wb below its
         stand-alone optimum).
      3. Repeat until the currents, TokaMaker's response and the plasma flux stop changing.
    Slices the seeding TokaMaker pass cannot solve get psi_p interpolated from their neighbours for the
    first convex program, and must solve in a later pass.

    Scenario keys of scen['ramp'] (required unless a default is given):
      time {'t_start', 't_end', 'dt'} [s]   slice times (t_start = stage 1's t_end)
      Ip                  [[t, I_p], ...] keyframes [s, A]
      shape               (default None) [[t, R0, Z0, a, kappa, delta], ...] [s, m, m, m, -, -]: targets from
                          create_isoflux, and TokaMaker's starting guess
      targets             (default None) boundary points per slice, a list or {slice: [[R, Z], ...]} [m];
                          override shape on the slices given
      n_iso               isoflux points per slice from shape (needed on slices without targets)
      iso_exact           (default False) isoflux targets as equalities: a bool, a per-slice list or {slice: bool}
      xpoints             (default None) x-point targets {slice: [[R, Z], ...]} [m]
      xpoint_exact        (default True) B_p = 0 there exactly; if False least squares with weight
                          1 / xpoint_tol [T] (xpoint_tol then required)
      ohm_scheme          (default 'trapezoid') 0-D Ohm's law: 'trapezoid', 'backward' (end-of-interval
                          R_p I_p) or 'current_weighted' (flux averaged with weights J_phi dA / I_p)
      Rp                  (default None) [[t, R_p], ...] [s, Ohm]; else Spitzer from
      Te, Zeff            [[t, T_e], ...] [s, eV] and Z_eff [-] (read only when Rp is None)
      ffp, pp             TokaMaker profiles create_power_flux_fun(40, *ffp) and (40, *pp)
      Ip_ratio            TokaMaker Ip_ratio target
      iso_weight          TokaMaker weight of the isoflux and saddle constraints
      seed_weight, pin_weight   TokaMaker coil-pin weights, seeding pass and later passes
      shape_tol [m]       boundary displacement that costs 1 in the convex program
      contact_tol [m]     boundary-flux / contact-point mismatch, as a displacement, that costs 1
      limiter_margin [m]  how far the rest of the wall must stay outside the boundary
      limiter_band [-]    wall points constrained: flux within limiter_band margins of the boundary
      contact_exclude [m] wall points this close to the target's innermost point are not constrained
      dV, V_reg [-]       weights of voltage smoothing (ramp, hand-off, stage-1 channel phase) and size
      relax [-]           under-relaxation of the plasma flux between Picard passes
      prox_weight [-]     weight of the proximal pull toward the previous Picard iterate
      joint {'flux_slack_Wb'} [Wb]   allowed drop of the hand-off flux below its stand-alone optimum
      final_flux          (default None) {'target' [Wb], 'scale' [Wb]}: least-squares target on 2 pi Fb_N
    The supply slew limit is stage 1's scen['slew'].
    """

    def __init__(self, machine, circuit, scen, stage1, log=print):
        """circuit: a plasma-free CircuitModel with stage 1's circuits; stage1: the SOLVED StartupProblem
        the ramp starts from (and re-optimises); log: function called with progress messages."""
        from OpenFUSIONToolkit.TokaMaker.util import create_isoflux
        self.m, self.c, self.scen, self.log, self.stage1 = machine, circuit, scen, log, stage1
        rp, tc = scen['ramp'], scen['ramp']['time']
        n = int(round((tc['t_end'] - tc['t_start']) / tc['dt']))
        self.t = np.linspace(tc['t_start'], tc['t_end'], n + 1)
        self.Ip = interp_keys(rp['Ip'], self.t).ravel()
        self.s = 1.0                                     # positive I_p only, as in stage 1
        if (abs(stage1.t[-1] - self.t[0]) > 1e-9
                or abs(stage1.s * stage1.x[-1, stage1.ip] - self.Ip[0]) > 0.02 * self.Ip[0]):
            raise ValueError('stage 1 must end where the ramp starts (same time and I_p)')
        self.sampler = UnitSampler(machine)
        gs, pu, om = machine.gs, machine.psi_u, circuit.output_matrix
        # boundary targets: explicit points per slice, else create_isoflux from the shape keyframes
        explicit = rp.get('targets') or {}
        explicit = dict(enumerate(explicit)) if isinstance(explicit, (list, tuple)) else explicit
        shape = interp_keys(rp['shape'], self.t) if rp.get('shape') is not None else None
        if shape is None and len(explicit) < n + 1:
            raise ValueError("give shape keyframes, or explicit targets for every slice")
        self.iso = [np.asarray(explicit[k], float) if k in explicit else create_isoflux(rp['n_iso'], *shape[k])
                    for k in range(n + 1)]
        # (R0, Z0, a, kappa, delta) per slice, for TokaMaker's starting guess only: keyframes, or from the targets
        self.shape = shape if shape is not None else np.array(
            [[i[:, 0].mean(), i[:, 1].mean(), np.ptp(i[:, 0]) / 2,
              max(np.ptp(i[:, 1]) / max(np.ptp(i[:, 0]), 1e-3), 1.0), 0.0] for i in self.iso])
        ex = rp.get('iso_exact', False)
        ex = dict(enumerate(ex)) if isinstance(ex, (list, tuple)) else ex
        self.iso_exact = [bool(ex.get(k, False)) if isinstance(ex, dict) else bool(ex) for k in range(n + 1)]
        # x-points per slice and the (B_R, B_Z) maps there
        self.xpts, self.G_xb = [None] * (n + 1), [None] * (n + 1)
        for k, v in (rp.get('xpoints') or {}).items():
            self.xpts[k] = np.atleast_2d(np.asarray(v, float))
            g = eval_fields(gs, pu, self.xpts[k])
            self.G_xb[k] = (om(g['Br'].T), om(g['Bz'].T))
        self.ohm_scheme = rp.get('ohm_scheme', 'trapezoid')
        if self.ohm_scheme not in ('trapezoid', 'backward', 'current_weighted'):
            raise ValueError(f'unknown ohm_scheme {self.ohm_scheme!r}')
        if self.ohm_scheme == 'current_weighted':        # plasma cells: TokaMaker region 1
            area, cen = cell_geometry(*machine.mesh[:2])
            cells = np.where(machine.mesh[2] == 1)[0]
            self.pl_pts, self.pl_area, self.wJ = cen[cells], area[cells], [None] * (n + 1)
            self.C_pl = om(eval_fields(gs, pu, self.pl_pts)['psi'].T)
        # flux maps (state -> psi) at the targets, the wall, and each slice's contact point: the wall point
        # nearest the target's innermost point
        self.G_iso = [om(eval_fields(gs, pu, iso)['psi'].T) for iso in self.iso]
        self.lim_pts = max(gs.lim_contours, key=len)[:-1, :2]
        self.contact = [self.lim_pts[np.argmin(np.hypot(*(self.lim_pts - iso[np.argmin(iso[:, 0])]).T))]
                        for iso in self.iso]
        self.G_contact = [om(eval_fields(gs, pu, ct[None, :])['psi'].T) for ct in self.contact]
        self.G_lim = om(eval_fields(gs, pu, self.lim_pts)['psi'].T)
        # joint optimisation: stage 1's program without its objective, which becomes a floor
        lin = stage1.P.linear_terms()
        L_star = float(sum((q @ stage1.P.value(v)).sum() for v, q in lin.items()))
        self.joint = {'P0': stage1.P.clone(objective=False), 'flux_star_Wb': stage1.flux_handoff(), 'lin': lin,
                      'L_cap': L_star + rp['joint']['flux_slack_Wb'] / (2 * np.pi),
                      'X0': {v: M[:circuit.nx] for v, M in stage1.handoff_terms().items()},
                      'h0': stage1.h[-1][:circuit.nx]}   # hand-off state x_0 = X0 . vars + h0
        self.x0 = stage1.x[-1][:circuit.nx].copy()

    # ----------------------------------------------------------- TokaMaker update
    def tm_update(self, I_units, w_pin, first=False):
        """Solve every slice with TokaMaker, coil + structure currents pinned to I_units [n_slices, nu] with
        weight w_pin; store the plasma's own flux psi_p = psi - psi_external (under-relaxed after the first
        pass). first: the seeding pass, marching in time (each slice pinned to the previous solved slice's
        currents). Returns the number of failed slices."""
        from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
        gs, m, rp = self.m.gs, self.m, self.scen['ramp']
        gs.set_profiles(ffp_prof=create_power_flux_fun(40, *rp['ffp']), pp_prof=create_power_flux_fun(40, *rp['pp']))
        gs.set_coil_bounds({n: [-1e9, 1e9] for n in m.unit_names})
        gs.settings.maxits = 300
        gs.update_settings()
        nfail, prev_I, dpsi = 0, None, []
        for k in range(len(self.t)):
            self.stale[k] = True
            Ik = I_units[k] if (not first or prev_I is None) else prev_I
            units = dict(zip(m.unit_names, Ik))
            gs.set_targets(Ip=float(abs(self.Ip[k])), Ip_ratio=rp['Ip_ratio'])
            gs.set_isoflux_constraints(self.iso[k], weights=np.full(len(self.iso[k]), rp['iso_weight']))
            gs.set_saddle_constraints(self.xpts[k], weights=None if self.xpts[k] is None
                                      else np.full(len(self.xpts[k]), rp['iso_weight']))
            gs.set_psi_dt(np.zeros(gs.np), 1.0e9)        # static equilibria (no eddy-current term)
            if self.psi_tot[k] is not None:
                gs.set_psi(self.psi_tot[k])
            elif first and k > 0 and self.psi_tot[k - 1] is not None:
                gs.set_psi(self.psi_tot[k - 1])
            else:
                gs.init_psi(*self.shape[k])
            for attempt, wfac in enumerate([1.0, 0.1, 0.01]):     # loosen the pin if TokaMaker fails
                gs.set_coil_reg(reg_terms=[gs.coil_reg_term({n: 1.0}, target=float(v), weight=w_pin * wfac)
                                           for n, v in units.items()])
                gs.set_coil_currents(units)
                if attempt > 0:
                    gs.init_psi(*self.shape[k])
                try:
                    gs.solve()
                    st = gs.get_stats(li_normalization='iter')
                    break
                except Exception as e:
                    err = e
            else:                                        # every attempt failed
                nfail += 1
                self.log(f'    TokaMaker slice {k} (t={self.t[k]:.2f} s) failed: {err}')
                continue
            if attempt > 0:
                self.log(f'    TokaMaker slice {k} (t={self.t[k]:.2f} s) solved only with the coil pin loosened {1 / wfac:.0f}x')
            cc, _ = gs.get_coil_currents()
            Itm = prev_I = np.array([cc[n] for n in m.unit_names])
            psi = self.psi_tot[k] = gs.get_psi(normalized=False)
            new = psi - m.psi_u.T @ Itm
            if self.psi_p[k] is not None:
                dpsi.append(np.abs(new - self.psi_p[k]).max() / max(np.abs(new).max(), 1e-12))
            alpha = rp['relax']                          # under-relax the plasma flux
            self.psi_p[k] = new if (first or self.psi_p[k] is None) else (1 - alpha) * self.psi_p[k] + alpha * new
            self.stale[k] = False
            xpts = gs.get_xpoints()[0]
            self.tm_info[k] = {'psi_bounds': [float(x) for x in gs.psi_bounds], 'Ip': float(st['Ip']),
                               'li3': float(st['l_i']), 'o_point': [float(v) for v in gs.o_point],
                               'L_i': float(MU0 * gs.o_point[0] * st['l_i'] / 2),
                               'ext_dev_rel': self.rel_change(Itm - Ik, Ik), 'pin_factor': wfac,
                               'diverted': bool(gs.diverted),
                               'x_points': np.zeros((0, 2)) if xpts is None else np.array(xpts),
                               'lim_point': [float(v) for v in gs.lim_point]}
            if self.ohm_scheme == 'current_weighted':    # weights J_phi dA / I_p over the plasma cells
                J = gs.get_field_eval('eval', values=np.ascontiguousarray(gs.calc_jtor_plasma())).eval(self.pl_pts)[:, 0]
                self.wJ[k] = J * self.pl_area / max(float((J * self.pl_area).sum()), 1e-30)
        self.dpsi_rel = float(max(dpsi)) if dpsi and not nfail else np.inf
        return nfail

    # ----------------------------------------------------------- joint convex program
    def solve_qp(self):
        """Build and solve the joint convex program (stage 1 + ramp), the plasma flux of the latest TokaMaker
        pass held fixed. On success re-load stage 1 and set self.V, self.Fb, self.x and the diagnostics
        shape_err_m [m] and flux_mismatch_Wb [Wb]. Returns the solver status."""
        c, m, joint, s1, rp = self.c, self.m, self.joint, self.stage1, self.scen['ramp']
        t, N, nv = self.t, len(self.t) - 1, c.ncirc
        Phi = np.array([self.sampler.linkage(m.gs, p) for p in self.psi_p])        # [N+1, nu]
        isoA = [eval_fields(m.gs, self.psi_p[k][None, :], self.iso[k])['psi'][0] for k in range(N + 1)]
        P = joint['P0'].clone()                          # stage 1's variables and constraints
        V = P.var('V_ramp', N * nv, np.tile(c.V_max, N))
        Fb = P.var('Fb', N + 1, 1.0)                     # boundary flux per slice [Wb/rad]
        eFb = -sp.eye(N + 1, format='csr')               # row k selects -Fb_k
        # x_k = Ph_k x_0 + Gv_k V + h_k, x_0 = X0 . (stage-1 variables) + h0; the plasma's changing flux
        # linkage enters as unit 'voltages' w = -dPhi/dt
        Gv, Ph, h = np.zeros((c.nx, N * nv)), np.eye(c.nx), np.zeros(c.nx)
        self.Gv, self.Ph, self.h = [Gv.copy()], [Ph.copy()], [h.copy()]
        for k in range(1, N + 1):
            A, B, Bw = c.zoh(t[k] - t[k - 1], forcing=c.U.T)
            Gv = A @ Gv
            Gv[:, (k - 1) * nv:k * nv] += B
            Ph, h = A @ Ph, A @ h + Bw @ (-(Phi[k] - Phi[k - 1]) / (t[k] - t[k - 1]))
            self.Gv.append(Gv.copy())
            self.Ph.append(Ph.copy())
            self.h.append(h.copy())

        def aff(C, k):                                   # output C x_k as ({Var: matrix}, constant)
            CP = C @ self.Ph[k]
            return {V: C @ self.Gv[k], **{v: CP @ M for v, M in joint['X0'].items()}}, C @ self.h[k] + CP @ joint['h0']

        # coil current and supply voltage limits on the ramp
        used = np.abs(c.T).sum(axis=1) > 0
        Tu, lo, hi = c.U[:m.nc][used], m.coil_lim[used, 0], m.coil_lim[used, 1]
        for k in range(1, N + 1):
            d, c0 = aff(Tu, k)
            P.le(d, hi - c0, tag='I_max')
            P.le({v: -M for v, M in d.items()}, c0 - lo, tag='I_min')
        P.le({V: sp.eye(N * nv)}, np.tile(c.V_hi, N), tag='V_max')
        P.le({V: -sp.eye(N * nv)}, -np.tile(c.V_lo, N), tag='V_min')
        # isoflux: limited at the contact point, the boundary is the flux surface through it, so psi
        # (external + plasma) at every target point = psi at the contact point. A flux error equal to a
        # displacement of shape_tol costs 1: dpsi = |grad psi| x displacement, |grad psi| = R B_p in the latest
        # TokaMaker solution, floored at a fifth of its mean
        psi_pc = [float(eval_fields(m.gs, self.psi_p[k][None, :], self.contact[k][None, :])['psi'][0, 0])
                  for k in range(N + 1)]                 # plasma flux at the contact points
        self.iso_grad = []
        for k in range(N + 1):
            f = eval_fields(m.gs, self.psi_tot[k][None, :], self.iso[k])
            gpsi = self.iso[k][:, 0] * np.hypot(f['Br'][0], f['Bz'][0])
            self.iso_grad.append(np.maximum(gpsi, 0.2 * gpsi.mean()))
            d, c0 = aff(self.G_iso[k] - np.ones((len(self.iso[k]), 1)) @ self.G_contact[k], k)
            if self.iso_exact[k]:
                P.eq(d, -(c0 + isoA[k] - psi_pc[k]), tag='isoflux_exact')
            else:
                W = sp.diags(1.0 / (self.iso_grad[k] * rp['shape_tol']))
                P.lsq({v: W @ sp.csr_matrix(M) for v, M in d.items()}, -(W @ (c0 + isoA[k] - psi_pc[k])), 1.0,
                      tag='isoflux')
        # x-points: B_p (external + plasma) = 0
        for k in range(N + 1):
            if self.xpts[k] is not None:
                fp = eval_fields(m.gs, self.psi_p[k][None, :], self.xpts[k])
                for C, b0, tag in ((self.G_xb[k][0], fp['Br'][0], 'xpt_Br'), (self.G_xb[k][1], fp['Bz'][0], 'xpt_Bz')):
                    d, c0 = aff(C, k)
                    if rp.get('xpoint_exact', True):
                        P.eq(d, -(c0 + b0), tag=tag)
                    else:
                        P.lsq(d, -(c0 + b0), 1.0 / rp['xpoint_tol'], tag=tag)
        # limited where intended: wall points farther than contact_exclude from the target's innermost point
        # keep their flux below Fb by the flux equivalent of limiter_margin. Only points within limiter_band
        # margins of the boundary flux in the latest TokaMaker solution are constrained (the others cannot
        # reach the boundary in one Picard step)
        for k in range(N + 1):
            margin = float(self.iso_grad[k].mean()) * rp['limiter_margin']     # |grad psi| x distance
            psi_wall = eval_fields(m.gs, self.psi_tot[k][None, :], self.lim_pts)['psi'][0]
            far = ((psi_wall > self.tm_info[k]['psi_bounds'][0] - rp['limiter_band'] * margin)
                   & (np.hypot(*(self.lim_pts - self.iso[k][np.argmin(self.iso[k][:, 0])]).T) > rp['contact_exclude']))
            if far.any():
                psi_pl = eval_fields(m.gs, self.psi_p[k][None, :], self.lim_pts[far])['psi'][0]
                d, c0 = aff(self.G_lim[far], k)
                P.le({**d, Fb: eFb[[k] * int(far.sum())]}, -margin - c0 - psi_pl, tag='limiter')
        # smooth voltages through the ramp, across the hand-off and over the stage-1 channel phase; slew limit
        Vs, V1, N1 = sp.diags(1.0 / c.V_max), s1.V, s1.V.n // nv
        D = sp.diags([-1.0, 1.0], [0, 1], shape=(N - 1, N))            # row j: x_{j+1} - x_j
        P.lsq({V: sp.kron(D, Vs)}, 0.0, rp['dV'], tag='reg_dV')
        E1 = sp.hstack([sp.eye(nv), sp.csr_matrix((nv, (N - 1) * nv))]).tocsr()     # first ramp interval
        EL = sp.hstack([sp.csr_matrix((nv, (N1 - 1) * nv)), sp.eye(nv)]).tocsr()    # last stage-1 interval
        P.lsq({V: Vs @ E1, V1: -Vs @ EL}, 0.0, rp['dV'], tag='smooth_handoff')
        slew = s1.scen.get('slew')
        if slew:
            cap = c.V_max * slew * (t[1] - t[0])
            Ds = sp.kron(D, sp.eye(nv))
            P.le({V: Ds}, np.tile(cap, N - 1), tag='slew_up')
            P.le({V: -Ds}, np.tile(cap, N - 1), tag='slew_dn')
            P.le({V: E1, V1: -EL}, cap, tag='slew_handoff_up')
            P.le({V: -E1, V1: EL}, cap, tag='slew_handoff_dn')
        # stage-1 steps V_{j+1} - V_j (rows j nv..(j+1) nv) whose earlier interval (t_j, t_{j+1}] has the channel
        rows = [j * nv + i for j in range(N1 - 1) if s1.t[j] >= s1.t_p0 - 1e-12 for i in range(nv)]
        D1 = sp.diags([-1.0, 1.0], [0, 1], shape=(N1 - 1, N1))
        P.lsq({V1: sp.kron(D1, Vs).tocsr()[rows]}, 0.0, rp['dV'], tag='smooth_channel')
        if joint['lin']:                                 # hand-off flux >= stand-alone optimum - slack
            P.le(joint['lin'], [joint['L_cap']], tag='objective_floor')
        P.lsq({V: sp.kron(sp.eye(N), Vs)}, 0.0, rp['V_reg'], tag='reg_V')
        if self.x_prev is not None:                      # proximal pull toward the previous Picard iterate
            Isc = np.abs(m.coil_lim[used]).max(axis=1)
            for k in range(1, N + 1):
                d, c0 = aff(Tu / Isc[:, None], k)
                P.lsq(d, (Tu @ self.x_prev[k]) / Isc - c0, rp['prox_weight'] * self.prox_factor, tag='prox')
        # Fb_k is the total flux at the contact point (Fb_0 that of the hand-off state): least squares, a
        # mismatch equal to a displacement of contact_tol costs 1 (equalities condition the solver badly)
        for k in range(N + 1):
            d, c0 = aff(self.G_contact[k], k)
            f = eval_fields(m.gs, self.psi_tot[k][None, :], self.contact[k][None, :])
            gpsi = max(self.contact[k][0] * float(np.hypot(f['Br'][0, 0], f['Bz'][0, 0])),
                       0.2 * float(self.iso_grad[k].mean()))
            P.lsq({**d, Fb: eFb[[k]]}, -(c0 + psi_pc[k]), 1.0 / (gpsi * rp['contact_tol']), tag='Fb_contact')
        # 0-D Ohm's law
        self.Rp_gs, self.Li_gs = self.plasma_RL()
        if self.ohm_scheme == 'current_weighted':
            # sum_x w_k(x) [psi_k(x) - psi_{k-1}(x)] = -R_p,k I_p,k dt / (2 pi), psi = external + plasma,
            # end-of-interval weights w_k = J_phi dA / I_p; Fb then follows from the contact-point term
            psi_pl = [eval_fields(m.gs, self.psi_p[k][None, :], self.pl_pts)['psi'][0] for k in range(N + 1)]
            for k in range(1, N + 1):
                w = self.wJ[k][None, :]
                (d1, c1), (d0, c0) = aff(w @ self.C_pl, k), aff(w @ self.C_pl, k - 1)
                rhs = (-self.s * self.Rp_gs[k] * abs(self.Ip[k]) * (t[k] - t[k - 1]) / (2 * np.pi)
                       - float((w @ (psi_pl[k] - psi_pl[k - 1]))[0]) - (c1[0] - c0[0]))
                P.eq({v: d1[v] - d0[v] for v in d1}, [rhs], tag='ohm_weighted')
        else:
            # energy form V_b I_p = R_p I_p^2 + d/dt(L_i I_p^2 / 2):
            #   2 pi (Fb_k - Fb_{k-1}) = -[R I dt + (L_i,k I_k^2 - L_i,k-1 I_k-1^2) / (I_k + I_k-1)],
            # R, I weighted th at the end of the interval (th = 1/2 trapezoid, 1 backward)
            th = 1.0 if self.ohm_scheme == 'backward' else 0.5
            Dfb = sp.diags([-1.0, 1.0], [0, 1], shape=(N, N + 1), format='csr')
            for k in range(1, N + 1):
                I1, I0 = abs(self.Ip[k]), abs(self.Ip[k - 1])
                res = ((th * self.Rp_gs[k] + (1 - th) * self.Rp_gs[k - 1]) * (th * I1 + (1 - th) * I0)
                       * (t[k] - t[k - 1]))
                ind = (self.Li_gs[k] * I1 ** 2 - self.Li_gs[k - 1] * I0 ** 2) / max(I1 + I0, 1e-9)
                P.eq({Fb: Dfb[[k - 1]]}, [-self.s * (res + ind) / (2 * np.pi)], tag='ohm_0d')
        ff = rp.get('final_flux')
        if ff:                                           # ((2 pi Fb_N - target) / scale)^2
            P.lsq({Fb: 2 * np.pi * np.eye(N + 1)[[N]]}, [ff['target']], 1.0 / ff['scale'], tag='final_flux')
        st = P.solve()
        self.P, self.status = P, st
        if P.accepted():
            s1.load(P)                                   # stage 1 re-optimised: new hand-off state
            self.x0 = s1.x[-1][:c.nx].copy()
            vv = P.value(V)
            self.V, self.Fb = vv.reshape(N, nv), P.value(Fb)
            self.x = np.array([self.Ph[k] @ self.x0 + self.h[k] + self.Gv[k] @ vv for k in range(N + 1)])
            psi_c = [float((self.G_contact[k] @ self.x[k])[0]) + psi_pc[k] for k in range(N + 1)]
            self.shape_err_m = [float((np.abs(self.G_iso[k] @ self.x[k] + isoA[k] - psi_c[k]) / self.iso_grad[k]).max())
                                for k in range(N + 1)]
            # boundary flux the coils deliver (contact point) minus Ohm's law [Wb] (0 with 'current_weighted')
            self.flux_mismatch_Wb = [float(2 * np.pi * (psi_c[k] - self.Fb[k])) for k in range(N + 1)]
        return st

    def plasma_RL(self):
        """Plasma resistance R_p [Ohm] per slice (R_p keyframes, else Spitzer from the T_e keyframes on an
        elliptical cross-section R0, a, kappa of the slice's shape) and internal inductance
        L_i = mu0 R l_i(3) / 2 [H] from each slice's TokaMaker equilibrium."""
        rp = self.scen['ramp']
        Li = np.array([i['L_i'] for i in self.tm_info])
        if rp.get('Rp') is not None:
            return interp_keys(rp['Rp'], self.t).ravel(), Li
        Te = interp_keys(rp['Te'], self.t).ravel()
        R0, a, kap = self.shape[:, 0], self.shape[:, 2], self.shape[:, 3]
        eta = 5.2e-5 * rp['Zeff'] * LN_LAMBDA / np.maximum(Te, 1.0) ** 1.5
        return eta * 2 * np.pi * R0 / (np.pi * a ** 2 * kap), Li

    def rel_change(self, dI, I_ref):
        """Largest relative change: coils relative to their limits, structure relative to the largest
        segment current (at least 1 A)."""
        nc = self.m.nc
        dI, I_ref = np.atleast_2d(dI), np.atleast_2d(I_ref)
        return max(float((np.abs(dI[:, :nc]) / np.abs(self.m.coil_lim).max(axis=1)).max()),
                   float(np.abs(dI[:, nc:]).max() / max(np.abs(I_ref[:, nc:]).max(), 1.0)))

    def unit_currents(self, x):
        """States [n, nx] -> unit (coil and segment) currents [n, nu]."""
        return (self.c.U @ x.T).T

    def voltage_steps(self):
        """Largest |V_j - V_j-1| / V_max per circuit from the end of the breakdown window to the end of the ramp."""
        s1 = self.stage1
        k0 = int(np.searchsorted(s1.t, s1.t_p0 - 1e-12))
        return np.abs(np.diff(np.vstack([s1.v[k0:], self.V]), axis=0)).max(axis=0) / self.c.V_max

    # ----------------------------------------------------------- Picard driver
    def run(self, max_iter=15, tol=2e-3, measure=None, shape_tol_m=None):
        """Picard iteration between TokaMaker and the joint convex program; returns True when converged:
        every TokaMaker slice solves, the currents change by < tol, TokaMaker's external currents stay
        within 5 tol of the convex program's, and the raw plasma flux changes by < 5 tol. measure: optional
        function (ramp, k) -> boundary error [m] of slice k on TokaMaker's solution (e.g. with
        startup_plots.lcfs_distance); its maximum is kept in self.history (one record per iteration) and,
        if shape_tol_m is given, must also be below it."""
        rp, n = self.scen['ramp'], len(self.t)
        self.x_prev, self.prox_factor = None, 1.0
        self.psi_p, self.psi_tot, self.tm_info = [None] * n, [None] * n, [None] * n
        self.stale, self.dpsi_rel = [True] * n, np.inf
        I_prev = self.unit_currents(np.repeat(self.x0[None, :], n, axis=0))
        self.history = []
        gs = self.m.gs
        maxits0 = gs.settings.maxits
        try:
            self.tm_update(I_prev, w_pin=rp['seed_weight'], first=True)
            # slices the seeding pass could not solve: interpolate from the solved neighbours as a starting
            # guess for the first convex program; they stay unsolved and must solve in a later pass
            ok = [k for k in range(n) if self.psi_p[k] is not None]
            for k in range(n):
                lo, hi = [j for j in ok if j < k], [j for j in ok if j > k]
                if self.psi_p[k] is not None or not lo or not hi:
                    continue                             # solved, or nothing to interpolate from (stops below)
                lo, hi = lo[-1], hi[0]
                w = (self.t[k] - self.t[lo]) / (self.t[hi] - self.t[lo])
                self.psi_p[k] = (1 - w) * self.psi_p[lo] + w * self.psi_p[hi]
                self.psi_tot[k] = (1 - w) * self.psi_tot[lo] + w * self.psi_tot[hi]
                a, b = self.tm_info[lo], self.tm_info[hi]
                self.tm_info[k] = {**a, 'L_i': (1 - w) * a['L_i'] + w * b['L_i'], 'ext_dev_rel': np.inf,
                                   'psi_bounds': [(1 - w) * x + w * y for x, y in zip(a['psi_bounds'], b['psi_bounds'])]}
                if self.ohm_scheme == 'current_weighted':
                    self.wJ[k] = (1 - w) * self.wJ[lo] + w * self.wJ[hi]
                self.log(f'    TokaMaker slice {k} (t={self.t[k]:.2f} s): starting guess interpolated from slices {lo} and {hi}')
            for it in range(max_iter):
                if any(p is None for p in self.psi_p):
                    self.log('  a TokaMaker slice is missing; stopping')
                    return False
                t_qp = time.time()
                self.prox_factor = 1.0
                st = self.solve_qp()
                # a failed convex step: retry with a 10x, then 100x stronger pull to the previous iterate
                while not self.P.accepted() and self.x_prev is not None and self.prox_factor < 100:
                    self.prox_factor *= 10
                    self.log(f'  Picard {it}: convex program {st}; retrying with the proximal weight x{self.prox_factor:.0f}')
                    st = self.solve_qp()
                t_qp = time.time() - t_qp
                if not self.P.accepted():
                    self.log(f'  Picard {it}: convex program {st} (solver attempts {self.P.solve_log}); stopping')
                    return False
                I_new = self.unit_currents(self.x)
                self.x_prev = self.x.copy()
                dI = self.rel_change(I_new - I_prev, I_new)
                t_tm = time.time()
                nf = self.tm_update(I_new, w_pin=rp['pin_weight'])
                t_tm = time.time() - t_tm
                devs = [i['ext_dev_rel'] for k, i in enumerate(self.tm_info) if i and k > 0 and not self.stale[k]]
                dev = max(devs) if devs else np.inf
                h = {'iter': it, 'status': st, 'dI_rel': dI, 'tm_dev': dev, 'tm_fail': nf, 'dpsi_rel': self.dpsi_rel,
                     'shape_err_max_cm': 100 * max(self.shape_err_m), 'flux_handoff_Wb': self.stage1.flux_handoff(),
                     'dV_max': float(self.voltage_steps().max()), 't_qp': t_qp, 't_tm': t_tm}
                if measure is not None and nf == 0:
                    h['shape_err_measured_cm'] = 100 * max(float(measure(self, k)) for k in range(n))
                self.history.append(h)
                measured = f" (measured {h['shape_err_measured_cm']:.2f} cm)" if 'shape_err_measured_cm' in h else ''
                self.log(f"  Picard {it}: {st:8s} current change {dI:.1e}, TokaMaker deviation {dev:.1e}, "
                         f"plasma flux change {self.dpsi_rel:.1e}, max shape error {h['shape_err_max_cm']:.2f} cm"
                         + measured + f" (convex {t_qp:.0f} s, TokaMaker {t_tm:.0f} s)")
                I_prev = I_new
                shape_ok = shape_tol_m is None or h.get('shape_err_measured_cm', np.inf) <= 100 * shape_tol_m
                if nf == 0 and dI < tol and dev < 5 * tol and self.dpsi_rel < 5 * tol and shape_ok:
                    return True
            return False
        finally:
            # undo the ramp's settings (the coil regularisation cannot be cleared through the API)
            gs.set_isoflux_constraints(None)
            gs.set_saddle_constraints(None)
            gs.set_targets()
            gs.set_coil_bounds(None)
            gs.set_psi_dt(np.zeros(gs.np), 1.0e9)
            gs.settings.maxits = maxits0
            gs.update_settings()

    def handoff_jump(self):
        """Change of the plasma's flux linkage with the coils and structure at the hand-off, where stage 1's
        uniform ring becomes the first TokaMaker equilibrium (same I_p); not modelled by the ramp. If sudden,
        flux conservation (M dx = -U^T dPhi) would kick the currents by dx; returned relative to the coil
        limits and to the largest structure current at the hand-off."""
        s1, nc = self.stage1, self.m.nc
        Phi_gs = self.sampler.linkage(self.m.gs, self.psi_p[0])
        dPhi = Phi_gs - s1.plasma.M_u * s1.x[-1, s1.ip]
        dI = self.unit_currents((-self.c.Minv @ (self.c.U.T @ dPhi))[None, :])[0]
        I0 = self.unit_currents(self.x0[None, :])[0]
        return {'dPhi_rel': float(np.abs(dPhi).max() / np.abs(Phi_gs).max()),
                'coil_kick_rel': float((np.abs(dI[:nc]) / np.abs(self.m.coil_lim).max(axis=1)).max()),
                'struct_kick_rel': float(np.abs(dI[nc:]).max() / max(np.abs(I0[nc:]).max(), 1e-30))}

    def summary(self):
        """Headline numbers of the converged (or last accepted) joint solution [Wb, cm, -]."""
        if not hasattr(self, 'Fb'):
            raise RuntimeError('the ramp has no solution: the convex program was never solved (see the log above)')
        return {'flux_handoff_standalone_Wb': self.joint['flux_star_Wb'],
                'flux_handoff_Wb': self.stage1.flux_handoff(),
                'flux_boundary_handoff_Wb': float(2 * np.pi * self.Fb[0]),
                'flux_boundary_end_Wb': float(2 * np.pi * self.Fb[-1]),
                'ramp_flux_cost_Wb': float(2 * np.pi * (self.Fb[0] - self.Fb[-1])),
                'shape_err_max_cm': float(100 * max(self.shape_err_m)),
                'flux_mismatch_max_Wb': float(max(np.abs(self.flux_mismatch_Wb))),
                'dV_max_after_breakdown': float(self.voltage_steps().max()),
                'picard_iterations': len(self.history)}
