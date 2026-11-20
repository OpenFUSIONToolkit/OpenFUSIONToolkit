#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
"""Plasma initiation (burn-through) model for the Breakdown_Startup machine and circuit models, after
DYON: H.-T. Kim, PhD thesis, Imperial College London (2013); Kim et al., Nucl. Fusion 60 126049
(2020), 62 126012 (2022).

    AtomicSpecies           ADAS iso-nuclear rates and powers of one element, per charge state
    fractional_abundance    charge-state balance with a transport loss (coronal for tau -> inf)
    WallModel               recycling and sputtering yields (thesis ch. 5)
    PlasmaModel             0D energy and particle balance of electrons, ions and every charge state,
                            with neutral screening, transport and the plasma-current circuit
    PlasmaCircuit           the vessel (segment) currents and I_p driven by given coil currents,
                            integrated with the plasma model (K22 eq. 5-7)

Atomic data: OPEN-ADAS adf11 files (H.P. Summers, The ADAS User Manual), downloaded on first use
into ADAS_DIR ($OFT_ADAS_DIR, default ~/.cache/oft_breakdown_startup/adas). They are not shipped:
the OPEN-ADAS terms forbid redistributing the files with a modelling code.

SI units (m, s, m^-3, m^3/s, W m^3) except temperatures [eV].
"""
import os
import urllib.request

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator

EV = 1.602176634e-19        # J per eV
AMU = 1.66053907e-27        # kg
MU0 = 4e-7 * np.pi
ME = 9.1093837e-31           # kg
EPS0 = 8.8541878128e-12     # F/m
C_LIGHT = 2.99792458e8      # m/s
ETA_SPITZER = 5.2649e-5     # parallel Spitzer resistivity coefficient [Ohm m eV^1.5] (K20 sec. 2.3.1)
EQUI_COEF = 7.75e-34        # electron-ion equilibration [W m^3 eV^0.5 amu] (thesis eq. 3.33)

ADAS_DIR = os.environ.get('OFT_ADAS_DIR', os.path.join(os.path.expanduser('~'), '.cache', 'oft_breakdown_startup', 'adas'))
ADAS_URL = 'https://open.adas.ac.uk/download/adf11/{cls}/{name}.dat'

# adf11 files per element: scd ionization, acd recombination, ccd charge-exchange recombination by
# neutral hydrogen, plt line power, prb recombination + bremsstrahlung power. 96 (the GCR project)
# everywhere, as DYON (K20 sec. 2.3); there is no ccd96 for O and Be, so ccd89 there.
ADAS_FILES = {
    'h': {'scd': 'scd96_h', 'acd': 'acd96_h', 'ccd': 'ccd96_h', 'plt': 'plt96_h', 'prb': 'prb96_h'},
    'd': {'scd': 'scd96_h', 'acd': 'acd96_h', 'ccd': 'ccd96_d', 'plt': 'plt96_h', 'prb': 'prb96_h'},
    'be': {'scd': 'scd96_be', 'acd': 'acd96_be', 'ccd': 'ccd89_be', 'plt': 'plt96_be', 'prb': 'prb96_be'},
    'c': {'scd': 'scd96_c', 'acd': 'acd96_c', 'ccd': 'ccd96_c', 'plt': 'plt96_c', 'prb': 'prb96_c'},
    'n': {'scd': 'scd96_n', 'acd': 'acd96_n', 'ccd': 'ccd96_n', 'plt': 'plt96_n', 'prb': 'prb96_n'},
    'o': {'scd': 'scd96_o', 'acd': 'acd96_o', 'ccd': 'ccd89_o', 'plt': 'plt96_o', 'prb': 'prb96_o'},
}
# Ionization energies of charge state z -> z+1 [eV]: NIST Atomic Spectra Database (retrieved 2026-11-19)
IONIZATION_EV = {
    'h': [13.598434599702], 'd': [13.602134],
    'be': [9.322699, 18.21115, 153.896205, 217.71858459],
    'c': [11.2602880, 24.383143, 47.88778, 64.49352, 392.09056, 489.99320779],
    'n': [14.53413, 29.60125, 47.4453, 77.4735, 97.8901, 552.06741, 667.0461377],
    'o': [13.618055, 35.12112, 54.93554, 77.41350, 113.8990, 138.1189, 739.32697, 871.4099138],
}
T_RATE_MIN = 0.02           # rates are extrapolated (log-log) below the ADAS tables down to this [eV] (~room temperature)
MASS_AMU = {'h': 1.007825, 'd': 2.014102, 'be': 9.01218, 'c': 12.0107, 'n': 14.0067, 'o': 15.9994}


def adas_file(name):
    """Local path of OPEN-ADAS adf11 file `name` (eg. 'scd96_c'), downloaded into ADAS_DIR if missing."""
    path = os.path.join(ADAS_DIR, name + '.dat')
    if not os.path.isfile(path):
        os.makedirs(ADAS_DIR, exist_ok=True)
        url = ADAS_URL.format(cls=name.split('_')[0], name=name)
        with urllib.request.urlopen(url, timeout=60) as r:
            text = r.read().decode('latin-1')
        if 'DOCTYPE' in text[:200] or '/' not in text.split('\n')[0]:     # OPEN-ADAS answers a missing file with a web page
            raise FileNotFoundError(f'{name}: not an adf11 file at {url}')
        with open(path + '.part', 'w') as f:                  # atomic: a broken download leaves no cache file
            f.write(text)
        os.replace(path + '.part', path)
    return path


def read_adf11(path):
    """An adf11 (iso-nuclear) file: {'z_max', 'logne' [log10 m^-3], 'logT' [log10 eV], 'blocks'
    [n_block, n_T, n_ne] log10 SI (m^3/s or W m^3), 'z1' [n_block] (the file's Z1 per block)}."""
    lines = open(path, encoding='latin-1').read().split('\n')
    head = lines[0].split()
    z_max, n_ne, n_T = int(head[0]), int(head[1]), int(head[2])
    nums, blocks, z1, k = [], [], [], 2                 # line 1 is a separator
    while len(nums) < n_ne + n_T:
        nums += [float(x) for x in lines[k].split()]
        k += 1
    logne, logT = np.array(nums[:n_ne]) + 6.0, np.array(nums[n_ne:n_ne + n_T])     # cm^-3 -> m^-3
    while k < len(lines):
        line = lines[k]
        if line.startswith('C') or line.strip() == '':
            break
        if 'Z1=' in line:
            z1.append(int(line.split('Z1=')[1].split('/')[0]))
            vals = []
            k += 1
            while len(vals) < n_ne * n_T:
                vals += [float(x) for x in lines[k].split()]
                k += 1
            blocks.append(np.array(vals).reshape(n_T, n_ne) - 6.0)   # cm^3/s, W cm^3 -> m^3/s, W m^3
            continue
        k += 1
    return {'z_max': z_max, 'logne': logne, 'logT': logT, 'blocks': np.array(blocks), 'z1': np.array(z1)}


class AtomicSpecies:
    """ADAS rates and powers of one element ('h', 'd', 'be', 'c', 'n', 'o'), per charge state z.

    Attributes: name, Z (nuclear charge), mass_amu, W [Z] ionization energy z -> z+1 [eV].
    rates(Te, ne) evaluates, at electron temperature Te [eV] and density ne [m^-3] (arrays of one shape):
        'S'    [Z, ...]  ionization z -> z+1, z = 0..Z-1                     [m^3/s]
        'alpha' [Z, ...] recombination z -> z-1, entry z-1 for z = 1..Z       [m^3/s]
        'cx'   [Z, ...]  charge-exchange recombination with neutral H/D, z -> z-1, entry z-1   [m^3/s]
        'plt'  [Z, ...]  line power of charge state z, z = 0..Z-1             [W m^3]
        'prb'  [Z, ...]  recombination + bremsstrahlung power of z, entry z-1 [W m^3]
    Interpolation is linear in (log T, log ne). Below the lowest tabulated temperature (0.2 eV for the
    96 files, 1 eV for ccd89) it extrapolates linearly in log-log down to T_RATE_MIN, and clamps there;
    it clamps in density and above the table. ccd tables are evaluated at the temperature given (which
    temperature belongs there is a model choice). The ccd89 tables have no threshold for endothermic
    capture (eg. Be1+ + H -> Be + H+, 4.3 eV: 6e-16 m^3/s at 1-3 eV, where the ccd96 C1+ channel gives
    ~1e-21); for those channels ccd89 is multiplied by exp(-dE/T), dE = W_H - W_(z-1) > 0."""

    def __init__(self, name):
        self.name = name.lower()
        files = ADAS_FILES[self.name]
        self.W = np.array(IONIZATION_EV[self.name])
        self.Z, self.mass_amu = len(self.W), MASS_AMU[self.name]
        self.tables, self._interp = {}, {}
        for kind, fname in files.items():
            t = self.tables[kind] = read_adf11(adas_file(fname))
            if t['z_max'] != self.Z or list(t['z1']) != list(range(1, self.Z + 1)):
                raise ValueError(f'{fname}: blocks Z1 = {list(t["z1"])} for Z = {self.Z}')
            self._interp[kind] = [RegularGridInterpolator((t['logT'], t['logne']), b, bounds_error=False, fill_value=None)
                                  for b in t['blocks']]
        dE = IONIZATION_EV['h'][0] - self.W                      # endothermicity of capture by z+1 -> z [eV]
        self.cx_threshold = np.where(dE > 0, dE, 0.0) if files['ccd'].startswith('ccd89') else np.zeros(self.Z)

    def _eval(self, kind, Te, ne):
        t = self.tables[kind]
        lT = np.clip(np.log10(np.maximum(Te, 1e-30)), min(np.log10(T_RATE_MIN), t['logT'][0]), t['logT'][-1])
        ln = np.clip(np.log10(np.maximum(ne, 1e-30)), t['logne'][0], t['logne'][-1])
        lT, ln = np.broadcast_arrays(lT, ln)
        pts = np.stack([lT.ravel(), ln.ravel()], axis=-1)
        return np.array([10.0 ** f(pts) for f in self._interp[kind]]).reshape((-1,) + lT.shape)

    def rates(self, Te, ne):
        Te, ne = np.asarray(Te, float), np.asarray(ne, float)
        cx = self._eval('ccd', Te, ne)
        if self.cx_threshold.any():
            T = np.maximum(Te, T_RATE_MIN)
            cx = cx * np.exp(-self.cx_threshold.reshape((-1,) + (1,) * T.ndim) / T)
        return {'S': self._eval('scd', Te, ne), 'alpha': self._eval('acd', Te, ne), 'cx': cx,
                'plt': self._eval('plt', Te, ne), 'prb': self._eval('prb', Te, ne)}


def fractional_abundance(sp, Te, ne, tau=np.inf):
    """Steady-state charge-state fractions [Z+1] of species sp at Te [eV], ne [m^-3] with ion
    confinement time tau [s], every lost ion returning as a neutral (thesis eq. 5.2; coronal for
    tau = inf). Scalar Te, ne."""
    r = sp.rates(Te, ne)
    S, al = np.ravel(r['S']), np.ravel(r['alpha'])      # [Z] at the one (Te, ne)
    Z = sp.Z
    A = np.zeros((Z + 1, Z + 1))                         # dn/dt / ne = A n
    loss = 0.0 if np.isinf(tau) else 1.0 / (ne * tau)
    for z in range(Z + 1):
        if z < Z:                                        # ionization z -> z+1
            A[z, z] -= S[z]
            A[z + 1, z] += S[z]
        if z > 0:                                        # recombination z -> z-1, transport loss back to neutral
            A[z, z] -= al[z - 1] + loss
            A[z - 1, z] += al[z - 1]
            A[0, z] += loss
    A[-1, :] = 1.0                                       # replace one equation by sum(n) = 1
    rhs = np.zeros(Z + 1)
    rhs[-1] = 1.0
    return np.linalg.solve(A, rhs)


# ============================================================================ wall
CARBON_WALL = {('d', 'c'): 0.03, ('h', 'c'): 0.03, ('o', 'o'): 1.0, ('o', 'c'): 1.0, ('c', 'c'): 0.0}
"""Yields Y[(incident ion, released atom)] of a carbon wall (thesis table 5.1; K22 table 1): chemical
sputtering of C by D ions 0.03, O ions release CO (one O and one C), C ions are not recycled."""


class WallModel:
    """Plasma-wall interaction. fuel_recycling c = (c1, c2, c3): Y_D^D(t) = c1 - c2 (1 - exp(-t/c3))
    (thesis eq. 5.8; (1, 0, 1) is Y = 1); yields {(incident, released): Y} for impurity release by ions
    of any species (constant yields, thesis sec. 5.4.1). Unlisted pairs are 0."""

    def __init__(self, fuel_recycling=(1.0, 0.0, 1.0), yields=None):
        self.c, self.yields = fuel_recycling, dict(yields or {})

    def fuel_yield(self, t):
        c1, c2, c3 = self.c
        return c1 - c2 * (1.0 - np.exp(-max(t, 0.0) / c3))

    def Y(self, incident, released):
        return self.yields.get((incident, released), 0.0)


# ============================================================================ circuit (K22 eq. 5-7)
class PlasmaCircuit:
    """Passive (segment) currents I_s and the plasma current I_p driven by given coil currents
    (control-room mode, K22 eq. 5):

        [M_ss  M_sp] d [I_s]   [R_s   0 ] [I_s]     [M_sc] dI_c     [    0    ]
        [M_ps  L_p ] - [I_p] + [ 0   R_p] [I_p] = - [M_pc] ----  -  [dL_p/dt I_p]
                     dt                                     dt

    Coil currents I_coil [nt, nc] (per turn, linear between the times t_coil [s]); L_units [nu, nu] and
    R_seg [ns] the machine's unit inductances (per turn) and segment resistances; M_p [nu] the mutual
    inductance of 1 A of plasma current with 1 A per turn in each unit. L_p, R_p come from the plasma model.

    The plasma's position is fixed (M_p constant) unless move() sets M_p(t) linear over an interval;
    then the d M_p/dt terms of K22 eq. 5 are added: -dM_sp/dt I_p to the structure rows and
    -(dM_ps/dt I_s + dM_pc/dt I_c) to the plasma row. The loop voltage induced at the plasma by the coils
    and the structure (K22 eq. 7) is V_ext = -(M_pc dI_c/dt + M_ps dI_s/dt + dM_pc/dt I_c + dM_ps/dt I_s).
    """

    def __init__(self, L_units, R_seg, nc, M_p, t_coil, I_coil):
        self.nc, self.ns = nc, len(R_seg)
        self.M_ss, self.M_sc = L_units[nc:, nc:], L_units[nc:, :nc]
        self.R_s = np.asarray(R_seg, float)
        self.t_coil, self.I_coil = np.asarray(t_coil, float), np.asarray(I_coil, float)
        self.dI_coil = np.diff(self.I_coil, axis=0) / np.diff(self.t_coil)[:, None]
        self.K_ss = np.linalg.inv(self.M_ss)
        self.move(M_p, M_p, 0.0, 1.0)

    def move(self, M_a, M_b, ta, dt):
        """M_p [nu] linear from M_a at ta to M_b at ta + dt, held at M_b after (M_a at and before ta)."""
        M_a, M_b = np.asarray(M_a, float), np.asarray(M_b, float)
        self._M = (M_a, M_b, ta, dt, not np.array_equal(M_a, M_b))
        self._set(M_a, np.zeros_like(M_a))

    def _set(self, M_p, dM_p):
        self.M_p, self.dM_p = M_p, dM_p
        self.M_ps, self.M_pc = M_p[self.nc:], M_p[:self.nc]
        self.K_sp, self.M_ps_K = self.K_ss @ self.M_ps, self.M_ps @ self.K_ss
        self.schur0 = float(self.M_ps @ self.K_sp)             # M_ps M_ss^-1 M_sp

    def at(self, t):
        """Set M_p and dM_p/dt for time t (move())."""
        M_a, M_b, ta, dt, moving = self._M
        if moving:
            s = (t - ta) / dt
            self._set(M_a + min(max(s, 0.0), 1.0) * (M_b - M_a), (M_b - M_a) / dt * (0.0 <= s <= 1.0))

    @classmethod
    def from_machine(cls, machine, t_coil, I_coil, R=None, Z=0.0, channel=None):
        """From a startup_model.MachineModel: the plasma as a filament at (R, Z) (K22 eq. 9:
        M_pu = 2 pi psi_u(R, Z)), or as a ChannelPlasma (its disk-averaged M_u)."""
        if channel is not None:
            M_p = channel.M_u
        else:
            from startup_model import eval_fields
            M_p = 2 * np.pi * eval_fields(machine.gs, machine.psi_u, np.array([[R, Z]]))['psi'][:, 0]
        return cls(machine.L, machine.R[machine.nc:], machine.nc, M_p, t_coil, I_coil)

    def coil_rate(self, t):
        """dI_c/dt [nc] on the interval holding t (constant between the given times)."""
        k = int(np.clip(np.searchsorted(self.t_coil, t, side='right') - 1, 0, len(self.dI_coil) - 1))
        return self.dI_coil[k]

    def derivs(self, t, I_s, I_p, Rp, Lp, dLp):
        """(dI_s/dt, dI_p/dt, V_ext) at time t. Lp is the coefficient of dI_p/dt in d(L_p I_p)/dt
        (L_p + I_p dL_p/dI_p when the geometry depends on I_p), dLp the explicit dL_p/dt."""
        self.at(t)
        dIc = self.coil_rate(t)
        dM, V_move = self.dM_p, 0.0
        if self._M[4]:
            I_c = np.array([np.interp(t, self.t_coil, self.I_coil[:, j]) for j in range(self.nc)])
            V_move = -(dM[:self.nc] @ I_c + dM[self.nc:] @ I_s)
        r_s = -self.M_sc @ dIc - self.R_s * I_s - dM[self.nc:] * I_p
        r_p = -self.M_pc @ dIc - Rp * I_p - dLp * I_p + V_move
        if Lp <= self.schur0:
            raise ValueError(f'L_p = {Lp:.3e} H <= M_ps M_ss^-1 M_sp = {self.schur0:.3e} H: the plasma is not '
                             'decoupled from the structure (filament too close to it, or L_p too small)')
        dIp = (r_p - self.M_ps_K @ r_s) / (Lp - self.schur0)
        dIs = self.K_ss @ r_s - self.K_sp * dIp
        V_ext = -(self.M_pc @ dIc + self.M_ps @ dIs) + V_move
        return dIs, dIp, V_ext

    def vacuum(self, t0, t1, I_s0=None, rtol=1e-8):
        """Segment currents [ns] at t1 with no plasma, starting from I_s0 (default 0) at t0."""
        def f(t, I_s):
            return self.K_ss @ (-self.M_sc @ self.coil_rate(t) - self.R_s * I_s)
        x0 = np.zeros(self.ns) if I_s0 is None else I_s0
        sol = solve_ivp(f, (t0, t1), x0, method='BDF', rtol=rtol, atol=1e-6, max_step=np.min(np.diff(self.t_coil)) * 50)
        return sol.y[:, -1]


# ============================================================================ 0D plasma model
def plasma_inductance(R, a, kappa, li, form='circular'):
    """Plasma self-inductance [H] of an ellipse (R, a, kappa): 'circular' mu0 R (ln(8R/a) - 2 + li/2)
    (thesis eq. 3.13); 'k22' K22 eq. 14 as printed, mu0 R (ln(8R/(a kappa^2)) + li/2 - 2) (negative for
    the elongated early plasmas of DIII-D, kappa ~ 4 at r ~ 0.25 m: O32); 'sqrt_kappa' the textbook
    elongated form mu0 R (ln(8R/(a sqrt(kappa))) + li/2 - 2) (a sqrt(kappa): the radius of the ellipse's
    area); 'k22_st'
    K22 eq. 15-16 (low aspect ratio, after Hirshman & Neilson 1986)."""
    if form == 'circular':
        return MU0 * R * (np.log(8 * R / a) - 2 + li / 2)
    if form == 'k22':
        return MU0 * R * (np.log(8 * R / (a * kappa ** 2)) + li / 2 - 2)
    if form == 'sqrt_kappa':
        return MU0 * R * (np.log(8 * R / (a * np.sqrt(kappa))) + li / 2 - 2)
    if form == 'k22_st':
        e = a / R
        fA = (1 + 1.81 * np.sqrt(e) + 2.05 * e) * np.log(8 / e) - (2 + 9.25 * np.sqrt(e) - 1.21 * e)
        fB = 0.73 * np.sqrt(e) * (1 + 2 * e ** 4 - 6 * e ** 5 + 3.7 * e ** 6)
        return MU0 * R * (fA * (1 - e) / (1 - e + fB * kappa) + li / 2)
    raise ValueError(f'unknown inductance form {form!r}')


class PlasmaModel:
    """0D plasma initiation model (DYON: thesis ch. 3 and 5; K20; K22 eq. 6-8, 17).

    State (inventories, so particles are conserved exactly and a changing plasma volume dilutes
    densities as K22 eq. 8 does): for the fuel and each impurity the neutral atoms in the vessel and the
    ions of every charge state in the plasma, then the electron and ion thermal energies [J] of the
    plasma volume, then I_p [A]:
        y = [N_0, N_1..N_Z (fuel), N_0..N_Z (each impurity), E_e, E_i, I_p]
    Neutral density of species A: N_0 / (gamma_A V_V), with gamma_A = 1 - (V_p - V_n^A)/V_V the
    neutral-volume coefficient (thesis eq. 3.66); ion densities N_z / V_p.

    Arguments:
        fuel             'd' or 'h'
        impurities       sequence of 'c', 'o', 'n', 'be'
        wall             WallModel
        V_V              vessel volume [m^3]
        T_n              neutral (and charge-exchange product) temperature [eV]
        lnL              Coulomb logarithm (DYON: 10, K20 sec. 2.3.1)
        li               internal inductance; inductance: plasma_inductance form
        screening        neutral screening on (thesis sec. 3.4.6)
        fixed_T          (Te, Ti) to hold the temperatures fixed (tests), else None
        paux             callable t -> absorbed auxiliary power [W] (thesis eq. 7.4), or None
        fuelling         callable t -> external fuel atoms per second (already times the efficiency), or None
        circuit          PlasmaCircuit: I_p and the segment currents I_s (appended to y) driven by the
                         coil currents (M3, K22 eq. 5); None: I_p driven by a prescribed loop voltage
        T_cx             energy [eV] of the fuel atoms in the plasma (charge exchange, penetration); None: T_n
                         (DYON). With a D2 fill the atoms come from electron-impact dissociation and carry
                         Franck-Condon energies (~2-3 eV per atom); option, a departure (prefill study)
        E_dis            electron energy [eV] paid per fuel atom ionized for its dissociation from D2 (binding
                         2.26 eV + the Franck-Condon energy); 0: DYON (fully dissociated fill)
        neutrals         'uniform' (DYON: one well-mixed vessel reservoir, n_0 = N_0 / (gamma V_V), screening by
                         the penetration depth) or 'kinetic' (option, prefill study): neutrals enter the plasma at
                         the one-way thermal flux n_out v/4 across its surface A_p and are ionized inside, so the
                         density inside is n_in = n_out k / (k + ne S_0), k = v A_p / (4 V_p) (quasi-steady: the
                         neutral transit time a/v ~ 0.4 ms); N_0 = n_out (V_V - V_p) + n_in V_p; the fuel enters as
                         D2 (H2) molecules at T_n (v of the molecule, 2 atoms each), impurities as atoms;
                         or 'sol' (two-zone core/edge option): as 'kinetic', but the neutrals first cross a
                         scrape-off layer of open lines outside the plasma region, of width
                         delta = sqrt(D_perp L_sol / C_s) (D_perp = T_e / (8 B), DYON's tau_perp convention;
                         L_sol from geometry 'L_sol', the open lines next to the region, else L_sol below).
                         A fraction 1 - exp(-delta / lambda) (lambda = v / (ne S_0)) is ionized there; those
                         ions reach the wall in L_sol / C_s and recycle, so they do not fuel the core, but
                         the electrons pay their ionization and line radiation (P_sol). SOL T_e, n_e: the
                         core values (0D)
        L_sol            connection length [m] of the scrape-off layer when the geometry gives none
        ec               electron cyclotron heating (option): dict(P=callable t -> injected power [W], freq [Hz]
                         (110e9), harmonic (2), R_wall (wall reflectivity for multiple passes, 0.9), c (pi)).
                         Absorbed power P f with f = (1 - exp(-tau)) / (1 - R_wall exp(-tau)) and the single-pass
                         optical depth tau = c (w_pe / w)^2 (T_e / m_e c^2) (w R / c) (low-density X2, Bornatici
                         et al. 1983 scaling; the prefactor c is order unity); f = 0 unless the resonance
                         R_res = harmonic e B_phi R / (m_e w) lies in [R - a, R + a]. Added to paux
        dis_per_cx       also pay E_dis per charge exchange of a fuel atom (each CX replaces the atom by a newly
                         dissociated one: strict energy conservation, the pessimistic bound); False: per ionization
        shrink           when V_p falls: 'compress' (K22 eq. 8: the ions and their energy stay, so the
                         densities rise) or 'lose' (O38, a departure: the ions of the volume leaving the
                         plasma go to the wall with their thermal energy at |dV_p/dt| / V_p, recycled
                         like the transport losses; densities and temperatures are unchanged by it)

    The geometry is a callable geometry(t, Ip) -> dict with R, a, kappa, Vp [m, m, -, m^3], Bphi [T at R],
    and either Lf [m] (connection length: tau_par = Lf/Cs, thesis eq. 3.49) or tau_par [s]; the loop
    voltage is a callable vloop(t) [V] at the plasma when there is no circuit.
    """

    def __init__(self, fuel='d', impurities=('c', 'o'), wall=None, V_V=35.0, T_n=0.026, lnL=10.0, li=0.5,
                 inductance='circular', screening=True, fixed_T=None, paux=None, fuelling=None, circuit=None,
                 shrink='compress', T_cx=None, E_dis=0.0, dis_per_cx=False,
                 neutrals='uniform', L_sol=30.0, ec=None):
        self.L_sol = L_sol
        self.ec = None if ec is None else dict({'freq': 110e9, 'harmonic': 2, 'R_wall': 0.9, 'c': np.pi}, **ec)
        if neutrals not in ('uniform', 'kinetic', 'sol'):
            raise ValueError(f'unknown neutrals {neutrals!r}')
        self.T_cx, self.E_dis, self.dis_per_cx, self.neutrals = (T_n if T_cx is None else T_cx), E_dis, dis_per_cx, neutrals
        if shrink not in ('compress', 'lose'):
            raise ValueError(f'unknown shrink {shrink!r}')
        self.shrink = shrink
        self.species = [AtomicSpecies(fuel)] + [AtomicSpecies(x) for x in impurities]
        self.names = [s.name for s in self.species]
        self.wall = wall or WallModel()
        self.V_V, self.T_n, self.lnL, self.li, self.inductance = V_V, T_n, lnL, li, inductance
        self.screening, self.fixed_T, self.paux, self.fuelling = screening, fixed_T, paux, fuelling
        self.m_fuel = self.species[0].mass_amu * AMU
        self.slices, k = [], 0                               # slice of y holding species i (charge states 0..Z)
        for s in self.species:
            self.slices.append(slice(k, k + s.Z + 1))
            k += s.Z + 1
        self.iEe, self.iEi, self.iIp = k, k + 1, k + 2
        self.circuit, self.iIs = circuit, k + 3
        self.n_state = k + 3 + (0 if circuit is None else circuit.ns)
        self.vloop = None

    # ------------------------------------------------------------------ initial state
    def initial_state(self, p0, geo, Te0=1.0, Ti0=None, gamma0=None, impurity_fraction=None, Ip0=0.0, I_s0=None):
        """y0 for a fill of D2 (or H2) at pressure p0 [Pa] and temperature T_n, fully dissociated:
        n_fill = 2 p0 / T_n [atoms/m^3] (K20 eq. 1; K22 eq. 18). A fraction gamma0 (default 2e-3 Te0^1.5,
        K20 sec. 2.1) of the atoms in the plasma volume is ionized, the rest of the fill stays neutral.
        impurity_fraction {name: f}: f n_fill V_V neutral impurity atoms (eg. {'o': 0.001}).
        I_s0: segment currents at the start (with a circuit; eg. PlasmaCircuit.vacuum), default 0."""
        Ti0 = self.T_n if Ti0 is None else Ti0
        gamma0 = 2e-3 * Te0 ** 1.5 if gamma0 is None else gamma0
        n_fill = 2 * p0 / (self.T_n * EV)
        y = np.zeros(self.n_state)
        N1 = gamma0 * n_fill * geo['Vp']
        y[self.slices[0].start] = n_fill * self.V_V - N1
        y[self.slices[0].start + 1] = N1
        for i, s in enumerate(self.species[1:], 1):
            y[self.slices[i].start] = (impurity_fraction or {}).get(s.name, 0.0) * n_fill * self.V_V
        y[self.iEe] = 1.5 * EV * Te0 * N1
        y[self.iEi] = 1.5 * EV * Ti0 * N1
        y[self.iIp] = Ip0
        if self.circuit is not None and I_s0 is not None:
            y[self.iIs:] = I_s0
        self.n_fill = n_fill
        return y

    # ------------------------------------------------------------------ physics
    def terms(self, t, y):
        """Every term of the model at (t, y): densities, temperatures, rates, powers, sources, and the
        derivative 'dy'. Powers are per unit plasma volume [W/m^3]."""
        Ip = y[self.iIp]
        geo = self.geometry(t, Ip)
        R, a, kappa, Vp, Bphi = geo['R'], geo['a'], geo['kappa'], geo['Vp'], geo['Bphi']
        N = [np.maximum(y[sl], 0.0) for sl in self.slices]
        Ne = sum(np.dot(np.arange(len(n)), n) for n in N)
        Ni = sum(n[1:].sum() for n in N)
        ne = max(Ne / Vp, 1.0)
        if self.fixed_T is None:
            Te = max(y[self.iEe] / (1.5 * EV * max(Ne, 1.0)), 1e-3)
            Ti = max(y[self.iEi] / (1.5 * EV * max(Ni, 1.0)), 1e-4)
        else:
            Te, Ti = self.fixed_T
        r = [s.rates(Te, ne) for s in self.species]
        cx_rate = [s.rates(Ti, ne)['cx'] for s in self.species]   # ion-neutral reaction: at Ti (O23)
        # neutral screening (thesis sec. 3.4.6): penetration depth lambda = v_n / (ne S_0)
        Vn, gam = [], []
        for s, rs in zip(self.species, r):
            Tv = self.T_cx if s is self.species[0] else self.T_n
            v_n = np.sqrt(8 * Tv * EV / (np.pi * s.mass_amu * AMU))             # mean thermal speed (O8)
            lam = v_n / (ne * rs['S'][0])
            if self.screening and lam < a:
                Vn.append(Vp * (1 - (1 - lam / a) ** 2))     # = 2 pi R pi kappa (a^2 - (a-lam)^2) for a torus
            else:
                Vn.append(Vp)
            gam.append(1 - (Vp - Vn[-1]) / self.V_V)
        n0 = [n[0] / (g * self.V_V) for n, g in zip(N, gam)]
        sol_iz = [0.0] * len(self.species)                    # ionizations per second in the SOL (recycled)
        delta_sol = 0.0
        if self.neutrals in ('kinetic', 'sol'):               # influx-limited neutrals inside the plasma
            A_p = 4 * np.pi ** 2 * R * a * np.sqrt((1 + kappa ** 2) / 2)
            Vo = max(self.V_V - Vp, 1e-3)
            if self.neutrals == 'sol':
                L_sol = geo.get('L_sol', np.nan)
                L_sol = L_sol if np.isfinite(L_sol) and L_sol > 0 else self.L_sol
                Cs0 = np.sqrt((Te + Ti) * EV / self.m_fuel)
                delta_sol = np.sqrt(Te / (8 * abs(Bphi)) * L_sol / Cs0)
            for i, (s, rs) in enumerate(zip(self.species, r)):
                m_n = (2 if i == 0 else 1) * s.mass_amu * AMU
                v_in = np.sqrt(8 * self.T_n * EV / (np.pi * m_n))
                k_in = v_in * A_p / (4 * Vp)
                f_in = k_in / (k_in + ne * rs['S'][0])
                att = np.exp(-delta_sol * ne * rs['S'][0] / v_in)            # transmitted through the SOL
                n_out = N[i][0] / (Vo + Vp * att * f_in)
                n0[i] = n_out * att * f_in                   # n_in
                sol_iz[i] = (1 - att) * k_in * Vp * n_out
                Vn[i], gam[i] = Vp, 1.0
        nz = [n / Vp for n in N]                              # entry 0 unused for ions
        # transport (thesis eq. 3.49-3.54; K22 eq. 17)
        Cs = np.sqrt((Te + Ti) * EV / self.m_fuel)
        tau_par = geo['tau_par'] if 'tau_par' in geo else geo['Lf'] / Cs
        tau_perp = 8 * a ** 2 * abs(Bphi) / Te
        tau = 1.0 / (1.0 / tau_par + 1.0 / tau_perp)
        # resistance and inductance
        Zeff = sum(np.dot(np.arange(len(n)) ** 2, n) for n in N) / max(Ne, 1.0)
        eta = ETA_SPITZER * Zeff * self.lnL / max(Te, 0.05) ** 1.5         # [Ohm m]
        Rp = eta * 2 * R / (kappa * a ** 2)
        Lp = plasma_inductance(R, a, kappa, self.li, self.inductance)
        h, hI = 1e-6, max(1.0, 1e-6 * abs(Ip))
        g2, g3 = self.geometry(t + h, Ip), self.geometry(t, Ip + hI)
        dLp = (plasma_inductance(g2['R'], g2['a'], g2['kappa'], self.li, self.inductance) - Lp) / h      # dL_p/dt at fixed I_p
        dLdI = (plasma_inductance(g3['R'], g3['a'], g3['kappa'], self.li, self.inductance) - Lp) / hI
        Lp_eff = Lp + Ip * dLdI          # d(L_p I_p)/dt = (L_p + I_p dL_p/dI_p) dI_p/dt + (dL_p/dt) I_p
        # sources and powers
        f = 0                                                 # fuel index
        dN = [np.zeros(len(n)) for n in N]
        P_iz = P_rad = P_equi = P_cx = P_fc = P_sol = 0.0
        cx_sum = 0.0                                          # sum over impurity ions of cx n  [1/s] per fuel neutral
        for i, (s, rs) in enumerate(zip(self.species, r)):
            Z, W = s.Z, s.W
            S, al, plt, prb = rs['S'], rs['alpha'], rs['plt'], rs['prb']
            cx = cx_rate[i]
            if sol_iz[i] > 0:                                 # SOL ionization of the influx: energy only
                P_sol += sol_iz[i] * ((W[0] + (self.E_dis if i == 0 else 0.0)) * EV + plt[0] / S[0]) / Vp
            # ionization and recombination between z and z+1 (z = 0: neutrals in V_n)
            for z in range(Z):
                V_src = Vn[i] if z == 0 else Vp
                dens = n0[i] if z == 0 else nz[i][z]
                ion = V_src * S[z] * ne * dens                    # ionizations per second, z -> z+1
                rec = Vp * al[z] * ne * nz[i][z + 1]              # recombinations per second, z+1 -> z
                dN[i][z] += rec - ion
                dN[i][z + 1] += ion - rec
                P_iz += V_src * S[z] * (W[z] + (self.E_dis if (i == 0 and z == 0) else 0.0)) * EV * ne * dens / Vp
                if i == 0 and z == 0:                             # new fuel ions bring the atoms' energy (0 in DYON)
                    P_fc += 1.5 * (self.T_cx - self.T_n) * EV * ion / Vp
                P_rad += V_src * plt[z] * ne * dens / Vp
                P_rad += (prb[z] - al[z] * W[z] * EV) * ne * nz[i][z + 1]
            # equilibration (thesis eq. 3.33) and ion transport loss
            zz = np.arange(1, Z + 1)
            P_equi += EQUI_COEF * (Te - Ti) * ne * self.lnL / Te ** 1.5 * np.dot(zz ** 2, nz[i][1:]) / s.mass_amu
            dN[i][1:] -= N[i][1:] / tau
            if i == f:
                P_cx += (Vn[f] / Vp) * 1.5 * (Ti - self.T_cx) * EV * n0[f] * cx[0] * nz[f][1]
                continue
            # charge exchange of impurity ions with fuel neutrals (in V_n of the fuel): z -> z-1, D0 -> D+
            cxe = Vn[f] * cx * n0[f] * nz[i][1:]                  # [Z] events per second from charge z
            dN[i][1:] -= cxe
            dN[i][:-1] += cxe
            cx_sum += cxe.sum()
            # ion energy: only 1+ -> 0 removes a hot ion (replaced by a cold D+); for z >= 2 the impurity
            # ion stays in the plasma (thesis eq. 3.38 sums the 1+ states only; O25)
            P_cx += 1.5 * (Ti - self.T_cx) * EV * cxe[0] / Vp
        dN[f][0] -= cx_sum
        dN[f][1] += cx_sum
        if self.dis_per_cx:                                   # fuel-fuel and impurity CX events
            P_iz += self.E_dis * EV * (Vn[f] * n0[f] * cx_rate[f][0] * nz[f][1] + cx_sum) / Vp
        # wall: fuel recycling (+ external fuelling) and impurity release by every ion species
        out = [N[i][1:].sum() / tau for i in range(len(self.species))]       # ions to the wall [1/s]
        shrink = max(-(g2['Vp'] - Vp) / h, 0.0) / Vp if self.shrink == 'lose' else 0.0     # [1/s] (O38)
        for i in range(len(self.species)):
            dN[i][1:] -= N[i][1:] * shrink
            out[i] += N[i][1:].sum() * shrink
        dN[f][0] += self.wall.fuel_yield(t) * out[f]
        if self.fuelling is not None:
            dN[f][0] += self.fuelling(t)
        for i, s in enumerate(self.species[1:], 1):
            dN[i][0] += sum(self.wall.Y(A, s.name) * o for A, o in zip(self.names, out))
        # energy
        P_oh = Ip ** 2 * Rp / Vp
        P_aux = 0.0 if self.paux is None else self.paux(t) / Vp
        f_ec = 0.0
        if self.ec is not None:
            ec = self.ec
            w = 2 * np.pi * ec['freq']
            R_res = ec['harmonic'] * EV * abs(Bphi) * R / (ME * w)
            if R - a <= R_res <= R + a:
                tau_ec = ec['c'] * ne * EV ** 2 / (EPS0 * ME * w ** 2) * Te * EV / (ME * C_LIGHT ** 2) * w * R / C_LIGHT
                f_ec = (1 - np.exp(-tau_ec)) / (1 - ec['R_wall'] * np.exp(-tau_ec))
            P_aux += f_ec * ec['P'](t) / Vp
        P_econ = 1.5 * ne * Te * EV / tau
        ni = Ni / Vp
        P_icon = 1.5 * ni * Ti * EV / tau
        dy = np.concatenate(dN + [np.zeros(self.n_state - self.iEe)])
        if self.fixed_T is None:
            dy[self.iEe] = Vp * (P_oh + P_aux - P_iz - P_rad - P_equi - P_econ - P_sol) - shrink * y[self.iEe]
            dy[self.iEi] = Vp * (P_equi + P_fc - P_cx - P_icon) - shrink * y[self.iEi]
        if self.circuit is None:
            V = self.vloop(t)
            dy[self.iIp] = (V - Rp * Ip - dLp * Ip) / Lp_eff
        else:
            dIs, dy[self.iIp], V = self.circuit.derivs(t, y[self.iIs:], Ip, Rp, Lp_eff, dLp)
            dy[self.iIs:] = dIs
        return {'t': t, 'dy': dy, 'ne': ne, 'Te': Te, 'Ti': Ti, 'Ip': Ip, 'n0': n0, 'nz': nz, 'N': N,
                'Vn': Vn, 'gamma': gam, 'tau': tau, 'tau_par': tau_par, 'tau_perp': tau_perp, 'Zeff': Zeff,
                'Rp': Rp, 'Lp': Lp, 'Lp_eff': Lp_eff, 'dLp_dt': dLp, 'Vloop': V, 'P_oh': P_oh, 'P_aux': P_aux, 'P_iz': P_iz, 'P_rad': P_rad,
                'P_equi': P_equi, 'P_econ': P_econ, 'P_cx': P_cx, 'P_fc': P_fc, 'P_sol': P_sol, 'delta_sol': delta_sol, 'f_ec': f_ec, 'P_icon': P_icon, 'Vp': Vp,
                'shrink_rate': shrink}

    def rhs(self, t, y):
        return self.terms(t, y)['dy']

    def run(self, geometry, vloop, y0, t_span, t_eval=None, rtol=1e-6, max_step=np.inf):
        """Integrate from y0 over t_span [s] (BDF, the counterpart of DYON's ode15s); vloop [V] callable,
        or None with a circuit. Returns a dict of arrays over the output times: every scalar of terms()
        ('Vloop' is the voltage induced at the plasma by coils and structure with a circuit), plus
        'y' [n_t, n_state] and 'status'."""
        self.geometry, self.vloop = geometry, vloop
        scale = np.abs(y0).copy()
        for sl in self.slices:
            scale[sl] = max(np.abs(y0[sl]).sum(), 1.0)
        scale[self.iEe] = scale[self.iEi] = max(abs(y0[self.iEe]), 1e-12)
        scale[self.iIp] = 1.0
        if self.circuit is not None:
            scale[self.iIs:] = 1.0
        sol = solve_ivp(self.rhs, t_span, y0, method='BDF', t_eval=t_eval, rtol=rtol, atol=1e-9 * scale,
                        max_step=max_step)
        out = {'t': np.asarray(sol.t), 'y': np.asarray(sol.y).reshape(len(y0), -1).T, 'status': sol.status,
               'message': sol.message}
        if len(out['t']) == 0:                                       # failed before the first output time
            return out
        rows = [self.terms(tk, yk) for tk, yk in zip(out['t'], out['y'])]
        for k, v in rows[0].items():
            if np.isscalar(v):
                out[k] = np.array([row[k] for row in rows])
        return out


# ============================================================================ field lines and plasma volume (K22 sec. 2.2)
TORR = 133.322368           # Pa
TOWNSEND_A, TOWNSEND_B = 510.0, 1.25e4     # 1/(m Torr), V/(m Torr): K22 eq. 13 (hydrogen isotopes)


def townsend_field(p, L):
    """E_Townsend [V/m] for pressure p [Pa] over an open line of length L [m] (K22 eq. 13):
    1.25e4 p[Torr] / ln(510 p[Torr] L) = 93.758 p / ln(3.825 p L); inf where 510 p L <= 1, and where L is
    not finite (closed lines, or no line found: no Townsend breakdown)."""
    pt = np.asarray(p, float) / TORR
    L = np.asarray(L, float)
    x = TOWNSEND_A * pt * np.where(np.isfinite(L), L, 0.0)
    return np.where(x > 1.0, TOWNSEND_B * pt / np.log(np.maximum(x, 1.0 + 1e-12)), np.inf)


def filament_psi(R, Z, Rf, Zf):
    """Flux [Wb/rad] at (R, Z) of 1 A in a toroidal filament at (Rf, Zf) (Maxwell; = M / 2 pi)."""
    from scipy.special import ellipe, ellipk
    k2 = np.clip(4 * R * Rf / ((R + Rf) ** 2 + (Z - Zf) ** 2), 1e-14, 1 - 1e-9)
    k = np.sqrt(k2)
    return MU0 * np.sqrt(R * Rf) * ((2 / k - k) * ellipk(k2) - 2 / k * ellipe(k2)) / (2 * np.pi)


class FieldLineMap:
    """Poloidal flux, fields and toroidal electric field on a regular (R, Z) grid inside the first wall,
    and DYON's field-line analysis of them (K22 sec. 2.2): open-line lengths, line-averaged parallel
    field, the Townsend test, the plasma volume (open lines that broke down plus the closed region),
    its shape and the parallel confinement length.

    Conventions (as startup_model): psi [Wb/rad], B_R = -(dpsi/dZ)/R, B_Z = (dpsi/dR)/R,
    B_phi = F0/R, E_phi = -(dpsi/dt)/R.

    Construction: FieldLineMap(R [nR], Z [nZ], wall [n, 2] (R, Z) polygon, F0 [T m], psi_units
    [nu, nZ, nR] optional unit flux maps (1 A per turn in each unit)); from_machine builds the unit
    maps of a startup_model.MachineModel by TokaMaker interpolation.

    Lines and flux tubes (the 'contour' method, used for the volume): field lines lie on psi contours.
    Open or closed is decided per cell by topology: the contour through a cell is closed when its
    superlevel or its sublevel connected set does not reach the wall (exact on the grid, any |grad psi|,
    either sign; _closed_cells). The band of the grid between psi_j and psi_j + dpsi, split into
    connected pieces (open and closed cells apart), is a set of flux tubes; dpsi is at least two cells'
    psi change over 90 % of the open cells. An open tube takes the length L = int |B|/B_p dl_pol and the
    length-weighted parallel field <E_par> (K22 eq. 12: the mean over equal steps along B) of the
    band-centre contour crossing it, extended to the wall polygon; a tube no contour crosses has L
    unknown (nan) and cannot break down. Volume sum(2 pi R dA). trace() follows single lines with K22's
    Euler steps (eq. 11) or RK4, as reference.
    """

    def __init__(self, R, Z, wall, F0, psi_units=None):
        from matplotlib.path import Path
        self.R, self.Z = np.asarray(R, float), np.asarray(Z, float)
        self.dR, self.dZ = self.R[1] - self.R[0], self.Z[1] - self.Z[0]
        if not (np.allclose(np.diff(self.R), self.dR) and np.allclose(np.diff(self.Z), self.dZ)):
            raise ValueError('FieldLineMap needs a uniform grid')
        self.RR, self.ZZ = np.meshgrid(self.R, self.Z)
        self.wall = np.asarray(wall, float)[:, :2]
        self.inside = Path(self.wall).contains_points(np.column_stack([self.RR.ravel(), self.ZZ.ravel()])).reshape(self.RR.shape)
        self.F0, self.psi_units = F0, psi_units
        self.dA = self.dR * self.dZ
        pad = np.pad(self.inside, 1, constant_values=False)          # cells next to the wall (or grid edge)
        self.edge = self.inside & ~(pad[:-2, 1:-1] & pad[2:, 1:-1] & pad[1:-1, :-2] & pad[1:-1, 2:])

    @classmethod
    def from_machine(cls, machine, wall, dx=0.01, margin=0.0):
        """Unit flux maps of a MachineModel on a grid of spacing dx [m] covering the wall polygon."""
        from startup_model import eval_fields
        wall = np.asarray(wall, float)[:, :2]
        R = np.arange(wall[:, 0].min() - margin + dx / 2, wall[:, 0].max() + margin, dx)
        Z = np.arange(wall[:, 1].min() - margin + dx / 2, wall[:, 1].max() + margin, dx)
        RR, ZZ = np.meshgrid(R, Z)
        psi = eval_fields(machine.gs, machine.psi_u, np.column_stack([RR.ravel(), ZZ.ravel()]))['psi']
        return cls(R, Z, wall, machine.F0, psi.reshape(len(psi), len(Z), len(R)))

    def fields(self, I_u, dI_u, Ip=0.0, dIp=0.0, Rp=None, Zp=0.0, dRp=0.0, dZp=0.0):
        """(psi, E_phi) grids [Wb/rad, V/m] of unit currents I_u [nu] (per turn) changing at dI_u [A/s],
        plus a plasma filament Ip [A] at (Rp, Zp) moving at (dRp, dZp) [m/s] (K22 eq. 9-10)."""
        psi = np.tensordot(I_u, self.psi_units, 1)
        dpsi = np.tensordot(dI_u, self.psi_units, 1)
        if Rp is not None and (Ip != 0.0 or dIp != 0.0):
            g = filament_psi(self.RR, self.ZZ, Rp, Zp)
            psi = psi + Ip * g
            h = 1e-4
            dg = ((filament_psi(self.RR, self.ZZ, Rp + h, Zp) - filament_psi(self.RR, self.ZZ, Rp - h, Zp)) / (2 * h) * dRp
                  + (filament_psi(self.RR, self.ZZ, Rp, Zp + h) - filament_psi(self.RR, self.ZZ, Rp, Zp - h)) / (2 * h) * dZp)
            dpsi = dpsi + dIp * g + Ip * dg
        return psi, -dpsi / self.RR

    def mutual(self, R, Z):
        """M_p [nu]: mutual inductance of a filament at (R, Z) with 1 A per turn in each unit,
        2 pi psi_u(R, Z) (K22 eq. 9), bilinear in the unit maps (inside the grid only)."""
        if not (self.R[0] <= R <= self.R[-1] and self.Z[0] <= Z <= self.Z[-1]):
            raise ValueError(f'({R}, {Z}) is outside the map grid')
        i = int(np.clip(np.searchsorted(self.R, R) - 1, 0, len(self.R) - 2))
        j = int(np.clip(np.searchsorted(self.Z, Z) - 1, 0, len(self.Z) - 2))
        u, v = (R - self.R[i]) / (self.R[i + 1] - self.R[i]), (Z - self.Z[j]) / (self.Z[j + 1] - self.Z[j])
        P = self.psi_units
        return 2 * np.pi * ((1 - u) * (1 - v) * P[:, j, i] + u * (1 - v) * P[:, j, i + 1]
                            + (1 - u) * v * P[:, j + 1, i] + u * v * P[:, j + 1, i + 1])

    def bfield(self, psi):
        """(B_R, B_Z, B_phi) grids [T] of a flux grid (second-order differences)."""
        dpsi_dZ, dpsi_dR = np.gradient(psi, self.dZ, self.dR)
        return -dpsi_dZ / self.RR, dpsi_dR / self.RR, self.F0 / self.RR

    def _contours(self, psi, Ephi, B, BR, BZ, Epar):
        """Contour machinery of one (psi, E_phi): L and <E_par> along psi contours (polylines, linear in
        each grid cell): L = int |B|/B_p dl_pol, <E_par> = int E_par |B|/B_p dl_pol / L (K22 eq. 12, equal
        steps along B). A contour that closes on itself is a closed line (L = inf); one that ends at the
        masked wall is open and is extended along its end segments to the wall polygon (the masked
        contour stops within a cell of it)."""
        from contourpy import contour_generator
        from scipy.interpolate import RegularGridInterpolator
        Bp = np.sqrt(BR ** 2 + BZ ** 2)
        f_w = RegularGridInterpolator((self.Z, self.R), B / np.maximum(Bp, 1e-12), bounds_error=False, fill_value=None)
        f_e = RegularGridInterpolator((self.Z, self.R), Epar, bounds_error=False, fill_value=None)
        f_psi = RegularGridInterpolator((self.Z, self.R), psi, bounds_error=False, fill_value=None)
        gen = contour_generator(self.R, self.Z, np.ma.masked_array(psi, ~self.inside))

        def stats(line):
            closed = len(line) > 3 and np.hypot(*(line[0] - line[-1])) < 1e-9
            if not closed:
                ends = []
                for e0, e1 in ((line[0], line[1]), (line[-1], line[-2])):
                    d = e0 - e1
                    d = d / max(np.hypot(*d), 1e-300)
                    ends.append(e0 + d * _ray_to_polygon(e0, d, self.wall))
                line = np.vstack([ends[0], line, ends[1]])
            mp, ds = 0.5 * (line[1:] + line[:-1]), np.hypot(*np.diff(line, axis=0).T)
            w = ds * f_w(mp[:, ::-1])
            Lr = w.sum()
            return (np.inf if closed else Lr), abs(np.sum(w * f_e(mp[:, ::-1]))) / max(Lr, 1e-300), closed

        def through(pts):
            pts = np.atleast_2d(pts)
            L, E, cl = np.full(len(pts), np.nan), np.zeros(len(pts)), np.zeros(len(pts), bool)
            for i, x in enumerate(pts):
                lines = [l for l in gen.lines(float(f_psi((x[1], x[0])))) if len(l) > 1]
                d = [np.min(np.hypot(*(l - x).T)) for l in lines]
                if lines and min(d) < 1.5 * max(self.dR, self.dZ):         # the contour through x, not another
                    L[i], E[i], cl[i] = stats(lines[int(np.argmin(d))])
            return L, E, cl
        return {'gen': gen, 'stats': stats, 'through': through}

    def lines_through(self, psi, Ephi, pts):
        """(L [m] (inf if closed, nan if no contour), |<E_par>| [V/m], closed) of the field lines through
        the points pts [n, 2], by the contour method."""
        BR, BZ, Bphi = self.bfield(psi)
        B = np.sqrt(BR ** 2 + BZ ** 2 + Bphi ** 2)
        return self._contours(psi, Ephi, B, BR, BZ, Bphi / B * Ephi)['through'](pts)

    def analyze(self, psi, Ephi, p, n_bands=400, n_mid=50, Z_mid=0.0, region='midplane', tau_form='loss_fraction',
                plasma_exists=False, sustained=None):
        """Field-line analysis of (psi, E_phi) at prefill pressure p [Pa]. Returns a dict:
            tube [nZ, nR]   flux-tube label per cell (-1 outside); tubes: per tube 'open', 'L' [m]
                            (inf when closed), 'Epar' [V/m] (|<E_par>|), 'ET' [V/m] (K22 eq. 13),
                            'bd' (open and <E_par> >= E_T), 'V' [m^3], 'A' [m^2]
            mid             the n_mid mid-plane lines (K22: 50 radial points at Z_mid): 'R', 'tube',
                            'open', 'L', 'Epar', 'ET', 'bd'
            breakdown       any mid-plane line breaks down (region='midplane') / any tube (region='all')
            plasma [nZ, nR] the plasma region: tubes that broke down (those through the mid-plane points,
                            or every tube with region='all'), tubes through the mid-plane points in
                            `sustained` [n_mid] (lines that broke down at an earlier map, the 'sustained'
                            option of run_with_field_lines; a held line that has closed counts as a
                            closed line), plus every closed tube (when something broke down here, or
                            plasma_exists: a plasma from an earlier breakdown); and its
                            Vp [m^3], A [m^2], R_vol, Z_vol (centre of mass of the cross-section),
                            r_vol (distance from it to the first wall), kappa = min(A/(pi r^2),
                            (Z_max - Z_min)/(2 r)) (K22 sec. 2.2), f_closed (volume fraction closed)
            L_open          mean length of the open mid-plane lines in the plasma [m] (K22 eq. 17)
            L_sol           mean length of the open mid-plane lines next to the plasma region, on either
                            side (the scrape-off layer; nan if none), for PlasmaModel(neutrals='sol')
            f_open          their fraction of the mid-plane lines in the plasma
            Lf              effective parallel length for tau_par = Lf / C_s: 'loss_fraction' L_open /
                            f_open (only the open fraction loses particles; closed lines have
                            L = inf, K22 sec. 2.2); 'k22_printed' L_open f_open (eq. 17 as printed;
                            O28); inf when nothing is open
        """
        BR, BZ, Bphi = self.bfield(psi)
        B = np.sqrt(BR ** 2 + BZ ** 2 + Bphi ** 2)
        Epar = Bphi / B * Ephi
        ins = self.inside
        if not np.isfinite(psi[ins]).all():
            raise ValueError('psi is not finite inside the wall')
        closed = _closed_cells(psi, ins, self.edge)                  # closed lines, by topology
        opn = ins & ~closed
        lo, hi = psi[ins].min(), psi[ins].max()
        gpsi = self.RR * np.sqrt(BR ** 2 + BZ ** 2) * max(self.dR, self.dZ)         # psi change across a cell
        dpsi = max((hi - lo) / n_bands, 2.0 * np.percentile(gpsi[opn if opn.any() else ins], 90), 1e-30)
        n_bands = max(int(np.ceil((hi - lo) / dpsi)), 1)              # open tubes >= ~2 cells wide (90 %)
        band = np.where(ins, np.minimum(((psi - lo) / dpsi).astype(int), n_bands - 1), -1)
        key = np.where(ins, 2 * band + closed, -1)                   # open and closed cells in separate tubes
        # connected pieces of each band: 4-neighbour links between inside cells of the same key
        import scipy.sparse as sps
        from scipy.sparse.csgraph import connected_components
        idx = np.arange(key.size).reshape(key.shape)
        same_h = (key[:, :-1] == key[:, 1:]) & (key[:, :-1] >= 0)
        same_v = (key[:-1, :] == key[1:, :]) & (key[:-1, :] >= 0)
        i0 = np.concatenate([idx[:, :-1][same_h], idx[:-1, :][same_v]])
        i1 = np.concatenate([idx[:, 1:][same_h], idx[1:, :][same_v]])
        G = sps.coo_matrix((np.ones(len(i0)), (i0, i1)), shape=(key.size, key.size))
        _, lab = connected_components(G, directed=False)
        cells = np.where(ins.ravel())[0]
        uniq, tube_c = np.unique(lab.ravel()[cells], return_inverse=True)
        tube = np.full(key.shape, -1)
        tube.ravel()[cells] = tube_c
        nt = len(uniq)
        V = np.bincount(tube_c, 2 * np.pi * self.RR.ravel()[cells] * self.dA, nt)
        A = np.bincount(tube_c, np.full(len(cells), self.dA), nt)
        is_open = np.bincount(tube_c, closed.ravel()[cells], nt) == 0
        ctr = self._contours(psi, Ephi, B, BR, BZ, Epar)
        stats = ctr['stats']
        cell_of = lambda x: (np.clip(np.rint((x[..., 1] - self.Z[0]) / self.dZ).astype(int), 0, len(self.Z) - 1),
                             np.clip(np.rint((x[..., 0] - self.R[0]) / self.dR).astype(int), 0, len(self.R) - 1))
        L, Ep = np.where(is_open, np.nan, np.inf), np.zeros(nt)      # open tubes no contour reaches: L unknown
        for jb in range(n_bands):                                   # one contour per band, at its centre
            for line in ctr['gen'].lines(lo + (jb + 0.5) * dpsi):
                if len(line) < 2:
                    continue
                kk, jj = cell_of(line)
                ok = (band[kk, jj] == jb) & opn[kk, jj]
                if not ok.any():
                    continue
                Lt, Et, cl = stats(line)
                if cl:
                    continue
                for t in np.unique(tube[kk[ok], jj[ok]]):           # every open piece of the band it crosses
                    if np.isnan(L[t]) or Lt > L[t]:
                        L[t], Ep[t] = Lt, Et
        ET = townsend_field(p, L)
        bd = is_open & np.isfinite(L) & (Ep >= ET)
        # mid-plane lines (K22: 50 radial points on the mid-plane): the contour through each point
        if n_mid < 2:
            raise ValueError('n_mid must be at least 2')
        k = np.argmin(np.abs(self.Z - Z_mid))
        row = np.where(ins[k])[0]
        if len(row) == 0:
            raise ValueError(f'no grid cell inside the wall at Z_mid = {Z_mid}')
        Rm = np.linspace(self.R[row[0]] + self.dR / 2, self.R[row[-1]] - self.dR / 2, n_mid)
        pts = np.column_stack([Rm, np.full(n_mid, self.Z[k])])
        km, jm = cell_of(pts)
        Lm, Em, _ = ctr['through'](pts)
        cm = closed[km, jm]                                         # open / closed by topology
        Lm = np.where(cm, np.inf, Lm)
        Lm[~ins[km, jm]] = np.nan                                   # outside the wall (a gap at the mid-plane)
        mid = {'R': Rm, 'tube': tube[km, jm].copy(), 'open': ~cm, 'L': Lm, 'Epar': Em}
        mid['tube'][np.isnan(Lm)] = -1                               # outside, or no contour found
        mid['ET'] = townsend_field(p, Lm)
        mid['bd'] = mid['open'] & np.isfinite(Lm) & (Em >= mid['ET'])
        mid['held'] = mid['bd'].copy()                               # in the plasma: broke down now or before (O29)
        if sustained is not None:
            mid['held'] |= np.asarray(sustained, bool) & mid['open'] & np.isfinite(Lm)
        tm = mid['tube']
        # each held mid-plane line stands for the open flux that crosses the mid-plane nearer to it than to
        # its neighbours (K22's 50 lines sample the region; a tube is narrower than the point spacing on a
        # fine grid, review cycle 2): on each side, psi half-way to the neighbour's, or psi at the wall for
        # the end lines and next to a gap in the mid-plane row; where psi along the mid-plane has an
        # extremum between two points, up to the extremum. Its cells: those whose centre psi lies in that
        # range, connected to the line through cells whose psi span (centre +- half the change across
        # them) overlaps it, 8-neighbour, so a strip narrower than a cell stays connected along its line
        # without taking a cell column from its neighbours (review cycles 3, 4).
        from scipy import ndimage
        pm = np.interp(Rm, self.R, psi[k])                          # psi at the points
        Rx = np.concatenate([[self.R[row[0]]], Rm, [self.R[row[-1]]]])
        px = np.concatenate([[psi[k, row[0]]], pm, [psi[k, row[-1]]]])
        lo_s, hi_s = np.minimum(px[:-1], px[1:]), np.maximum(px[:-1], px[1:])       # the n_mid + 1 gaps
        e_up = (px[:-1] + px[1:]) / 2                                # bound for the point below each gap
        e_dn = e_up.copy()                                           # bound for the point above it
        e_dn[0], e_up[-1] = px[0], px[-1]                            # the end lines reach the wall
        ext = []
        for g in range(n_mid + 1):
            jj = np.where((self.R > Rx[g]) & (self.R < Rx[g + 1]))[0]
            inn = ins[k, jj]
            if g not in (0, n_mid) and len(jj) and not inn.all():    # the wall crosses the mid-plane here
                e_up[g] = psi[k, jj[np.argmin(inn)] - 1] if inn[0] else px[g]
                e_dn[g] = psi[k, jj[len(inn) - 1 - np.argmin(inn[::-1])] + 1] if inn[-1] else px[g + 1]
            seg = psi[k, jj[inn]]                                    # psi extremum overshooting the gap ends
            ext.append((min(seg.min(), lo_s[g]), max(seg.max(), hi_s[g])) if len(seg) else (lo_s[g], hi_s[g]))
        span = 0.5 * gpsi
        strips = np.zeros_like(opn)
        for i in np.where(mid['held'] & (tm >= 0))[0]:
            r = [pm[i], psi[km[i], jm[i]]]
            for g, e in ((i, e_dn[i]), (i + 1, e_up[i + 1])):        # the gaps below and above point i
                over = ext[g][0] < lo_s[g] or ext[g][1] > hi_s[g]
                r += [x for x in ext[g] if (x < lo_s[g] or x > hi_s[g])] if over else [e]
            a, b = min(r), max(r)
            lab, _ = ndimage.label(opn & (psi + span >= a) & (psi - span <= b), structure=np.ones((3, 3)))
            if lab[km[i], jm[i]]:
                strips |= (lab == lab[km[i], jm[i]]) & (psi >= a) & (psi <= b)
        take = np.zeros(nt, bool)
        if region == 'midplane':
            breakdown = bool(mid['bd'].any())
        else:
            take |= bd
            breakdown = bool(bd.any())
        if breakdown or plasma_exists:
            take |= ~is_open                                         # closed lines count once there is plasma
        plasma = np.isin(tube, np.where(take)[0]) | strips
        out = {'tube': tube, 'tubes': {'open': is_open, 'L': L, 'Epar': Ep, 'ET': ET, 'bd': bd, 'V': V, 'A': A},
               'mid': mid, 'breakdown': breakdown, 'plasma': plasma, 'dpsi': dpsi}
        out.update(self.shape(plasma))
        w = np.where(plasma, self.RR, 0.0)
        out['f_closed'] = float(w[closed].sum() / max(w.sum(), 1e-300))
        in_pl = plasma[km, jm] & (tm >= 0)                           # mid-plane lines inside the plasma
        if region == 'midplane':
            in_pl &= mid['held'] | ~mid['open']
        op = in_pl & mid['open'] & np.isfinite(mid['L'])
        out['f_open'] = float(op.sum() / max(in_pl.sum(), 1))
        ip = np.where(plasma[km, jm] & (tm >= 0))[0]                 # SOL: open mid-plane lines next to the region
        Ls = [mid['L'][j] for j in ((ip.min() - 1, ip.max() + 1) if ip.size else ())
              if 0 <= j < n_mid and mid['open'][j] and np.isfinite(mid['L'][j])]
        out['L_sol'] = float(np.mean(Ls)) if Ls else np.nan
        out['L_open'] = float(mid['L'][op].mean()) if op.any() else np.inf
        if not op.any():
            out['Lf'] = np.inf
        elif tau_form == 'loss_fraction':
            out['Lf'] = out['L_open'] / out['f_open']
        elif tau_form == 'k22_printed':
            out['Lf'] = out['L_open'] * out['f_open']
        else:
            raise ValueError(f'unknown tau_form {tau_form!r}')
        return out

    def shape(self, mask):
        """K22 sec. 2.2 shape of a region [nZ, nR]: Vp, A, R_vol, Z_vol, r_vol, kappa (zeros if empty)."""
        if not mask.any():
            return {'Vp': 0.0, 'A': 0.0, 'R_vol': np.nan, 'Z_vol': np.nan, 'r_vol': np.nan, 'kappa': np.nan}
        A = mask.sum() * self.dA
        Rv, Zv = self.RR[mask].mean(), self.ZZ[mask].mean()
        Vp = float(2 * np.pi * self.RR[mask].sum() * self.dA)
        r = _dist_to_polygon(np.array([Rv, Zv]), self.wall)
        kappa = min(A / (np.pi * r ** 2), (self.ZZ[mask].max() - self.ZZ[mask].min() + self.dZ) / (2 * r))
        return {'Vp': Vp, 'A': float(A), 'R_vol': float(Rv), 'Z_vol': float(Zv), 'r_vol': float(r), 'kappa': float(kappa)}

    def trace(self, psi, Ephi, x0, dl=0.05, scheme='euler', L_max=1e5):
        """Follow the line through x0 = (R, Z) both ways with steps dl [m] along B (K22 eq. 11 'euler',
        or 'rk4'), fields interpolated bilinearly on the grid. Returns (L [m] (inf if it closes on
        itself), |<E_par>| [V/m] (mean over the step points, K22 eq. 12), closed)."""
        from scipy.interpolate import RegularGridInterpolator
        from matplotlib.path import Path
        BR, BZ, Bphi = self.bfield(psi)
        B = np.sqrt(BR ** 2 + BZ ** 2 + Bphi ** 2)
        fR = RegularGridInterpolator((self.Z, self.R), BR / B, bounds_error=False, fill_value=None)
        fZ = RegularGridInterpolator((self.Z, self.R), BZ / B, bounds_error=False, fill_value=None)
        fE = RegularGridInterpolator((self.Z, self.R), Bphi / B * Ephi, bounds_error=False, fill_value=None)
        wall = Path(self.wall)
        v = lambda x: np.array([fR((x[1], x[0])), fZ((x[1], x[0]))])
        x0 = np.asarray(x0, float)
        n_tot, E_sum, L_end = 1, float(fE((x0[1], x0[0]))), 0.0
        for sgn in (1.0, -1.0):
            x, s_pol = x0.copy(), 0.0
            for _ in range(int(L_max / dl)):
                x_old = x
                if scheme == 'euler':
                    step = dl * v(x)
                else:
                    k1 = v(x); k2 = v(x + 0.5 * sgn * dl * k1); k3 = v(x + 0.5 * sgn * dl * k2); k4 = v(x + sgn * dl * k3)
                    step = dl / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
                x = x + sgn * step
                s_pol += np.hypot(*step)
                if not wall.contains_point(x):
                    s = np.hypot(*step)                             # the part of the last step inside the wall
                    L_end += dl * min(_ray_to_polygon(x_old, sgn * step / max(s, 1e-300), self.wall) / max(s, 1e-300), 1.0)
                    break
                n_tot += 1
                E_sum += float(fE((x[1], x[0])))
                if s_pol > 4 * max(self.dR, np.hypot(*step)) and np.hypot(*(x - x0)) < max(4 * np.hypot(*step), 0.5 * self.dR):
                    return np.inf, abs(E_sum / n_tot), True
            else:
                return np.inf, abs(E_sum / n_tot), False
        return (n_tot - 1) * dl + L_end, abs(E_sum / n_tot), False


def _closed_cells(psi, inside, edge):
    """[nZ, nR] True where the psi contour through the cell is closed: its superlevel or its sublevel
    connected set (4-neighbour, cells inside the wall) does not reach a cell next to the wall
    (union-find with cells added in psi order, both ways; cells of equal psi are added as one group
    and judged after it is merged)."""
    nZ, nR = psi.shape
    cells = np.flatnonzero(inside.ravel())
    p, e = psi.ravel(), edge.ravel()
    out = np.zeros(psi.size, bool)
    for sgn in (1.0, -1.0):
        parent, touch, done = np.arange(psi.size), e.copy(), np.zeros(psi.size, bool)

        def find(n):
            r = n
            while parent[r] != r:
                r = parent[r]
            while parent[n] != r:                                  # path compression
                parent[n], n = r, parent[n]
            return r
        order = cells[np.argsort(-sgn * p[cells], kind='stable')]
        new = np.flatnonzero(np.r_[True, p[order][1:] != p[order][:-1]])
        for g0, g1 in zip(new, np.r_[new[1:], len(order)]):
            for c in order[g0:g1]:
                done[c] = True
                i, j = divmod(c, nR)
                for n in ((c - nR) if i > 0 else -1, (c + nR) if i < nZ - 1 else -1, (c - 1) if j > 0 else -1,
                          (c + 1) if j < nR - 1 else -1):
                    if n >= 0 and done[n]:
                        r, rc = find(n), find(c)
                        if r != rc:
                            parent[r] = rc
                            touch[rc] |= touch[r]
            for c in order[g0:g1]:
                out[c] |= not touch[find(c)]
    return out.reshape(psi.shape)


def _ray_to_polygon(x, d, poly):
    """Distance [m] from x along the unit direction d to the first crossing of the closed polygon poly
    [n, 2] (0 if none ahead)."""
    a, b = poly, np.roll(poly, -1, axis=0)
    e = b - a
    den = d[0] * e[:, 1] - d[1] * e[:, 0]
    ok = np.abs(den) > 1e-300
    w = a - x
    t = np.where(ok, (w[:, 0] * e[:, 1] - w[:, 1] * e[:, 0]) / np.where(ok, den, 1), -1)      # along the ray
    u = np.where(ok, (w[:, 0] * d[1] - w[:, 1] * d[0]) / np.where(ok, den, 1), -1)            # along the edge
    hit = ok & (t >= 0) & (u >= 0) & (u <= 1)
    return float(t[hit].min()) if hit.any() else 0.0


def _dist_to_polygon(x, poly):
    """Shortest distance [m] from point x to the closed polygon poly [n, 2]."""
    a, b = poly, np.roll(poly, -1, axis=0)
    ab = b - a
    t = np.clip(np.einsum('ij,ij->i', x - a, ab) / np.maximum(np.einsum('ij,ij->i', ab, ab), 1e-300), 0, 1)
    return float(np.min(np.hypot(*(a + t[:, None] * ab - x).T)))


# ============================================================================ the DYON workflow (K22)
def geometry_between(ga, gb, ta, dt):
    """geometry(t, Ip) for PlasmaModel.run from map geometry ga to gb over [ta, ta + dt]: linear in each
    quantity, harmonic in Lf (1/Lf, the loss rate, is linear); held at gb after ta + dt."""
    def geometry(tt, Ip):
        s = min(max((tt - ta) / dt, 0.0), 1.0)
        g = {k: ga[k] + s * (gb[k] - ga[k]) for k in ga}
        g['Lf'] = 1.0 / (1.0 / ga['Lf'] + s * (1.0 / gb['Lf'] - 1.0 / ga['Lf']))
        return g
    return geometry


def run_with_field_lines(model, fmap, p_fill, t_start, t_end, y0=None, dt_map=1e-3, Rp=None, Zp=0.0, Te0=1.0,
                         gamma0=None, impurity_fraction=None, region='midplane', tau_form='loss_fraction',
                         n_mid=50, Z_mid=0.0, hold='k22', pressure='neutral', position='fixed', tau_pos=5e-3, log=None):
    """DYON's electromagnetic plasma initiation (K22): a PlasmaModel with a PlasmaCircuit, its geometry
    from the field-line analysis of FieldLineMap every dt_map (K22 table 1: 1 ms).

    Before breakdown the structure currents follow the coils (vacuum). From t_start the map is
    analysed every dt_map at the effective prefill pressure p_fill [Pa]; at the first time a line meets
    the Townsend criterion (K22 sec. 2.2; O6) the plasma model starts from the K22 initial state (fill,
    Te0, gamma0, impurity_fraction) unless y0 is given. Thereafter each map uses the coils, the
    structure and the plasma current, with its own -M dI_p/dt (K22 eq. 10).

    Plasma current position (O27): position='fixed' (default) a filament held at (Rp, Zp), which must be
    the circuit's filament position (K22 prescribes the position); 'region' (a departure: K22 prescribes
    R, Z, the magnetic axis once surfaces close) a filament that relaxes toward the area centre of the
    plasma region (R_vol, Z_vol) with time constant tau_pos [s] (starting there at breakdown), moving
    linearly between maps with the geometry's lag (the circuit's M_p from FieldLineMap.mutual, with its
    d M_p/dt terms: PlasmaCircuit.move; Rp, Zp are then not used, and the circuit's own M_p is restored on
    return, also on an exception). The relaxation keeps the filament's speed independent of dt_map: the
    region can flicker between maps (lines on and off the Townsend criterion), and following it map by
    map gave a motion voltage ~1/dt_map (review cycle 4). The returned 'out' holds the terms at the
    moving position; model.terms after the run uses the restored filament.
    No force balance: the mechanical work I_p dM_p/dt . I_u is not in the plasma energy.

    Geometry handed to the plasma model: R = R_vol, a = r_vol, kappa, Vp, Bphi = F0/R_vol, Lf (field-line
    parallel length, capped at 1e7 m). Between maps it is interpolated linearly from the map before
    the last to the last (1/Lf, the loss rate, for Lf; geometry_between), a lag of one interval, so it is
    continuous for the BDF solver. The output is at the map times and half-way between them.

    Plasma region after breakdown (O29): hold='k22' (default, K22 fig. 4 caption: particles are on the
    lines that are closed or meet the Townsend criterion, re-evaluated at every map; ITPA 2026 sec. 5:
    as the neutrals deplete fewer open lines meet it and the closed volume replaces them). 'sustained'
    (a departure, not in K22): a mid-plane line that met the criterion at any map since breakdown stays
    in the plasma while it is open, with no release. If the region is empty the last geometry is held.
    Townsend pressure after breakdown: pressure='neutral' (default) the fuel neutrals of the model,
    n_0 T_n / 2 (molecules, the inverse of K22 eq. 18): depletion turns short lines off, but it lowers
    E_T on lines long enough that A p L > e (left of the Paschen minimum); 'fill' p_fill.

    Convergence in dt_map (review cycle 3): with the geometry one map behind, results converge slowly in
    dt_map (toy machine, I_p at 2 ms: 1922, 2430, 2693, 2855, 3098 A at 2, 1, 0.5, 0.25, 0.125 ms; at
    20 ms ~92 A for dt_map <= 0.5 ms, 102 at 1 ms). K22's 1 ms is not converged (O37).

    A given y0 gets the structure currents of the vacuum phase at t_bd.

    Returns a dict: 't_bd' (None if no breakdown before t_end), 'out' (the plasma model's output arrays,
    concatenated; None without breakdown), 'status' (0, or the solver's status where it stopped),
    'message', 't_stop' (the last time reached; status -2: the model raised ValueError, eg. L_p below
    the structure coupling M_ps M_ss^-1 M_sp when the region's R_vol is far from Rp, O27), 'maps' (per map time: t, p (the Townsend pressure), Rp, Zp (the filament), Vp, R_vol, Z_vol, r_vol,
    kappa, Lf, L_open, f_open, f_closed, n_bd (mid-plane lines broken down at that map), n_held (lines
    in the plasma), breakdown).
    """
    if Rp is None and position == 'fixed':
        raise ValueError('Rp (the plasma current position, the circuit filament radius) must be given')
    if hold not in ('k22', 'sustained') or pressure not in ('neutral', 'fill') or position not in ('fixed', 'region'):
        raise ValueError(f'unknown hold {hold!r}, pressure {pressure!r} or position {position!r}')
    circ = model.circuit
    M_saved = circ._M
    nc = circ.nc
    keys = ('Vp', 'R_vol', 'Z_vol', 'r_vol', 'kappa', 'Lf', 'L_open', 'f_open', 'f_closed', 'L_sol')
    maps = {k: [] for k in ('t', 'p', 'n_bd', 'n_held', 'breakdown', 'Rp', 'Zp') + keys}
    held = np.zeros(n_mid, bool)
    I_coil = lambda t: np.array([np.interp(t, circ.t_coil, circ.I_coil[:, j]) for j in range(nc)])
    log = log or (lambda *a: None)

    def analyze(t, I_s, dI_s, Ip, dIp, plasma_exists=False, p=p_fill, pos=(Rp, Zp), vel=(0.0, 0.0)):
        nonlocal held
        I_u = np.concatenate([I_coil(t), I_s])
        dI_u = np.concatenate([circ.coil_rate(t), dI_s])
        psi, Ephi = fmap.fields(I_u, dI_u, Ip, dIp, pos[0], pos[1], vel[0], vel[1])
        sus = held if plasma_exists and hold == 'sustained' else None
        o = fmap.analyze(psi, Ephi, p, n_mid=n_mid, Z_mid=Z_mid, region=region, tau_form=tau_form,
                         plasma_exists=plasma_exists, sustained=sus)
        if hold == 'sustained' and (plasma_exists or o['breakdown']):
            held = held | o['mid']['held']
        maps['t'].append(t)
        maps['p'].append(p)
        maps['Rp'].append(pos[0] if pos[0] is not None else np.nan)
        maps['Zp'].append(pos[1])
        maps['n_bd'].append(int(o['mid']['bd'].sum()))
        maps['n_held'].append(int(o['mid']['held'].sum()))
        maps['breakdown'].append(o['breakdown'])
        for k in keys:
            maps[k].append(o[k] if k in o else np.nan)
        return o

    def geo_of(o):
        return {'R': o['R_vol'], 'a': o['r_vol'], 'kappa': o['kappa'], 'Vp': o['Vp'], 'Bphi': abs(fmap.F0) / o['R_vol'],
                'Lf': min(o['Lf'], 1e7), 'L_sol': o['L_sol']}

    # vacuum until the first Townsend breakdown
    f = lambda t, I: circ.K_ss @ (-circ.M_sc @ circ.coil_rate(t) - circ.R_s * I)
    I_s = circ.vacuum(circ.t_coil[0], t_start)
    t, o, t_bd = t_start, None, None
    while t <= t_end + 1e-12:
        o = analyze(t, I_s, f(t, I_s), 0.0, 0.0)
        if o['breakdown'] and o['Vp'] > 0:
            t_bd = t
            break
        I_s = solve_ivp(f, (t, t + dt_map), I_s, method='BDF', rtol=1e-8, atol=1e-6).y[:, -1]
        t += dt_map
    if t_bd is None:
        log(f'no Townsend breakdown in {t_start * 1e3:.0f}..{t_end * 1e3:.0f} ms at p = {p_fill * 1e3:.3f} mPa')
        return {'t_bd': None, 'out': None, 'status': 0, 'message': 'no breakdown', 't_stop': t_end,
                'maps': {k: np.array(v) for k, v in maps.items()}}
    g_prev = g_last = geo_of(o)
    P_prev = P_last = (Rp, Zp) if position == 'fixed' else (o['R_vol'], o['Z_vol'])
    if y0 is None:
        y0 = model.initial_state(p_fill, g_last, Te0=Te0, gamma0=gamma0, impurity_fraction=impurity_fraction)
    y0 = np.array(y0, float)
    y0[model.iIs:] = I_s                                             # the structure currents of the vacuum phase
    log(f'breakdown at {t_bd * 1e3:.1f} ms: {maps["n_bd"][-1]} mid-plane lines, V_p = {o["Vp"]:.2f} m^3, '
        f'R_vol = {o["R_vol"]:.2f} m, L_f = {o["Lf"]:.0f} m')
    outs, y, t, status, message = [], y0, t_bd, 0, 'reached t_end'
    try:
        while t < t_end - 1e-12:
            t1 = min(t + dt_map, t_end)
            geometry = geometry_between(dict(g_prev), dict(g_last), t, dt_map)
            if position == 'region':
                circ.move(fmap.mutual(*P_prev), fmap.mutual(*P_last), t, t1 - t)
            try:
                o_run = model.run(geometry, None, y, (t, t1), t_eval=[t, 0.5 * (t + t1), t1])
            except ValueError as e:                                      # eg. L_p below the structure coupling (O27)
                o_run = {'status': -2, 'message': str(e)}
            if o_run['status'] != 0:
                status, message = o_run['status'], o_run['message']
                log(f'solver stopped at {t * 1e3:.1f} ms: {message}')
                break
            outs.append({k: v[1:] if len(outs) else v for k, v in o_run.items() if isinstance(v, np.ndarray)})
            y, t0, t = o_run['y'][-1], t, t1
            tt = model.terms(t, y)
            p_now = tt['n0'][0] * model.T_n * EV / 2 if pressure == 'neutral' else p_fill
            vel = ((P_last[0] - P_prev[0]) / (t1 - t0), (P_last[1] - P_prev[1]) / (t1 - t0)) if position == 'region' else (0.0, 0.0)
            o = analyze(t, y[model.iIs:], tt['dy'][model.iIs:], y[model.iIp], tt['dy'][model.iIp], plasma_exists=True,
                        p=p_now, pos=P_last, vel=vel)
            g_prev, P_prev = g_last, P_last
            if o['Vp'] > 0:
                g_last = geo_of(o)
                if position == 'region':
                    w = min(dt_map / tau_pos, 1.0)
                    P_last = (P_last[0] + w * (o['R_vol'] - P_last[0]), P_last[1] + w * (o['Z_vol'] - P_last[1]))
    finally:
        circ.move(M_saved[0], M_saved[1], M_saved[2], M_saved[3])
    out = {k: np.concatenate([d[k] for d in outs]) for k in outs[0]} if outs else None
    return {'t_bd': t_bd, 'out': out, 'status': status, 'message': message, 't_stop': t,
            'maps': {k: np.array(v) for k, v in maps.items()}}
