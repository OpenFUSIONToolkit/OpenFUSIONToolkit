"""Convert TRANSP CDF/MDS output into TokaMaker_TORAX-compatible input dicts.

This reads a TRANSP output file (the CDF archive of the MDSplus tree — same
data, same variable names) and produces the same {time: value} /
{time: {rho: value}} dict formats that kf_pulse_ex_tokamaker_torax.ipynb feeds
into TokaMaker_TORAX's tmtx.set_* methods, so real DIII-D TRANSP data can
replace that notebook's illustrative arrays.

Coordinate note
----------------
TRANSP's X (zone-centered) and XB (zone-boundary) grids are, by construction,
sqrt(toroidal flux / toroidal flux at the LCFS) — i.e. already TORAX's
rho_toroidal coordinate. This was verified against TRFLX for this file to
machine precision (X and XB are also time-independent). So profiles only need
to be *resampled* onto TORAX's rho grid, not coordinate-transformed.

Scope
-----
Kinetic profiles, composition, heating, fueling, and Ip only — geometry/seed
eqdsk construction is intentionally out of scope (for a real device you'd get
the equilibrium from EFIT/MDS, not from TRANSP's internal moment expansion).

Known limitations (see each method's docstring for detail):
- TORAX's "generic" heat/particle sources are Gaussian-only (P_total,
  location, width). Real TRANSP deposition profiles are not Gaussian, so
  get_nbi_heating/get_particle_source *fit* a Gaussian to the actual profile
  (weighted centroid/std) rather than reproducing its shape exactly.
- Edge/boundary values come from TRANSP's outermost analyzed flux surface
  (rho ~ 0.99); TRANSP does not resolve the SOL, so these are not true
  measured separatrix values.
- Impurity species is inferred from TRANSP's average-Z/A profile by rounding
  to the nearest element TORAX's Mavrin model supports. TRANSP only tracks an
  effective average impurity, not real species mixtures.
"""

import numpy as np
import scipy.io as sio

# Elements TORAX's Mavrin impurity-radiation model supports, by atomic number.
_MAVRIN_Z_TO_SYMBOL = {
    2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
    10: 'Ne', 18: 'Ar', 36: 'Kr', 54: 'Xe', 74: 'W',
}

# Candidate main-ion density variables TRANSP may carry, keyed by species symbol.
_MAIN_ION_CANDIDATES = {'H': 'NH', 'D': 'ND', 'T': 'NT'}


def _weighted_centroid_width(rho, profile):
    """Weighted centroid/std of `profile` over `rho`, for fitting a Gaussian
    deposition shape (location, width) to an arbitrary TRANSP profile."""
    w = np.clip(profile, 0.0, None)
    wsum = w.sum()
    if wsum <= 0:
        return float('nan'), float('nan')
    loc = float((w * rho).sum() / wsum)
    width = float(np.sqrt((w * (rho - loc) ** 2).sum() / wsum))
    return loc, width


def _linear_extrapolate_to_axis(x_sorted, y_sorted, x_axis=0.0):
    """Linearly extrapolate y(x) to x_axis using the two points nearest
    x_axis (x_sorted/y_sorted ordered by ascending *distance* from x_axis -
    see _fix_axis_point, which builds that ordering correctly regardless of
    whether x increases or decreases moving away from the axis), instead of
    holding the innermost native value flat out to the axis. A flat hold
    creates a zero-slope segment from the axis to the first native grid
    point in TokaMaker's piecewise-linear profile, which measurably
    under-drives on-axis current density: confirmed on a real solve, this
    was consistent with TokaMaker's q_0 coming out ~15% high vs. TRANSP's
    own q_0 even after the axis-dense native grid fix (see
    get_ffprime_pprime_profiles's "Output grid" note) - the remaining gap,
    not the grid density."""
    x0, x1 = x_sorted[0], x_sorted[1]
    y0, y1 = y_sorted[0], y_sorted[1]
    slope = (y1 - y0) / (x1 - x0)
    return y0 + slope * (x_axis - x0)


def _fix_axis_point(x, y, x_axis=0.0):
    """Replace y's point nearest x_axis with _linear_extrapolate_to_axis's
    2-point extrapolation, rather than leaving TRANSP's own native
    axis-adjacent sample as-is. Orders by *distance* from x_axis (not raw
    ascending x) so this works whether x increases or decreases moving away
    from the axis - e.g. psi here is often negative and decreases outward,
    where a plain ascending sort would pick the far edge, not the axis."""
    order = np.argsort(np.abs(x - x_axis))
    out = y.copy()
    out[order[0]] = _linear_extrapolate_to_axis(x[order], y[order], x_axis=x_axis)
    return out


def _monotonic_gradient(y, x):
    """np.gradient(y, x), safe against a non-monotonic `x`: sorts by `x`
    before differentiating and restores the original ordering afterward.
    A no-op safety net when `x` is already monotonic (increasing or
    decreasing) - see get_ffprime_pprime_profiles."""
    d = np.diff(x)
    if np.all(d > 0) or np.all(d < 0):
        return np.gradient(y, x)
    order = np.argsort(x)
    dy_sorted = np.gradient(y[order], x[order])
    unsort = np.empty_like(order)
    unsort[order] = np.arange(order.size)
    return dy_sorted[unsort]


def _fit_exponential_decay_length(rho, profile, decay_start=1.0):
    """Least-squares decay length for TORAX's gas_puff source model,
    profile ~ C * exp(-(decay_start - rho) / decay_length) (see torax's
    sources/gas_puff_source.py + formulas.exponential_profile) - fit by
    linear regression of log(profile) against (decay_start - rho), which is
    exact for this model (unlike the Gaussian-moment fit used elsewhere in
    this module, which assumes a different functional form entirely and
    fits gas-puff-like edge-peaked profiles poorly)."""
    valid = profile > 0
    if valid.sum() < 2:
        return float('nan')
    x = decay_start - rho[valid]
    y = np.log(profile[valid])
    slope, _ = np.polyfit(x, y, 1)
    if slope >= 0:
        return float('nan')
    return float(-1.0 / slope)


class TranspRun:
    """Wraps one TRANSP CDF file and extracts TokaMaker_TORAX-ready inputs."""

    def __init__(self, cdf_path):
        self._f = sio.netcdf_file(cdf_path, mmap=False)
        self.time = self._var('TIME')          # s
        self.rho = self._var('X')[0]            # zone-centered rho_tor, shape (nx,)
        self.rho_b = self._var('XB')[0]          # zone-boundary rho_tor, shape (nxb,)

        # X/XB are supposed to be time-independent rho_tor grids (see module
        # docstring). Fail loudly if a given file violates that assumption
        # instead of silently producing wrong profiles.
        x_all = self._var('X')
        xb_all = self._var('XB')
        if not np.allclose(x_all, x_all[0]) or not np.allclose(xb_all, xb_all[0]):
            raise ValueError(
                "X/XB vary with time in this file — the rho_tor shortcut this "
                "module relies on doesn't hold; profiles need a real per-time "
                "coordinate transform via TRFLX before resampling."
            )

    def has(self, name):
        return name in self._f.variables

    def _var(self, name):
        return self._f.variables[name].data.astype(float).copy()

    def _units(self, name):
        v = self._f.variables[name]
        return v.units.decode().strip() if hasattr(v, 'units') else None

    def _long_name(self, name):
        v = self._f.variables[name]
        return v.long_name.decode().strip() if hasattr(v, 'long_name') else None

    # -- resampling -----------------------------------------------------

    def profile(self, name, rho_grid, source_rho='X'):
        """`name`(time, native_rho) linearly resampled onto `rho_grid`.

        Beyond the native grid's range, np.interp holds the edge value flat —
        appropriate here since TRANSP's grid doesn't reach exactly rho=0 or 1.
        """
        raw = self._var(name)
        native_rho = self.rho if source_rho == 'X' else self.rho_b
        rho_grid = np.asarray(rho_grid, dtype=float)
        out = np.empty((raw.shape[0], rho_grid.size))
        for it in range(raw.shape[0]):
            out[it] = np.interp(rho_grid, native_rho, raw[it])
        return out

    def to_time_rho_dict(self, arr, rho_grid, time=None):
        time = self.time if time is None else time
        rho_grid = np.asarray(rho_grid, dtype=float)
        return {
            float(t): {float(r): float(v) for r, v in zip(rho_grid, row)}
            for t, row in zip(time, arr)
        }

    def to_time_dict(self, arr, time=None):
        time = self.time if time is None else time
        return {float(t): float(v) for t, v in zip(time, arr)}

    # -- kinetic profiles (for set_ne / set_Te / set_Ti) -----------------

    def get_ne(self, rho_grid):
        """Electron density in m^-3, {time: {rho: value}}."""
        arr = self.profile('NE', rho_grid) * 1e6   # N/cm^3 -> N/m^3
        return self.to_time_rho_dict(arr, rho_grid)

    def get_Te(self, rho_grid):
        """Electron temperature in keV, {time: {rho: value}}."""
        arr = self.profile('TE', rho_grid) / 1000.0  # eV -> keV
        return self.to_time_rho_dict(arr, rho_grid)

    def get_Ti(self, rho_grid):
        """Ion temperature in keV, {time: {rho: value}}."""
        arr = self.profile('TI', rho_grid) / 1000.0
        return self.to_time_rho_dict(arr, rho_grid)

    def get_fast_ion_pressure_axis(self):
        """On-axis fast-ion (NBI beam, nonthermal) pressure in Pa, {time: value}.

        For use as TokaMaker_TORAX.set_fast_ion_pressure()'s external correction
        added directly onto TORAX's own thermal-only pax - TORAX has no NBI
        fast-ion physics model in this version (numerics.enable_fast_ions /
        profile_conditions.fast_ions only takes effect for species an active
        source declares via zero_fast_ions(), and currently only
        ion_cyclotron_source/ICRH does that - see set_fast_ion_pressure's
        docstring), so this is a plain measured input, not something TORAX
        evolves or models.

        PMHDF_IN ("nonthermal pressure to MHD solver") is TRANSP's own on-axis
        beam-ion pressure, already isotropized for its own MHD solver -
        PMHD_IN = PMHDT_IN (thermal) + PMHDF_IN (this) + PMHDR_IN (rotation,
        not modeled by either code here), verified additive on this file to
        ~1e-7 relative precision.
        """
        arr = self._var('PMHDF_IN')[:, 0]
        return self.to_time_dict(arr)

    def get_Zeff_profile(self, rho_grid):
        """Zeff profile (dimensionless), {time: {rho: value}}."""
        arr = self.profile('ZEFFP', rho_grid)
        return self.to_time_rho_dict(arr, rho_grid)

    def get_right_bc(self, quantity, rho_edge=None):
        """Edge value time trace for use as set_ne/set_Te/set_Ti's `right_bc`.

        `quantity` is one of 'ne' (m^-3), 'Te'/'Ti' (keV). Default rho_edge is
        the outermost point of TRANSP's own grid (~0.99) since TRANSP doesn't
        resolve the SOL — this is an outer-core value, not a true separatrix
        measurement.
        """
        var, scale = {'ne': ('NE', 1e6), 'Te': ('TE', 1e-3), 'Ti': ('TI', 1e-3)}[quantity]
        rho_edge = self.rho[-1] if rho_edge is None else rho_edge
        arr = self.profile(var, [rho_edge]) * scale
        return self.to_time_dict(arr[:, 0])

    def get_transport_coefficients(self, rho_grid):
        """TRANSP's interpretive heat diffusivities, in m^2/s, {time: {rho: value}}.

        CONDE/CONDI are TRANSP's power-balance-derived *total* effective
        electron/ion heat diffusivities (whatever combination of neoclassical
        and turbulent transport is needed to match the measured profile
        evolution) — not a first-principles prediction. This is the natural
        quantity to compare against a predictive code's own total chi_e/chi_i
        (eg. TORAX's chi_neo_e + chi_turb_e), as a check on whether the
        predictive transport model reproduces the real discharge's effective
        transport level, not as a component-by-component (neoclassical vs.
        turbulent) comparison — TRANSP doesn't separate those the way a
        predictive turbulent-transport model does.
        """
        # CONDE/CONDI are defined on TRANSP's zone-boundary (XB) grid, not the
        # zone-centered X grid profile() defaults to.
        chi_e = self.profile('CONDE', rho_grid, source_rho='XB') * 1e-4  # cm^2/s -> m^2/s
        chi_i = self.profile('CONDI', rho_grid, source_rho='XB') * 1e-4
        return {
            'chi_e': self.to_time_rho_dict(chi_e, rho_grid),
            'chi_i': self.to_time_rho_dict(chi_i, rho_grid),
        }

    def get_pedestal_values(self, ped_top=0.9):
        """ne (m^-3), Te (keV), Ti (keV) sampled at rho=ped_top, for use as
        n_e_ped/T_e_ped/T_i_ped in set_pedestal. TRANSP profiles are transport
        solutions, not a fitted pedestal — this just reads the value there."""
        ne = self.profile('NE', [ped_top])[:, 0] * 1e6
        Te = self.profile('TE', [ped_top])[:, 0] / 1000.0
        Ti = self.profile('TI', [ped_top])[:, 0] / 1000.0
        return {
            'n_e_ped': self.to_time_dict(ne),
            'T_e_ped': self.to_time_dict(Te),
            'T_i_ped': self.to_time_dict(Ti),
        }

    # -- Ip ---------------------------------------------------------------

    def get_Ip(self):
        """Measured plasma current in A, {time: value}, for set_Ip."""
        return self.to_time_dict(self._var('PCUR'))

    # -- plasma composition (for set_plasma_composition) ------------------

    def get_plasma_composition(self):
        """Main-ion fractions, Zeff time trace, and an inferred impurity
        species — for set_plasma_composition(main_ion=..., Zeff=..., impurity=...).

        Main-ion fractions come from whichever of NH/ND/NT are present in this
        file (checked by variable existence, since composition is run-specific).
        Impurity species is TRANSP's average Z/A (XZIMP/AIMP) rounded to the
        nearest element TORAX's Mavrin model supports — TRANSP tracks one
        effective impurity, not a real species mixture, so this is a
        best-effort label, not a measurement of actual impurity content.
        """
        species_density = {}
        for symbol, varname in _MAIN_ION_CANDIDATES.items():
            if self.has(varname):
                species_density[symbol] = self._var(varname)  # (time, X), N/cm^3

        if not species_density:
            raise ValueError("No main-ion density variables (NH/ND/NT) found in this file.")

        total = sum(species_density.values())
        # Flux-surface-volume-averaged fraction per species, one number per
        # species for the whole discharge (TORAX's main_ion fractions are
        # static, not time/rho-dependent).
        main_ion_fractions = {
            symbol: float(np.mean(dens / total))
            for symbol, dens in species_density.items()
        }

        zeff_trace = self.to_time_dict(np.mean(self._var('ZEFFP'), axis=1))

        impurity = None
        if self.has('XZIMP') and self.has('AIMP'):
            z_avg = float(np.mean(self._var('XZIMP')))
            a_avg = float(np.mean(self._var('AIMP')))
            nearest_z = min(_MAVRIN_Z_TO_SYMBOL, key=lambda z: abs(z - z_avg))
            impurity = {
                'symbol': _MAVRIN_Z_TO_SYMBOL[nearest_z],
                'measured_avg_Z': z_avg,
                'measured_avg_A': a_avg,
            }

        return {
            'main_ion_fractions': main_ion_fractions,
            'Zeff_time_avg': float(np.mean(list(zeff_trace.values()))),
            'Zeff_trace': zeff_trace,
            'impurity': impurity,
        }

    # -- geometry / seed-equilibrium shape (for TokaMaker) ------------------
    #
    # NOTE on scope: this does NOT reconstruct a real EFIT-quality equilibrium
    # (that needs a real Grad-Shafranov solve against magnetics, which TRANSP
    # itself consumes rather than produces independently). What follows are
    # TRANSP's own internal flux-surface geometry quantities — same numbers
    # the notebook would otherwise pull from an eqdsk's LCFS boundary (rcentr,
    # bcentr, zaxis, a, kappa, delta upper/lower) — useful as seed-equilibrium
    # shape targets for TokaMaker, not as a substitute for a measured
    # equilibrium.

    _NUM_BOUNDARY_MOMENTS = 8  # RMC00..RMC08 / RMS01..RMS08 etc. in this file
    _MU0 = 4 * np.pi * 1e-7    # vacuum permeability, for get_reconstructed_Ip's Grad-Shafranov relation

    def get_geometry(self):
        """Scalar/LCFS geometry time traces for seeding TokaMaker equilibria.

        - F0 = B0*R0 [T*m] comes directly from TRANSP's BZXR (vacuum "Bz*R"),
          which is exactly the product TokaMaker's setup(F0=...) wants — no
          need to separately pick an R0/B0 convention for that call.
        - R0 is reported two ways: R_axis (TRANSP's magnetic axis major
          radius) and R0_geo (LCFS geometric center, (Rmax+Rmin)/2 from
          RMJMP) — these differ by a few cm; the notebook's eqdsk-based
          r0_flattop is the RMJMP-style geometric center, so prefer R0_geo
          for shape-target calculations and R_axis only if you specifically
          need the magnetic axis location.
        - a, kappa (LCFS elongation from ELONG), delta_upper/delta_lower
          (from TRIANGU/TRIANGL) are read directly from TRANSP rather than
          derived from a boundary contour. kappa_upper/kappa_lower are
          derived from the boundary moments (see get_boundary_contour) since
          TRANSP doesn't carry a single-sided elongation variable. Cross-
          checked at one time slice against a full moment-expansion contour
          reconstruction (get_boundary_contour): a/R0/kappa/delta agreed to
          within ~1-2% (residual from truncating the moment expansion at
          n=8 and finite theta resolution).
        - pax is on-axis total pressure (PMHD_IN, includes fast-ion/rotation
          pressure — TRANSP does not separately expose a thermal-only
          on-axis pressure), for eqdsk['pres'][0]'s role in mygs.set_targets.
        - psi_lcfs is TRANSP's own enclosed poloidal flux (PLFLXA), sign-
          flipped to match TokaMaker's/TORAX's COCOS convention (see
          "Sign convention" note below), for an eqdsk's psibry's role as a
          psi_lcfs constraint target.

        Sign convention: verified 2026-08-19 against a real coupled run —
        TRANSP's raw PLFLXA is positive and *increases* through an Ip
        ramp-up (PCUR and BZXR/F0 are both positive too, so this isn't an
        Ip-sign issue), while TokaMaker's/TORAX's own psi *decreases*
        (normal flux consumption during a current ramp; TORAX's own raw psi
        needs its own -1/(2*pi) flip+scale to match TokaMaker's convention,
        see pulse_design.py's plot_scalars). Negating PLFLXA here makes both
        conventions decrease in the same direction. This is still only a
        sign fix, not a gauge/offset fix — psi is meaningful only up to an
        additive constant, so absolute values need not (and generally will
        not) match TokaMaker/TORAX exactly; treat this as a soft constraint
        for that reason, not a precise, convention-matched flux value.
        """
        R_axis = self._var('RAXIS') / 100.0       # cm -> m
        Z_axis = self._var('YAXIS') / 100.0
        F0 = self._var('BZXR') / 100.0             # T*cm -> T*m
        R0_geo = self._var('RMJMP')[:, -1] / 100.0  # LCFS geometric center
        a = self._var('RMNMP')[:, -1] / 100.0       # LCFS minor radius
        kappa = self._var('ELONG')[:, -1]
        delta_upper = self._var('TRIANGU')[:, -1]
        delta_lower = self._var('TRIANGL')[:, -1]
        kappa_upper, kappa_lower = self._kappa_updown_from_moments()
        pax = self._var('PMHD_IN')[:, 0]
        psi_lcfs = -self._var('PLFLXA')  # Wb/rad, sign-flipped — see "Sign convention" above

        return {
            'time': self.time,
            'R_axis': self.to_time_dict(R_axis),
            'Z_axis': self.to_time_dict(Z_axis),
            'F0': self.to_time_dict(F0),
            'R0_geo': self.to_time_dict(R0_geo),
            'a': self.to_time_dict(a),
            'kappa': self.to_time_dict(kappa),
            'kappa_upper': self.to_time_dict(kappa_upper),
            'kappa_lower': self.to_time_dict(kappa_lower),
            'delta_upper': self.to_time_dict(delta_upper),
            'delta_lower': self.to_time_dict(delta_lower),
            'pax': self.to_time_dict(pax),
            'psi_lcfs': self.to_time_dict(psi_lcfs),
            'Ip': self.get_Ip(),
        }

    def _kappa_updown_from_moments(self):
        """kappa_upper/kappa_lower at each time, from the LCFS boundary
        moment expansion (TRANSP has no direct single-sided-elongation
        variable — only combined ELONG)."""
        n = self._NUM_BOUNDARY_MOMENTS
        theta = np.linspace(0.0, 2 * np.pi, 200, endpoint=False)
        rmc0 = self._var('RMC00')[:, -1]
        ymc0 = self._var('YMC00')[:, -1]
        R = rmc0[:, None] * np.ones_like(theta)
        Z = ymc0[:, None] * np.ones_like(theta)
        for k in range(1, n + 1):
            rmc = self._var(f'RMC0{k}')[:, -1][:, None]
            rms = self._var(f'RMS0{k}')[:, -1][:, None]
            ymc = self._var(f'YMC0{k}')[:, -1][:, None]
            yms = self._var(f'YMS0{k}')[:, -1][:, None]
            R = R + rmc * np.cos(k * theta) + rms * np.sin(k * theta)
            Z = Z + ymc * np.cos(k * theta) + yms * np.sin(k * theta)
        a = (R.max(axis=1) - R.min(axis=1)) / 2.0
        z_axis = self._var('YAXIS')
        kappa_upper = (Z.max(axis=1) - z_axis) / a
        kappa_lower = (z_axis - Z.min(axis=1)) / a
        return kappa_upper, kappa_lower

    def get_boundary_contour(self, time_index, ntheta=400):
        """Full LCFS (R, Z) contour [m] at a given time index, reconstructed
        from TRANSP's boundary Fourier moments (RMC0n/RMS0n/YMC0n/YMS0n) —
        the same role as an eqdsk's `rzout`, for feeding TokaMaker's
        create_isoflux / xpoints_from_moments with real shot geometry instead
        of the notebook's illustrative shape.
        """
        n = self._NUM_BOUNDARY_MOMENTS
        theta = np.linspace(0.0, 2 * np.pi, ntheta, endpoint=False)
        R = self._var('RMC00')[time_index, -1] * np.ones_like(theta)
        Z = self._var('YMC00')[time_index, -1] * np.ones_like(theta)
        for k in range(1, n + 1):
            R = R + (self._var(f'RMC0{k}')[time_index, -1] * np.cos(k * theta)
                      + self._var(f'RMS0{k}')[time_index, -1] * np.sin(k * theta))
            Z = Z + (self._var(f'YMC0{k}')[time_index, -1] * np.cos(k * theta)
                      + self._var(f'YMS0{k}')[time_index, -1] * np.sin(k * theta))
        return np.stack([R / 100.0, Z / 100.0], axis=1)  # cm -> m, shape (ntheta, 2)

    def get_flux_surface_grid(self, time_index, ntheta=200):
        """R(rho, theta), Z(rho, theta) [m] - the full poloidal cross
        section (not just the LCFS) at one time index, reconstructed from
        TRANSP's boundary Fourier moments (RMC0n/RMS0n/YMC0n/YMS0n). These
        carry this same radial dependence at every XB grid index -
        get_boundary_contour only reads their outermost ([:, -1]) index;
        this method is the same expansion evaluated at every radial index,
        for get_reconstructed_Ip's Grad-Shafranov area integral, which
        needs (R, Z) over the whole cross-section, not just its boundary.

        Returns {'R': (nxb+1, ntheta), 'Z': (nxb+1, ntheta), 'rho':
        (nxb+1,), 'theta': (ntheta,)} [m, m, dimensionless rho_tor,
        radians]. The XB grid's innermost point is rho~0.02, not exactly
        the magnetic axis - a degenerate rho=0 "surface" (R=RAXIS, Z=YAXIS
        at every theta) is prepended so the grid spans rho in [0, 1] with
        no core-region gap that an area integral over it would otherwise
        miss.

        Same moment-expansion truncation caveat as get_boundary_contour /
        get_geometry's kappa_upper/kappa_lower (n=8, finite ntheta)
        applies at every radial index here, not just the LCFS - this is
        TRANSP's own internal flux-surface representation, not a
        substitute for a real equilibrium reconstruction.
        """
        n = self._NUM_BOUNDARY_MOMENTS
        theta = np.linspace(0.0, 2 * np.pi, ntheta, endpoint=False)
        R = self._var('RMC00')[time_index][:, None] * np.ones_like(theta)
        Z = self._var('YMC00')[time_index][:, None] * np.ones_like(theta)
        for k in range(1, n + 1):
            rmc = self._var(f'RMC0{k}')[time_index][:, None]
            rms = self._var(f'RMS0{k}')[time_index][:, None]
            ymc = self._var(f'YMC0{k}')[time_index][:, None]
            yms = self._var(f'YMS0{k}')[time_index][:, None]
            R = R + rmc * np.cos(k * theta) + rms * np.sin(k * theta)
            Z = Z + ymc * np.cos(k * theta) + yms * np.sin(k * theta)
        R = R / 100.0  # cm -> m, (nxb, ntheta)
        Z = Z / 100.0

        R_axis = self._var('RAXIS')[time_index] / 100.0
        Z_axis = self._var('YAXIS')[time_index] / 100.0
        rho_full = np.concatenate([[0.0], self.rho_b])
        R_full = np.vstack([np.full(ntheta, R_axis), R])
        Z_full = np.vstack([np.full(ntheta, Z_axis), Z])

        return {'R': R_full, 'Z': Z_full, 'rho': rho_full, 'theta': theta}

    # -- heating (for set_heating / load_TORAX_config) ---------------------

    def get_nbi_heating(self, rho_grid):
        """NBI heating power, plus a Gaussian fit to its deposition profile.

        Returns P_total(t) [W] from PINJ (injected power — includes
        shine-through/CX losses that TORAX's Gaussian model doesn't represent,
        so absorbed power is systematically lower than this in reality; the
        volume integral of PBI+PBE, i.e. absorbed power, is also returned for
        comparison), plus deposition_location(t)/deposition_width(t) fit from
        the actual PBI+PBE profile shape.

        LIMITATION: TORAX's generic heat source is a single time-dependent
        Gaussian (location, width) — it cannot reproduce an arbitrary
        deposition profile. This fits the closest Gaussian (weighted
        centroid/std) to the real profile; multi-lobed or hollow deposition
        shapes will be represented poorly. Use get_absorbed_power_profile if
        you need the full shape (e.g. to build a 'scaled_profile'-style
        source instead of a Gaussian one).
        """
        if not (self.has('PBI') and self.has('PBE') and self.has('PINJ')):
            return None

        pbi = self._var('PBI')
        pbe = self._var('PBE')
        dep = pbi + pbe  # W/cm^3, on the X grid
        loc = np.empty(dep.shape[0])
        width = np.empty(dep.shape[0])
        for it in range(dep.shape[0]):
            loc[it], width[it] = _weighted_centroid_width(self.rho, dep[it])

        dvol = self._var('DVOL')  # cm^3
        absorbed_power = (dep * dvol).sum(axis=1)  # W

        return {
            'P_total': self.to_time_dict(self._var('PINJ')),
            'P_absorbed': self.to_time_dict(absorbed_power),
            'deposition_location': self.to_time_dict(loc),
            'deposition_width': self.to_time_dict(width),
            # nanmean: loc/width are nan at timesteps where the source is momentarily
            # zero (eg. a chopped/modulated beam) - a plain mean would propagate those
            # nans into the scalar average used by set_heating.
            'deposition_location_avg': float(np.nanmean(loc)),
            'deposition_width_avg': float(np.nanmean(width)),
        }

    def get_ech_heating(self, rho_grid):
        """ECRH power and deposition shape, for set_heating(ecrh=..., ecrh_loc=...,
        ecrh_width=...). Uses PEECH (electron deposition, W/cm^3) - ECRH is ~100%
        absorbed by electrons, so P_total is that profile's volume integral, not
        a separate injected-power variable like NBI's PINJ.

        Same Gaussian-fit limitation as get_nbi_heating. Variable name not
        verified against a shot that actually has ECRH - this file doesn't.
        """
        if not self.has('PEECH'):
            return None
        dep = self._var('PEECH')  # W/cm^3
        dvol = self._var('DVOL')  # cm^3
        p_total = (dep * dvol).sum(axis=1)  # W

        loc = np.empty(dep.shape[0])
        width = np.empty(dep.shape[0])
        for it in range(dep.shape[0]):
            loc[it], width[it] = _weighted_centroid_width(self.rho, dep[it])

        return {
            'P_total': self.to_time_dict(p_total),
            'deposition_location_avg': float(np.nanmean(loc)),
            'deposition_width_avg': float(np.nanmean(width)),
        }

    def get_ich_heating(self, rho_grid):
        """ICRH power and deposition shape. set_heating has no dedicated ICRH
        slot (only generic_heat/ecrh) - add via tmtx.load_TORAX_config's
        'scaled_profile' source model instead. PIICRH/PEICRH (ion/electron
        deposition, W/cm^3) and PORICRH (total injected power) are TRANSP's
        ICRH variables.

        Same Gaussian-fit limitation as get_nbi_heating. Variable names not
        verified against a shot that actually has ICRH - this file doesn't.
        """
        if not (self.has('PIICRH') and self.has('PEICRH')):
            return None
        dep = self._var('PIICRH') + self._var('PEICRH')  # W/cm^3
        dvol = self._var('DVOL')
        absorbed = (dep * dvol).sum(axis=1)  # W

        loc = np.empty(dep.shape[0])
        width = np.empty(dep.shape[0])
        for it in range(dep.shape[0]):
            loc[it], width[it] = _weighted_centroid_width(self.rho, dep[it])

        p_total = self._var('PORICRH') if self.has('PORICRH') else absorbed
        return {
            'P_total': self.to_time_dict(p_total),
            'P_absorbed': self.to_time_dict(absorbed),
            'deposition_location_avg': float(np.nanmean(loc)),
            'deposition_width_avg': float(np.nanmean(width)),
        }

    def get_absorbed_power_profile(self, rho_grid):
        """Full NBI absorbed-power deposition profile [W/m^3], {time: {rho: value}}.

        For building an explicit-profile source (TORAX's 'scaled_profile'
        model, as used for ICRH in the example notebook) instead of fitting a
        Gaussian — use this if your TORAX/TokaMaker_TORAX version exposes an
        explicit-profile option for NBI heating.
        """
        pbi = self.profile('PBI', rho_grid) * 1e6  # W/cm^3 -> W/m^3
        pbe = self.profile('PBE', rho_grid) * 1e6
        return {
            'ion': self.to_time_rho_dict(pbi, rho_grid),
            'electron': self.to_time_rho_dict(pbe, rho_grid),
        }

    # -- particle sources (for set_fueling) --------------------------------

    def get_particle_source(self, varname):
        """Volume-integrate a TRANSP particle-source-rate profile (N/CM3/SEC)
        into S_total(t) [particles/s], plus a Gaussian fit (location, width)
        of its deposition profile — for set_fueling's gas_puff/pellet/
        generic_particle *_S_total / *_location / *_width kwargs.

        Same Gaussian-only limitation as get_nbi_heating applies.
        """
        if not self.has(varname):
            return None
        src = self._var(varname)          # (time, X), N/cm^3/s
        dvol = self._var('DVOL')          # (time, X), cm^3
        s_total = (src * dvol).sum(axis=1)  # particles/s

        loc = np.empty(src.shape[0])
        width = np.empty(src.shape[0])
        for it in range(src.shape[0]):
            loc[it], width[it] = _weighted_centroid_width(self.rho, src[it])

        return {
            'S_total': self.to_time_dict(s_total),
            'location': self.to_time_dict(loc),
            'width': self.to_time_dict(width),
            # nanmean - see get_nbi_heating's comment on the same pattern.
            'location_avg': float(np.nanmean(loc)),
            'width_avg': float(np.nanmean(width)),
        }

    def get_gas_puff_decay_length(self, varname, decay_start=1.0):
        """Exponential decay-length fit (per-timestep + nanmean) of a TRANSP
        edge particle-source profile, for set_fueling's gas_puff_decay_length.

        TORAX's gas_puff source is a single exponential decaying inward from
        `decay_start` (the edge), not a Gaussian — get_particle_source's
        Gaussian-moment fit (used for pellet/generic_particle, which really
        are Gaussian in TORAX) fits this functional form poorly, especially
        for edge-peaked profiles like recycling that keep rising all the way
        to the domain boundary.
        """
        if not self.has(varname):
            return None
        src = self._var(varname)  # (time, X), N/cm^3/s
        decay = np.array([
            _fit_exponential_decay_length(self.rho, src[it], decay_start)
            for it in range(src.shape[0])
        ])
        return {
            'decay_length': self.to_time_dict(decay),
            'decay_length_avg': float(np.nanmean(decay)),
        }

    def get_particle_source_profile(self, varname, rho_grid):
        """Full particle-source-rate profile [particles/m^3/s], {time: {rho: value}}.

        For checking `get_particle_source`'s Gaussian fit (location/width)
        against the actual source shape, the same role
        `get_absorbed_power_profile` plays for `get_nbi_heating`.
        """
        if not self.has(varname):
            return None
        src = self.profile(varname, rho_grid) * 1e6  # N/cm^3/s -> N/m^3/s
        return self.to_time_rho_dict(src, rho_grid)

    def get_fueling(self):
        """Best-effort split of edge/beam particle sources into gas-puff-like
        and beam-fueling-like sources, for set_fueling.

        TRANSP's nominal external-gas-puff variable (SESGF) is often
        negligible compared to the recycling source (SERC_D) — as it is in
        this file (~1 particle/s vs ~1e22/s) — because most DIII-D edge
        fueling comes from wall recycling, not the puff valve, once the
        plasma is loaded. This returns *both* SESGF and SERC_D so you can
        see which actually dominates for your shot, rather than silently
        picking one.
        """
        return {
            'gas_puff_SESGF': self.get_particle_source('SESGF'),
            'recycling_SERC_D': self.get_particle_source('SERC_D'),
            'recycling_SERC_D_decay_length': self.get_gas_puff_decay_length('SERC_D'),
            'beam_fueling_SBE': self.get_particle_source('SBE'),
        }

    # -- derived profiles (for TORAX-vs-TRANSP profile comparison) ----------

    def get_current_density(self, rho_grid):
        """Current density components [A/m^2], {name: {time: {rho: value}}},
        for comparison against TORAX's own j_total/j_ohmic/j_bootstrap/
        j_generic_current/j_non_inductive profiles.

        - 'total': CUR, TRANSP's total toroidal plasma current density, for
          TORAX's j_total.
        - 'ohmic': CUROH, for TORAX's j_ohmic.
        - 'bootstrap': CURBS, for TORAX's j_bootstrap. This file carries
          several bootstrap-current models (CURBSSAU/CURBSNEO/CURBSWNC/...);
          CURBS is the one this run actually summed into CUR - close to but
          not exactly CURBSSAU (the Sauter formula "as used"), which is a
          slightly different variant.
        - 'beam_driven': CURB (NBI-driven current), for TORAX's
          j_generic_current - this notebook's set_heating(nbi_current=True)
          routes NBI-driven current through TORAX's generic-current channel,
          not j_external.
        - 'non_inductive': CURXT ("driven plasma current, smoothed") =
          CURBS + CURB to ~1e-5 relative precision in this file, for
          TORAX's j_non_inductive.

        All are read from TRANSP's zone-centered X grid and resampled onto
        `rho_grid`, same as every other profile getter in this module.
        Missing channels return None instead of a fabricated zero profile.
        """
        scale = 1e4  # A/cm^2 -> A/m^2
        out = {}
        for key, varname in [('total', 'CUR'), ('ohmic', 'CUROH'),
                              ('bootstrap', 'CURBS'), ('beam_driven', 'CURB'),
                              ('non_inductive', 'CURXT')]:
            if not self.has(varname):
                out[key] = None
                continue
            arr = self.profile(varname, rho_grid) * scale
            out[key] = self.to_time_rho_dict(arr, rho_grid)
        return out

    def get_q_profile(self, rho_grid):
        """Safety factor q(rho, t) [dimensionless], {time: {rho: value}] -
        TRANSP's Q profile, for comparison against TORAX's own q profile.

        Q is defined on TRANSP's zone-boundary (XB) grid, not the
        zone-centered X grid `profile()` defaults to (same situation as
        get_transport_coefficients's CONDE/CONDI).
        """
        arr = self.profile('Q', rho_grid, source_rho='XB')
        return self.to_time_rho_dict(arr, rho_grid)

    def get_ion_density(self, rho_grid):
        """Total ion density n_i [m^-3], {time: {rho: value}} - TRANSP's NI,
        for comparison against TORAX's own n_i profile."""
        arr = self.profile('NI', rho_grid) * 1e6  # N/cm^3 -> N/m^3
        return self.to_time_rho_dict(arr, rho_grid)

    def get_pressure_profile(self, rho_grid):
        """Thermal plasma pressure [Pa], {time: {rho: value}} - TRANSP's
        PMHDT_IN (thermal-only component of the pressure fed to its MHD
        solver; identical to PPLAS in this file), for TORAX's
        pressure_thermal_total. Not the same variable as get_geometry's
        'pax' (PMHD_IN, on-axis *total* pressure incl. fast-ion/rotation) -
        this is the full radial profile, thermal component only, matching
        TORAX's quantity in kind as well as location.
        """
        arr = self.profile('PMHDT_IN', rho_grid)  # already Pa
        return self.to_time_rho_dict(arr, rho_grid)

    def get_psi_profile(self, rho_grid):
        """Poloidal flux profile psi(rho, t) [Wb/rad], {time: {rho: value}} -
        TRANSP's PLFLX, sign-flipped, for comparison against TORAX's own psi
        profile.

        PLFLX is defined on TRANSP's zone-boundary (XB) grid. Same sign
        convention as get_geometry's psi_lcfs (PLFLXA, PLFLX's own edge
        value) - see that method's "Sign convention" docstring note. Still
        only a sign fix, not a gauge/offset fix - compare shape/relative
        evolution here, not absolute normalization.
        """
        arr = -self.profile('PLFLX', rho_grid, source_rho='XB')
        return self.to_time_rho_dict(arr, rho_grid)

    def get_ffprime_pprime_profiles(self, psi_norm_grid=None, smooth_ffprime=True):
        """p' = dP/dpsi [Pa/Wb] and ff' = F*dF/dpsi [T^2*m^2/Wb], {time:
        {'psi_norm': array, 'pprime': array, 'ffprime': array}} - for
        TokaMaker's set_profiles(pp_prof=..., ffp_prof=...), replacing the
        notebook's illustrative placeholder (psi_sample = linspace(0,1,50);
        ffp_prof_ref = pp_prof_ref = 1 - psi_sample) with profiles derived
        directly from this shot's TRANSP data. TRANSP carries no native
        ff'/p' output, so both are built here from more basic quantities:

        - P: PMHD_IN, the *total* pressure fed to TRANSP's own MHD/
          equilibrium solver (thermal + fast-ion + rotation - see
          get_fast_ion_pressure_axis's PMHD_IN = PMHDT_IN + PMHDF_IN +
          PMHDR_IN breakdown), **not** get_pressure_profile's PMHDT_IN
          (thermal-only, matching TORAX's pressure_thermal_total). This
          looks like it should be the thermal-only quantity by analogy
          with get_pressure_profile, but p' here feeds the Grad-Shafranov
          force balance TokaMaker's solve() actually enforces, which
          total pressure sets (same reasoning as get_geometry's 'pax',
          also PMHD_IN, already feeding this same equilibrium fit's
          mygs.set_targets(pax=...)) - thermal-only would be internally
          inconsistent with that. Empirically confirmed by
          get_reconstructed_Ip: using PMHDT_IN reproduced measured Ip to
          <1% only when fast-ion pressure was negligible, but under-shot
          by 8-24% (worse, the larger the fast-ion pressure fraction) at
          every other time slice on this file; switching to PMHD_IN
          brought every slice to ~1-2% agreement, matching the pre-
          existing few-percent floor from get_geometry's own boundary-
          moment-truncation cross-check. On the X grid.
        - F = R*Bt: **not** GFUN directly, despite its long_name ("G:
          PARA/DIAMAGNETISM") reading like the physical quantity itself -
          GFUN is TRANSP's internal boundary-*normalized* diamagnetic
          function (g = F/F_lcfs); verified on this file to be exactly 1.0
          at the outermost XB point at every time slice, not the physical
          F~3.5 T*m a DIII-D-scale F0 would suggest. The physical F is
          recovered as F = GFUN * BZXR (BZXR = vacuum B0*R0 [T*cm], the
          same variable get_geometry's F0 reads) - GFUN's boundary value of
          1.0 times BZXR reproduces F0 exactly, confirming the
          normalization. On the XB grid. (GFUNC is a different variable -
          a Grad-Shafranov equilibrium check, not this diamagnetism
          function - do not substitute it here.)
        - psi: -PLFLX, sign-flipped, same convention as get_psi_profile /
          get_geometry's psi_lcfs (see get_geometry's "Sign convention"
          docstring note). On the XB grid.

        Grid handling: P lives on X, F and psi on XB - both 50 points but
        at different rho_tor locations (see module docstring), so P is
        resampled from X onto XB by linear interpolation (the same
        np.interp pattern profile() uses) before differentiating, giving
        p' and ff' a shared native grid to differentiate against psi on.

        Differentiation is np.gradient(y, psi) at each time, guarded
        against a non-monotonic psi (sorted before differentiating,
        restored to native order after - see _monotonic_gradient); psi
        here is monotonic at every time slice checked on this file, so
        that guard is a no-op safety net rather than a correction this
        file actually exercises.

        psi normalization: psi_norm = (psi - psi_axis)/(psi_lcfs -
        psi_axis). psi_lcfs is this same array's own outermost value (XB's
        last point is exactly rho=1.0, and -PLFLX there equals
        get_geometry's psi_lcfs (-PLFLXA) to machine precision on this
        file - no separate read needed). psi_axis is taken as exactly 0.0,
        not from PSI0 (an unpopulated scalar placeholder in this file -
        shape (), see get_scalars) - PLFLX is referenced from the magnetic
        axis by construction (confirmed here: PLFLX at the innermost XB
        point, rho~0.02, is ~1e-3 of its edge value, consistent with ->0
        as rho->0).

        Output grid: `psi_norm_grid=None` (default) returns each time's own
        *native* psi_norm values (XB's 50 points, sorted, plus an explicit
        psi_norm=0 axis anchor - 51 points total) rather than resampling
        onto a fixed uniform grid. This matters: psi grows very slowly
        with rho near the magnetic axis (a normal, expected equilibrium
        property), so XB's inner ~6-8 points (rho=0.02 to ~0.15) - real,
        distinct TRANSP data - land in only the first few percent of
        psi_norm's [0, 1] range. A uniform output grid (e.g. the notebook's
        old linspace(0, 1, 50) placeholder convention) samples that region
        at only 1-2 points, discarding most of the near-axis structure the
        native XB grid actually resolves; everywhere outside the core the
        native and a uniform grid are comparably dense, so this matters
        specifically for on-axis quantities. CONFIRMED costly in practice:
        feeding a downstream TokaMaker solve a uniform-resampled version of
        these profiles under-resolved the near-axis current density enough
        to visibly bias the solved q_0 high relative to TRANSP's own q0 (an
        independent axis-dense-grid fix closed roughly half that gap on
        this file - see the notebook's seed-equilibrium markdown for the
        measured before/after). Pass an explicit array only if you
        specifically need one fixed, shared grid across every time (e.g.
        for a simple elementwise time-series comparison) - this resamples
        (with the same near-axis information loss described above) onto
        that grid instead. Because the native grid's exact psi_norm values
        differ time to time (each time has its own psi_lcfs), don't reuse
        one time's returned 'psi_norm' array as a fixed x-grid for another
        time's 'pprime'/'ffprime' - each time's three returned arrays are a
        matched set.

        SIGN CAVEAT - read before wiring into set_profiles: with this
        file's -PLFLX convention, psi *decreases* with rho (this is about
        the *radial* behavior at fixed time - a different axis from
        get_geometry's note that TokaMaker's/TORAX's psi decreases through
        the current ramp in *time*). Since pressure also decreases
        outward, dP/dpsi comes out POSITIVE under this convention
        (verified: ~1e4-1e6 Pa/Wb on this file, not negative) - the
        opposite sign from the "textbook"/eqdsk convention where psi
        increases from axis to edge and p'(psi) is negative for a normal
        peaked pressure profile. This is a real consequence of which
        radial direction *this file's* psi convention points, not a bug -
        and it means the placeholder's positive, monotonically-decreasing
        shape (1 - psi_norm) is actually the sign that matches what comes
        out here, while the generic "p' should be negative" intuition
        assumes the other convention. If TokaMaker's solve() comes out
        with the wrong-sign pressure gradient, negate both returned arrays
        (equivalent to differentiating against +PLFLX instead) rather than
        re-deriving from scratch.

        NOISE CAVEAT - ff' is far noisier than p': F only varies ~1-3.5%
        across the whole profile (see GFUN's range), so dF/dpsi amplifies
        GFUN's point-to-point differencing noise much more than PMHD_IN's
        smooth, order-unity pressure drop does for p' - confirmed here:
        the raw derivative swings sign multiple times within a single
        time slice and reaches |ff'| ~ 1 T^2*m^2/Wb between adjacent grid
        points, an order of magnitude too jagged to be physical, worst
        right near the axis (a handful of points span a sign change).
        Tried and rejected as fixes: a Savitzky-Golay pre-filter on GFUN,
        a global low-order polynomial fit to F(psi), and a local polynomial
        fit through just the innermost few points - all three *increased*
        get_reconstructed_Ip's PCUR disagreement (to varying degrees, worst
        for the global fit) rather than reducing it, because each one
        perturbs the bulk of the profile - which dominates the Ip
        integral - to chase a noise source concentrated where the area
        element is smallest. `smooth_ffprime=True` (default) therefore
        does nothing to the bulk array; it only replaces the single
        innermost-psi point via _fix_axis_point's 2-point extrapolation
        (see that function's docstring) rather than leaving TRANSP's own
        possibly-noisy axis-adjacent sample in place - confirmed neutral
        for PCUR agreement either way (Ip barely sees the axis point), so
        this is a bet on physical reasoning (removing the zero-slope flat
        spot right at the axis) rather than something this check can
        itself confirm helps; pass smooth_ffprime=False to skip even that
        and get the fully raw, untouched derivative - plot both before
        trusting either.

        VALIDATION: get_reconstructed_Ip cross-checks these profiles
        end-to-end by reconstructing Ip via the Grad-Shafranov relation
        and comparing against measured PCUR - see that method's docstring
        for the result and what it does/doesn't rule out.
        """
        psi_norm_xb, pprime_xb, ffprime_xb = self._ffprime_pprime_native(smooth_ffprime)

        out = {}
        for it, t in enumerate(self.time):
            order = np.argsort(psi_norm_xb[it])
            psi_sorted = psi_norm_xb[it][order]
            pp_sorted = pprime_xb[it][order]
            ffp_sorted = ffprime_xb[it][order]
            # native grid's innermost point is close to but not exactly psi_norm=0 -
            # prepend an explicit axis anchor (see docstring's "Output grid" note).
            # pp_sorted[0]/ffp_sorted[0] are already axis-corrected by
            # _ffprime_pprime_native's own _fix_axis_point call, so just reuse them.
            psi_native = np.concatenate([[0.0], psi_sorted])
            pp_native = np.concatenate([[pp_sorted[0]], pp_sorted])
            ffp_native = np.concatenate([[ffp_sorted[0]], ffp_sorted])

            if psi_norm_grid is None:
                out[float(t)] = {'psi_norm': psi_native, 'pprime': pp_native, 'ffprime': ffp_native}
            else:
                grid = np.asarray(psi_norm_grid, dtype=float)
                out[float(t)] = {
                    'psi_norm': grid,
                    'pprime': np.interp(grid, psi_native, pp_native),
                    'ffprime': np.interp(grid, psi_native, ffp_native),
                }
        return out

    def get_transp2imas_profiles_1d(self, smooth_ffprime=True):
        """equilibrium.time_slice.profiles_1d fields - psi, phi, q, pressure,
        dpressure_dpsi, f, f_df_dpsi - built directly from this shot's TRANSP
        CDF the same way transp/transp-imas-translator's transp2imas.f90
        does (github.com/transp/transp-imas-translator, transp2imas.f90's
        "Calculate p and p'"/"calculate dF/dPLFLX" blocks; cross-checked
        against imas2transp/imas-eq-q.py, which validates exactly these two
        fields - profiles_1d.dpressure_dpsi and .f_df_dpsi - against TRANSP's
        own fort.313/fort.314 dumps). See the IMAS data dictionary's
        equilibrium.time_slice.profiles_1d.psi entry for the target field
        this feeds.

        phi (toroidal flux, Wb) is -TRFLX, and q (safety factor) is read
        directly off TRANSP's own 'Q' - both on XB, both validated together
        on a sample CDF: raw d(TRFLX)/d(PLFLX2PI) (no sign flips at all)
        matches the CDF's own 'Q' to <0.2% at every interior point, so
        phi's "-" sign has to match whatever sign psi uses (here, psi's
        validated "-PLFLX2PI") for q = dphi/dpsi to keep the right sign -
        it's a companion to psi's sign choice, not an independent one.

        Differs from get_ffprime_pprime_profiles() in one load-bearing way:
        **psi here is the CDF's PLFLX2PI ("TOTAL POLOIDAL FLUX", WEBERS),
        not PLFLX ("POLOIDAL FLUX", Wb/rad)** that
        get_ffprime_pprime_profiles/get_geometry read. Confirmed empirically
        on this file: PLFLX2PI == 2*pi*PLFLX to 6 significant figures at
        every point, every time. transp2imas.f90 reads PLFLX2PI (into a
        confusingly-named local array it calls "PLFLX") for
        profiles_1d%psi - see its "Calculate dXB/dPLFLX, it will be used
        later to calculate p' and FF'" block - i.e. the official translator
        uses the genuine total poloidal flux, not flux-per-radian, for the
        IMAS dictionary's Wb-unit psi field.

        Practical effect: because psi is 2*pi times larger here,
        dpressure_dpsi and f_df_dpsi both come out a factor of 2*pi
        *smaller* than get_ffprime_pprime_profiles's 'pprime'/'ffprime' at
        the same physical point (psi_norm itself is unaffected - the 2*pi
        cancels in that ratio). Same PMHD_IN/GFUN/BZXR quantities, sign
        flip, GFUN smoothing and native-XB-grid reasoning as
        get_ffprime_pprime_profiles - see that method's docstring for the
        full derivation and the sign/noise caveats, which apply unchanged
        here.

        CHOSEN as the TokaMaker-feeding cell's replacement despite that gap
        (see the notebook markdown above that cell): get_reconstructed_Ip's
        validation of the Wb/rad convention does not carry over here
        automatically, since dpressure_dpsi/f_df_dpsi both shift by 2*pi -
        re-run an equivalent PCUR/q0 check against this method's output
        before trusting the solved equilibria the same way.

        AXIS ANCHOR: pressure/dpressure_dpsi/q/f/f_df_dpsi are all linearly
        extrapolated to the axis from the two innermost native points (see
        _linear_extrapolate_to_axis and _fix_axis_point) rather than held
        flat - a flat hold creates a zero-slope segment right at the axis
        that under-drives on-axis current density in a downstream
        Grad-Shafranov solve. f/f_df_dpsi get this same treatment despite
        being noisy near the axis (see get_ffprime_pprime_profiles's NOISE
        CAVEAT); whether that's actually better than a flat hold for THIS
        field couldn't be settled by get_reconstructed_Ip's PCUR check - Ip
        integrates current density over the whole cross-section, weighted
        by area that's tiny near the axis and large near the edge, so it's
        structurally almost insensitive to what happens at a single axis
        point (confirmed: extrapolating vs. holding flat changed PCUR
        agreement by <0.01 percentage points on a sample CDF - noise-level,
        not a real signal either way). q_0 specifically, by contrast, is a
        purely local axis quantity that this check can't stand in for -
        only a direct TokaMaker/G-S re-solve can confirm whether this
        actually narrows a q_0 gap.

        Returns {time: {'psi_norm': array, 'psi': array, 'phi': array, 'q':
        array, 'pressure': array, 'dpressure_dpsi': array, 'f': array,
        'f_df_dpsi': array}}. Each array is sorted by ascending psi_norm
        and axis-anchored (see above) - 51 points total, matching
        get_ffprime_pprime_profiles's native-grid output - see that
        method's "Output grid" docstring note for why (a uniform grid
        under-resolves the near-axis region). psi_norm = (psi -
        psi_axis)/(psi_lcfs - psi_axis) with psi_axis=0.0 and
        psi_lcfs=psi's own outermost native value - identical
        normalization to get_ffprime_pprime_profiles, and identical
        numerically too (the 2*pi cancels in the ratio), so psi_norm is
        interchangeable between the two methods even though psi itself
        isn't.
        """
        pressure_x = self._var('PMHD_IN')     # (time, nx), Pa, on X - total pressure
        F0 = self._var('BZXR') / 100.0        # T*cm -> T*m, (time,)
        F_xb = self._var('GFUN') * F0[:, None]  # T*m, on XB - raw, noisy near axis (see docstring)
        psi_xb = -self._var('PLFLX2PI')       # Wb (real, 2*pi-included) - see docstring
        phi_xb = -self._var('TRFLX')          # Wb, toroidal flux - same sign as psi, see docstring
        q_xb = self._var('Q')                 # dimensionless safety factor, on XB - no sign flip

        out = {}
        for it, t in enumerate(self.time):
            p_xb = np.interp(self.rho_b, self.rho, pressure_x[it])  # X -> XB
            psi = psi_xb[it]
            F = F_xb[it]
            dpressure_dpsi = _fix_axis_point(psi, _monotonic_gradient(p_xb, psi))
            ffp = F * _monotonic_gradient(F, psi)
            f_df_dpsi = _fix_axis_point(psi, ffp) if smooth_ffprime else ffp

            psi_axis = 0.0
            psi_lcfs = psi[-1]  # XB's outermost point is exactly rho=1.0
            psi_norm = (psi - psi_axis) / (psi_lcfs - psi_axis)
            q = _fix_axis_point(psi_norm, q_xb[it])
            p_fixed = _fix_axis_point(psi_norm, p_xb)
            F_fixed = _fix_axis_point(psi_norm, F) if smooth_ffprime else F

            order = np.argsort(psi_norm)
            fields = {
                'psi_norm': psi_norm[order],
                'psi': psi[order],
                'phi': phi_xb[it][order],
                'q': q[order],
                'pressure': p_fixed[order],
                'dpressure_dpsi': dpressure_dpsi[order],
                'f': F_fixed[order],
                'f_df_dpsi': f_df_dpsi[order],
            }
            # native grid's innermost point is close to but not exactly the
            # axis - prepend an explicit anchor (see docstring's AXIS ANCHOR
            # note): psi_norm/psi/phi vanish at the axis by construction
            # (exact 0.0); everything else already has its own innermost
            # point axis-corrected above via _fix_axis_point, so just reuse it.
            exact_axis_value = {'psi_norm': 0.0, 'psi': psi_axis, 'phi': 0.0}
            axis_row = {
                key: exact_axis_value.get(key, arr[0])
                for key, arr in fields.items()
            }
            out[float(t)] = {
                key: np.concatenate([[axis_row[key]], arr])
                for key, arr in fields.items()
            }
        return out

    def get_transp2imas_current_profiles(self):
        """core_profiles.time_slice.profiles_1d current-density fields -
        j_tor, j_ohmic, j_bootstrap, j_non_inductive [A/m^2] - built
        directly from this shot's TRANSP CDF the same way
        transp/transp-imas-translator's transp2imas.f90 does (its "CUR"/
        "CUROH"/"CURBS" blocks), for future TORAX current-source inputs.
        On the X (zone-center) grid, matching TRANSP's own CUR/CUROH/CURBS
        grid - not interpolated onto XB like get_transp2imas_profiles_1d's
        fields, since there's no psi array on this grid to differentiate
        against here (TORAX's current sources are keyed by rho, not psi).

        SIGN CORRECTION vs. the official translator: transp2imas.f90 negates
        all three (CUR, CUROH, CURBS) under the same "Hardcode minus sign
        for ITER for now" comment used throughout that file. Verified WRONG
        for this dataset - on a sample CDF, sum(CUR * DAREA) == PCUR exactly
        (ratio 1.0000) with no sign flip at all - and, since this codebase's
        own PCUR is already used unflipped everywhere else (get_geometry,
        get_reconstructed_Ip), almost certainly wrong for DIII-D too, not
        just this sample. j_ohmic/j_bootstrap aren't independently
        re-verified the same way, but they're components of the same TRANSP
        current budget as CUR (CUR ~= CUROH + CURBS + any auxiliary-driven
        current), so they're read here with the matching, un-flipped sign
        rather than transp2imas.f90's ITER-specific one.

        Auxiliary-driven current (NBI/EC/IC, TRANSP's CURB/ECCUR/ICCUR) maps
        to per-source core_sources.source[].profiles_1d.j_parallel in the
        official translator, not here - which of those exist depends on
        which heating systems this shot actually used. Add them the same
        way (read the CDF variable, do NOT apply transp2imas.f90's sign
        flip without re-deriving it) if/when needed.

        Returns {time: {'j_tor': array, 'j_ohmic': array, 'j_bootstrap':
        array, 'j_non_inductive': array}}, each (nx,).
        j_non_inductive = j_tor - j_ohmic, matching transp2imas.f90's own
        definition.
        """
        j_tor_x = self._var('CUR') * 1.0e4          # A/cm^2 -> A/m^2, on X
        j_ohmic_x = self._var('CUROH') * 1.0e4      # A/cm^2 -> A/m^2, on X
        j_bootstrap_x = self._var('CURBS') * 1.0e4  # A/cm^2 -> A/m^2, on X

        out = {}
        for it, t in enumerate(self.time):
            j_tor = j_tor_x[it]
            j_ohmic = j_ohmic_x[it]
            out[float(t)] = {
                'j_tor': j_tor,
                'j_ohmic': j_ohmic,
                'j_bootstrap': j_bootstrap_x[it],
                'j_non_inductive': j_tor - j_ohmic,
            }
        return out

    def get_geqdsk_dict(self, time_index):
        """Raw G-EQDSK dict for the equilibrium at run.time[time_index],
        built entirely from this shot's TRANSP CDF - PSIRZ for the 2D flux
        map, get_ffprime_pprime_profiles's pipeline for the 1D F/p/q/p'/ff'
        profiles, get_boundary_contour for the plasma boundary - for
        OpenFUSIONToolkit.TokaMaker.eqdsk.GEQDSKEquilibrium.from_raw()/
        write_geqdsk(), bypassing TokaMaker's own Grad-Shafranov solve
        entirely. Use this to seed TokaMaker/TORAX equilibria straight from
        TRANSP's own already-self-consistent reconstruction, rather than an
        independently re-solved one.

        COCOS: cocos=2 (matching this notebook's own existing
        mygs.save_eqdsk(cocos=2) convention), with psi in the CDF's native
        Wb/rad convention (-PLFLX/-PSIRZ) - **not** the 2*pi-scaled
        -PLFLX2PI convention get_transp2imas_profiles_1d uses for the IMAS
        equilibrium IDS's psi field. These are genuinely different targets
        requiring opposite conventions - confirmed by an end-to-end check
        independent of anything assumed elsewhere:
        GEQDSKEquilibrium's own flux-surface-averaged q_profile, computed
        purely from this dict's 2D psi(R,Z) map and FPOL via a completely
        different code path than TRANSP's own 'Q', reproduces 'Q' to <1%
        across the *entire* profile with this exact combination (cocos=2 +
        Wb/rad psi). Using cocos=2 with the 2*pi-scaled convention instead
        reproduces the right shape but a uniform factor of 2*pi off in
        magnitude - get_transp2imas_profiles_1d's convention is correct for
        the IMAS dictionary's own psi field, but not for this GEQDSK
        target; don't conflate the two.

        2D psi(R,Z): TRANSP's PSIRZ (Wb/rad, on the RGRID x ZGRID
        rectangular grid) reshaped as PSIRZ[time_index].reshape(len(ZGRID),
        len(RGRID)) - verified this ordering's minimum lands exactly at
        (RAXIS, YAXIS); the transposed ordering doesn't - then sign-flipped
        (-PSIRZ) to match this module's -PLFLX convention everywhere else -
        confirmed empirically that raw PSIRZ shares PLFLX's exact values
        along the midplane (same sign, same magnitude, same rho
        dependence), so the same "-" flip applies. RGRID/ZGRID are
        themselves time-independent despite being stored per-timestep
        (confirmed via np.allclose across all times) and already regularly
        spaced, so no further grid construction is needed.

        1D profiles (FPOL, PRES, FFPRIM, PPRIME, QPSI): built from GFUN*F0
        (F), PMHD_IN (pressure), TRANSP's own 'Q', and
        _ffprime_pprime_native's already axis-corrected p'/ff' - all on
        their native XB grid - resampled onto a uniform psi_norm grid of
        NW=len(RGRID) points (the standard GEQDSK convention: profile
        arrays are always length NW, matching the R-grid point count).

        Boundary: get_boundary_contour(time_index) (TRANSP's own LCFS
        Fourier moments) - reused for both RBBBS/ZBBBS and, as a
        placeholder (TRANSP carries no separate wall/limiter contour in
        this file), RLIM/ZLIM.

        SIMAG/SIBRY: SIMAG is psi_2D's own value at the (RAXIS, YAXIS) grid
        cell (the genuine, slightly nonzero, grid-quantized value - not
        forced to exactly 0.0 like get_ffprime_pprime_profiles's psi_axis
        convenience convention - for exact self-consistency with the
        PSIRZ array being written); SIBRY is -PLFLXA (this module's
        existing, validated LCFS flux convention).

        RCENTR/BCENTR: RCENTR is taken as RAXIS (this file has no separate
        machine-reference-radius scalar); BCENTR = F0/RCENTR with
        F0=BZXR/100 (T*m) - their product is what matters physically
        (R*Bt), not either value individually.

        Returns a dict with exactly the keys write_geqdsk/
        GEQDSKEquilibrium.from_raw expect: NW, NH, RDIM, ZDIM, RCENTR,
        RLEFT, ZMID, RMAXIS, ZMAXIS, SIMAG, SIBRY, BCENTR, CURRENT, FPOL,
        PRES, FFPRIM, PPRIME, PSIRZ, QPSI, NBBBS, RBBBS, ZBBBS, LIMITR,
        RLIM, ZLIM. Pass to
        OpenFUSIONToolkit.TokaMaker.eqdsk.GEQDSKEquilibrium.from_raw(d,
        cocos=2) or write_geqdsk(d, filename) directly.
        """
        rgrid = self._var('RGRID')[time_index] / 100.0   # cm -> m
        zgrid = self._var('ZGRID')[time_index] / 100.0
        nR, nZ = len(rgrid), len(zgrid)
        psi2d = -self._var('PSIRZ')[time_index].reshape(nZ, nR)  # Wb/rad, -PLFLX convention

        raxis = self._var('RAXIS')[time_index] / 100.0
        zaxis = self._var('YAXIS')[time_index] / 100.0
        F0 = self._var('BZXR')[time_index] / 100.0  # T*m
        pcur = self._var('PCUR')[time_index]

        psi_norm_native, pprime_native, ffprime_native = self._ffprime_pprime_native()
        psi_norm_it = psi_norm_native[time_index]
        order = np.argsort(psi_norm_it)
        psi_norm_sorted = psi_norm_it[order]

        F_xb = self._var('GFUN')[time_index] * F0
        p_xb = np.interp(self.rho_b, self.rho, self._var('PMHD_IN')[time_index])
        q_xb = self._var('Q')[time_index]

        NW = nR
        psi_uniform = np.linspace(0.0, 1.0, NW)
        fpol = np.interp(psi_uniform, psi_norm_sorted, F_xb[order])
        pres = np.interp(psi_uniform, psi_norm_sorted, p_xb[order])
        qpsi = np.interp(psi_uniform, psi_norm_sorted, q_xb[order])
        pprime = np.interp(psi_uniform, psi_norm_sorted, pprime_native[time_index][order])
        ffprim = np.interp(psi_uniform, psi_norm_sorted, ffprime_native[time_index][order])

        ir_axis = np.argmin(np.abs(rgrid - raxis))
        iz_axis = np.argmin(np.abs(zgrid - zaxis))
        simag = psi2d[iz_axis, ir_axis]
        sibry = -self._var('PLFLXA')[time_index]

        boundary = self.get_boundary_contour(time_index, ntheta=200)

        return {
            'NW': NW, 'NH': nZ,
            'RDIM': rgrid[-1] - rgrid[0], 'ZDIM': zgrid[-1] - zgrid[0],
            'RCENTR': raxis, 'RLEFT': rgrid[0], 'ZMID': 0.5 * (zgrid[0] + zgrid[-1]),
            'RMAXIS': raxis, 'ZMAXIS': zaxis,
            'SIMAG': simag, 'SIBRY': sibry,
            'BCENTR': F0 / raxis, 'CURRENT': pcur,
            'FPOL': fpol, 'PRES': pres, 'FFPRIM': ffprim, 'PPRIME': pprime,
            'PSIRZ': psi2d, 'QPSI': qpsi,
            'NBBBS': len(boundary), 'RBBBS': boundary[:, 0], 'ZBBBS': boundary[:, 1],
            'LIMITR': len(boundary), 'RLIM': boundary[:, 0], 'ZLIM': boundary[:, 1],
        }

    def _ffprime_pprime_native(self, smooth_ffprime=True):
        """psi_norm, pprime, ffprime on the native XB grid (each (ntime,
        nxb)), before resampling onto a uniform psi_norm grid - the shared
        computation behind get_ffprime_pprime_profiles (which resamples
        this onto a fixed output grid) and get_reconstructed_Ip (which
        needs p'/ff' evaluated at the specific psi_norm values its flux-
        surface geometry grid lands on). See get_ffprime_pprime_profiles's
        docstring for the full derivation and caveats (GFUN normalization,
        sign convention, differencing noise/smoothing).

        Both arrays' innermost-psi point is corrected via
        _fix_axis_point (p' unconditionally, ff' when smooth_ffprime=True)
        rather than left as TRANSP's own native axis-adjacent sample -
        see _linear_extrapolate_to_axis's docstring for why, including for
        ff' despite its own near-axis noise. This happens here (not just in
        get_ffprime_pprime_profiles's wrapper) so get_reconstructed_Ip
        benefits too - it uses this method's output directly, and
        implicitly takes native_ffprime[0] as its own axis value wherever
        it interpolates below the native grid's innermost point."""
        pressure_x = self._var('PMHD_IN')    # (time, nx), Pa, on X - total pressure, see docstring
        F0 = self._var('BZXR') / 100.0       # T*cm -> T*m, (time,)
        F_xb = self._var('GFUN') * F0[:, None]  # T*m, on XB - raw, noisy near axis, see docstring
        psi_xb = -self._var('PLFLX')         # Wb/rad, sign-flipped, on XB

        psi_norm_xb = np.empty_like(psi_xb)
        pprime_xb = np.empty_like(psi_xb)
        ffprime_xb = np.empty_like(psi_xb)
        for it in range(psi_xb.shape[0]):
            p_xb = np.interp(self.rho_b, self.rho, pressure_x[it])  # X -> XB
            psi = psi_xb[it]
            F = F_xb[it]

            pprime_xb[it] = _fix_axis_point(psi, _monotonic_gradient(p_xb, psi))
            ffp = F * _monotonic_gradient(F, psi)
            ffprime_xb[it] = _fix_axis_point(psi, ffp) if smooth_ffprime else ffp

            psi_lcfs = psi[-1]  # XB's outermost point is exactly rho=1.0
            psi_axis = 0.0      # see get_ffprime_pprime_profiles - PLFLX is referenced from the axis
            psi_norm_xb[it] = (psi - psi_axis) / (psi_lcfs - psi_axis)

        return psi_norm_xb, pprime_xb, ffprime_xb

    def get_reconstructed_Ip(self, ntheta=200, smooth_ffprime=True):
        """Plasma current Ip(t) [A] reconstructed purely from this module's
        own p'/ff' profiles (get_ffprime_pprime_profiles) via the
        Grad-Shafranov toroidal current relation

            j_tor(R, Z) = R*p'(psi) + ff'(psi) / (mu0*R)

        integrated over the full poloidal cross-section
        (get_flux_surface_grid) at each time, vs. the shot's own measured
        PCUR (get_Ip()) - {'Ip_reconstructed': {time: value},
        'Ip_measured': {time: value}, 'percent_diff': {time: value}}.

        This is a self-consistency check, not an independent measurement:
        it asks whether get_ffprime_pprime_profiles's derived p'/ff' and
        get_flux_surface_grid's moment-expansion geometry together
        reproduce the current TRANSP itself measured (PCUR), a much
        stronger end-to-end test than checking either profile's shape in
        isolation - useful precisely because it doesn't depend on ever
        running a full transp2imas/IMAS-based comparison.

        Method: p'(psi_norm)/ff'(psi_norm) are flux-surface quantities
        (from _ffprime_pprime_native, i.e. on their *native* XB-derived
        psi_norm grid, not resampled onto get_ffprime_pprime_profiles's
        uniform output grid - avoids a second, unnecessary resampling
        here), interpolated onto get_flux_surface_grid's rho grid via its
        psi_norm at each radial index; j_tor then varies over the
        cross-section only through the explicit R dependence above. The
        area element at each (rho, theta) grid point is
        |dR/drho * dZ/dtheta - dR/dtheta * dZ/drho| drho dtheta (the
        standard Jacobian determinant for the (rho, theta) -> (R, Z) map);
        theta derivatives use a periodic centered difference (np.roll,
        since theta is a uniform full-period grid), rho derivatives use
        np.gradient (non-uniform-spacing-aware, though rho here is close
        to uniform). Ip is then a Riemann sum over theta (exact for smooth
        periodic data at fixed ntheta) followed by np.trapezoid over rho.

        RESULT on this file (176858Z11.CDF, all 6 ramp-up slices): ~1-2%
        agreement with PCUR - consistent with the pre-existing few-percent
        floor from get_geometry's own boundary-moment-truncation
        cross-check (n=8, finite ntheta), i.e. no further unexplained
        discrepancy beyond that floor.

        This check is also what caught a real bug during development:
        an earlier version of get_ffprime_pprime_profiles used
        get_pressure_profile's thermal-only PMHDT_IN for p' (a reasonable-
        looking default, matching TORAX's pressure_thermal_total). That
        version reproduced PCUR to <1% only at the one time slice with
        ~zero fast-ion pressure, but under-shot by 8-24% - worse, the
        larger the fast-ion pressure fraction - at every other slice,
        while get_flux_surface_grid's own area cross-check stayed flat
        at ~2-5% throughout (i.e. not correlated with the error pattern,
        ruling geometry out as the cause) and GFUNC (TRANSP's internal
        Grad-Shafranov residual check) confirmed TRANSP's own equilibrium
        was self-consistent to <0.1% throughout (ruling out bad source
        data). That combination pointed squarely at the pressure variable
        choice - fixed by switching to PMHD_IN (total pressure), which is
        what actually enters the force balance TokaMaker's solve()
        enforces (see get_ffprime_pprime_profiles's docstring). If you
        change that variable again, rerun this check before trusting the
        result - a plausible-looking profile choice does not guarantee
        force-balance consistency.

        CAVEAT: this doesn't mean a mismatch always localizes this
        cleanly - a large or systematic mismatch against PCUR implicates
        something upstream (p'/ff' construction or get_flux_surface_grid's
        geometry), but not always which; cross-checking against a
        component that isolates one piece (as get_flux_surface_grid's own
        area check and GFUNC did here) is what makes it diagnosable rather
        than just "wrong somewhere."
        """
        psi_norm_xb, pprime_xb, ffprime_xb = self._ffprime_pprime_native(smooth_ffprime)

        Ip_reconstructed = np.empty(len(self.time))
        for it in range(len(self.time)):
            grid = self.get_flux_surface_grid(it, ntheta=ntheta)
            R, Z, rho, theta = grid['R'], grid['Z'], grid['rho'], grid['theta']
            dtheta = theta[1] - theta[0]

            order = np.argsort(psi_norm_xb[it])
            psi_norm_sorted = psi_norm_xb[it][order]
            psi_norm_full = np.concatenate([[0.0], psi_norm_xb[it]])  # rho=0 axis point, see get_flux_surface_grid
            pp_at_grid = np.interp(psi_norm_full, psi_norm_sorted, pprime_xb[it][order])
            ffp_at_grid = np.interp(psi_norm_full, psi_norm_sorted, ffprime_xb[it][order])

            dR_dtheta = (np.roll(R, -1, axis=1) - np.roll(R, 1, axis=1)) / (2 * dtheta)
            dZ_dtheta = (np.roll(Z, -1, axis=1) - np.roll(Z, 1, axis=1)) / (2 * dtheta)
            dR_drho = np.gradient(R, rho, axis=0)
            dZ_drho = np.gradient(Z, rho, axis=0)
            jacobian = np.abs(dR_drho * dZ_dtheta - dR_dtheta * dZ_drho)

            j_tor = R * pp_at_grid[:, None] + ffp_at_grid[:, None] / (self._MU0 * R)
            theta_integral = (j_tor * jacobian).sum(axis=1) * dtheta
            Ip_reconstructed[it] = np.trapezoid(theta_integral, rho)

        Ip_measured = self._var('PCUR')
        percent_diff = 100.0 * (Ip_reconstructed - Ip_measured) / Ip_measured

        return {
            'Ip_reconstructed': self.to_time_dict(Ip_reconstructed),
            'Ip_measured': self.to_time_dict(Ip_measured),
            'percent_diff': self.to_time_dict(percent_diff),
        }

    def get_conductivity(self, rho_grid):
        """Parallel electrical conductivity sigma_parallel [S/m],
        {time: {rho: value}} - 1/ETA_USE, where ETA_USE is TRANSP's
        "resistivity used or inferred" [Ohm*cm] (the resistivity actually
        driving the ohmic current in this run, as opposed to the several
        alternative Spitzer/neoclassical resistivity models this file also
        carries - ETA_SP/ETA_SPS/ETA_NC/ETA_TSC), for comparison against
        TORAX's own sigma_parallel profile.
        """
        eta_ohm_m = self.profile('ETA_USE', rho_grid) * 1e-2  # Ohm*cm -> Ohm*m
        sigma = 1.0 / eta_ohm_m
        return self.to_time_rho_dict(sigma, rho_grid)

    def get_particle_diffusivity(self, rho_grid):
        """TRANSP's interpretive electron particle diffusivity D_e [m^2/s],
        {time: {rho: value}} - DIFFE, for comparison against TORAX's own
        predicted D_neo_e/D_turb_e.

        CAVEAT (stronger than get_transport_coefficients's heat-diffusivity
        one): particle-balance-derived diffusivity is ill-conditioned
        wherever the interpretive particle flux is small (division by a
        near-zero flux/gradient), which is common - unlike CONDE/CONDI
        (heat), this file's DIFFE swings both very large and negative
        (roughly -5e6 to +4e7 cm^2/s at points in this shot). Treat this as
        a rough order-of-magnitude sanity check at best, not a shape
        comparison; TRANSP's DIFFNE ("total" variant) is present in this
        file's variable list but is uniformly zero here (unpopulated for
        this run), which is itself a sign of how unreliable this channel is
        for this discharge.
        """
        arr = self.profile('DIFFE', rho_grid, source_rho='XB') * 1e-4  # cm^2/s -> m^2/s
        return self.to_time_rho_dict(arr, rho_grid)

    # -- global scalars (for TORAX-vs-TRANSP scalar comparison) -------------

    def get_scalars(self):
        """Global 0-D time traces, for comparison against TORAX's own scalar
        outputs (`tmtx._data_tree.scalars`) - the TRANSP-side counterpart to
        `plot_scalars()`'s TokaMaker-vs-TORAX panels, used as TORAX-vs-TRANSP
        instead. Everything here is read directly off the CDF's own scalar
        (or per-zone, volume-integrated) diagnostics, not derived from any
        set_*() input, so none of it depends on modeling choices made
        elsewhere in this notebook.

        - q95/q0/li3: TRANSP's own Q95/Q0/LI_3 - same definitions TORAX's
          q95/q0/li3 post-processing targets, so directly comparable.
        - v_loop: VSURC ("calc. avg. surface voltage") - a different
          specific definition than TORAX's v_loop_lcfs (d(psi_lcfs)/dt), so
          treat this as a sanity check, not an exact match.
        - P_ohmic/P_radiation/P_alpha: volume integrals of POH/PRAD/
          (PALE+PALI) [W/cm^3] * DVOL [cm^3] -> W. Matches TORAX's
          P_ohmic_e/P_radiation_e/P_alpha_total in kind (both are
          volume-integrated powers) but not necessarily in exact accounting
          - eg. PRAD is TRANSP's "net radiated power used", which may bundle
          line/brem/cyclotron differently than TORAX's Mavrin-model
          radiation.
        - P_aux: PINJ (+ PEECH/PORICRH if this file has them) - total
          auxiliary heating power, for TORAX's P_aux_total.
        - tau_E/H98: TAUEA (thermal energy confinement time) and H98Y2 (the
          TauE98y,2 H-factor) - TRANSP's own, for TORAX's post-processed
          tau_E/H98.
        - psi_lcfs: PLFLXA, sign-flipped - see get_geometry's "Sign
          convention" docstring note (still only a sign fix, not a gauge/
          offset fix).
        - psi_axis: PSI0, sign-flipped for the same reason, only if this
          file actually populated it as a time trace - many TRANSP runs
          leave it an unused scalar placeholder (shape () instead of
          (ntime,)); None in that case rather than a fabricated per-time
          value.
        - beta_N: not returned - this file carries no direct normalized-beta
          variable (only poloidal-beta channels, a different definition), so
          there is no TRANSP-side value to compare TORAX's beta_N against.
        - pax: PMHD_IN on-axis value, same quantity/caveat as get_geometry's
          'pax' (includes fast-ion/rotation pressure).
        """
        dvol = self._var('DVOL')  # cm^3
        p_ohmic = (self._var('POH') * dvol).sum(axis=1)  # W/cm^3 -> W
        p_radiation = (self._var('PRAD') * dvol).sum(axis=1)
        p_alpha = ((self._var('PALE') + self._var('PALI')) * dvol).sum(axis=1)

        p_aux = self._var('PINJ').copy() if self.has('PINJ') else np.zeros_like(self.time)
        if self.has('PEECH'):
            p_aux = p_aux + (self._var('PEECH') * dvol).sum(axis=1)
        if self.has('PORICRH'):
            p_aux = p_aux + self._var('PORICRH')

        psi0_flipped = -self._var('PSI0') if self.has('PSI0') else None
        psi_axis = (self.to_time_dict(psi0_flipped)
                    if psi0_flipped is not None and psi0_flipped.shape == self.time.shape
                    else None)

        return {
            'Ip': self.get_Ip(),
            'q95': self.to_time_dict(self._var('Q95')),
            'q0': self.to_time_dict(self._var('Q0')),
            'li3': self.to_time_dict(self._var('LI_3')),
            'v_loop': self.to_time_dict(self._var('VSURC')),
            'P_ohmic': self.to_time_dict(p_ohmic),
            'P_radiation': self.to_time_dict(p_radiation),
            'P_alpha': self.to_time_dict(p_alpha),
            'P_aux': self.to_time_dict(p_aux),
            'tau_E': self.to_time_dict(self._var('TAUEA')),
            'H98': self.to_time_dict(self._var('H98Y2')),
            'psi_lcfs': self.to_time_dict(-self._var('PLFLXA')),
            'psi_axis': psi_axis,
            'pax': self.to_time_dict(self._var('PMHD_IN')[:, 0]),
        }

    # -- diagnostics --------------------------------------------------------

    def summary(self):
        """Print what this file actually contains, so nothing gets silently
        assumed. Run this before trusting any of the above."""
        print(f"time range: {self.time[0]:.4f} - {self.time[-1]:.4f} s "
              f"({len(self.time)} points)")
        print(f"rho grid: X has {len(self.rho)} zone centers "
              f"({self.rho[0]:.3f}-{self.rho[-1]:.3f}), "
              f"XB has {len(self.rho_b)} boundaries "
              f"({self.rho_b[0]:.3f}-{self.rho_b[-1]:.3f})")

        comp = self.get_plasma_composition()
        print(f"main ion fractions: {comp['main_ion_fractions']}")
        print(f"Zeff (time-avg, vol-avg): {comp['Zeff_time_avg']:.3f}")
        if comp['impurity']:
            imp = comp['impurity']
            print(f"impurity: measured avg Z={imp['measured_avg_Z']:.2f}, "
                  f"A={imp['measured_avg_A']:.2f} -> nearest Mavrin species "
                  f"'{imp['symbol']}'")

        ip = np.array(list(self.get_Ip().values()))
        print(f"Ip: {ip.min()/1e6:.3f} - {ip.max()/1e6:.3f} MA")

        geom = self.get_geometry()
        it = len(self.time) // 2
        t_mid = self.time[it]
        print(f"geometry @ t={t_mid:.3f}s: R0_geo={geom['R0_geo'][t_mid]:.3f} m, "
              f"a={geom['a'][t_mid]:.3f} m, F0={geom['F0'][t_mid]:.3f} T*m, "
              f"kappa={geom['kappa'][t_mid]:.2f} (up={geom['kappa_upper'][t_mid]:.2f}, "
              f"lo={geom['kappa_lower'][t_mid]:.2f}), "
              f"delta_up={geom['delta_upper'][t_mid]:.2f}, "
              f"delta_lo={geom['delta_lower'][t_mid]:.2f}, "
              f"pax={geom['pax'][t_mid]:.0f} Pa")

        has_ech = any(self.has(v) for v in ('ECHE', 'ECHI', 'PEECH'))
        has_ich = any(self.has(v) for v in ('PORICRH', 'PIICRH', 'PEICRH'))
        print(f"heating sources present: NBI={self.has('PINJ')}, "
              f"ECH={has_ech}, ICH={has_ich}")
        if self.has('PINJ'):
            pinj = self._var('PINJ')
            print(f"  PINJ: {pinj.min()/1e6:.2f} - {pinj.max()/1e6:.2f} MW")

        fueling = self.get_fueling()
        for name, result in fueling.items():
            if result is None:
                print(f"fueling[{name}]: not present in this file")
                continue
            if 'decay_length' in result:
                print(f"fueling[{name}]: decay_length_avg={result['decay_length_avg']:.3f}")
                continue
            s_total = np.array(list(result['S_total'].values()))
            print(f"fueling[{name}]: S_total {s_total.min():.2e} - "
                  f"{s_total.max():.2e} particles/s, "
                  f"loc_avg={result['location_avg']:.3f}, "
                  f"width_avg={result['width_avg']:.3f}")


if __name__ == '__main__':
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else '176858Z11.CDF'
    run = TranspRun(path)
    run.summary()
