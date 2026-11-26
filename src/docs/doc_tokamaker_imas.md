IMAS input and output for TokaMaker     {#doc_tokamaker_imas}
================

TokaMaker can write its equilibria to, and read machine descriptions, profiles, targets and
magnetic measurements from, [IMAS](https://imas-data-dictionary.readthedocs.io) Interface Data
Structures (IDSs) through [imas-python](https://imas-python.readthedocs.io). The interface is in
\ref OpenFUSIONToolkit.TokaMaker.imas_io "TokaMaker.imas_io".

imas-python is optional: install it with `pip install imas-python` (or the `imas` extra of the
OpenFUSIONToolkit package). It is only imported when the IMAS functions are used.

## Data Dictionary versions, time modes and conventions
Only Data Dictionary (DD) 4.x is supported. Every IDS read is checked, and older ones are rejected
(convert them first, eg. with `imas.convert_ids(ids, "4.1.1")`). New IDSs use imas-python's default
DD version unless `dd_version` is given; `save_imas` writes the entry with the version of its IDSs,
which must all match. DD 4.0 lacks `boundary.phi` and `boundary.rho_tor` (not written there) and
coil-level `pf_active.coil.geometry` (coils need elements there).

Time modes:
- Output IDSs use homogeneous time (`homogeneous_time = 1`). When currents are written into an
  existing IDS, a static one (eg. a machine description) becomes homogeneous, a homogeneous one must
  already use the same times (or none, otherwise it is rejected unchanged), and a heterogeneous one
  gets the times on the new currents only, so other signals keep their own times. The previous
  current signal (including its errors) is replaced.
- Writing equilibria into an existing `equilibrium` IDS keeps its slices (which must be at the same
  times: the root `time` for homogeneous IDSs, the slice times otherwise) and their reconstruction
  `constraints`; the rest of each slice is replaced. The new slices are computed first, so a
  failure leaves the IDS unchanged.
- Inputs accept homogeneous (1) or heterogeneous (0) time; slices and signals are taken at the time
  nearest the one requested, from the root `time` or their own `time` respectively. Machine
  descriptions (`wall`, coil and loop geometry) may be static (2).

DD 4 uses COCOS 17. TokaMaker
uses COCOS 7, which has the same signs but poloidal flux per radian, so:
- \f$ \psi \f$-like quantities (`psi`, `psi_axis`, `psi_boundary`, flux loop signals) are multiplied by \f$ 2\pi \f$
- \f$ d/d\psi \f$-like quantities (`dpressure_dpsi`, `f_df_dpsi`, `dvolume_dpsi`, `darea_dpsi`) are divided by \f$ 2\pi \f$
- everything else (\f$ I_p \f$, \f$ q \f$, \f$ F \f$, \f$ B \f$, currents) is unchanged

In DD 4, \f$ B_R = -\partial_Z\psi/(2\pi R) \f$ and \f$ B_Z = \partial_R\psi/(2\pi R) \f$, and for positive \f$ I_p \f$
the flux decreases from the magnetic axis.

Data are written and read as in-memory IDS objects, so any imas-python data entry can be used. A
path ending in `.nc` is a netCDF file, which needs no IMAS Access Layer:
~~~{.py}
from OpenFUSIONToolkit.TokaMaker import imas_io
eq_ids = imas_io.equilibrium_to_imas([eq1, eq2], [0.0, 1.0])   # One time slice per equilibrium
pf_ids = imas_io.pf_active_to_imas([eq1, eq2], [0.0, 1.0])
imas_io.save_imas('pulse.nc', [eq_ids, pf_ids])
# Or, for a single equilibrium
mygs.save_imas('equilibrium.nc', time=0.0)
~~~

## Output
`equilibrium_to_imas` fills, for each time slice:
- `vacuum_toroidal_field`: `r0` (default: geometric axis) and \f$ b_0 = F_0/r_0 \f$
- `global_quantities`: `ip`, `psi_axis`, `psi_boundary` (at the selected boundary), `magnetic_axis`,
  `current_centre`, `q_axis`, `q_95`, `q_min`, `volume` (whole plasma), `length_pol` (of the boundary),
  `energy_mhd`, `beta_pol`, `beta_tor`, `beta_tor_norm`, `li_3`, `rho_tor_boundary`
- `boundary`: `type`, `outline`, `psi`, `psi_norm`, `geometric_axis`, `minor_radius`, `elongation`,
  `triangularity` (upper and lower), `phi`, `rho_tor` (DD >= 4.1)
- `contour_tree.node`: the O-point and X-points
- `profiles_1d`: `psi`, `psi_norm`, `phi`, `rho_tor`, `rho_tor_norm`, `dpsi_drho_tor`, `f`, `f_df_dpsi`,
  `pressure`, `dpressure_dpsi`, `q`, `magnetic_shear`, `j_phi`, `j_parallel`, `volume`, `dvolume_dpsi`,
  `dvolume_drho_tor`, `rho_volume_norm`, `area`, `darea_dpsi`, `surface`, `gm1`-`gm9`, `b_field_min`,
  `b_field_max`, `trapped_fraction`, `r_inboard`, `r_outboard`, `geometric_axis`, `elongation`,
  `triangularity_upper`, `triangularity_lower`
- `profiles_2d[0]`: `psi` and `b_field_r/phi/z` on a rectangular grid spanning the plasma region

1D profiles are sampled uniformly in \f$ \sqrt{\hat{\psi}} \f$ from the axis to \f$ \hat{\psi} = 1 - \f$ `lcfs_pad`,
the surface also used for `boundary` (as in `save_eqdsk`); `profiles_1d.psi_norm` is normalized to the
LCFS, and `boundary.psi_norm` is \f$ 1 - \f$ `lcfs_pad`. Axis values that cannot be traced are filled
with their limits or by extrapolation. The flux-surface averages come from
\ref OpenFUSIONToolkit.TokaMaker.TokaMaker_equilibrium.get_fsa "get_fsa", so they share one quadrature.

Some definitions differ from \ref OpenFUSIONToolkit.TokaMaker.TokaMaker_equilibrium.get_stats "get_stats":
- `beta_pol` \f$ = 4 \int p\, dV / (R_0 \mu_0 I_p^2) \f$ (IMAS), not the circumference-based value of `get_stats`
- `beta_tor` is a fraction (not %), with \f$ B_0 \f$ at `r0`
- `li_3` is \f$ l_i \f$ with `li_normalization='iter'`
- `j_phi` is \f$ \langle j_\phi/R \rangle / \langle 1/R \rangle \f$ and `j_parallel` is
  \f$ \langle \mathbf{J}\cdot\mathbf{B} \rangle / B_0 \f$ (see \ref doc_tokamaker_current_conventions)

Not written: LCFS `area` and `surface` in `global_quantities` (only the padded boundary is traced),
squareness, gaps, inverse (flux-coordinate) grids, `coordinate_system`, and GGD representations. GGD
would represent TokaMaker's high-order solution only approximately: GGD `coefficients` have no
defined basis. Only tokamak equilibria (`psi_convention=0`, not mirror or dipole mode) with a vacuum
toroidal field (\f$ F_0 \neq 0 \f$) can be written. The vacuum field is read from the solved field
(the `F0` property is not updated by a reconstruction with `fitF0`). Equilibria with a fitted flux
scale (`psiscale`, from fixed-boundary reconstructions) are rejected: TokaMaker's field evaluators
include it but its flux bounds, profiles and averages (and `save_eqdsk`) do not. The check compares
flux differences, so it does not depend on the zero of the flux.

`pf_active_to_imas` writes the current per turn of each coil set, with one `element` per coil
region. `pf_passive_to_imas` writes the current in each conducting region. Either can fill an
existing IDS (eg. the machine description the mesh was built from), matching objects by name (exact,
then as `domain_from_imas` forms region names). A multi-element passive loop gets the mean of its
element currents per turn (`turns_with_sign`), with a warning if they differ. All equilibria
written together must come from the same TokaMaker object.

## Input
- \ref OpenFUSIONToolkit.TokaMaker.imas_io.domain_from_imas "domain_from_imas" builds a
  \ref OpenFUSIONToolkit.TokaMaker.meshing.gs_Domain "gs_Domain" from `wall` (limiter outline and
  vessel units), `pf_active` (one coil set per coil, one region per element with `turns_with_sign`
  turns) and, optionally, `pf_passive`. Element geometries of all six DD types are converted to
  polygons; an annulus with a hole gets a vacuum region `<name>_HOLE`. Resistivities must be given
  (or, for a single-element passive loop, the loop resistance, converted with
  \f$ R_{loop} = 2\pi\eta N^2 / \int dA/R \f$ for its `turns_with_sign` N). Mesh sizes are not part of IMAS and
  are given per region type. Circuits are not mapped. TokaMaker cannot constrain the elements of a
  multi-element passive loop to carry the same current, so they become separate conductors.
- \ref OpenFUSIONToolkit.TokaMaker.imas_io.profiles_from_imas "profiles_from_imas" returns, from an
  `equilibrium` IDS, `ffp_prof`/`pp_prof` for `set_profiles`, `Ip`/`pax` targets, \f$ F_0 \f$, and
  the boundary as isoflux points. From a `core_profiles` IDS it returns kinetic profiles for
  `set_kinetic_profiles` on normalized toroidal flux (\f$ \hat{\Phi} = \rho_{tor,norm}^2 \f$), with
  temperatures converted to keV. Densities are thermal (`density_thermal`, else `density` less
  `density_fast`). TokaMaker uses one effective ion species, so `ni_prof` is the total thermal ion
  density and `ti_prof` the thermal-density-weighted ion temperature (exact thermal ion pressure;
  `t_i_average` where the ion density vanishes).
  Without `profiles_1d.psi_norm`, the coordinate is recovered from `psi`, the global flux bounds and
  `boundary.psi_norm`.
- \ref OpenFUSIONToolkit.TokaMaker.imas_io.targets_from_pulse_schedule "targets_from_pulse_schedule"
  interpolates `position_control` (boundary outline, strike points, X-points, magnetic axis) and
  `flux_control.ip` references to a given time.
- \ref OpenFUSIONToolkit.TokaMaker.imas_io.load_magnetics "load_magnetics" adds reconstruction
  constraints from `magnetics` flux loops, poloidal field probes and \f$ I_p \f$, and from `pf_active`
  coil currents. Only full poloidal flux loops given by one complete (R,Z) position are used (saddle,
  diamagnetic, differential and multi-segment loops are skipped with a warning), as are probes with a position and poloidal angle (differential probes are
  skipped), samples flagged invalid by `validity`/`validity_timed` are skipped, every
  measurement needs a positive uncertainty, and the pf_active currents must cover all coil sets.
  Probe normals follow DD 4.1:
  \f$ n = (\cos\alpha\cos\beta, \cos\alpha\sin\beta, -\sin\alpha) \f$ in \f$ (R,\phi,Z) \f$.
  \ref OpenFUSIONToolkit.TokaMaker.imas_io.recon_constraints_to_imas "recon_constraints_to_imas"
  writes the measured and reconstructed values, uncertainties, \f$ \chi^2 \f$ and source paths to
  `equilibrium.time_slice[].constraints`, replacing the slice's constraints so they describe this
  fit only. Imported constraints go at the index of their source sensor or coil (skipped ones stay
  empty, and two from the same source are an error); others are appended. `chi_squared_reduced` is
  written only if the degrees of freedom are given.
- TokaMaker requires \f$ I_p > 0 \f$; inputs with a negative plasma current are rejected.

## TORAX
The equilibrium IDS contains the quantities read by TORAX's IMAS geometry loader (`geometry_type='imas'`,
TORAX 2.0), and is tested with it. For coupling TokaMaker to TORAX in one process, TORAX's `tokamaker`
geometry, which takes the output of `get_fsa` directly, needs no IMAS.
