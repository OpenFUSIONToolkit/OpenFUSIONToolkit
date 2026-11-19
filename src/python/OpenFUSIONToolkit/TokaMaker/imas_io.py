#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
r'''! IMAS (Data Dictionary 4.x) input and output for TokaMaker, using imas-python

Output: `equilibrium`, `pf_active` and `pf_passive` IDSs from TokaMaker equilibria, and the
`equilibrium.time_slice.constraints` of a reconstruction. Input: a TokaMaker mesh domain from a
machine description (`wall`, `pf_active`, `pf_passive`), source and kinetic profiles from
`equilibrium` or `core_profiles`, shape targets from `pulse_schedule`, and reconstruction
constraints from `magnetics` and `pf_active`.

Only Data Dictionary 4.x (COCOS 17) is supported; IDSs from older versions are rejected. TokaMaker
uses COCOS 7 (\f$ \psi \f$ per radian, same signs as COCOS 17), so \f$ \psi \f$-like quantities are
multiplied by \f$ 2\pi \f$ and \f$ d/d\psi \f$-like ones divided by \f$ 2\pi \f$; all other
quantities are unchanged. See @ref doc_tokamaker_imas for the fields read and written, supported
time modes and limitations.

imas-python is an optional dependency (`pip install imas-python`) and is imported on first use.

@authors Matthew Pharr
@date November 2026
@ingroup doxy_oft_python
'''
import warnings
import numpy
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from .._interface import mu0
from .meshing import gs_Domain

## @cond
_TWO_PI = 2.0*numpy.pi
## @endcond


def _imas():
    '''! Import imas-python, with an explanatory error if it is missing'''
    try:
        import imas
    except ImportError as err:
        raise ImportError('IMAS support requires imas-python (pip install imas-python)') from err
    return imas


def _dd_version(ids):
    '''! Data Dictionary version string of an IDS'''
    return _imas().util.get_data_dictionary_version(ids)


def check_dd4(ids):
    r'''! Check that an IDS uses Data Dictionary 4.x

    @param ids IDS object
    @result `ids`, unchanged
    '''
    version = _dd_version(ids)
    if int(version.split('.')[0]) != 4:
        raise ValueError('TokaMaker IMAS support requires Data Dictionary 4.x, got version {0} for "{1}"; '
                         'convert the data to DD 4 first (eg. `imas.convert_ids(ids, "4.1.1")`)'.format(
                             version, ids.metadata.name))
    return ids


def _has_node(struct, name):
    '''! Whether a structure has child `name` in its Data Dictionary version'''
    try:
        struct.metadata[name]
    except KeyError:
        return False
    return True


def _new_ids(name, dd_version=None):
    '''! Create an empty IDS with homogeneous time and TokaMaker code metadata'''
    imas = _imas()
    ids = imas.IDSFactory(dd_version).new(name)
    check_dd4(ids)
    ids.ids_properties.homogeneous_time = 1
    ids.ids_properties.comment = 'Produced by TokaMaker (Open FUSION Toolkit)'
    ids.code.name = 'TokaMaker'
    try:
        from .. import __version__
        ids.code.version = str(__version__)
    except ImportError:
        pass
    return ids


def save_imas(uri, ids_list, mode='w'):
    r'''! Write IDSs to an IMAS data entry

    A path ending in `.nc` writes a netCDF file, which needs no IMAS Access Layer; any other
    imas-python URI (eg. `imas:hdf5?path=...`) needs `imas_core`. The entry is opened with the Data
    Dictionary version of the IDSs, which must all use the same version.

    @param uri Data entry URI or `.nc` file path
    @param ids_list IDS object or list of IDS objects
    @param mode Data entry open mode ('w' to create or overwrite, 'a' to append)
    '''
    imas = _imas()
    ids_list = _as_list(ids_list)
    versions = sorted(set(_dd_version(check_dd4(ids)) for ids in ids_list))
    if len(versions) != 1:
        raise ValueError('All IDSs must use the same Data Dictionary version (got {0})'.format(', '.join(versions)))
    with imas.DBEntry(uri, mode, dd_version=versions[0]) as entry:
        for ids in ids_list:
            entry.put(ids)


def load_imas(uri, ids_name, occurrence=0):
    r'''! Read one IDS from an IMAS data entry (Data Dictionary 4.x only)

    @param uri Data entry URI or `.nc` file path
    @param ids_name Name of the IDS (eg. 'equilibrium')
    @param occurrence IDS occurrence
    @result IDS object
    '''
    imas = _imas()
    with imas.DBEntry(uri, 'r') as entry:
        ids = entry.get(ids_name, occurrence, autoconvert=False)
    return check_dd4(ids)


def _as_equil(obj):
    '''! TokaMaker_equilibrium for a TokaMaker or TokaMaker_equilibrium object'''
    equil = getattr(obj, '_tMaker_equil', None)
    if equil is not None:
        return equil
    if hasattr(obj, '_equil_ptr'):
        return obj
    raise TypeError('Expected a TokaMaker or TokaMaker_equilibrium object')


def _as_list(vals):
    '''! Wrap a single object in a list (lists and tuples are returned as lists)'''
    if isinstance(vals, (list, tuple)):
        return list(vals)
    return [vals]


def _equilibria_and_times(eqs, times):
    '''! Equilibrium objects and their times, checking they match in number and share one mesh'''
    eqs = [_as_equil(eq) for eq in _as_list(eqs)]
    times = numpy.atleast_1d(numpy.asarray(times, dtype=numpy.float64))
    if times.shape[0] != len(eqs):
        raise ValueError('Number of times ({0}) does not match number of equilibria ({1})'.format(times.shape[0], len(eqs)))
    if any(eq._tMaker is not eqs[0]._tMaker for eq in eqs):
        raise ValueError('All equilibria must come from the same TokaMaker object (mesh and coils)')
    return eqs, times


def _vacuum_F(equil):
    r'''! Vacuum \f$ F \f$ of the solved equilibrium, from the field outside the plasma

    This is the live value, which a reconstruction with `fitF0` changes without updating the
    `F0` property.
    '''
    tmaker = equil._tMaker
    outer = tmaker.r[numpy.argmax(tmaker.r[:, 0]), :2]
    return float(numpy.asarray(equil.get_field_eval('F').eval(outer)).reshape(-1)[0])


def _check_supported(equil):
    '''! Reject equilibria the IMAS writer does not support'''
    settings = equil._tMaker.settings
    if getattr(settings, 'mirror_mode', False) or getattr(settings, 'dipole_mode', False):
        raise ValueError('IMAS equilibrium output supports tokamak equilibria only (not mirror or dipole mode)')
    if equil.psi_convention != 0:
        raise ValueError('IMAS equilibrium output supports tokamak equilibria only (psi_convention=0)')
    if equil.psi_bounds[0] < -1.E98:
        raise ValueError('IMAS equilibrium output needs an equilibrium with plasma (no LCFS found)')
    # A reconstruction can fit an overall flux scale (psiscale), which TokaMaker's field evaluators
    # include but its flux bounds, profiles and flux-surface averages do not
    # Compared on the flux difference between the O-point and the LCFS point, which does not depend on
    # the arbitrary zero of the flux
    psi_lcfs, psi_axis = equil.psi_bounds
    psi_eval = equil.get_field_eval('PSI')
    dpsi_field = float(numpy.asarray(psi_eval.eval(numpy.asarray(equil.o_point))).reshape(-1)[0]) \
        - float(numpy.asarray(psi_eval.eval(numpy.asarray(equil.lim_point))).reshape(-1)[0])
    scale = dpsi_field/(psi_axis-psi_lcfs)
    if abs(scale-1.0) > 1.E-6:
        raise ValueError('IMAS equilibrium output does not support a fitted flux scale (psiscale = {0:.6g}); '
                         'TokaMaker\'s profiles and averages would not match its field'.format(scale))
    if _vacuum_F(equil) == 0.0:
        raise ValueError('IMAS equilibrium output needs a vacuum toroidal field (F0 != 0)')


def _extrapolate_axis(x, y):
    '''! Quadratic extrapolation of y(x) to x=0 from its first three points'''
    good = numpy.isfinite(y)
    if numpy.count_nonzero(good) < 3:
        return numpy.nan
    xs, ys = x[good][:3], y[good][:3]
    return numpy.polyval(numpy.polyfit(xs, ys, 2), 0.0)


def _outline_shape(pts):
    '''! Geometric axis, minor radius, elongation and triangularities of a closed (R,Z) outline'''
    rmin, rmax = pts[:, 0].min(), pts[:, 0].max()
    a = (rmax-rmin)/2.0
    r_geo = (rmax+rmin)/2.0
    i_top, i_bot = numpy.argmax(pts[:, 1]), numpy.argmin(pts[:, 1])
    return {
        'r_geo': r_geo,
        'z_geo': (pts[:, 1].max()+pts[:, 1].min())/2.0,
        'a': a,
        'kappa': (pts[i_top, 1]-pts[i_bot, 1])/(2.0*a),
        'delta_upper': (r_geo-pts[i_top, 0])/a,
        'delta_lower': (r_geo-pts[i_bot, 0])/a,
    }


def _midplane_radii(equil, psi_n):
    '''! Inboard and outboard major radius of each surface at the height of the magnetic axis'''
    psi_eval = equil.get_field_eval('PSI')
    psi_lcfs, psi_axis = equil.psi_bounds
    r_axis, z_axis = equil.o_point
    tmaker = equil._tMaker
    r_lo, r_hi = tmaker.r[:, 0].min(), tmaker.r[:, 0].max()

    def psi_norm(r):
        '''! Normalized flux along the axis height'''
        pts = numpy.column_stack((numpy.atleast_1d(r), z_axis*numpy.ones_like(numpy.atleast_1d(r))))
        vals = numpy.asarray(psi_eval.eval(pts)).reshape(-1)
        return (vals-psi_axis)/(psi_lcfs-psi_axis)

    out = []
    for r_end in (max(r_lo, 1.E-3), r_hi):
        rs = numpy.linspace(r_axis, r_end, 2001)
        vals = psi_norm(rs)
        # Keep the monotonic part from the axis to the first point at or past the LCFS
        stop = numpy.argmax(vals >= 1.0) if numpy.any(vals >= 1.0) else len(vals)-1
        rs, vals = rs[:stop+1], vals[:stop+1]
        radii = numpy.full(psi_n.shape, numpy.nan)
        for i, target in enumerate(psi_n):
            if target <= 0.0:
                radii[i] = r_axis
                continue
            j = numpy.searchsorted(vals, target)
            if j <= 0 or j >= len(vals):
                continue
            radii[i] = brentq(lambda r, t=target: psi_norm(r)[0]-t, rs[j-1], rs[j])
        out.append(radii)
    return out[0], out[1]


def _fill_equilibrium_slice(ts, equil, b0, r0, npsi, nr, nz, lcfs_pad):
    '''! Fill one `equilibrium.time_slice` from a TokaMaker equilibrium'''
    psi_lcfs, psi_axis = equil.psi_bounds
    r_axis, z_axis = equil.o_point
    Ip, centroid, vol, pvol, _, _, _ = equil.get_globals()
    stats = equil.get_stats(li_normalization='iter')
    psi_b = 1.0-lcfs_pad
    # Uniform in sqrt(psi_N) (~rho), which resolves the axis region for transport codes
    psi_n = numpy.linspace(0.0, numpy.sqrt(psi_b), npsi)**2
    dpsi_dpsin = psi_lcfs-psi_axis  # d(psi_TM)/d(psi_N)
    # Flux-surface averages (the axis is filled by its limit or by extrapolation)
    fsa = equil.get_fsa(psi=psi_n[1:])
    _, f, fp, p, pp = equil.get_profiles(psi=psi_n)
    _, fc, _, _ = equil.calc_sauter_fc(psi=psi_n[1:])

    def with_axis(key, axis_val=None):
        '''! get_fsa profile with the axis value prepended (given, or extrapolated)'''
        vals = numpy.asarray(fsa[key], dtype=numpy.float64)
        if axis_val is None:
            axis_val = _extrapolate_axis(psi_n[1:], vals)
        return numpy.r_[axis_val, vals]

    q = with_axis('q')
    r_avg = with_axis('<R>', r_axis)
    inv_r = with_axis('<1/R>', 1.0/r_axis)
    inv_r2 = with_axis('<1/R^2>', 1.0/r_axis**2)
    dvdpsi = with_axis('dV/dPsi')
    gpsi = with_axis('<|grad psi|>', 0.0)
    gpsi2 = with_axis('<|grad psi|^2>', 0.0)
    bp2 = with_axis('<Bp^2>', 0.0)
    inv_b2 = with_axis('<1/B^2>', r_axis**2/f[0]**2)
    gpsi2_b2 = with_axis('<|grad psi|^2/B^2>', 0.0)
    area = with_axis('area', 0.0)
    surface = with_axis('surface', 0.0)
    b_axis = abs(f[0])/r_axis
    b_min = with_axis('B_min', b_axis)
    b_max = with_axis('B_max', b_axis)
    b2 = bp2+f**2*inv_r2
    # Toroidal flux and volume by spline integration in psi_N
    dpsi_full = _TWO_PI*dpsi_dpsin  # d(psi_IMAS)/d(psi_N)
    phi = numpy.sign(b0)*numpy.abs(CubicSpline(psi_n, q*dpsi_full).antiderivative()(psi_n))
    volume = CubicSpline(psi_n, dvdpsi*dpsi_dpsin).antiderivative()(psi_n)
    rho_tor = numpy.sqrt(numpy.abs(phi)/(numpy.pi*abs(b0)))
    # d(rho)/d(psi_N) from d(rho^2)/d(psi_N) = |q dpsi_IMAS/dpsi_N|/(pi |B0|), finite at the axis
    drho2 = numpy.abs(q*dpsi_full)/(numpy.pi*abs(b0))
    with numpy.errstate(divide='ignore', invalid='ignore'):
        drho_dpsi = drho2/(2.0*rho_tor)/dpsi_dpsin  # d(rho)/d(psi_TM)
    drho_dpsi[0] = numpy.nan
    gm2 = drho_dpsi**2*bp2
    gm3 = drho_dpsi**2*gpsi2
    gm6 = drho_dpsi**2*gpsi2_b2
    gm7 = numpy.abs(drho_dpsi)*gpsi
    for vals in (gm2, gm3, gm6, gm7):
        vals[0] = _extrapolate_axis(psi_n[1:], vals[1:])
    # Per-surface shape
    r_min = with_axis('R_min', r_axis)
    r_max = with_axis('R_max', r_axis)
    z_min = with_axis('Z_min', z_axis)
    z_max = with_axis('Z_max', z_axis)
    r_zmin = with_axis('R_at_Zmin', r_axis)
    r_zmax = with_axis('R_at_Zmax', r_axis)
    with numpy.errstate(divide='ignore', invalid='ignore'):
        a_surf = (r_max-r_min)/2.0
        r_geo_surf = (r_max+r_min)/2.0
        kappa = (z_max-z_min)/(r_max-r_min)
        delta_u = (r_geo_surf-r_zmax)/a_surf
        delta_l = (r_geo_surf-r_zmin)/a_surf
    for vals in (kappa, delta_u, delta_l):
        vals[0] = _extrapolate_axis(psi_n[1:], vals[1:])
    r_in, r_out = _midplane_radii(equil, psi_n)
    # Current densities: IMAS j_phi = <j_phi/R>/<1/R> and j_parallel = <J.B>/B0
    j_phi = (pp+inv_r2*f*fp/mu0)/inv_r
    j_par = (f*pp+fp*b2/mu0)/b0

    # Boundary: the surface at psi_N = psi_b (the separatrix itself if lcfs_pad = 0)
    outline = equil.trace_surf(psi_b)
    shape = _outline_shape(outline)
    dl = numpy.linalg.norm(numpy.diff(numpy.vstack((outline, outline[:1])), axis=0), axis=1).sum()
    psi_bnd = _TWO_PI*(psi_axis+psi_b*dpsi_dpsin)
    bnd = ts.boundary
    bnd.type = 1 if equil.diverted else 0
    bnd.outline.r = outline[:, 0]
    bnd.outline.z = outline[:, 1]
    bnd.psi_norm = psi_b
    bnd.psi = psi_bnd
    bnd.geometric_axis.r = shape['r_geo']
    bnd.geometric_axis.z = shape['z_geo']
    bnd.minor_radius = shape['a']
    bnd.elongation = shape['kappa']
    bnd.triangularity_upper = shape['delta_upper']
    bnd.triangularity_lower = shape['delta_lower']
    bnd.triangularity = (shape['delta_upper']+shape['delta_lower'])/2.0
    if _has_node(bnd, 'phi'):  # DD >= 4.1
        bnd.phi = phi[-1]
        bnd.rho_tor = rho_tor[-1]

    # Global quantities; `psi_boundary` is at the selected boundary above, `volume` is the whole
    # plasma (inside the LCFS). LCFS area and surface are not written as only the padded surface is traced.
    gq = ts.global_quantities
    gq.ip = Ip
    gq.psi_axis = _TWO_PI*psi_axis
    gq.psi_boundary = psi_bnd
    gq.magnetic_axis.r = r_axis
    gq.magnetic_axis.z = z_axis
    gq.magnetic_axis.b_field_phi = f[0]/r_axis
    gq.current_centre.r = centroid[0]
    gq.current_centre.z = centroid[1]
    gq.q_axis = q[0]
    gq.q_95 = numpy.interp(0.95, psi_n, q) if psi_b >= 0.95 else stats['q_95']
    i_qmin = numpy.nanargmin(numpy.abs(q))
    gq.q_min.value = q[i_qmin]
    gq.q_min.psi_norm = psi_n[i_qmin]
    gq.q_min.psi = _TWO_PI*(psi_axis+psi_n[i_qmin]*dpsi_dpsin)
    gq.q_min.rho_tor_norm = rho_tor[i_qmin]/rho_tor[-1]
    gq.volume = vol
    gq.length_pol = dl
    gq.energy_mhd = 1.5*pvol
    gq.beta_pol = 4.0*pvol/(r0*mu0*Ip**2)
    gq.beta_tor = 2.0*mu0*pvol/vol/b0**2
    gq.beta_tor_norm = 100.0*gq.beta_tor*shape['a']*abs(b0)/(abs(Ip)/1.E6)
    gq.li_3 = stats['l_i']
    gq.rho_tor_boundary = rho_tor[-1]

    # Critical points of the flux map: O-point and X-points
    xpts, _ = equil.get_xpoints()
    nodes = [(2 if psi_axis > psi_lcfs else 0, r_axis, z_axis, psi_axis)]
    if xpts is not None:
        psi_eval = equil.get_field_eval('PSI')
        for pt in xpts:
            nodes.append((1, pt[0], pt[1], float(numpy.asarray(psi_eval.eval(pt)).reshape(-1)[0])))
    ts.contour_tree.node.resize(len(nodes))
    for node, (ctype, rn, zn, psin_abs) in zip(ts.contour_tree.node, nodes):
        node.critical_type = ctype
        node.r = rn
        node.z = zn
        node.psi = _TWO_PI*psin_abs

    # 1D profiles (psi_norm is normalized to the LCFS, as in the Data Dictionary)
    p1 = ts.profiles_1d
    p1.psi = _TWO_PI*(psi_axis+psi_n*dpsi_dpsin)
    p1.psi_norm = psi_n
    p1.phi = phi
    p1.rho_tor = rho_tor
    p1.rho_tor_norm = rho_tor/rho_tor[-1]
    with numpy.errstate(divide='ignore', invalid='ignore'):
        dpsi_drho = _TWO_PI/drho_dpsi  # d(psi_IMAS)/d(rho) = 2 pi / (d rho/d psi_TM)
    dpsi_drho[0] = 0.0
    p1.dpsi_drho_tor = dpsi_drho
    p1.f = f
    p1.f_df_dpsi = f*fp/_TWO_PI
    p1.pressure = p
    p1.dpressure_dpsi = pp/_TWO_PI
    p1.q = q
    with numpy.errstate(divide='ignore', invalid='ignore'):
        shear = rho_tor/q*numpy.gradient(q, rho_tor)
    shear[0] = 0.0
    p1.magnetic_shear = shear
    p1.j_phi = j_phi
    p1.j_parallel = j_par
    p1.volume = volume
    p1.dvolume_dpsi = dvdpsi/_TWO_PI
    p1.rho_volume_norm = numpy.sqrt(volume/volume[-1])
    p1.dvolume_drho_tor = dvdpsi*dpsi_drho/_TWO_PI
    p1.area = area
    p1.darea_dpsi = CubicSpline(psi_n, area).derivative()(psi_n)/dpsi_full
    p1.surface = surface
    p1.gm1 = inv_r2
    p1.gm2 = gm2
    p1.gm3 = gm3
    p1.gm4 = inv_b2
    p1.gm5 = b2
    p1.gm6 = gm6
    p1.gm7 = gm7
    p1.gm8 = r_avg
    p1.gm9 = inv_r
    p1.b_field_min = b_min
    p1.b_field_max = b_max
    p1.trapped_fraction = numpy.r_[0.0, 1.0-fc]
    p1.r_inboard = r_in
    p1.r_outboard = r_out
    p1.geometric_axis.r = r_geo_surf
    p1.geometric_axis.z = (z_max+z_min)/2.0
    p1.elongation = kappa
    p1.triangularity_upper = delta_u
    p1.triangularity_lower = delta_l

    # 2D profiles on a rectangular grid spanning the plasma region
    tmaker = equil._tMaker
    plasma_pts = tmaker.r[tmaker.lc[tmaker.reg == 1].reshape(-1), :2]
    rg = numpy.linspace(plasma_pts[:, 0].min(), plasma_pts[:, 0].max(), nr)
    zg = numpy.linspace(plasma_pts[:, 1].min(), plasma_pts[:, 1].max(), nz)
    r2d, z2d = numpy.meshgrid(rg, zg, indexing='ij')
    pts = numpy.column_stack((r2d.reshape(-1), z2d.reshape(-1)))
    psi2d = numpy.asarray(equil.get_field_eval('PSI').eval(pts)).reshape(nr, nz)
    b2d = numpy.asarray(equil.get_field_eval('B').eval(pts)).reshape(nr, nz, 3)
    ts.profiles_2d.resize(1)
    p2 = ts.profiles_2d[0]
    p2.type.index = 0
    p2.type.name = 'total'
    p2.grid_type.index = 1
    p2.grid_type.name = 'rectangular'
    p2.grid_type.description = 'Cylindrical R,Z ala eqdsk'
    p2.grid.dim1 = rg
    p2.grid.dim2 = zg
    p2.r = r2d
    p2.z = z2d
    p2.psi = _TWO_PI*psi2d
    p2.b_field_r = b2d[:, :, 0]
    p2.b_field_phi = b2d[:, :, 1]
    p2.b_field_z = b2d[:, :, 2]


def equilibrium_to_imas(eqs, times, ids=None, npsi=129, nr=129, nz=129, lcfs_pad=0.01, r0=None, dd_version=None):
    r'''! Write TokaMaker equilibria to an `equilibrium` IDS (Data Dictionary 4.x, COCOS 17)

    One time slice is written per equilibrium (tokamak equilibria with a vacuum field only). The
    boundary is the surface at \f$ \hat{\psi} = 1 - \f$ `lcfs_pad` (as in `save_eqdsk`), and 1D
    profiles are sampled uniformly in \f$ \sqrt{\hat{\psi}} \f$ from the axis to that surface. See
    @ref doc_tokamaker_imas for the fields written and definitions that differ from `get_stats`.

    @param eqs TokaMaker or TokaMaker_equilibrium object, or a list of them (same TokaMaker object)
    @param times Time of each equilibrium [s] (scalar or list)
    @param ids Existing `equilibrium` IDS to fill (default: create one)
    @param npsi Number of 1D profile points (at least 4)
    @param nr Number of radial points in the 2D grid
    @param nz Number of vertical points in the 2D grid
    @param lcfs_pad Padding in normalized flux at the LCFS for 1D profiles and the boundary
    @param r0 Reference major radius for `vacuum_toroidal_field` (default: geometric axis of the first equilibrium)
    @param dd_version Data Dictionary version for a new IDS (default: imas-python default, must be 4.x)
    @result `equilibrium` IDS
    '''
    if npsi < 4:
        raise ValueError('npsi must be at least 4 (the axis values are extrapolated from the first three surfaces)')
    if (nr < 2) or (nz < 2):
        raise ValueError('nr and nz must be at least 2')
    eqs, times = _equilibria_and_times(eqs, times)
    for eq in eqs:
        _check_supported(eq)
    target = None
    if ids is not None:
        # Slices of an existing IDS are kept, with their reconstruction constraints; the rest of each slice
        # is replaced. Everything is computed in a new IDS first, so a failure leaves `ids` unchanged.
        check_dd4(ids)
        if len(ids.time_slice) > 0:
            mode = int(ids.ids_properties.homogeneous_time.value) if ids.ids_properties.homogeneous_time.has_value else 0
            if mode == 1:
                old_times = numpy.asarray(ids.time, dtype=numpy.float64)
            else:
                old_times = numpy.array([float(ts.time) if ts.time.has_value else numpy.nan for ts in ids.time_slice])
            if (len(old_times) != len(times)) or not numpy.allclose(old_times, times, rtol=0.0, atol=1.E-12):
                raise ValueError('equilibrium IDS already has time slices at other times; write to an empty IDS '
                                 'or one with matching slices')
        target = ids
        dd_version = _dd_version(ids)
    ids = _new_ids('equilibrium', dd_version)
    if r0 is None:
        r0 = eqs[0].get_stats()['R_geo']
    b0 = numpy.array([_vacuum_F(eq)/r0 for eq in eqs])
    ids.time = times
    ids.vacuum_toroidal_field.r0 = r0
    ids.vacuum_toroidal_field.b0 = b0
    ids.time_slice.resize(len(eqs))
    for i, eq in enumerate(eqs):
        ts = ids.time_slice[i]
        ts.time = times[i]
        _fill_equilibrium_slice(ts, eq, b0[i], r0, npsi, nr, nz, lcfs_pad)
    if target is None:
        return ids
    # Transfer into the existing IDS (assignments only)
    target.ids_properties.homogeneous_time = 1
    target.time = times
    target.vacuum_toroidal_field.r0 = r0
    target.vacuum_toroidal_field.b0 = b0
    target.time_slice.resize(len(eqs), keep=True)
    for old, new in zip(target.time_slice, ids.time_slice):
        for child in list(old.iter_nonempty_()):
            if child.metadata.name != 'constraints':
                _clear_node(old, child)
        for child in new.iter_nonempty_():
            if child.metadata.name != 'constraints':
                _copy_node(old, child)
    return target


def _copy_node(dst_parent, src):
    '''! Copy a filled IDS node (recursively) into the child of the same name of `dst_parent`'''
    imas = _imas()
    name = src.metadata.name
    dst = getattr(dst_parent, name)
    if isinstance(src, imas.ids_struct_array.IDSStructArray):
        dst.resize(len(src))
        for dst_item, src_item in zip(dst, src):
            for child in src_item.iter_nonempty_():
                _copy_node(dst_item, child)
    elif isinstance(src, imas.ids_structure.IDSStructure):
        for child in src.iter_nonempty_():
            _copy_node(dst, child)
    else:
        setattr(dst_parent, name, src.value)


def _set_times(ids, times):
    r'''! Prepare an existing IDS for signals at `times`, without retiming its other signals

    @result True if the IDS uses heterogeneous time (each signal carries its own `time`)
    '''
    mode = int(ids.ids_properties.homogeneous_time.value) if ids.ids_properties.homogeneous_time.has_value else 2
    if mode == 0:
        return True
    if mode == 1:
        if (len(ids.time) > 0) and ((len(ids.time) != len(times)) or not numpy.allclose(ids.time, times, rtol=0.0, atol=1.E-12)):
            raise ValueError('"{0}" already has a different time base; write to an IDS with matching times, an empty '
                             'or static IDS, or one with heterogeneous time'.format(ids.metadata.name))
    else:
        # Static machine description: becomes homogeneous in time
        ids.ids_properties.homogeneous_time = 1
    ids.time = times
    return False


def _match_coil_set(name, index, coil_sets):
    '''! TokaMaker coil set for an IMAS coil: exact name first, then the name as @ref domain_from_imas forms it'''
    name = str(name)
    if name in coil_sets:
        return name
    key = _region_name(name, 'COIL', index)
    return key if key in coil_sets else None


def pf_active_to_imas(eqs, times, ids=None, dd_version=None):
    r'''! Write coil currents to a `pf_active` IDS

    Each TokaMaker coil set is one `coil` (current per turn, as in IMAS), with one `element` per
    coil region. If `ids` is given (eg. the machine description the mesh was built from with
    @ref domain_from_imas), coils are matched by name and only their currents are written: a static
    IDS becomes homogeneous in time, a homogeneous one must already use `times` (or none), and a
    heterogeneous one gets the times on each current.

    @param eqs TokaMaker or TokaMaker_equilibrium object, or a list of them (same TokaMaker object)
    @param times Time of each equilibrium [s] (scalar or list)
    @param ids Existing `pf_active` IDS to fill (default: create one)
    @param dd_version Data Dictionary version for a new IDS (default: imas-python default, must be 4.x)
    @result `pf_active` IDS
    '''
    eqs, times = _equilibria_and_times(eqs, times)
    tmaker = eqs[0]._tMaker
    currents = [eq.get_coil_currents()[0] for eq in eqs]
    names = sorted(tmaker.coil_sets, key=lambda name: tmaker.coil_sets[name]['id'])
    if ids is None:
        ids = _new_ids('pf_active', dd_version)
        ids.coil.resize(len(names))
        for coil, name in zip(ids.coil, names):
            coil.name = name
            sub_coils = tmaker.coil_sets[name]['sub_coils']
            coil.element.resize(len(sub_coils))
            for element, sub_coil in zip(coil.element, sub_coils):
                element.name = sub_coil['name']
                element.turns_with_sign = sub_coil.get('nturns', 1.0)
        coil_map = {name: coil for name, coil in zip(names, ids.coil)}
    else:
        check_dd4(ids)
        coil_map = {}
        for i, coil in enumerate(ids.coil):
            key = _match_coil_set(coil.name, i, tmaker.coil_sets)
            if key is not None:
                coil_map[key] = coil
        missing = [name for name in names if name not in coil_map]
        if missing:
            raise KeyError('Coil sets not found in pf_active: {0}'.format(', '.join(missing)))
    hetero = _set_times(ids, times)
    for name in names:
        coil = coil_map[name]
        _clear_struct(coil.current)  # Including the errors of a previous signal
        coil.current.data = numpy.array([curr[name] for curr in currents])
        if hetero:
            coil.current.time = times
    return ids


def _loop_current(element_currents, turns, name):
    r'''! Loop current from the currents of its element regions, each carrying `turns` times the loop current

    @param element_currents Region currents [nelem, ntimes]
    @param turns `turns_with_sign` of each element [nelem]
    @param name Loop name (for warnings)
    @result Loop current [ntimes]
    '''
    per_turn = numpy.asarray(element_currents, dtype=numpy.float64)/numpy.asarray(turns, dtype=numpy.float64)[:, None]
    current = per_turn.mean(axis=0)
    spread = numpy.max(numpy.abs(per_turn-current[None, :])) if per_turn.shape[0] > 1 else 0.0
    if spread > 1.E-3*max(numpy.max(numpy.abs(current)), 1.E-12):
        warnings.warn('Elements of pf_passive loop "{0}" carry different currents per turn (TokaMaker treats them as '
                      'separate conductors); writing their mean'.format(name))
    return current


def pf_passive_to_imas(eqs, times, psis=None, ids=None, dd_version=None):
    r'''! Write the induced current in each conducting region to a `pf_passive` IDS

    Each TokaMaker conducting region is one `loop`. If `ids` is given, loops are matched by name
    (time handling as in @ref pf_active_to_imas); the current of a loop is the mean of its element
    currents per turn (`turns_with_sign`), with a warning if they differ.

    @param eqs TokaMaker or TokaMaker_equilibrium object, or a list of them (same TokaMaker object)
    @param times Time of each equilibrium [s] (scalar or list)
    @param psis \f$ \psi \f$ [np] with the conductor currents for each time (default: each equilibrium's own \f$ \psi \f$)
    @param ids Existing `pf_passive` IDS to fill (default: create one)
    @param dd_version Data Dictionary version for a new IDS (default: imas-python default, must be 4.x)
    @result `pf_passive` IDS
    '''
    eqs, times = _equilibria_and_times(eqs, times)
    tmaker = eqs[0]._tMaker
    if psis is None:
        psis = [eq.get_psi(normalized=False) for eq in eqs]
    psis = _as_list(psis)
    if len(psis) != len(eqs):
        raise ValueError('Number of psi arrays ({0}) does not match number of equilibria ({1})'.format(len(psis), len(eqs)))
    psis = [numpy.ascontiguousarray(psi, dtype=numpy.float64) for psi in psis]
    for psi in psis:
        if psi.shape != (tmaker.np,):
            raise ValueError('psi must have shape ({0},), got {1}'.format(tmaker.np, psi.shape))
    # Tessellated cell areas
    tri = tmaker.r[tmaker.lc, :2]
    cell_area = 0.5*numpy.abs((tri[:, 1, 0]-tri[:, 0, 0])*(tri[:, 2, 1]-tri[:, 0, 1])
                              - (tri[:, 2, 0]-tri[:, 0, 0])*(tri[:, 1, 1]-tri[:, 0, 1]))
    names = [name for name, cond in tmaker._cond_dict.items() if cond.get('eta', -1.0) > 0.0]
    masks = {name: (tmaker.reg == tmaker._cond_dict[name]['reg_id']) for name in names}
    region_currents = {name: numpy.zeros((len(eqs),)) for name in names}
    for k, (eq, psi) in enumerate(zip(eqs, psis)):
        _, j_cell = eq.calc_conductor_currents(psi, cell_centered=True)
        for name in names:
            region_currents[name][k] = numpy.sum(j_cell[masks[name]]*cell_area[masks[name]])
    if ids is None:
        ids = _new_ids('pf_passive', dd_version)
        ids.loop.resize(len(names))
        for loop, name in zip(ids.loop, names):
            loop.name = name
            loop.resistivity = tmaker._cond_dict[name]['eta']
            loop.current = region_currents[name]
        ids.time = times
    else:
        check_dd4(ids)
        hetero = _set_times(ids, times)
        for i, loop in enumerate(ids.loop):
            nelem = len(loop.element)
            # Exact names first (eg. native TokaMaker conductors), then the names domain_from_imas forms
            exact = [str(loop.name) if nelem <= 1 else '{0}_E{1}'.format(loop.name, j+1) for j in range(max(nelem, 1))]
            formed = [_region_name(str(loop.name), 'LOOP', i, j if nelem > 1 else None) for j in range(max(nelem, 1))]
            keys = exact if any(key in region_currents for key in exact) else formed
            found = [j for j, key in enumerate(keys) if key in region_currents]
            if not found:
                continue
            if len(found) != len(keys):
                raise KeyError('pf_passive loop "{0}" has elements without a TokaMaker region'.format(loop.name))
            turns = [float(loop.element[j].turns_with_sign) if (nelem > 0 and loop.element[j].turns_with_sign.has_value)
                     else 1.0 for j in found]
            for field in ('current', 'current_error_upper', 'current_error_lower', 'time'):
                if _has_node(loop, field):
                    setattr(loop, field, numpy.zeros((0,)))
            loop.current = _loop_current([region_currents[keys[j]] for j in found], turns, loop.name)
            if hetero:
                loop.time = times
    return ids


#------------------------------------------------------------------------------
# Input: machine description
#------------------------------------------------------------------------------
def _region_name(name, prefix, index, element=None):
    '''! TokaMaker region name for an IMAS object (upper case, unique per object and element)'''
    name = str(name).strip().upper().replace(' ', '_')
    if (name == '') or name.startswith('#'):
        name = '{0}{1}'.format(prefix, index+1)
    if element is not None:
        name = '{0}_E{1}'.format(name, element+1)
    return name


def _circle(rc, zc, radius, npts=32):
    '''! Polygon (N,2) approximating a circle'''
    theta = numpy.linspace(0.0, _TWO_PI, npts, endpoint=False)
    return numpy.column_stack((rc+radius*numpy.cos(theta), zc+radius*numpy.sin(theta)))


def _arcs_polygon(rs, zs, radii, npts=8):
    '''! Polygon through points joined by circular arcs (zero radius for straight segments)'''
    out = []
    n = len(rs)
    for i in range(n):
        p1 = numpy.array([rs[i], zs[i]])
        p2 = numpy.array([rs[(i+1) % n], zs[(i+1) % n]])
        out.append(p1)
        rad = radii[i] if i < len(radii) else 0.0
        chord = numpy.linalg.norm(p2-p1)
        if (rad == 0.0) or (not numpy.isfinite(rad)) or (abs(rad) < chord/2.0):
            continue
        mid = (p1+p2)/2.0
        normal = numpy.array([-(p2-p1)[1], (p2-p1)[0]])/chord
        center = mid-numpy.sign(rad)*normal*numpy.sqrt(rad**2-(chord/2.0)**2)
        a1 = numpy.arctan2(*(p1-center)[::-1])
        a2 = numpy.arctan2(*(p2-center)[::-1])
        da = numpy.angle(numpy.exp(1j*(a2-a1)))
        for t in numpy.linspace(0.0, 1.0, npts+1)[1:-1]:
            ang = a1+t*da
            out.append(center+abs(rad)*numpy.array([numpy.cos(ang), numpy.sin(ang)]))
    return numpy.array(out)


def _geometry_shape(geom):
    r'''! Outer polygon (N,2) and inner hole polygon (or None) for a DD4 `pf_active`/`pf_passive` element geometry'''
    gtype = int(geom.geometry_type.value) if geom.geometry_type.has_value else 0
    if gtype <= 0:
        # Infer the type from the filled structure
        if len(geom.outline.r) > 0:
            gtype = 1
        elif geom.rectangle.width.has_value:
            gtype = 2
        elif geom.oblique.length_alpha.has_value:
            gtype = 3
        elif len(geom.arcs_of_circle.r) > 0:
            gtype = 4
        elif geom.annulus.radius_outer.has_value:
            gtype = 5
        elif geom.thick_line.thickness.has_value:
            gtype = 6
    if gtype == 1:
        return numpy.column_stack((geom.outline.r, geom.outline.z)), None
    if gtype == 2:
        rect = geom.rectangle
        w, h = rect.width/2.0, rect.height/2.0
        return numpy.array([[rect.r-w, rect.z-h], [rect.r+w, rect.z-h], [rect.r+w, rect.z+h], [rect.r-w, rect.z+h]]), None
    if gtype == 3:
        obl = geom.oblique
        p0 = numpy.array([obl.r, obl.z])
        va = obl.length_alpha*numpy.array([numpy.cos(obl.alpha), numpy.sin(obl.alpha)])
        vb = obl.length_beta*numpy.array([-numpy.sin(obl.beta), numpy.cos(obl.beta)])
        return numpy.array([p0, p0+va, p0+va+vb, p0+vb]), None
    if gtype == 4:
        arcs = geom.arcs_of_circle
        return _arcs_polygon(numpy.asarray(arcs.r), numpy.asarray(arcs.z), numpy.asarray(arcs.curvature_radii)), None
    if gtype == 5:
        ann = geom.annulus
        hole = None
        if ann.radius_inner.has_value and ann.radius_inner > 0.0:
            hole = _circle(ann.r, ann.z, ann.radius_inner)
        return _circle(ann.r, ann.z, ann.radius_outer), hole
    if gtype == 6:
        tl = geom.thick_line
        p1 = numpy.array([tl.first_point.r, tl.first_point.z])
        p2 = numpy.array([tl.second_point.r, tl.second_point.z])
        normal = numpy.array([-(p2-p1)[1], (p2-p1)[0]])
        normal *= tl.thickness/2.0/numpy.linalg.norm(normal)
        return numpy.array([p1-normal, p2-normal, p2+normal, p1+normal]), None
    raise ValueError('Element geometry is empty or has an unknown geometry_type ({0})'.format(gtype))


def _polygon_area(pts):
    '''! Area of a closed polygon (N,2)'''
    return 0.5*abs(numpy.dot(pts[:, 0], numpy.roll(pts[:, 1], -1))-numpy.dot(pts[:, 1], numpy.roll(pts[:, 0], -1)))


def _clear_node(parent, child):
    '''! Empty one filled child node of an IDS structure (recursively)'''
    imas = _imas()
    name = child.metadata.name
    if isinstance(child, imas.ids_struct_array.IDSStructArray):
        child.resize(0)
    elif isinstance(child, imas.ids_structure.IDSStructure):
        _clear_struct(child)
    elif child.metadata.ndim > 0:
        setattr(parent, name, numpy.zeros((0,)*child.metadata.ndim))
    elif str(child.metadata.data_type).endswith('STR'):
        setattr(parent, name, '')
    elif str(child.metadata.data_type).endswith('INT'):
        setattr(parent, name, imas.ids_defs.EMPTY_INT)
    else:
        setattr(parent, name, imas.ids_defs.EMPTY_FLOAT)


def _clear_struct(struct):
    '''! Empty every filled field of an IDS structure (recursively)'''
    for child in list(struct.iter_nonempty_()):
        _clear_node(struct, child)


def _polygon_inv_r(pts):
    r'''! \f$ \int dA/R \f$ over a closed polygon (N,2), exactly (Green's theorem: \f$ \oint \ln R\, dZ \f$)'''
    r1, z1 = pts[:, 0], pts[:, 1]
    r2, z2 = numpy.roll(r1, -1), numpy.roll(z1, -1)
    with numpy.errstate(divide='ignore', invalid='ignore'):
        mean_ln = numpy.where(numpy.abs(r2-r1) > 1.E-14*numpy.maximum(r1, r2),
                              (r2*numpy.log(r2)-r1*numpy.log(r1))/(r2-r1)-1.0, numpy.log(r1))
    return abs(numpy.sum(mean_ln*(z2-z1)))


def _clean_outline(r, z):
    '''! Closed (R,Z) outline as an (N,2) array without a repeated last point'''
    pts = numpy.column_stack((r, z))
    if (len(pts) > 1) and numpy.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    return pts


def _positive(node, what):
    '''! Value of a filled, positive IDS float node, or an error naming it'''
    if (not node.has_value) or (not numpy.isfinite(float(node))) or (float(node) <= 0.0):
        raise ValueError('{0} must be given and positive'.format(what))
    return float(node)


def _add_shape(domain, name, shape, dx_hole):
    '''! Add a polygon region, or an annulus around a vacuum hole region `<name>_HOLE`'''
    outer, hole = shape
    if hole is None:
        domain.add_polygon(outer, name)
    else:
        domain.define_region(name+'_HOLE', dx_hole, 'vacuum')
        domain.add_annulus(hole, name+'_HOLE', outer, name)


def domain_from_imas(wall, pf_active, dx, pf_passive=None, description_index=0, rextent=None, zextents=(None, None)):
    r'''! Build a TokaMaker mesh domain from an IMAS machine description

    Regions created:
    - `PLASMA`: inside `wall.description_2d[description_index].limiter.unit[0].outline`
    - Each annular vessel unit: a conductor (resistivity from the unit) and a vacuum region
      (`VACUUM<k>`) inside it; vessel units given as elements become one conductor per element
    - Each `pf_passive.loop[].element`: a conductor (resistivity, or for a single-element loop from the
      loop resistance)
    - Each `pf_active.coil`: a coil set, with one coil region per `element` and
      `nTurns = turns_with_sign`
    - A vacuum region `<name>_HOLE` inside each element with an annulus geometry

    Region names are the IMAS names (upper case, spaces replaced), with `_E<j>` appended for
    multi-element objects, so @ref pf_active_to_imas and @ref pf_passive_to_imas can match them.
    Circuits (`pf_active.circuit`) are not mapped. TokaMaker cannot force the elements of a
    multi-element passive loop to carry the same current; they become separate conductors.

    @param wall `wall` IDS
    @param pf_active `pf_active` IDS
    @param dx Target mesh size for each region type: dictionary with keys 'plasma', 'vacuum',
    'boundary', 'coil' and 'conductor'
    @param pf_passive `pf_passive` IDS (optional)
    @param description_index Index of `wall.description_2d` to use
    @param rextent Radial extent of the domain (default: padded from the geometry, see @ref TokaMaker.meshing.gs_Domain "gs_Domain")
    @param zextents Vertical extents of the domain (default: padded from the geometry)
    @result @ref TokaMaker.meshing.gs_Domain "gs_Domain" with all regions defined
    '''
    check_dd4(wall)
    check_dd4(pf_active)
    if pf_passive is not None:
        check_dd4(pf_passive)
    for key in ('plasma', 'vacuum', 'boundary', 'coil', 'conductor'):
        if key not in dx:
            raise KeyError('Missing mesh size "{0}" in dx'.format(key))
    domain = gs_Domain(rextent=rextent, zextents=zextents)
    domain.define_region('AIR', dx['boundary'], 'boundary')
    desc = wall.description_2d[description_index]
    if len(desc.limiter.unit) == 0:
        raise ValueError('wall.description_2d[{0}] has no limiter outline'.format(description_index))
    limiter = _clean_outline(desc.limiter.unit[0].outline.r, desc.limiter.unit[0].outline.z)
    domain.define_region('PLASMA', dx['plasma'], 'plasma')
    domain.add_polygon(limiter, 'PLASMA')
    # Vessel: annular units (inner region becomes a vacuum region), then element units
    annular = []
    for i, unit in enumerate(desc.vessel.unit):
        if len(unit.annular.outline_inner.r) > 0:
            annular.append((i, unit))
        else:
            for j, element in enumerate(unit.element):
                name = _region_name(str(element.name), _region_name(str(unit.name), 'VESSEL', i)+'_E', j)
                eta = _positive(element.resistivity, 'Resistivity of vessel unit {0} element {1}'.format(i, j))
                domain.define_region(name, dx['conductor'], 'conductor', eta=eta)
                domain.add_polygon(_clean_outline(element.outline.r, element.outline.z), name)
    annular.sort(key=lambda item: _polygon_area(_clean_outline(item[1].annular.outline_inner.r, item[1].annular.outline_inner.z)))
    for k, (i, unit) in enumerate(annular):
        inner = _clean_outline(unit.annular.outline_inner.r, unit.annular.outline_inner.z)
        outer = _clean_outline(unit.annular.outline_outer.r, unit.annular.outline_outer.z)
        vac_name = 'VACUUM{0}'.format(k+1)
        name = _region_name(str(unit.name), 'VESSEL', i)
        eta = _positive(unit.annular.resistivity, 'Resistivity of annular vessel unit {0}'.format(i))
        domain.define_region(vac_name, dx['vacuum'], 'vacuum')
        domain.define_region(name, dx['conductor'], 'conductor', eta=eta)
        domain.add_annulus(inner, vac_name, outer, name)
    # Passive loops
    if pf_passive is not None:
        for i, loop in enumerate(pf_passive.loop):
            nelem = len(loop.element)
            if nelem > 1:
                warnings.warn('pf_passive loop "{0}" has {1} elements; each becomes a separate conductor'.format(loop.name, nelem))
            for j, element in enumerate(loop.element):
                shape = _geometry_shape(element.geometry)
                if loop.resistivity.has_value and loop.resistivity > 0.0:
                    eta = float(loop.resistivity)
                elif (nelem == 1) and loop.resistance.has_value and loop.resistance > 0.0:
                    # Loop resistance -> resistivity of a solid ring with this cross-section, for which
                    # R = 2 pi eta / int(dA/R). The region carries N times the loop current, so the same
                    # ohmic power needs R_loop/N^2.
                    inv_r = _polygon_inv_r(shape[0])-(_polygon_inv_r(shape[1]) if shape[1] is not None else 0.0)
                    nturns = float(element.turns_with_sign) if element.turns_with_sign.has_value else 1.0
                    eta = float(loop.resistance)*inv_r/(_TWO_PI*nturns**2)
                else:
                    raise ValueError('pf_passive loop "{0}" needs a positive resistivity (or, for a single element, '
                                     'resistance)'.format(loop.name))
                name = _region_name(str(loop.name), 'LOOP', i, j if nelem > 1 else None)
                domain.define_region(name, dx['conductor'], 'conductor', eta=eta)
                _add_shape(domain, name, shape, dx['vacuum'])
    # Coils
    for i, coil in enumerate(pf_active.coil):
        coil_set = _region_name(str(coil.name), 'COIL', i)
        elements = list(coil.element)
        if len(elements) == 0:
            if not _has_node(coil, 'geometry'):  # Coil-level geometry is DD >= 4.1
                raise ValueError('pf_active coil "{0}" has no elements'.format(coil.name))
            domain.define_region(coil_set, dx['coil'], 'coil', coil_set=coil_set, nTurns=1.0)
            _add_shape(domain, coil_set, _geometry_shape(coil.geometry), dx['vacuum'])
            continue
        for j, element in enumerate(elements):
            name = _region_name(str(coil.name), 'COIL', i, j if len(elements) > 1 else None)
            nturns = float(element.turns_with_sign) if element.turns_with_sign.has_value else 1.0
            domain.define_region(name, dx['coil'], 'coil', coil_set=coil_set, nTurns=nturns)
            _add_shape(domain, name, _geometry_shape(element.geometry), dx['vacuum'])
    return domain


#------------------------------------------------------------------------------
# Input: profiles and targets
#------------------------------------------------------------------------------
def _nearest(times, time):
    '''! Index of the nearest time (0 if `time` is None)'''
    times = numpy.atleast_1d(numpy.asarray(times, dtype=numpy.float64))
    if (time is None) or (times.shape[0] == 0):
        return 0
    return int(numpy.argmin(numpy.abs(times-time)))


def _slice_index(ids, slices, time):
    '''! Index of the time slice of an array of structures nearest `time`, for homogeneous or heterogeneous time'''
    if len(slices) == 0:
        raise ValueError('"{0}" has no time slices'.format(ids.metadata.name))
    mode = int(ids.ids_properties.homogeneous_time.value)
    if mode == 1:
        times = ids.time
    elif mode == 0:
        times = [float(sl.time) for sl in slices]
    else:
        raise ValueError('"{0}" has no time-dependent data (homogeneous_time={1})'.format(ids.metadata.name, mode))
    return _nearest(times, time)


def _signal_index(ids, signal, time):
    '''! Index of the sample of a `data`/`time` signal nearest `time`, for homogeneous or heterogeneous time'''
    mode = int(ids.ids_properties.homogeneous_time.value)
    if mode == 1:
        return _nearest(ids.time, time)
    if mode == 0:
        return _nearest(signal.time, time)
    raise ValueError('"{0}" has no time-dependent data (homogeneous_time={1})'.format(ids.metadata.name, mode))


def _normalized(y):
    '''! Profile normalized to its axis value (or to its largest magnitude if that is zero; zero profiles stay zero)'''
    y = numpy.asarray(y, dtype=numpy.float64)
    if not numpy.any(y != 0.0):
        return y.copy()
    scale = y[0] if y[0] != 0.0 else y[numpy.argmax(numpy.abs(y))]
    return y/scale


def _positive_ip(ip, source):
    '''! Plasma current, checking it is positive (TokaMaker's targets and constraints require Ip > 0)'''
    if ip <= 0.0:
        raise ValueError('Non-positive plasma current in {0}; TokaMaker requires Ip > 0 (mirror the discharge to a '
                         'positive-Ip convention first)'.format(source))
    return ip


def profiles_from_imas(source, time=None):
    r'''! Source profiles and targets for TokaMaker from an `equilibrium` IDS, or kinetic profiles from a `core_profiles` IDS

    From `equilibrium` (nearest time slice), the result has:
    - 'profiles': `ffp_prof` and `pp_prof` for `set_profiles` (normalized to their axis values)
    - 'targets': `Ip` and `pax` for `set_targets`
    - 'F0': \f$ R B_\phi \f$ in vacuum, from `vacuum_toroidal_field`
    - 'isoflux': boundary outline points (N,2) for `set_isoflux_constraints`

    From `core_profiles` (nearest time), the result has keyword arguments for
    `set_kinetic_profiles` on normalized toroidal flux (`coord='phi_n'`, \f$ \hat{\Phi} = \rho_{tor,norm}^2 \f$):
    `te_prof` and `ti_prof` [keV], `ne_prof` and `ni_prof` [m^-3], and `Zeff` (each only if present).
    Densities are thermal (fast populations are excluded). TokaMaker uses one effective ion species,
    so `ni_prof` is the total thermal ion density and `ti_prof` the thermal-density-weighted ion
    temperature, which keeps the thermal ion pressure exact.

    @param source `equilibrium` or `core_profiles` IDS (Data Dictionary 4.x, homogeneous or heterogeneous time)
    @param time Time [s] (default: first slice)
    @result Dictionary of profiles and targets
    '''
    check_dd4(source)
    name = source.metadata.name
    if name == 'equilibrium':
        it = _slice_index(source, source.time_slice, time)
        ts = source.time_slice[it]
        p1 = ts.profiles_1d
        if len(p1.psi_norm) > 0:
            psi_n = numpy.asarray(p1.psi_norm)
        else:
            # psi_boundary is at the selected boundary, which lies at boundary.psi_norm of the LCFS flux
            psi = numpy.asarray(p1.psi)
            bnd_norm = float(ts.boundary.psi_norm) if ts.boundary.psi_norm.has_value else 1.0
            psi_n = bnd_norm*(psi-ts.global_quantities.psi_axis)/(ts.global_quantities.psi_boundary-ts.global_quantities.psi_axis)
        for label, vals in (('f_df_dpsi', p1.f_df_dpsi), ('dpressure_dpsi', p1.dpressure_dpsi), ('pressure', p1.pressure)):
            vals = numpy.asarray(vals)
            if (vals.shape != psi_n.shape) or (not numpy.all(numpy.isfinite(vals))):
                raise ValueError('equilibrium.time_slice[{0}].profiles_1d.{1} must be finite and match the profile '
                                 'coordinate ({2} points)'.format(it, label, psi_n.shape[0]))
        if (psi_n.shape[0] < 2) or (not numpy.all(numpy.isfinite(psi_n))):
            raise ValueError('equilibrium.time_slice[{0}].profiles_1d has no valid flux coordinate'.format(it))
        if not ts.global_quantities.ip.has_value:
            raise ValueError('equilibrium.time_slice[{0}].global_quantities.ip is not set'.format(it))
        r0 = _positive(source.vacuum_toroidal_field.r0, 'equilibrium.vacuum_toroidal_field.r0')
        # B0 is on the root time base: use the same index for homogeneous time, else the time nearest the slice
        b0 = numpy.asarray(source.vacuum_toroidal_field.b0)
        if (b0.shape[0] == 0) or (not numpy.all(numpy.isfinite(b0))):
            raise ValueError('equilibrium.vacuum_toroidal_field.b0 is empty or not finite')
        if int(source.ids_properties.homogeneous_time.value) == 1:
            ib = it
        elif b0.shape[0] == 1:
            ib = 0
        else:
            ib = _nearest(source.time, float(ts.time))
        # COCOS 17 -> 7: d/dpsi gains a factor 2 pi (shape normalization keeps the relative sign)
        out = {
            'profiles': {
                'ffp_prof': {'type': 'linterp', 'x': psi_n, 'y': _normalized(_TWO_PI*numpy.asarray(p1.f_df_dpsi))},
                'pp_prof': {'type': 'linterp', 'x': psi_n, 'y': _normalized(_TWO_PI*numpy.asarray(p1.dpressure_dpsi))},
            },
            'targets': {
                'Ip': _positive_ip(float(ts.global_quantities.ip), 'equilibrium.time_slice[{0}].global_quantities.ip'.format(it)),
                'pax': float(p1.pressure[0]),
            },
            'F0': r0*float(b0[ib]),
        }
        if len(ts.boundary.outline.r) > 0:
            out['isoflux'] = _clean_outline(ts.boundary.outline.r, ts.boundary.outline.z)
        return out
    if name == 'core_profiles':
        it = _slice_index(source, source.profiles_1d, time)
        p1 = source.profiles_1d[it]
        x = numpy.asarray(p1.grid.rho_tor_norm)**2
        if len(x) == 0:
            raise ValueError('core_profiles.profiles_1d[{0}].grid.rho_tor_norm is empty'.format(it))

        def prof(vals, scale=1.0):
            '''! Piecewise-linear profile on normalized toroidal flux'''
            return {'type': 'linterp', 'x': x, 'y': scale*numpy.asarray(vals), 'coord': 'phi_n'}

        out = {}
        if len(p1.electrons.temperature) > 0:
            out['te_prof'] = prof(p1.electrons.temperature, 1.E-3)
        ne = _thermal_density(p1.electrons)
        if ne is not None:
            out['ne_prof'] = prof(ne)
        ni, ti = _effective_ions(p1)
        if ti is not None:
            out['ti_prof'] = prof(ti, 1.E-3)
        if ni is not None:
            out['ni_prof'] = prof(ni)
        if len(p1.zeff) > 0:
            out['Zeff'] = prof(p1.zeff)
        return out
    raise ValueError('Unsupported IDS "{0}" (expected "equilibrium" or "core_profiles")'.format(name))


def _thermal_density(species):
    '''! Thermal density of an electron or ion population: `density_thermal`, else `density` less `density_fast`, else `density`'''
    if len(species.density_thermal) > 0:
        return numpy.asarray(species.density_thermal)
    if len(species.density) == 0:
        return None
    dens = numpy.asarray(species.density)
    if len(species.density_fast) > 0:
        dens = dens-numpy.asarray(species.density_fast)
    return dens


def _effective_ions(p1):
    r'''! Total thermal ion density and thermal-density-weighted ion temperature of a `core_profiles.profiles_1d` slice

    Where the thermal ion density vanishes, the temperature is `t_i_average` if given, else the mean
    of the ion temperatures.

    @result Density [m^-3] and temperature [eV] (each None if they cannot be formed)
    '''
    dens = [_thermal_density(ion) for ion in p1.ion]
    temps = [numpy.asarray(ion.temperature) for ion in p1.ion]
    have_dens = (len(dens) > 0) and all((d is not None) and d.shape[0] > 0 for d in dens)
    have_temps = (len(temps) > 0) and all(t.shape[0] > 0 for t in temps)
    ni = None
    if len(p1.n_i_thermal_total) > 0:
        ni = numpy.asarray(p1.n_i_thermal_total)
    elif have_dens:
        ni = numpy.sum(dens, axis=0)
    t_avg = numpy.asarray(p1.t_i_average) if len(p1.t_i_average) > 0 else None
    ti = None
    if have_dens and have_temps:
        total = numpy.sum(dens, axis=0)
        fallback = t_avg if t_avg is not None else numpy.mean(temps, axis=0)
        with numpy.errstate(divide='ignore', invalid='ignore'):
            ti = numpy.where(total > 0.0, numpy.sum([d*t for d, t in zip(dens, temps)], axis=0)/total, fallback)
    elif t_avg is not None:
        ti = t_avg
    elif have_temps and (len(temps) == 1):
        ti = temps[0]
    return ni, ti


def F0_from_tf(tf, time=None):
    r'''! Vacuum \f$ F = R B_\phi \f$ from a `tf` IDS

    @param tf `tf` IDS (Data Dictionary 4.x)
    @param time Time [s] (default: first time)
    @result F0 [T m]
    '''
    check_dd4(tf)
    vals = numpy.asarray(tf.b_field_phi_vacuum_r.data)
    if vals.shape[0] == 0:
        raise ValueError('tf.b_field_phi_vacuum_r is empty')
    return float(vals[_signal_index(tf, tf.b_field_phi_vacuum_r, time)])


def _reference_at(signal, times, time):
    '''! Value of a pulse_schedule reference at `time` (None if empty)'''
    ref = numpy.asarray(signal.reference)
    if ref.shape[0] == 0:
        return None
    if (time is None) or (ref.shape[0] == 1):
        return float(ref[0])
    return float(numpy.interp(time, times, ref))


def targets_from_pulse_schedule(ps, time=None):
    r'''! Shape and current targets for TokaMaker from a `pulse_schedule` IDS

    @param ps `pulse_schedule` IDS (Data Dictionary 4.x, homogeneous or heterogeneous time)
    @param time Time [s] (default: first time)
    @result Dictionary with 'isoflux' (boundary outline and strike points, (N,2)),
    'saddles' (X-points, (N,2)), and 'R0', 'Z0' and 'Ip' targets, each only if present
    '''
    check_dd4(ps)
    pc = ps.position_control
    mode = int(ps.ids_properties.homogeneous_time.value)
    if mode not in (0, 1):
        raise ValueError('pulse_schedule has no time-dependent data (homogeneous_time={0})'.format(mode))
    pc_time = numpy.asarray(ps.time if mode == 1 else pc.time)
    fc_time = numpy.asarray(ps.time if mode == 1 else ps.flux_control.time)

    def points(arr):
        '''! (N,2) reference points at `time`'''
        pts = []
        for pt in arr:
            r, z = _reference_at(pt.r, pc_time, time), _reference_at(pt.z, pc_time, time)
            if (r is not None) and (z is not None):
                pts.append([r, z])
        return numpy.array(pts).reshape(-1, 2)

    out = {}
    isoflux = numpy.vstack((points(pc.boundary_outline), points(pc.strike_point)))
    if isoflux.shape[0] > 0:
        out['isoflux'] = isoflux
    saddles = points(pc.x_point)
    if saddles.shape[0] > 0:
        out['saddles'] = saddles
    r0 = _reference_at(pc.magnetic_axis.r, pc_time, time)
    z0 = _reference_at(pc.magnetic_axis.z, pc_time, time)
    if r0 is not None:
        out['R0'] = r0
    if z0 is not None:
        out['Z0'] = z0
    ip = _reference_at(ps.flux_control.ip, fc_time, time)
    if ip is not None:
        out['Ip'] = _positive_ip(ip, 'pulse_schedule.flux_control.ip')
    return out


def coil_currents_from_imas(pf_active, mygs, time=None):
    r'''! Coil currents per turn for each TokaMaker coil set from a `pf_active` IDS

    Coils are matched to coil sets by exact name, then by the name @ref domain_from_imas forms.

    @param pf_active `pf_active` IDS (Data Dictionary 4.x)
    @param mygs TokaMaker object
    @param time Time [s] (default: first time)
    @result Dictionary of coil currents [A] (per turn), dictionary of their upper errors (None if
    absent), and dictionary of their (source path, index, number of coils)
    '''
    check_dd4(pf_active)
    currents = {}
    errors = {}
    sources = {}
    for i, coil in enumerate(pf_active.coil):
        name = _match_coil_set(coil.name, i, mygs.coil_sets)
        if name is None:
            continue
        data = numpy.asarray(coil.current.data)
        if data.shape[0] == 0:
            continue
        it = _signal_index(pf_active, coil.current, time)
        currents[name] = float(data[it])
        err = numpy.asarray(coil.current.data_error_upper)
        errors[name] = float(err[it]) if err.shape[0] > it else None
        sources[name] = ('pf_active/coil({0})/current'.format(i+1), i, len(pf_active.coil))
    return currents, errors, sources


def _sample(ids, signal, time, min_validity):
    '''! Value and upper error (None if absent) of a `data` signal at the nearest time, or None if absent or invalid'''
    data = numpy.asarray(signal.data)
    if data.shape[0] == 0:
        return None, None
    it = _signal_index(ids, signal, time)
    if _has_node(signal, 'validity') and signal.validity.has_value and int(signal.validity.value) < min_validity:
        return None, None
    if _has_node(signal, 'validity_timed'):
        validity_timed = numpy.asarray(signal.validity_timed)
        if (validity_timed.shape[0] > it) and (validity_timed[it] < min_validity):
            return None, None
    err = numpy.asarray(signal.data_error_upper)
    val = float(data[it])
    err = float(err[it]) if err.shape[0] > it else None
    if (not numpy.isfinite(val)) or ((err is not None) and (not numpy.isfinite(err))):
        raise ValueError('Non-finite value or error at {0} (sample {1})'.format(_imas().util.get_full_path(signal), it))
    return val, err


def load_magnetics(recon, magnetics, pf_active=None, time=None, rel_err=0.01, abs_err=None, min_validity=0):
    r'''! Add reconstruction constraints from `magnetics` (and `pf_active` coil currents)

    - `flux_loop` of type `poloidal_flux` (or unset) with one complete (R,Z) position: poloidal flux
      [Wb] there (COCOS 17 \f$ \psi \f$ is TokaMaker's flux loop signal). Saddle, diamagnetic,
      differential and multi-segment loops are skipped with a warning.
    - `b_field_pol_probe` with a position and poloidal angle (differential probes are skipped): field along the probe normal
      \f$ n = (\cos\alpha\cos\beta, \cos\alpha\sin\beta, -\sin\alpha) \f$ in \f$ (R,\phi,Z) \f$, with
      \f$ \alpha \f$ the poloidal and \f$ \beta \f$ the toroidal angle (DD 4.1)
    - `ip`: plasma current (first method; must be positive)
    - `pf_active` coil currents, which must cover every coil set (checked before anything is added)

    Samples whose `validity` or `validity_timed` is below `min_validity` are skipped (the DD uses
    negative values for problems: -1 needs verification, -2 and below are invalid). Errors are
    `data_error_upper` where present, otherwise `rel_err` times the value with a floor of
    `abs_err[kind]` (keys 'flux_loop', 'b_field_pol_probe', 'ip', 'coil'); every error must end up
    positive, so give a floor for signals that can be zero.

    The source path and array index of each constraint are kept with it for @ref recon_constraints_to_imas.

    @param recon @ref TokaMaker.reconstruction.reconstruction "reconstruction" object
    @param magnetics `magnetics` IDS (Data Dictionary 4.x)
    @param pf_active `pf_active` IDS for coil current constraints (optional)
    @param time Time [s] (default: first time)
    @param rel_err Relative error used where no error is given
    @param abs_err Dictionary of absolute error floors by kind (optional)
    @param min_validity Lowest accepted validity flag
    @result Number of constraints added
    '''
    check_dd4(magnetics)
    abs_err = {} if abs_err is None else abs_err

    def error(val, err, kind, source):
        '''! Given error if positive, else the relative error with its floor (must be positive)'''
        if (err is not None) and (err > 0.0):
            return err
        out = max(rel_err*abs(val), abs_err.get(kind, 0.0))
        if out <= 0.0:
            raise ValueError('No positive uncertainty for {0}; give data_error_upper or abs_err["{1}"]'.format(source, kind))
        return out

    # Collect everything first, so nothing is added to `recon` if an input is incomplete
    flux_loops, probes, skipped = [], [], []
    for i, loop in enumerate(magnetics.flux_loop):
        source = 'magnetics/flux_loop({0})/flux'.format(i+1)
        ltype = int(loop.type.index.value) if loop.type.index.has_value else 0
        # Only full toroidal loops given by one complete (R,Z) position (multi-segment loops are not supported)
        if (ltype not in (0, 1)) or (len(loop.position) != 1) or (not loop.position[0].r.has_value) \
                or (not loop.position[0].z.has_value):
            skipped.append(source)
            continue
        val, err = _sample(magnetics, loop.flux, time, min_validity)
        if val is None:
            continue
        loc = numpy.array([loop.position[0].r, loop.position[0].z])
        flux_loops.append((loc, val, error(val, err, 'flux_loop', source), source, i, len(magnetics.flux_loop)))
    if skipped:
        warnings.warn('Skipped flux loops that are not single-position poloidal flux loops: {0}'.format(', '.join(skipped)))
    skipped_probes = []
    for i, probe in enumerate(magnetics.b_field_pol_probe):
        source = 'magnetics/b_field_pol_probe({0})/field'.format(i+1)
        ptype = int(probe.type.index.value) if probe.type.index.has_value else 0
        if (ptype == 6) or (len(probe.indices_differential) > 0) or (not probe.position.r.has_value) \
                or (not probe.position.z.has_value) or (not probe.poloidal_angle.has_value):
            skipped_probes.append(source)
            continue
        val, err = _sample(magnetics, probe.field, time, min_validity)
        if val is None:
            continue
        alpha = float(probe.poloidal_angle)
        beta = float(probe.toroidal_angle) if probe.toroidal_angle.has_value else 0.0
        norm = numpy.array([numpy.cos(alpha)*numpy.cos(beta), numpy.cos(alpha)*numpy.sin(beta), -numpy.sin(alpha)])
        loc = numpy.array([probe.position.r, probe.position.z])
        probes.append((loc, norm, val, error(val, err, 'b_field_pol_probe', source), source, i, len(magnetics.b_field_pol_probe)))
    if skipped_probes:
        warnings.warn('Skipped differential probes and probes without position or poloidal angle: {0}'.format(
            ', '.join(skipped_probes)))
    ip = None
    if len(magnetics.ip) > 0:
        val, err = _sample(magnetics, magnetics.ip[0], time, min_validity)
        if val is not None:
            ip = (_positive_ip(val, 'magnetics.ip'), error(val, err, 'ip', 'magnetics/ip(1)'), 'magnetics/ip(1)')
    coils = None
    if pf_active is not None:
        currents, errors, coil_sources = coil_currents_from_imas(pf_active, recon._tMaker_obj, time)
        missing = [name for name in recon._tMaker_obj.coil_sets if name not in currents]
        if missing:
            raise KeyError('No pf_active current for coil sets: {0}'.format(', '.join(missing)))
        errs = {name: error(currents[name], errors[name], 'coil', coil_sources[name][0]) for name in currents}
        coils = (currents, errs, coil_sources)
    # Add constraints, each tagged with its source path and slot (index, size) in the source array
    for loc, val, err, source, index, size in flux_loops:
        recon.add_flux_loop(loc, val, err)
        recon._flux_loops[-1].imas_source = (source, index, size)
    for loc, norm, val, err, source, index, size in probes:
        recon.add_Mirnov(loc, norm, val, err)
        recon._mirnovs[-1].imas_source = (source, index, size)
    ncons = len(flux_loops)+len(probes)
    if ip is not None:
        recon.set_Ip(ip[0], ip[1])
        recon._Ip_con.imas_source = (ip[2], 0, 1)
        ncons += 1
    if coils is not None:
        recon.set_coil_currents(coils[0], coils[1])
        coil_names = {cset['id']: name for name, cset in recon._tMaker_obj.coil_sets.items()}
        for con in recon._coil_current_cons:
            con.imas_source = coils[2][coil_names[con.ind]]
        ncons += len(coils[0])
    return ncons




def recon_constraints_to_imas(recon, eq_ids, slice_index=0, error_mat=None, freedom_degrees_n=None):
    r'''! Write reconstruction constraints to `equilibrium.time_slice[].constraints`

    The constraints of the slice are replaced by those currently set on `recon`, so the structure
    describes exactly this fit. Uses the error matrix of
    @ref TokaMaker.reconstruction.reconstruction.eval_error "eval_error" (evaluated here if not given;
    constraints must be set up). Flux loop, poloidal probe, plasma current and coil current
    constraints are written, each with its measured and reconstructed value, its uncertainty
    (`measured_error_upper`), unit `weight`, and its \f$ \chi^2 \f$. Constraints from
    @ref load_magnetics are placed at the index of their source (`magnetics/flux_loop`,
    `magnetics/b_field_pol_probe` and `pf_active/coil`, as these arrays are defined) with their
    `source` path, leaving skipped sensors empty; two constraints from the same source are an error.
    Other constraints are appended after them. `constraints_n` counts all constraints of the fit
    (including saddle loop, pressure and diamagnetic flux constraints, which have no slot here), and
    `chi_squared_reduced` is written only if the number of degrees of freedom of the fit is given.

    @param recon @ref TokaMaker.reconstruction.reconstruction "reconstruction" object
    @param eq_ids `equilibrium` IDS to fill
    @param slice_index Time slice to fill
    @param error_mat Error matrix from `eval_error` (optional)
    @param freedom_degrees_n Number of degrees of freedom of the fit (constraints minus fitted parameters, optional)
    @result `eq_ids`
    '''
    check_dd4(eq_ids)
    if error_mat is None:
        recon.setup_constraints()
        error_mat = recon.eval_error()
    cons = eq_ids.time_slice[slice_index].constraints

    def fill(node, con, row):
        '''! One constraint node from an error-matrix row'''
        node.measured = error_mat[row, 2]
        node.measured_error_upper = con.err
        node.reconstructed = error_mat[row, 1]
        node.weight = 1.0
        node.exact = 0
        node.chi_squared = error_mat[row, 0]**2
        source = getattr(con, 'imas_source', None)
        if source is not None:
            node.source = source[0]

    def slots(group, name):
        '''! Array index of each constraint (source slot if from IMAS, else appended), and the array size'''
        tags = [getattr(con, 'imas_source', None) for con in group]
        size = max([tag[2] for tag in tags if tag is not None], default=0)
        index = []
        used = set()
        for tag in tags:
            if tag is None:
                index.append(size)
                size += 1
            else:
                if tag[1] in used:
                    raise ValueError('Two {0} constraints come from the same source ({1})'.format(name, tag[0]))
                used.add(tag[1])
                index.append(tag[1])
        return index, size

    groups = [
        ('flux_loop', recon._flux_loops),
        ('b_field_pol_probe', recon._mirnovs),
        (None, recon._saddles),
        (None, recon._pressure_cons),
        ('pf_current', recon._coil_current_cons),
    ]
    # Validate everything before the destination is modified
    layout = {name: slots(group, name) for name, group in groups if name is not None}
    nrows = sum(len(group) for _, group in groups)+(recon._Ip_con is not None)+(recon._Dflux_con is not None)
    error_mat = numpy.asarray(error_mat)
    if (error_mat.ndim != 2) or (error_mat.shape[0] != nrows) or (error_mat.shape[1] < 3):
        raise ValueError('error_mat has shape {0}, expected ({1}, 4) for the constraints set on recon'.format(error_mat.shape, nrows))
    _clear_struct(cons)
    row = 0
    chi2 = 0.0
    for name, group in groups:
        if name is not None:
            getattr(cons, name).resize(layout[name][1])
        for k, con in enumerate(group):
            if name is not None:
                fill(getattr(cons, name)[layout[name][0][k]], con, row)
            chi2 += error_mat[row, 0]**2
            row += 1
    if recon._Ip_con is not None:
        fill(cons.ip, recon._Ip_con, row)
        chi2 += error_mat[row, 0]**2
        row += 1
    if recon._Dflux_con is not None:
        chi2 += error_mat[row, 0]**2
        row += 1
    cons.constraints_n = row
    if freedom_degrees_n is not None:
        cons.freedom_degrees_n = int(freedom_degrees_n)
        if freedom_degrees_n > 0:
            cons.chi_squared_reduced = chi2/freedom_degrees_n
    return eq_ids
