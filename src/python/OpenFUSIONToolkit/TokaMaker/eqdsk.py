#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! GEQDSK (g-file) and Osborne p-file (kinetic profile) readers and writers

Reads standard GEQDSK equilibrium files and computes flux-surface-averaged
quantities, and reads/writes Osborne p-files for kinetic profiles.  Adapted
from OMFIT `omfit_eqdsk` / `omfit_osborne` with no OMFIT runtime dependency.

Dependencies: numpy, scipy, contourpy (ships with matplotlib).

@authors Daniel Burgess
@date April 2026
@ingroup doxy_oft_python
'''

import io
import re
import tempfile
import warnings
from collections import OrderedDict

import contourpy
import numpy as np
from scipy import constants, integrate, interpolate


# ===========================================================================
# GEQDSK (g-file) reader/writer and equilibrium analysis
# ===========================================================================


# ---------------------------------------------------------------------------
# COCOS parameter table
# ---------------------------------------------------------------------------

def _cocos_params(cocos_index):
    r'''! Return COCOS sign/exponent parameters for a given COCOS index

    @param cocos_index COCOS convention index (1-8 or 11-18)
    @result Dict with keys: `sigma_Bp`, `sigma_RpZ`, `sigma_rhotp`, `exp_Bp`
    '''
    if cocos_index < 1 or cocos_index > 18 or cocos_index in (9, 10):
        raise ValueError(f"Invalid COCOS index: {cocos_index}")

    exp_Bp = 0 if cocos_index < 10 else 1
    base = cocos_index if cocos_index < 10 else cocos_index - 10

    sigma_Bp = +1 if base in (1, 2, 5, 6) else -1
    sigma_RpZ = +1 if base in (1, 2, 7, 8) else -1
    sigma_rhotp = +1 if base in (1, 3, 5, 7) else -1

    return {
        "sigma_Bp": sigma_Bp,
        "sigma_RpZ": sigma_RpZ,
        "sigma_rhotp": sigma_rhotp,
        "exp_Bp": exp_Bp,
    }


# ---------------------------------------------------------------------------
# GEQDSK parser
# ---------------------------------------------------------------------------

def _read_geqdsk(filename):
    r'''! Parse a GEQDSK (g-file) into a plain dict

    @param filename Path to the g-file
    @result Dict whose keys match the standard GEQDSK field names (`NW`, `NH`,
            `RDIM`, ..., `FPOL`, `PRES`, `FFPRIM`, `PPRIME`, `PSIRZ`, `QPSI`,
            `RBBBS`, `ZBBBS`, etc.)
    '''

    def splitter(text, step=16):
        return [text[step * k : step * (k + 1)] for k in range(len(text) // step)]

    def merge(lines):
        if not lines:
            return ""
        if len(lines[0]) > 80:
            return "".join(lines).replace(" ", "")
        return "".join(lines)

    with open(filename, "r") as f:
        EQDSK = f.read().splitlines()

    g = {}

    # --- Header line ---
    g["CASE"] = np.array(splitter(EQDSK[0][:48], 8))
    try:
        tmp = [x for x in EQDSK[0][48:].split() if x]
        _idum, g["NW"], g["NH"] = (int(x) for x in tmp[:3])
    except ValueError:
        _idum = int(EQDSK[0][48:52])
        g["NW"] = int(EQDSK[0][52:56])
        g["NH"] = int(EQDSK[0][56:60])
    offset = 1

    # --- 20 scalar values (4 lines x 5 values) ---
    scalars = list(map(float, splitter(merge(EQDSK[offset : offset + 4]))))
    (
        g["RDIM"], g["ZDIM"], g["RCENTR"], g["RLEFT"], g["ZMID"],
        g["RMAXIS"], g["ZMAXIS"], g["SIMAG"], g["SIBRY"], g["BCENTR"],
        g["CURRENT"], g["SIMAG"], _, g["RMAXIS"], _,
        g["ZMAXIS"], _, g["SIBRY"], _, _,
    ) = scalars
    offset += 4

    NW, NH = int(g["NW"]), int(g["NH"])
    nlNW = int(np.ceil(NW / 5.0))

    # --- 1-D profiles (each NW long) ---
    for name in ("FPOL", "PRES", "FFPRIM", "PPRIME"):
        g[name] = np.array(list(map(float, splitter(merge(EQDSK[offset : offset + nlNW])))))
        offset += nlNW

    # --- 2-D poloidal flux: PSIRZ (NH x NW) ---
    try:
        nlNWNH = int(np.ceil(NW * NH / 5.0))
        flat = np.fromiter(splitter(merge(EQDSK[offset : offset + nlNWNH])),
                           dtype=np.float64)[: NH * NW]
        g["PSIRZ"] = flat.reshape((NH, NW))
        offset += nlNWNH
    except ValueError:
        nlNWNH = NH * nlNW
        flat = np.fromiter(splitter(merge(EQDSK[offset : offset + nlNWNH])),
                           dtype=np.float64)[: NH * NW]
        g["PSIRZ"] = flat.reshape((NH, NW))
        offset += nlNWNH

    # --- Safety factor ---
    g["QPSI"] = np.array(list(map(float, splitter(merge(EQDSK[offset : offset + nlNW])))))
    offset += nlNW

    # --- Boundary and limiter ---
    if len(EQDSK) > offset + 1:
        parts = [x for x in EQDSK[offset].split() if x]
        g["NBBBS"] = int(parts[0])
        g["LIMITR"] = int(parts[1])
        offset += 1

        nlNBBBS = int(np.ceil(g["NBBBS"] * 2 / 5.0))
        bnd_vals = list(map(float, splitter(merge(EQDSK[offset : offset + nlNBBBS]))))
        g["RBBBS"] = np.array(bnd_vals[0::2])[: g["NBBBS"]]
        g["ZBBBS"] = np.array(bnd_vals[1::2])[: g["NBBBS"]]
        offset += max(nlNBBBS, 1)

        try:
            nlLIMITR = int(np.ceil(g["LIMITR"] * 2 / 5.0))
            lim_vals = list(map(float, splitter(merge(EQDSK[offset : offset + nlLIMITR]))))
            g["RLIM"] = np.array(lim_vals[0::2])[: g["LIMITR"]]
            g["ZLIM"] = np.array(lim_vals[1::2])[: g["LIMITR"]]
            offset += nlLIMITR
        except ValueError:
            g["LIMITR"] = 5
            dd = g["RDIM"] / 10.0
            R = np.linspace(0, g["RDIM"], 2) + g["RLEFT"]
            Z = np.linspace(0, g["ZDIM"], 2) - g["ZDIM"] / 2.0 + g["ZMID"]
            rmin = max(R[0], np.min(g["RBBBS"]) - dd)
            rmax = min(R[1], np.max(g["RBBBS"]) + dd)
            zmin = max(Z[0], np.min(g["ZBBBS"]) - dd)
            zmax = min(Z[1], np.max(g["ZBBBS"]) + dd)
            g["RLIM"] = np.array([rmin, rmax, rmax, rmin, rmin])
            g["ZLIM"] = np.array([zmin, zmin, zmax, zmax, zmin])
    else:
        g["NBBBS"] = 0
        g["LIMITR"] = 0
        g["RBBBS"] = np.array([])
        g["ZBBBS"] = np.array([])
        g["RLIM"] = np.array([])
        g["ZLIM"] = np.array([])

    # --- Optional extended data (RHOVN, PCURRT, etc.) ---
    try:
        parts = [float(x) for x in EQDSK[offset].split() if x]
        g["KVTOR"], g["RVTOR"], g["NMASS"] = parts[:3]
        offset += 1

        if g["KVTOR"] > 0:
            g["PRESSW"] = np.array(list(map(float, splitter(merge(EQDSK[offset : offset + nlNW])))))
            offset += nlNW
            g["PWPRIM"] = np.array(list(map(float, splitter(merge(EQDSK[offset : offset + nlNW])))))
            offset += nlNW

        if g["NMASS"] > 0:
            g["DMION"] = np.array(list(map(float, splitter(merge(EQDSK[offset : offset + nlNW])))))
            offset += nlNW

        g["RHOVN"] = np.array(list(map(float, splitter(merge(EQDSK[offset : offset + nlNW])))))
        offset += nlNW
    except Exception:
        pass

    if "RHOVN" not in g or len(g.get("RHOVN", [])) == 0:
        g["RHOVN"] = np.sqrt(np.linspace(0, 1, NW))

    if not np.any(g["PRES"]):
        pres = integrate.cumulative_trapezoid(
            g["PPRIME"],
            np.linspace(g["SIMAG"], g["SIBRY"], len(g["PPRIME"])),
            initial=0,
        )
        g["PRES"] = pres - pres[-1]

    return g


# ---------------------------------------------------------------------------
# GEQDSK writer
# ---------------------------------------------------------------------------

def _write_fortran_block(stream, values):
    r'''! Write *values* in standard GEQDSK format (5 per line, 16 chars each)

    @param stream Open text stream to write to
    @param values 1-D sequence of floats to serialise
    '''
    for i, v in enumerate(values):
        stream.write(f"{v:16.9E}")
        if (i + 1) % 5 == 0:
            stream.write("\n")
    if len(values) % 5 != 0:
        stream.write("\n")


def _write_geqdsk_to_stream(g, stream):
    r'''! Serialise a raw g-file dict *g* to an open text *stream*

    @param g Raw g-file dict (as produced by `_read_geqdsk`)
    @param stream Open text stream to write to
    '''
    NW = int(g["NW"])
    NH = int(g["NH"])

    case_arr = g.get("CASE", np.array([" " * 8] * 6))
    case_str = "".join(c.ljust(8)[:8] for c in case_arr)[:48]
    stream.write(f"{case_str}{0:4d}{NW:4d}{NH:4d}\n")

    scalars = [
        g["RDIM"], g["ZDIM"], g["RCENTR"], g["RLEFT"], g["ZMID"],
        g["RMAXIS"], g["ZMAXIS"], g["SIMAG"], g["SIBRY"], g["BCENTR"],
        g["CURRENT"], g["SIMAG"], 0.0, g["RMAXIS"], 0.0,
        g["ZMAXIS"], 0.0, g["SIBRY"], 0.0, 0.0,
    ]
    _write_fortran_block(stream, scalars)

    for name in ("FPOL", "PRES", "FFPRIM", "PPRIME"):
        _write_fortran_block(stream, g[name])

    _write_fortran_block(stream, g["PSIRZ"].ravel())
    _write_fortran_block(stream, g["QPSI"])

    nbbbs = int(g.get("NBBBS", len(g.get("RBBBS", []))))
    nlim = int(g.get("LIMITR", len(g.get("RLIM", []))))
    stream.write(f" {nbbbs:5d} {nlim:5d}\n")

    if nbbbs > 0:
        bnd = np.empty(2 * nbbbs)
        bnd[0::2] = g["RBBBS"][:nbbbs]
        bnd[1::2] = g["ZBBBS"][:nbbbs]
        _write_fortran_block(stream, bnd)

    if nlim > 0:
        lim = np.empty(2 * nlim)
        lim[0::2] = g["RLIM"][:nlim]
        lim[1::2] = g["ZLIM"][:nlim]
        _write_fortran_block(stream, lim)

    # --- Optional extended blocks (KVTOR/RVTOR/NMASS, PRESSW, PWPRIM,
    # DMION, RHOVN).  Emitted only when the dict contains them, so
    # round-tripping an extended g-file produced by `_read_geqdsk` is
    # lossless.  Standard minimal g-files are unaffected.
    has_kvtor = "KVTOR" in g and "RVTOR" in g and "NMASS" in g
    if has_kvtor:
        kvtor = int(g["KVTOR"])
        rvtor = float(g["RVTOR"])
        nmass = int(g["NMASS"])
        stream.write(f" {kvtor:5d} {rvtor:16.9E} {nmass:5d}\n")

        if kvtor > 0 and "PRESSW" in g and "PWPRIM" in g:
            _write_fortran_block(stream, g["PRESSW"])
            _write_fortran_block(stream, g["PWPRIM"])

        if nmass > 0 and "DMION" in g:
            _write_fortran_block(stream, g["DMION"])

        if "RHOVN" in g and len(g["RHOVN"]) > 0:
            _write_fortran_block(stream, g["RHOVN"])


def write_geqdsk(g, filename):
    r'''! Write a raw g-file dict to *filename*

    @param g Raw g-file dict (as produced by `_read_geqdsk`)
    @param filename Output path for the g-file
    '''
    with open(filename, "w") as f:
        _write_geqdsk_to_stream(g, f)


# ---------------------------------------------------------------------------
# Contour tracing
# ---------------------------------------------------------------------------

def _trace_contours(R, Z, PSI, levels):
    r'''! Extract contour lines of \f$\psi\f$ at given levels using contourpy

    @param R 1-D array of R grid coordinates [m]
    @param Z 1-D array of Z grid coordinates [m]
    @param PSI 2-D array (`len(Z)` x `len(R)`) of poloidal flux on the grid [Wb]
    @param levels 1-D array of \f$\psi\f$ values at which to extract contours
    @result List (one per level) of lists of `(N,2)` arrays of `(R,Z)` points
    '''
    cg = contourpy.contour_generator(R, Z, PSI, name="serial", line_type="Separate")
    all_contours = []
    for lev in levels:
        segments = cg.lines(lev)
        all_contours.append(segments)
    return all_contours


def _select_main_contour(segments, R0, Z0, sigma_RpZ, sigma_rhotp):
    r'''! Select the main closed contour encircling the magnetic axis

    Uses the winding / angular-coverage criterion: the contour whose
    double integral of angle vs. arc-fraction has the largest amplitude
    is most likely the one that wraps around the axis.

    @param segments List of `(N,2)` candidate contour segments at a single \f$\psi\f$ level
    @param R0 Magnetic axis R [m]
    @param Z0 Magnetic axis Z [m]
    @param sigma_RpZ COCOS right-handedness sign
    @param sigma_rhotp COCOS \f$\theta_{\mathrm{pol}}\cdot\phi_{\mathrm{tor}}\f$ sign
    @result `(N,2)` ndarray of the selected contour with correct orientation, or `None`
    '''
    if not segments:
        return None

    sign_theta = sigma_RpZ * sigma_rhotp
    best = None
    best_score = -1.0

    for seg in segments:
        r, z = seg[:, 0], seg[:, 1]
        if len(r) < 4 or np.any(np.isnan(r * z)):
            continue

        r = r.copy()
        z = z.copy()
        r[0] = r[-1] = 0.5 * (r[0] + r[-1])
        z[0] = z[-1] = 0.5 * (z[0] + z[-1])

        theta = np.unwrap(np.arctan2(z - Z0, r - R0))
        theta -= np.mean(theta)
        s = np.linspace(0, 1, len(theta))

        score = np.max(np.abs(
            integrate.cumulative_trapezoid(
                integrate.cumulative_trapezoid(theta, s, initial=0), s
            )
        ))

        if score > best_score:
            best_score = score
            orient = int(np.sign(
                (z[0] - Z0) * (r[1] - r[0]) - (r[0] - R0) * (z[1] - z[0])
            ))
            if orient != 0:
                best = seg[::sign_theta * orient, :]
            else:
                best = seg

    return best


def _crop_at_xpoint(seg, R0, Z0):
    r'''! Crop an open separatrix contour at the X-point to produce a closed LCFS

    @param seg `(N,2)` open contour segment with columns `(R, Z)`
    @param R0 Magnetic axis R [m]
    @param Z0 Magnetic axis Z [m]
    @result `(M,2)` closed contour (first point == last point)
    '''
    r, z = seg[:, 0], seg[:, 1]
    n = len(r)

    gap = np.sqrt((r[0] - r[-1])**2 + (z[0] - z[-1])**2)
    perimeter = np.sum(np.sqrt(np.diff(r)**2 + np.diff(z)**2))
    if gap < 0.01 * perimeter:
        out = seg.copy()
        out[-1] = out[0]
        return out

    r_xpt = 0.5 * (r[0] + r[-1])
    z_xpt = 0.5 * (z[0] + z[-1])

    dist = np.sqrt((r - r_xpt)**2 + (z - z_xpt)**2)

    from scipy.ndimage import uniform_filter1d
    dist_smooth = uniform_filter1d(dist, size=max(5, n // 50))

    local_min = np.zeros(n, dtype=bool)
    for i in range(1, n - 1):
        if dist_smooth[i] < dist_smooth[i-1] and dist_smooth[i] < dist_smooth[i+1]:
            local_min[i] = True

    min_indices = np.where(local_min)[0]

    if len(min_indices) < 2:
        out = seg.copy()
        out[-1] = out[0]
        return out

    sorted_by_dist = min_indices[np.argsort(dist[min_indices])]

    idx1 = sorted_by_dist[0]
    idx2 = None
    for candidate in sorted_by_dist[1:]:
        separation = abs(candidate - idx1)
        if separation > n // 10:
            idx2 = candidate
            break

    if idx2 is None:
        out = seg.copy()
        out[-1] = out[0]
        return out

    if idx1 > idx2:
        idx1, idx2 = idx2, idx1

    plasma_seg = seg[idx1:idx2+1]

    theta = np.unwrap(np.arctan2(plasma_seg[:, 1] - Z0,
                                  plasma_seg[:, 0] - R0))
    winding = abs(theta[-1] - theta[0])

    if winding < np.pi:
        plasma_seg = np.vstack([seg[idx2:], seg[:idx1+1]])

    closed = np.vstack([plasma_seg, plasma_seg[:1]])
    return closed


# ---------------------------------------------------------------------------
# Flux-surface geometry helper
# ---------------------------------------------------------------------------

def _flux_geometry(R, Z):
    r'''! Compute geometric properties of a single flux-surface contour

    @param R 1-D array of contour R coordinates (closed: first == last, or nearly so)
    @param Z 1-D array of contour Z coordinates (closed: first == last, or nearly so)
    @result Dict with keys `R` (geometric center), `Z`, `a` (minor radius),
            `kappa`, `kapu`, `kapl`, `delta`, `delu`, `dell`, `perimeter`,
            `surfArea`, `eps`
    '''
    geo = {}

    def _parabola_extremum(idx, main, other):
        n = len(main)
        if n < 3:
            return other[idx], main[idx]

        im = (idx - 1) % (n - 1) if main[0] == main[-1] else max(idx - 1, 0)
        ip = (idx + 1) % (n - 1) if main[0] == main[-1] else min(idx + 1, n - 1)

        ym1, y0, yp1 = main[im], main[idx], main[ip]
        xm1, x0, xp1 = other[im], other[idx], other[ip]

        denom = 2.0 * (ym1 - 2.0 * y0 + yp1)
        if abs(denom) < 1e-30:
            return x0, y0
        frac = (ym1 - yp1) / denom
        frac = np.clip(frac, -0.5, 0.5)
        x_ext = x0 + frac * 0.5 * (xp1 - xm1)
        y_ext = y0 + frac * 0.5 * (yp1 - ym1)
        return x_ext, y_ext

    imaxr = np.argmax(R)
    z_at_max_r, max_r = _parabola_extremum(imaxr, R, Z)

    iminr = np.argmin(R)
    z_at_min_r, min_r = _parabola_extremum(iminr, R, Z)

    imaxz = np.argmax(Z)
    r_at_max_z, max_z = _parabola_extremum(imaxz, Z, R)

    iminz = np.argmin(Z)
    r_at_min_z, min_z = _parabola_extremum(iminz, Z, R)

    dR = np.diff(R, append=R[0])
    dZ = np.diff(Z, append=Z[0])
    dl_segs = np.sqrt(dR**2 + dZ**2)
    dl = 0.5 * (dl_segs + np.roll(dl_segs, 1))

    geo["R"] = 0.5 * (max_r + min_r)
    geo["Z"] = 0.5 * (max_z + min_z)
    geo["a"] = 0.5 * (max_r - min_r)
    if geo["a"] > 0:
        geo["kappa"] = 0.5 * (max_z - min_z) / geo["a"]
        geo["kapu"] = (max_z - z_at_max_r) / geo["a"]
        geo["kapl"] = (z_at_max_r - min_z) / geo["a"]
        geo["delu"] = (geo["R"] - r_at_max_z) / geo["a"]
        geo["dell"] = (geo["R"] - r_at_min_z) / geo["a"]
        geo["delta"] = 0.5 * (geo["delu"] + geo["dell"])
    else:
        geo["kappa"] = 1.0
        geo["kapu"] = 1.0
        geo["kapl"] = 1.0
        geo["delu"] = 0.0
        geo["dell"] = 0.0
        geo["delta"] = 0.0
    geo["perimeter"] = np.sum(dl)
    geo["surfArea"] = 2 * np.pi * np.sum(R * dl)
    geo["eps"] = geo["a"] / geo["R"] if geo["R"] > 0 else 0.0

    return geo


def _resample_contour(R, Z, npts=257, periodic=True):
    r'''! Resample a closed contour to *npts* equally-spaced-in-arc-length points

    @param R 1-D array of contour R coordinates (nearly closed: first \f$\approx\f$ last)
    @param Z 1-D array of contour Z coordinates (nearly closed: first \f$\approx\f$ last)
    @param npts Number of output points (including the repeated closing point)
    @param periodic If `True` (default), use a periodic cubic spline. Set to `False`
                    for the separatrix surface where the X-point cusp produces a
                    derivative discontinuity.
    @result Tuple `(R_new, Z_new)` of 1-D arrays of length *npts*
    '''
    R = np.asarray(R, dtype=float).copy()
    Z = np.asarray(Z, dtype=float).copy()
    R[-1], Z[-1] = R[0], Z[0]

    ds = np.sqrt(np.diff(R) ** 2 + np.diff(Z) ** 2)
    s = np.empty(len(R))
    s[0] = 0.0
    s[1:] = np.cumsum(ds)
    s_total = s[-1]

    if s_total < 1e-14:
        return R, Z

    tck_R = interpolate.splrep(s, R, k=3, per=periodic)
    tck_Z = interpolate.splrep(s, Z, k=3, per=periodic)

    s_new = np.linspace(0, s_total, npts)
    R_new = interpolate.splev(s_new, tck_R)
    Z_new = interpolate.splev(s_new, tck_Z)

    if not periodic:
        R_new[-1], Z_new[-1] = R_new[0], Z_new[0]

    return R_new, Z_new


def _detect_xpoint_angle(R, Z, R0, Z0):
    r'''! Detect the X-point angle on a separatrix contour

    Finds the location of the sharpest bend (minimum cosine of angle between
    adjacent tangent vectors) and returns \f$\theta = \mathrm{atan2}(Z-Z_0, R-R_0)\f$
    at that point.

    @param R 1-D array of separatrix R coordinates (nearly closed)
    @param Z 1-D array of separatrix Z coordinates (nearly closed)
    @param R0 Magnetic axis R [m]
    @param Z0 Magnetic axis Z [m]
    @result Poloidal angle \f$\theta\f$ of the X-point, or `None` if detection fails
    '''
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)
    n = len(R)
    if n < 6:
        return None

    n1 = n - 1
    dR_a = np.gradient(R[:n1])
    dZ_a = np.gradient(Z[:n1])
    dR_b = np.gradient(np.concatenate([R[1:n1], R[0:1]]))
    dZ_b = np.gradient(np.concatenate([Z[1:n1], Z[0:1]]))
    dot = dR_a * dR_b + dZ_a * dZ_b
    norm = np.sqrt(dR_a**2 + dZ_a**2) * np.sqrt(dR_b**2 + dZ_b**2)
    norm = np.maximum(norm, 1e-30)
    cos_angle = dot / norm
    idx_xpt = np.argmin(cos_angle)

    return float(np.arctan2(Z[idx_xpt] - Z0, R[idx_xpt] - R0))


def _resample_contour_theta(R, Z, R0, Z0, npts=257, theta_xpt=None):
    r'''! Resample a closed contour to *npts* equally-spaced-in-angle points

    Fits periodic cubic splines \f$R(\theta)\f$ and \f$Z(\theta)\f$ where
    \f$\theta = \mathrm{atan2}(Z-Z_0, R-R_0)\f$, sorted starting from the
    X-point angle.

    @param R 1-D array of contour R coordinates (nearly closed)
    @param Z 1-D array of contour Z coordinates (nearly closed)
    @param R0 Magnetic axis R [m]
    @param Z0 Magnetic axis Z [m]
    @param npts Number of output points (including the repeated closing point)
    @param theta_xpt Pre-computed X-point angle; if `None`, detected from this contour
    @result Tuple `(R_new, Z_new)` of 1-D arrays of length *npts*
    '''
    R = np.asarray(R, dtype=float).copy()
    Z = np.asarray(Z, dtype=float).copy()
    R[-1], Z[-1] = R[0], Z[0]
    n = len(R)

    if n < 6:
        return R, Z

    n1 = n - 1

    if theta_xpt is None:
        theta_xpt = _detect_xpoint_angle(R, Z, R0, Z0)
        if theta_xpt is None:
            return R, Z

    t_raw = np.arctan2(Z - Z0, R - R0)

    t_rel = (t_raw[:n1] - theta_xpt) % (2.0 * np.pi)
    idx = np.argsort(t_rel)
    idx = np.concatenate([idx, idx[0:1]])

    t_sorted = np.unwrap(t_raw[idx])
    R_sorted = R[idx]
    Z_sorted = Z[idx]

    if t_sorted[0] > t_sorted[1]:
        t_sorted = -t_sorted

    t_span = t_sorted[-1] - t_sorted[0]
    if abs(t_span) < 1e-6:
        return R, Z

    theta0 = np.linspace(0, 2 * np.pi, npts) + theta_xpt
    idx_eval = np.argsort((theta0[:-1] + theta_xpt) % (2 * np.pi) - theta_xpt)
    idx_eval = np.concatenate([idx_eval, idx_eval[0:1]])
    theta_eval = np.unwrap(theta0[idx_eval])
    if t_sorted[0] > 0 and theta_eval[0] < 0:
        theta_eval = -theta_eval
    elif t_sorted[0] < 0 and theta_eval[0] > 0:
        theta_eval = -theta_eval

    tck_R = interpolate.splrep(t_sorted, R_sorted, k=3, per=True)
    tck_Z = interpolate.splrep(t_sorted, Z_sorted, k=3, per=True)

    t_range = t_sorted[-1] - t_sorted[0]
    theta_wrapped = (theta_eval - t_sorted[0]) % t_range + t_sorted[0]
    R_new = interpolate.splev(theta_wrapped, tck_R, ext=0)
    Z_new = interpolate.splev(theta_wrapped, tck_Z, ext=0)

    R_new[-1], Z_new[-1] = R_new[0], Z_new[0]
    return R_new, Z_new


# ---------------------------------------------------------------------------
# GEQDSKEquilibrium class
# ---------------------------------------------------------------------------

class GEQDSKEquilibrium:
    r'''! GEQDSK equilibrium with flux-surface-averaged quantities

    Reads a standard g-file and lazily computes magnetic fields, flux-surface
    contours, averages, geometry, and inductance.
    '''

    def __init__(self, filename, cocos=1, nlevels=None, resample="theta",
                 extrapolate_edge=True):
        r'''! Initialize GEQDSKEquilibrium by reading a g-file from disk

        @param filename Path to the g-file
        @param cocos COCOS convention index (default 1, standard EFIT)
        @param nlevels Number of normalised-\f$\psi\f$ levels for flux-surface
                       analysis (defaults to `NW` from the g-file)
        @param resample Contour resampling method for near-separatrix surfaces
                        (\f$\hat{\psi} \geq 0.99\f$): `"theta"` or `"arc_length"`
        @param extrapolate_edge If `True` (default), extrapolate \f$p'\f$ and
                                \f$FF'\f$ at the separatrix when the g-file has
                                them forced to zero
        '''
        self._raw = _read_geqdsk(filename)
        self._cocos_index = int(cocos)
        self._cocos = _cocos_params(cocos)
        self._nlevels = nlevels if nlevels is not None else int(self._raw["NW"])
        if resample not in ("theta", "arc_length"):
            raise ValueError(
                f"resample must be 'theta' or 'arc_length', got {resample!r}"
            )
        self._resample_method = resample
        self._extrapolate_edge = bool(extrapolate_edge)
        self._cache = {}

    @classmethod
    def from_bytes(cls, raw_bytes, cocos=1, nlevels=None, resample="theta",
                   extrapolate_edge=True):
        r'''! Construct from in-memory bytes (e.g. from HDF5 storage)

        @param raw_bytes Raw content of a g-file
        @param cocos COCOS convention index
        @param nlevels Number of \f$\hat{\psi}\f$ levels
        @param resample Contour resampling method (`"theta"` or `"arc_length"`)
        @param extrapolate_edge Extrapolate \f$p'\f$ and \f$FF'\f$ at separatrix
                                when forced to zero
        @result New `GEQDSKEquilibrium` instance
        '''
        import os as _os
        with tempfile.NamedTemporaryFile(mode="wb", suffix=".geqdsk",
                                         delete=False) as tmp:
            tmp.write(raw_bytes)
            tmp_path = tmp.name
        try:
            return cls(tmp_path, cocos=cocos, nlevels=nlevels, resample=resample,
                       extrapolate_edge=extrapolate_edge)
        finally:
            try:
                _os.remove(tmp_path)
            except OSError:
                pass

    @classmethod
    def from_raw(cls, raw_dict, cocos=1, nlevels=None, resample="theta",
                 extrapolate_edge=True):
        r'''! Construct directly from a raw g-file dict (no file I/O)

        @param raw_dict Dict with standard GEQDSK keys (as returned by `_read_geqdsk`)
        @param cocos COCOS convention index
        @param nlevels Number of \f$\hat{\psi}\f$ levels
        @param resample Contour resampling method (`"theta"` or `"arc_length"`)
        @param extrapolate_edge Extrapolate \f$p'\f$ and \f$FF'\f$ at separatrix
                                when forced to zero
        @result New `GEQDSKEquilibrium` instance
        '''
        obj = object.__new__(cls)
        obj._raw = {k: (v.copy() if isinstance(v, np.ndarray) else v)
                     for k, v in raw_dict.items()}
        obj._cocos_index = int(cocos)
        obj._cocos = _cocos_params(cocos)
        obj._nlevels = nlevels if nlevels is not None else int(obj._raw["NW"])
        if resample not in ("theta", "arc_length"):
            raise ValueError(
                f"resample must be 'theta' or 'arc_length', got {resample!r}"
            )
        obj._resample_method = resample
        obj._extrapolate_edge = bool(extrapolate_edge)
        obj._cache = {}
        return obj

    # --- Raw data properties ---

    @property
    def R_grid(self):
        r'''! 1-D R grid [m]'''
        if "R_grid" not in self._cache:
            NW = int(self._raw["NW"])
            self._cache["R_grid"] = np.linspace(
                self._raw["RLEFT"],
                self._raw["RLEFT"] + self._raw["RDIM"],
                NW,
            )
        return self._cache["R_grid"]

    @property
    def Z_grid(self):
        r'''! 1-D Z grid [m]'''
        if "Z_grid" not in self._cache:
            NH = int(self._raw["NH"])
            self._cache["Z_grid"] = np.linspace(
                self._raw["ZMID"] - self._raw["ZDIM"] / 2.0,
                self._raw["ZMID"] + self._raw["ZDIM"] / 2.0,
                NH,
            )
        return self._cache["Z_grid"]

    @property
    def psi_RZ(self):
        r'''! 2-D poloidal flux array (`NH` x `NW`) [Wb]'''
        return self._raw["PSIRZ"]

    @property
    def psi_N_RZ(self):
        r'''! 2-D normalised poloidal flux \f$\hat{\psi}(R,Z)\f$ on the (R,Z) grid'''
        return (self._raw["PSIRZ"] - self.psi_axis) / (self.psi_boundary - self.psi_axis)

    @property
    def psi_axis(self):
        r'''! Poloidal flux at the magnetic axis [Wb]'''
        return self._raw["SIMAG"]

    @property
    def psi_boundary(self):
        r'''! Poloidal flux at the last closed flux surface [Wb]'''
        return self._raw["SIBRY"]

    @property
    def psi_N(self):
        r'''! Normalised raw g-file \f$\hat{\psi}\f$ grid for 1-D profiles (`NW` points)

        Matches the length and spacing of `fpol`, `pres`, `pprime`, `ffprim`,
        `qpsi`, and `rhovn`.  For the flux-surface-analysis grid (used by the
        derived FSA quantities) see `psi_N_levels`.
        '''
        return np.linspace(0, 1, int(self._raw["NW"]))

    @property
    def psi_N_levels(self):
        r'''! Normalised \f$\hat{\psi}\f$ levels used for flux-surface analysis (`nlevels` points)

        This is the grid on which lazy FSA quantities (`q_profile`,
        `j_tor_averaged`, `geometry[...]`, `averages[...]`, `midplane[...]`,
        etc.) are evaluated.  Defaults to `NW` (identical to `psi_N`) unless
        a different `nlevels` was passed at construction.
        '''
        return np.linspace(0, 1, self._nlevels)

    @property
    def fpol(self):
        r'''! \f$F = R B_t\f$ poloidal current function, on the raw `psi_N` grid (`NW` points) [T m]'''
        return self._raw["FPOL"]

    @property
    def pres(self):
        r'''! Pressure profile on the raw `psi_N` grid (`NW` points) [Pa]'''
        return self._raw["PRES"]

    @property
    def pprime(self):
        r'''! \f$dP/d\psi\f$ on the raw `psi_N` grid (`NW` points) [Pa/Wb]'''
        return self._raw["PPRIME"]

    @property
    def ffprim(self):
        r'''! \f$F\,dF/d\psi\f$ on the raw `psi_N` grid (`NW` points) [T\f$^2\f$ m\f$^2\f$/Wb]'''
        return self._raw["FFPRIM"]

    @property
    def qpsi(self):
        r'''! Safety factor from the g-file, on the raw `psi_N` grid (`NW` points)'''
        return self._raw["QPSI"]

    @property
    def Ip(self):
        r'''! Plasma current \f$I_p\f$ [A]'''
        return self._raw["CURRENT"]

    @property
    def R_mag(self):
        r'''! R of the magnetic axis [m]'''
        return self._raw["RMAXIS"]

    @property
    def Z_mag(self):
        r'''! Z of the magnetic axis [m]'''
        return self._raw["ZMAXIS"]

    @property
    def R_center(self):
        r'''! Reference geometric center R [m]'''
        return self._raw["RCENTR"]

    @property
    def B_center(self):
        r'''! Vacuum toroidal field at `R_center` [T]'''
        return self._raw["BCENTR"]

    @property
    def boundary_R(self):
        r'''! R coordinates of the plasma boundary [m]'''
        return self._raw["RBBBS"]

    @property
    def boundary_Z(self):
        r'''! Z coordinates of the plasma boundary [m]'''
        return self._raw["ZBBBS"]

    @property
    def limiter_R(self):
        r'''! R coordinates of the limiter [m]'''
        return self._raw["RLIM"]

    @property
    def limiter_Z(self):
        r'''! Z coordinates of the limiter [m]'''
        return self._raw["ZLIM"]

    @property
    def rhovn(self):
        r'''! Normalised toroidal flux coordinate \f$\rho = \sqrt{\Phi_{\mathrm{tor}}/\Phi_{\mathrm{tor,edge}}}\f$

        Always computed from the safety factor profile.
        '''
        NW = int(self._raw["NW"])
        psi_grid = np.linspace(self.psi_axis, self.psi_boundary, NW)
        phi = integrate.cumulative_trapezoid(self._raw["QPSI"], psi_grid, initial=0)
        phi_edge = phi[-1]
        if abs(phi_edge) < 1e-30:
            return np.sqrt(np.linspace(0, 1, NW))
        return np.sqrt(np.abs(phi / phi_edge))

    @property
    def cocos(self):
        r'''! Current COCOS convention index'''
        return self._cocos_index

    # --- COCOS conversion and sign-flip methods ---

    def cocosify(self, cocos_out, copy=False):
        r'''! Convert the raw g-file data from the current COCOS to *cocos_out*

        Applies the multiplicative sign / \f$2\pi\f$ factors to every relevant
        field following Sauter & Medvedev, Comput. Phys. Commun. 184 (2013) 293,
        Eq. 14/23.

        @param cocos_out Target COCOS convention index (1-8 or 11-18)
        @param copy If `True`, return a new `GEQDSKEquilibrium` with converted
                    data, leaving this object unchanged.  If `False` (default),
                    convert in place and return `self`.
        @result Converted object (`self` when *copy=False*)
        '''
        cc_in = _cocos_params(self._cocos_index)
        cc_out = _cocos_params(cocos_out)

        sigma_Bp_eff = cc_out["sigma_Bp"] * cc_in["sigma_Bp"]
        sigma_RpZ_eff = cc_out["sigma_RpZ"] * cc_in["sigma_RpZ"]
        sigma_rhotp_eff = cc_out["sigma_rhotp"] * cc_in["sigma_rhotp"]
        exp_Bp_eff = cc_out["exp_Bp"] - cc_in["exp_Bp"]

        twopi_exp = (2.0 * np.pi) ** exp_Bp_eff

        psi_factor = sigma_RpZ_eff * sigma_Bp_eff * twopi_exp
        dpsi_factor = sigma_RpZ_eff * sigma_Bp_eff / twopi_exp
        bt_factor = sigma_RpZ_eff
        ip_factor = sigma_RpZ_eff
        q_factor = sigma_rhotp_eff

        target = self._copy_for_mutation() if copy else self

        g = target._raw
        for key in ("SIMAG", "SIBRY"):
            g[key] *= psi_factor
        g["PSIRZ"] = g["PSIRZ"] * psi_factor
        g["PPRIME"] = g["PPRIME"] * dpsi_factor
        g["FFPRIM"] = g["FFPRIM"] * dpsi_factor
        g["FPOL"] = g["FPOL"] * bt_factor
        g["BCENTR"] *= bt_factor
        g["CURRENT"] *= ip_factor
        g["QPSI"] = g["QPSI"] * q_factor

        target._cocos_index = int(cocos_out)
        target._cocos = cc_out
        target._cache.clear()
        return target

    def flip_Bt_Ip(self, copy=False):
        r'''! Reverse the signs of \f$B_t\f$ and \f$I_p\f$ in the raw g-file data

        Negates `BCENTR`, `FPOL`, `CURRENT`, `PSIRZ`, `SIMAG`, `SIBRY`, `PPRIME`,
        and `FFPRIM` -- equivalent to flipping the direction of both the toroidal
        field and the plasma current while keeping the COCOS index unchanged.

        @param copy If `True`, return a new object; otherwise modify in place
        @result Modified `GEQDSKEquilibrium` object
        '''
        target = self._copy_for_mutation() if copy else self
        g = target._raw

        g["BCENTR"] *= -1
        g["FPOL"] = g["FPOL"] * -1
        g["CURRENT"] *= -1
        g["SIMAG"] *= -1
        g["SIBRY"] *= -1
        g["PSIRZ"] = g["PSIRZ"] * -1
        g["PPRIME"] = g["PPRIME"] * -1
        g["FFPRIM"] = g["FFPRIM"] * -1

        target._cache.clear()
        return target

    def _copy_for_mutation(self):
        return GEQDSKEquilibrium.from_raw(
            self._raw, cocos=self._cocos_index,
            nlevels=self._nlevels, resample=self._resample_method,
            extrapolate_edge=self._extrapolate_edge,
        )

    # --- Save / serialise ---

    def save(self, filename):
        r'''! Write the (possibly modified) g-file data to *filename*

        @param filename Output path for the g-file
        '''
        write_geqdsk(self._raw, filename)

    def to_bytes(self):
        r'''! Serialise to in-memory bytes (round-trips with `from_bytes`)

        @result Raw g-file content as bytes
        '''
        buf = io.StringIO()
        _write_geqdsk_to_stream(self._raw, buf)
        return buf.getvalue().encode("ascii")

    # --- Lazy field computation ---

    def _compute_fields(self):
        r'''! Compute \f$B_R\f$, \f$B_Z\f$, \f$J_t\f$ on the full (R,Z) grid'''
        if "Br" in self._cache:
            return

        cc = self._cocos
        R = self.R_grid
        Z = self.Z_grid
        PSI = self.psi_RZ
        RR, _ZZ = np.meshgrid(R, Z)

        dR = R[1] - R[0]
        dZ = Z[1] - Z[0]
        dPSIdZ, dPSIdR = np.gradient(PSI, dZ, dR)

        twopi_exp = (2.0 * np.pi) ** cc["exp_Bp"]
        self._cache["Br"] = cc["sigma_RpZ"] * cc["sigma_Bp"] * dPSIdZ / (RR * twopi_exp)
        self._cache["Bz"] = -cc["sigma_RpZ"] * cc["sigma_Bp"] * dPSIdR / (RR * twopi_exp)

        Br = self._cache["Br"]
        Bz = self._cache["Bz"]
        dBrdZ, dBrdR = np.gradient(Br, dZ, dR)
        _dBzdZ, dBzdR = np.gradient(Bz, dZ, dR)
        self._cache["Jt"] = cc["sigma_RpZ"] * (dBrdZ - dBzdR) / constants.mu_0

    # --- Lazy flux-surface tracing and averaging ---

    def _trace_surfaces(self):
        r'''! Trace flux surfaces and compute all averaged quantities'''
        if "avg" in self._cache:
            return

        self._compute_fields()

        cc = self._cocos
        R = self.R_grid
        Z = self.Z_grid

        psi_N_levels = self.psi_N_levels
        dpsi = self.psi_boundary - self.psi_axis
        psi_levels = psi_N_levels * dpsi + self.psi_axis

        Br_interp = interpolate.RectBivariateSpline(Z, R, self._cache["Br"])
        Bz_interp = interpolate.RectBivariateSpline(Z, R, self._cache["Bz"])
        Jt_interp = interpolate.RectBivariateSpline(Z, R, self._cache["Jt"])

        NW = int(self._raw["NW"])
        psi_N_raw = np.linspace(0, 1, NW)
        F_interp = interpolate.InterpolatedUnivariateSpline(psi_N_raw, self._raw["FPOL"])

        pprime_raw = self._raw["PPRIME"].copy()
        ffprim_raw = self._raw["FFPRIM"].copy()
        if self._extrapolate_edge:
            for prof in (pprime_raw, ffprim_raw):
                if (
                    NW >= 4
                    and prof[-1] == 0.0
                    and abs(prof[-2]) > 1e-30
                ):
                    prof[-1] = 3.0 * prof[-2] - 3.0 * prof[-3] + prof[-4]
                    if np.sign(prof[-1]) != np.sign(prof[-2]):
                        prof[-1] = 0.0
        pprime_interp = interpolate.InterpolatedUnivariateSpline(psi_N_raw, pprime_raw)
        ffprim_interp = interpolate.InterpolatedUnivariateSpline(psi_N_raw, ffprim_raw)

        contours = _trace_contours(R, Z, self.psi_RZ, psi_levels)

        nc = len(psi_N_levels)
        R0 = self.R_mag
        Z0 = self.Z_mag

        avg = {key: np.zeros(nc) for key in [
            "R", "1/R", "1/R**2", "R**2",
            "Bp", "Bp**2", "Bt", "Bt**2", "Btot**2",
            "Jt", "Jt/R", "Jt/R_num",
            "vp", "q", "ip",
            "F", "PPRIME", "FFPRIM",
        ]}
        geo_arrays = {key: np.zeros(nc) for key in [
            "R", "Z", "a", "kappa", "kapu", "kapl",
            "delta", "delu", "dell", "perimeter", "surfArea", "eps",
        ]}
        contour_data = []

        dpsi_arr = np.abs(np.gradient(psi_levels))

        Bp2_vol = 0.0
        for k in range(nc):
            pn = psi_N_levels[k]
            F_k = float(F_interp(pn))
            pprime_k = float(pprime_interp(pn))
            ffprim_k = float(ffprim_interp(pn))

            avg["F"][k] = F_k
            avg["PPRIME"][k] = pprime_k
            avg["FFPRIM"][k] = ffprim_k

            if pn == 0:
                contour_data.append(np.array([[R0, Z0]]))
                avg["R"][k] = R0
                avg["1/R"][k] = 1.0 / R0
                avg["1/R**2"][k] = 1.0 / R0**2
                avg["R**2"][k] = R0**2
                Bt_axis = F_k / R0
                avg["Bp"][k] = 0.0
                avg["Bp**2"][k] = 0.0
                avg["Bt"][k] = Bt_axis
                avg["Bt**2"][k] = Bt_axis**2
                avg["Btot**2"][k] = Bt_axis**2
                Jt0 = float(Jt_interp.ev(Z0, R0))
                avg["Jt"][k] = Jt0
                avg["Jt/R_num"][k] = Jt0 / R0
                avg["vp"][k] = 0.0
                avg["ip"][k] = 0.0
                for gk in geo_arrays:
                    geo_arrays[gk][k] = 0.0
                geo_arrays["R"][k] = R0
                geo_arrays["Z"][k] = Z0
                continue

            if pn == 1.0 and len(self.boundary_R) >= 4:
                seg = np.column_stack([self.boundary_R, self.boundary_Z])
                if np.sqrt((seg[0,0]-seg[-1,0])**2 + (seg[0,1]-seg[-1,1])**2) > 1e-6:
                    seg = np.vstack([seg, seg[:1]])
            else:
                seg = _select_main_contour(
                    contours[k], R0, Z0,
                    cc["sigma_RpZ"], cc["sigma_rhotp"],
                )
                if seg is not None and len(seg) >= 10:
                    seg = _crop_at_xpoint(seg, R0, Z0)

            if seg is None or len(seg) < 4:
                contour_data.append(np.array([[R0, Z0]]))
                avg["R"][k] = R0
                avg["1/R"][k] = 1.0 / R0
                avg["1/R**2"][k] = 1.0 / R0**2
                avg["R**2"][k] = R0**2
                continue

            if pn == 1.0:
                r_s, z_s = seg[:, 0].copy(), seg[:, 1].copy()
                r_s[-1], z_s[-1] = r_s[0], z_s[0]
            elif pn >= 0.99 and len(seg) >= 20:
                # Near-separatrix: choose theta-based vs arc-length
                # resampling per self._resample_method.  Theta-based
                # resampling places points uniformly in angle around the
                # magnetic axis, which preserves the X-point cusp better
                # than arc-length for separatrix-adjacent surfaces.
                if self._resample_method == "theta":
                    r_s, z_s = _resample_contour_theta(
                        seg[:, 0], seg[:, 1], R0, Z0, npts=257,
                    )
                else:
                    r_s, z_s = _resample_contour(
                        seg[:, 0], seg[:, 1], npts=257, periodic=False,
                    )
            elif len(seg) >= 20:
                r_s, z_s = _resample_contour(seg[:, 0], seg[:, 1], npts=257)
            else:
                r_s, z_s = seg[:, 0].copy(), seg[:, 1].copy()
                r_s[-1], z_s[-1] = r_s[0], z_s[0]
            contour_data.append(np.column_stack([r_s, z_s]))

            dR = np.diff(r_s, append=r_s[0])
            dZ = np.diff(z_s, append=z_s[0])
            dl_segs = np.sqrt(dR**2 + dZ**2)
            dl = 0.5 * (dl_segs + np.roll(dl_segs, 1))

            Br_s = Br_interp.ev(z_s, r_s)
            Bz_s = Bz_interp.ev(z_s, r_s)
            Jt_s = Jt_interp.ev(z_s, r_s)

            Bp2_s = Br_s**2 + Bz_s**2
            Bp_mod = np.sqrt(Bp2_s)

            signBp = (
                cc["sigma_rhotp"] * cc["sigma_RpZ"]
                * np.sign((z_s - Z0) * Br_s - (r_s - R0) * Bz_s)
            )
            Bp_signed = signBp * Bp_mod

            Bt_s = F_k / r_s
            B2_s = Bp2_s + Bt_s**2

            Bp_floor = max(1e-6 * np.max(Bp_mod), 1e-14)
            Bp_safe = np.maximum(Bp_mod, Bp_floor)
            fe_dl = dl / Bp_safe
            int_fe_dl = np.sum(fe_dl)

            def flx_avg(quantity):
                return np.sum(fe_dl * quantity) / int_fe_dl

            avg["R"][k] = flx_avg(r_s)
            avg["1/R"][k] = flx_avg(1.0 / r_s)
            avg["1/R**2"][k] = flx_avg(1.0 / r_s**2)
            avg["R**2"][k] = flx_avg(r_s**2)
            avg["Bp"][k] = flx_avg(Bp_signed)
            avg["Bp**2"][k] = flx_avg(Bp2_s)
            avg["Bt"][k] = flx_avg(Bt_s)
            avg["Bt**2"][k] = flx_avg(Bt_s**2)
            avg["Btot**2"][k] = flx_avg(B2_s)
            # Numerical <Jt> from curl(B)/mu_0; stored for cross-checks.
            # Note: avg["Jt/R"] is written analytically below (after the
            # per-surface loop), so we intentionally do NOT populate a
            # numerical `Jt/R` here -- save it under `Jt/R_num` for
            # downstream consumers that want the numerical version.
            avg["Jt"][k] = flx_avg(Jt_s)
            avg["Jt/R_num"][k] = flx_avg(Jt_s / r_s)

            avg["vp"][k] = (
                cc["sigma_rhotp"] * cc["sigma_Bp"]
                * np.sign(avg["Bp"][k])
                * int_fe_dl
                * (2.0 * np.pi) ** (1.0 - cc["exp_Bp"])
            )

            avg["q"][k] = (
                cc["sigma_rhotp"] * cc["sigma_Bp"]
                * avg["vp"][k] * F_k * avg["1/R**2"][k]
                / ((2 * np.pi) ** (2.0 - cc["exp_Bp"]))
            )

            avg["ip"][k] = (
                cc["sigma_rhotp"]
                * np.sum(dl * Bp_signed)
                / (constants.mu_0)
            )

            geo_k = _flux_geometry(r_s, z_s)
            for gk in geo_arrays:
                if gk in geo_k:
                    geo_arrays[gk][k] = geo_k[gk]

            Bpl = np.sum(Bp_mod * dl * 2 * np.pi)
            Bp2_vol += Bpl * dpsi_arr[k]

        # Fix q on axis by quadratic extrapolation
        if psi_N_levels[0] == 0 and nc > 3:
            x = psi_N_levels[1:4]
            y = avg["q"][1:4]
            coeffs = np.polyfit(x, y, 2)
            avg["q"][0] = coeffs[2]
        elif psi_N_levels[0] == 0 and nc > 2:
            x = psi_N_levels[1:3]
            y = avg["q"][1:3]
            avg["q"][0] = y[0] - (y[1] - y[0]) / (x[1] - x[0]) * x[0]

        # Fix near-axis geometry
        if psi_N_levels[0] == 0 and nc > 2:
            x = psi_N_levels[1:]
            for gk in ["kapu", "kapl"]:
                y = geo_arrays[gk][1:]
                geo_arrays[gk][0] = y[1] - ((y[1] - y[0]) / (x[1] - x[0])) * x[1]
            geo_arrays["kappa"][0] = 0.5 * (geo_arrays["kapu"][0] + geo_arrays["kapl"][0])
            for gk in ["delta", "delu", "dell"]:
                geo_arrays[gk][0] = 0.0

        # Analytic Grad-Shafranov current densities
        avg["Jt/R"] = (
            -cc["sigma_Bp"]
            * (avg["PPRIME"] + avg["FFPRIM"] * avg["1/R**2"] / constants.mu_0)
            * (2.0 * np.pi) ** cc["exp_Bp"]
        ) + 0.0
        avg["Jt_GS"] = (
            -cc["sigma_Bp"]
            * (avg["PPRIME"] * avg["R"] + avg["FFPRIM"] * avg["1/R"] / constants.mu_0)
            * (2.0 * np.pi) ** cc["exp_Bp"]
        ) + 0.0

        jt_check = avg["Jt_GS"][:-1] if self._extrapolate_edge else avg["Jt_GS"]
        jt_nz = jt_check[jt_check != 0]
        if len(jt_nz) > 1 and np.any(jt_nz > 0) and np.any(jt_nz < 0):
            warnings.warn(
                "The direct Grad-Shafranov Jt profile changes sign across "
                "flux surfaces; this may indicate unusual equilibrium data.",
                stacklevel=2,
            )

        # Geometry: volume and cross-section area
        psi_arr = psi_N_levels * dpsi + self.psi_axis
        geo_arrays["vol"] = np.abs(
            integrate.cumulative_trapezoid(avg["vp"], psi_arr, initial=0)
        )
        geo_arrays["cxArea"] = np.abs(
            integrate.cumulative_trapezoid(
                avg["vp"] * avg["1/R"], psi_arr, initial=0
            ) / (2.0 * np.pi)
        )

        # Internal inductance
        ip = self.Ip
        r_0 = self.R_center if self.R_center else R0

        Rb = self._raw["RBBBS"]
        Zb = self._raw["ZBBBS"]
        if len(Rb) > 3:
            dRb = np.diff(Rb, append=Rb[0])
            dZb = np.diff(Zb, append=Zb[0])
            circum = np.sum(np.sqrt(dRb**2 + dZb**2))
            vol = np.pi * np.abs(np.sum(Rb**2 * dZb))
            a_bdry = 0.5 * (np.max(Rb) - np.min(Rb))
            kappa_bdry = 0.5 * (np.max(Zb) - np.min(Zb)) / a_bdry if a_bdry > 0 else 1.0
        else:
            circum = geo_arrays["perimeter"][-1] if geo_arrays["perimeter"][-1] > 0 else 1.0
            vol = geo_arrays["vol"][-1] if geo_arrays["vol"][-1] > 0 else 1.0
            a_bdry = geo_arrays["a"][-1] if geo_arrays["a"][-1] > 0 else 1.0
            kappa_bdry = geo_arrays["kappa"][-1] if geo_arrays["kappa"][-1] > 0 else 1.0

        kappa_a = vol / (2.0 * np.pi * r_0 * np.pi * a_bdry * a_bdry) if a_bdry > 0 else 1.0
        correction_factor = (1 + kappa_bdry**2) / (2.0 * kappa_a) if kappa_a > 0 else 1.0

        if abs(ip) > 0:
            li_def = Bp2_vol / vol / constants.mu_0**2 / ip**2 * circum**2
        else:
            li_def = 0.0

        li_info = {
            "li_from_definition": li_def,
            "li(1)_EFIT": li_def,
            "li(1)_TLUCE": li_def / circum**2 * 2 * vol / r_0 * correction_factor if circum > 0 else 0.0,
            "li(1)": li_def,
            "li(2)": li_def / circum**2 * 2 * vol / R0 if circum > 0 else 0.0,
            "li(3)": 2 * Bp2_vol / r_0 / ip**2 / constants.mu_0**2 if abs(ip) > 0 else 0.0,
        }

        # Betas
        betas = {}
        if np.any(self._raw["PRES"]):
            P_interp = interpolate.InterpolatedUnivariateSpline(psi_N_raw, self._raw["PRES"])
            P_on_levels = np.array([float(P_interp(pn)) for pn in psi_N_levels])
            Btvac = self.B_center * self.R_center / geo_arrays["R"][-1]
            P_vol = integrate.cumulative_trapezoid(avg["vp"] * P_on_levels, psi_arr, initial=0)
            if vol > 0 and abs(Btvac) > 0:
                betas["beta_t"] = abs(P_vol[-1] / (Btvac**2 / 2.0 / constants.mu_0) / vol)
                i_MA = ip / 1e6
                betas["beta_n"] = betas["beta_t"] / abs(i_MA / a_bdry / Btvac) * 100 if abs(i_MA * a_bdry * Btvac) > 0 else 0.0
            Bpave = ip * constants.mu_0 / circum if circum > 0 else 1.0
            if vol > 0 and abs(Bpave) > 0:
                betas["beta_p"] = abs(P_vol[-1] / (Bpave**2 / 2.0 / constants.mu_0) / vol)

        # Outboard midplane profiles
        mid = {}
        mid["R"] = geo_arrays["R"] + geo_arrays["a"]
        mid["Z"] = np.full(nc, Z0)

        mid["Br"] = Br_interp.ev(mid["Z"], mid["R"])
        mid["Bz"] = Bz_interp.ev(mid["Z"], mid["R"])

        signBp = (
            -cc["sigma_rhotp"] * cc["sigma_RpZ"]
            * np.sign(mid["Bz"])
        )
        mid["Bp"] = signBp * np.sqrt(mid["Br"] ** 2 + mid["Bz"] ** 2)
        mid["Bp"][0] = 0.0

        mid["Bt"] = np.array([
            float(F_interp(psi_N_levels[k])) / mid["R"][k]
            for k in range(nc)
        ])
        mid["Btot"] = np.sqrt(mid["Bp"] ** 2 + mid["Bt"] ** 2)
        self._cache["midplane"] = mid

        self._cache["avg"] = avg
        self._cache["geo"] = geo_arrays
        self._cache["contours"] = contour_data
        self._cache["li"] = li_info
        self._cache["betas"] = betas

    # --- Public analysis properties ---

    @property
    def j_tor_averaged(self):
        r'''! Flux-surface-averaged toroidal current density \f$\langle J_t/R\rangle / \langle 1/R\rangle\f$ [A/m\f$^2\f$]'''
        self._trace_surfaces()
        avg = self._cache["avg"]
        return avg["Jt/R"] / avg["1/R"]

    @property
    def j_tor_averaged_direct(self):
        r'''! Direct flux-surface average \f$\langle J_t\rangle\f$ from the Grad-Shafranov equation [A/m\f$^2\f$]'''
        self._trace_surfaces()
        return self._cache["avg"]["Jt_GS"]

    @property
    def j_tor_averaged_numerical(self):
        r'''! Numerically flux-surface-averaged \f$\langle J_t\rangle\f$ via \f$\nabla\times B / \mu_0\f$ [A/m\f$^2\f$]'''
        self._trace_surfaces()
        return self._cache["avg"]["Jt"]

    @property
    def j_tor_over_R(self):
        r'''! \f$\langle J_t/R\rangle\f$ from the Grad-Shafranov equation [A/m\f$^3\f$]'''
        self._trace_surfaces()
        return self._cache["avg"]["Jt/R"]

    @property
    def q_profile(self):
        r'''! Safety factor from flux-surface averaging'''
        self._trace_surfaces()
        return self._cache["avg"]["q"]

    @property
    def li(self):
        r'''! Internal inductance dict (keys: `li_from_definition`, `li(1)`, `li(1)_EFIT`, `li(1)_TLUCE`, `li(2)`, `li(3)`)'''
        self._trace_surfaces()
        return self._cache["li"]

    @property
    def geometry(self):
        r'''! Per-surface geometric quantities dict

        Keys: `R`, `Z`, `a`, `kappa`, `kapu`, `kapl`, `delta`, `delu`, `dell`,
        `perimeter`, `surfArea`, `eps`, `vol`, `cxArea`
        '''
        self._trace_surfaces()
        return self._cache["geo"]

    @property
    def averages(self):
        r'''! All flux-surface-averaged quantities dict

        Keys: `R`, `1/R`, `1/R**2`, `R**2`, `Bp`, `Bp**2`, `Bt`, `Bt**2`,
        `Btot**2`, `Jt`, `Jt/R`, `Jt/R_num`, `Jt_GS`, `vp`, `q`, `ip`, `F`,
        `PPRIME`, `FFPRIM`.

        `Jt/R` is the analytic Grad-Shafranov form
        \f$-\sigma_{Bp}(p' + FF'\langle 1/R^2\rangle/\mu_0)(2\pi)^{e_{Bp}}\f$;
        `Jt_GS` is the companion \f$\langle J_t\rangle\f$ analytic form used by
        `j_tor_averaged_direct`.  `Jt` and `Jt/R_num` are the numerical
        counterparts computed from \f$\nabla\times B/\mu_0\f$ on the contour,
        retained for cross-checks.
        '''
        self._trace_surfaces()
        return self._cache["avg"]

    @property
    def midplane(self):
        r'''! Outboard midplane quantities on the \f$\hat{\psi}\f$ grid (keys: `R`, `Z`, `Br`, `Bz`, `Bp`, `Bt`, `Btot`)'''
        self._trace_surfaces()
        return self._cache["midplane"]

    @property
    def betas(self):
        r'''! Plasma beta values (keys: `beta_t`, `beta_p`, `beta_n`)'''
        self._trace_surfaces()
        return self._cache["betas"]

    @property
    def contours(self):
        r'''! List of `(N,2)` contour arrays for each \f$\hat{\psi}\f$ level'''
        self._trace_surfaces()
        return self._cache["contours"]

    # --- Integration methods ---

    def volume_integral(self, what):
        r'''! Volume integral of a quantity on the FSA grid

        @param what Array-like of length `nlevels`, sampled at `self.psi_N_levels`
        @result 1-D ndarray (`nlevels`,) cumulative integral from core to each surface
        '''
        self._trace_surfaces()
        dpsi = self.psi_boundary - self.psi_axis
        psi_arr = self.psi_N_levels * dpsi + self.psi_axis
        return integrate.cumulative_trapezoid(
            self._cache["avg"]["vp"] * np.asarray(what), psi_arr, initial=0
        )

    def surface_integral(self, what):
        r'''! Cross-section (area) integral of a quantity on the FSA grid

        @param what Array-like of length `nlevels`, sampled at `self.psi_N_levels`
        @result 1-D ndarray (`nlevels`,) cumulative area integral
        '''
        self._trace_surfaces()
        dpsi = self.psi_boundary - self.psi_axis
        psi_arr = self.psi_N_levels * dpsi + self.psi_axis
        return integrate.cumulative_trapezoid(
            self._cache["avg"]["vp"] * self._cache["avg"]["1/R"] * np.asarray(what),
            psi_arr, initial=0,
        ) / (2.0 * np.pi)

    def flux_integral(self, psi_N_val, profile):
        r'''! Total flux integral at a given \f$\hat{\psi}\f$ value (scalar)

        @param psi_N_val Normalised poloidal flux location (0 = axis, 1 = boundary)
        @param profile Array-like of length `nlevels` to integrate (sampled at
                       `self.psi_N_levels`)
        @result Value of the volume integral at *psi_N_val*
        '''
        cum = self.volume_integral(profile)
        return float(np.interp(psi_N_val, self.psi_N_levels, cum))

    # --- Introspection / discovery ---

    # Catalogue of all accessible quantities, grouped by category.
    # Each entry: (attribute_name, short_description, units_or_shape_hint)
    _CATALOGUE = [
        ("Global scalars", [
            ("Ip",           "Plasma current",                       "A"),
            ("R_mag",        "Magnetic axis R",                      "m"),
            ("Z_mag",        "Magnetic axis Z",                      "m"),
            ("R_center",     "Reference geometric center R",         "m"),
            ("B_center",     "Vacuum toroidal field at R_center",    "T"),
            ("psi_axis",     "Poloidal flux at magnetic axis",       "Wb"),
            ("psi_boundary", "Poloidal flux at LCFS",                "Wb"),
            ("cocos",        "COCOS convention index",               ""),
        ]),
        ("1-D profiles (raw, on psi_N grid from g-file, length NW)", [
            ("fpol",   "F = R * B_t",                "T m"),
            ("pres",   "Pressure",                   "Pa"),
            ("pprime", "dP/dpsi",                    "Pa/Wb"),
            ("ffprim", "F * dF/dpsi",                "T^2 m^2/Wb"),
            ("qpsi",   "Safety factor (from file)",  ""),
            ("rhovn",  "Normalised toroidal flux",   ""),
            ("psi_N",  "Normalised psi grid (NW)",   ""),
        ]),
        ("Boundary and limiter", [
            ("boundary_R", "LCFS R coordinates", "m"),
            ("boundary_Z", "LCFS Z coordinates", "m"),
            ("limiter_R",  "Limiter R coordinates", "m"),
            ("limiter_Z",  "Limiter Z coordinates", "m"),
        ]),
        ("2-D grids", [
            ("R_grid",    "1-D R grid",                        "m"),
            ("Z_grid",    "1-D Z grid",                        "m"),
            ("psi_RZ",    "Poloidal flux on (R,Z) grid",       "Wb"),
            ("psi_N_RZ",  "Normalised flux on (R,Z) grid",     ""),
        ]),
        ("Derived 1-D profiles (lazy, FSA, on psi_N_levels, length nlevels)", [
            ("psi_N_levels",               "Normalised psi FSA grid (nlevels)", ""),
            ("q_profile",                  "Safety factor from FSA",            ""),
            ("j_tor_averaged",             "<Jt/R>/<1/R> (standard)",          "A/m^2"),
            ("j_tor_averaged_direct",      "<Jt> from GS",                      "A/m^2"),
            ("j_tor_averaged_numerical",   "<Jt> via curl(B)/mu0",              "A/m^2"),
            ("j_tor_over_R",               "<Jt/R> from GS",                    "A/m^3"),
        ]),
        ("Derived containers (lazy, dicts of arrays)", [
            ("geometry",  "R,Z,a,kappa,delta,vol,cxArea,perimeter,...", "dict"),
            ("averages",  "All FSA quantities (R,1/R,Bp,Bt,vp,q,ip,...)", "dict"),
            ("midplane",  "Outboard midplane R,Z,Br,Bz,Bp,Bt,Btot",     "dict"),
            ("li",        "Internal inductance li(1),li(2),li(3),...",  "dict"),
            ("betas",     "Plasma betas beta_t, beta_p, beta_n",         "dict"),
            ("contours",  "List of (N,2) flux-surface (R,Z) contours",   "list"),
        ]),
    ]

    @property
    def keys(self):
        r'''! Flat list of all available attribute names (scalars + profiles + derived)'''
        return [name for _, items in self._CATALOGUE for name, *_ in items]

    def describe(self, verbose=False):
        r'''! Print a categorised summary of all available quantities

        Access any listed name as an attribute, e.g. `eq.q_profile` or
        `eq.geometry['kappa']`.  Dict-valued entries (`geometry`, `li`, etc.)
        can be inspected with their `.keys()`.

        @param verbose If `True`, also triggers lazy evaluation of derived
                       quantities so actual shapes and sample values are printed.
                       If `False` (default), shape info is only printed for
                       already-cached quantities to avoid an unexpected
                       flux-surface trace.
        '''
        def _shape_str(val):
            if isinstance(val, np.ndarray):
                if val.ndim == 0:
                    return f"scalar={val.item():.4g}"
                return f"array {val.shape}"
            if isinstance(val, dict):
                return f"dict[{len(val)}] keys={list(val.keys())}"
            if isinstance(val, list):
                return f"list[{len(val)}]"
            if isinstance(val, (int, float, np.floating, np.integer)):
                return f"{float(val):.6g}"
            return str(type(val).__name__)

        print(f"GEQDSKEquilibrium (COCOS {self._cocos_index}, "
              f"NW={int(self._raw['NW'])}, NH={int(self._raw['NH'])}, "
              f"nlevels={self._nlevels})")
        print("=" * 78)
        for section, items in self._CATALOGUE:
            print(f"\n{section}")
            print("-" * len(section))
            for name, desc, units in items:
                # Only evaluate non-lazy attrs unless verbose is requested
                is_lazy = name in ("q_profile", "j_tor_averaged",
                                    "j_tor_averaged_direct",
                                    "j_tor_averaged_numerical",
                                    "j_tor_over_R", "geometry", "averages",
                                    "midplane", "li", "betas", "contours")
                if is_lazy and not verbose and "avg" not in self._cache:
                    shape = "(lazy; call describe(verbose=True) or access to compute)"
                else:
                    try:
                        val = getattr(self, name)
                        shape = _shape_str(val)
                    except Exception as exc:
                        shape = f"<error: {exc}>"
                unit_str = f" [{units}]" if units else ""
                print(f"  {name:30s} {desc:48s}{unit_str}")
                if shape:
                    print(f"  {'':30s}   -> {shape}")


# ===========================================================================
# Osborne p-file (kinetic profile) reader/writer
# ===========================================================================

# Unit conversion: n [10^20/m^3] * T [keV] -> p [kPa]
_NT_TO_KPA = 1.602176634e-19 * 1e20  # exactly 16.02176634

# ---------------------------------------------------------------------------
# Known profile metadata
# ---------------------------------------------------------------------------

PFILE_DESCRIPTIONS = OrderedDict([
    ("ne", "Electron density"),
    ("te", "Electron temperature"),
    ("ni", "Ion density"),
    ("ti", "Ion temperature"),
    ("nb", "Fast ion density"),
    ("pb", "Fast ion pressure"),
    ("ptot", "Total pressure"),
    ("omeg", "Toroidal rotation: VTOR/R"),
    ("omegp", "Poloidal rotation: Bt * VPOL / (RBp)"),
    ("omgvb", "VxB rotation term in the ExB rotation frequency"),
    ("omgpp", "Diamagnetic term in the ExB rotation frequency"),
    ("omgeb", "ExB rotation frequency"),
    ("er", "Radial electric field from force balance"),
    ("ommvb", "Main ion VxB term of Er/RBp"),
    ("ommpp", "Main ion pressure term of Er/RBp"),
    ("omevb", "Electron VxB term of Er/RBp"),
    ("omepp", "Electron pressure term of Er/RBp"),
    ("kpol", "KPOL = VPOL/Bp"),
    ("omghb", "Hahm-Burrell ExB velocity shearing rate"),
    ("nz1", "Density of the 1st impurity species"),
    ("vtor1", "Toroidal velocity of the 1st impurity species"),
    ("vpol1", "Poloidal velocity of the 1st impurity species"),
    ("nz2", "Density of the 2nd impurity species"),
    ("vtor2", "Toroidal velocity of the 2nd impurity species"),
    ("vpol2", "Poloidal velocity of the 2nd impurity species"),
])

PFILE_UNITS = OrderedDict([
    ("ne", "10^20/m^3"),
    ("te", "KeV"),
    ("ni", "10^20/m^3"),
    ("ti", "KeV"),
    ("nb", "10^20/m^3"),
    ("pb", "KPa"),
    ("ptot", "KPa"),
    ("omeg", "kRad/s"),
    ("omegp", "kRad/s"),
    ("omgvb", "kRad/s"),
    ("omgpp", "kRad/s"),
    ("omgeb", "kRad/s"),
    ("er", "kV/m"),
    ("ommvb", ""),
    ("ommpp", ""),
    ("omevb", ""),
    ("omepp", ""),
    ("kpol", "km/s/T"),
    ("omghb", ""),
    ("nz1", "10^20/m^3"),
    ("vtor1", "km/s"),
    ("vpol1", "km/s"),
    ("nz2", "10^20/m^3"),
    ("vtor2", "km/s"),
    ("vpol2", "km/s"),
])

_PFILE_HEADER_RE = re.compile(
    r"^(\d+)\s+(\S+)\s+(\S+)\(([^)]*)\)\s+(.*?)\s*$"
)


# ---------------------------------------------------------------------------
# P-file parser / writer
# ---------------------------------------------------------------------------

def _read_pfile(filename):
    r'''! Parse an Osborne p-file into an OrderedDict

    @param filename Path to the p-file
    @result OrderedDict keyed by profile name (`"ne"`, `"te"`, ...). Each value
            is a dict with keys `"psinorm"`, `"data"`, `"derivative"`, `"units"`,
            and `"deriv_label"`.  The special key `"N Z A"` (if present) maps
            to a dict with `"N"`, `"Z"`, `"A"` arrays.
    '''
    with open(filename, "r") as f:
        lines = f.read().strip().splitlines()

    profiles = OrderedDict()
    idx = 0
    while idx < len(lines):
        header = lines[idx]
        tokens = header.split()
        if len(tokens) < 2:
            idx += 1
            continue

        count = int(tokens[0])

        if "N Z A of ION SPECIES" in header:
            rows = []
            for i in range(idx + 1, idx + 1 + count):
                rows.append(list(map(float, lines[i].split())))
            cols = list(zip(*rows))
            profiles["N Z A"] = {
                "N": np.array(cols[0]),
                "Z": np.array(cols[1]),
                "A": np.array(cols[2]),
            }
            idx += 1 + count
            continue

        m = _PFILE_HEADER_RE.match(header)
        if m is None:
            idx += 1
            continue

        _xkey = m.group(2)
        key = m.group(3)
        units = m.group(4)
        deriv_label = m.group(5)

        rows = []
        for i in range(idx + 1, idx + 1 + count):
            rows.append(list(map(float, lines[i].split())))
        cols = list(zip(*rows))

        profiles[key] = {
            "psinorm": np.array(cols[0]),
            "data": np.array(cols[1]),
            "derivative": np.array(cols[2]),
            "units": units,
            "deriv_label": deriv_label,
        }

        idx += 1 + count

    return profiles


def _write_pfile(profiles, filename):
    r'''! Write an OrderedDict of profiles to an Osborne p-file

    @param profiles OrderedDict with the same structure as returned by `_read_pfile`
    @param filename Output path for the p-file
    '''
    buf = []
    for key, val in profiles.items():
        if key == "N Z A":
            n = len(val["A"])
            buf.append(f"{n} N Z A of ION SPECIES\n")
            for i in range(n):
                buf.append(
                    f" {val['N'][i]:f}   {val['Z'][i]:f}   {val['A'][i]:f}\n"
                )
        else:
            n = len(val["data"])
            if n <= 1:
                continue
            units = val.get("units", "")
            deriv_label = val.get("deriv_label", f"d{key}/dpsiN")
            buf.append(f"{n} psinorm {key}({units}) {deriv_label}\n")
            for i in range(n):
                buf.append(
                    f" {val['psinorm'][i]:f}   {val['data'][i]:f}"
                    f"   {val['derivative'][i]:f}\n"
                )

    with open(filename, "w") as f:
        f.writelines(buf)


# ---------------------------------------------------------------------------
# PFile class
# ---------------------------------------------------------------------------

class PFile:
    r'''! Interface for Osborne p-file kinetic profiles

    Example usage:
    ```
    pf = PFile("p123456.01234")
    pf.ne                    # electron density array
    pf.te                    # electron temperature array
    pf.psinorm_for("ne")     # psinorm grid for ne
    "omgeb" in pf            # check if profile exists
    ```
    '''

    def __init__(self, filename):
        r'''! Initialize PFile by reading an Osborne p-file from disk

        @param filename Path to the p-file
        '''
        self._raw = _read_pfile(filename)

    @classmethod
    def from_bytes(cls, raw_bytes):
        r'''! Construct from in-memory bytes

        @param raw_bytes Raw p-file content
        @result New `PFile` instance
        '''
        import os as _os
        with tempfile.NamedTemporaryFile(
            mode="wb", suffix=".pfile", delete=False
        ) as tmp:
            tmp.write(raw_bytes)
            tmp_path = tmp.name
        try:
            return cls(tmp_path)
        finally:
            try:
                _os.remove(tmp_path)
            except OSError:
                pass

    def save(self, filename):
        r'''! Write the profiles to *filename* in p-file format

        @param filename Output path for the p-file
        '''
        _write_pfile(self._raw, filename)

    def to_bytes(self):
        r'''! Serialise to in-memory bytes (round-trips with `from_bytes`)

        @result Raw p-file content as bytes
        '''
        import os as _os
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pfile", delete=False
        ) as tmp:
            tmp_path = tmp.name
        try:
            _write_pfile(self._raw, tmp_path)
            with open(tmp_path, "rb") as fh:
                data = fh.read()
        finally:
            try:
                _os.remove(tmp_path)
            except OSError:
                pass
        return data

    # --- Dict-like access ---

    @property
    def keys(self):
        r'''! Profile names in file order (list of str)'''
        return list(self._raw.keys())

    def __contains__(self, key):
        return key in self._raw

    def __getitem__(self, key):
        return self._raw[key]

    def __iter__(self):
        return iter(self._raw)

    def __len__(self):
        return len(self._raw)

    # --- Per-profile accessors ---

    def psinorm_for(self, key):
        r'''! Return the psinorm grid for profile *key*

        @param key Profile name
        @result 1-D ndarray of psinorm values, or `None` if not present
        '''
        entry = self._raw.get(key)
        if entry is None or key == "N Z A":
            return None
        return entry["psinorm"]

    def derivative_for(self, key):
        r'''! Return the derivative array for profile *key*

        @param key Profile name
        @result 1-D ndarray of d(data)/dpsinorm, or `None` if not present
        '''
        entry = self._raw.get(key)
        if entry is None or key == "N Z A":
            return None
        return entry["derivative"]

    def units_for(self, key):
        r'''! Return the units string for profile *key*

        @param key Profile name
        @result Units string, or `None` if the key is not a profile
        '''
        entry = self._raw.get(key)
        if entry is None or key == "N Z A":
            return None
        return entry.get("units", "")

    # --- Named properties for common profiles ---

    def _get_data(self, key):
        entry = self._raw.get(key)
        if entry is None:
            return None
        return entry["data"]

    @property
    def ne(self):
        r'''! Electron density \f$n_e\f$ [10\f$^{20}\f$/m\f$^3\f$]'''
        return self._get_data("ne")

    @property
    def te(self):
        r'''! Electron temperature \f$T_e\f$ [keV]'''
        return self._get_data("te")

    @property
    def ni(self):
        r'''! Ion density \f$n_i\f$ [10\f$^{20}\f$/m\f$^3\f$]'''
        return self._get_data("ni")

    @property
    def ti(self):
        r'''! Ion temperature \f$T_i\f$ [keV]'''
        return self._get_data("ti")

    @property
    def nb(self):
        r'''! Fast ion density \f$n_b\f$ [10\f$^{20}\f$/m\f$^3\f$]'''
        return self._get_data("nb")

    @property
    def pb(self):
        r'''! Fast ion pressure \f$p_b\f$ [kPa]'''
        return self._get_data("pb")

    @property
    def ptot(self):
        r'''! Total pressure \f$p_{\mathrm{tot}}\f$ [kPa]'''
        return self._get_data("ptot")

    @property
    def omeg(self):
        r'''! Toroidal rotation \f$V_{\mathrm{TOR}}/R\f$ [kRad/s]'''
        return self._get_data("omeg")

    @property
    def omegp(self):
        r'''! Poloidal rotation \f$B_t V_{\mathrm{POL}}/(R B_p)\f$ [kRad/s]'''
        return self._get_data("omegp")

    @property
    def omgvb(self):
        r'''! \f$V\times B\f$ rotation term [kRad/s]'''
        return self._get_data("omgvb")

    @property
    def omgpp(self):
        r'''! Diamagnetic rotation term [kRad/s]'''
        return self._get_data("omgpp")

    @property
    def omgeb(self):
        r'''! \f$E\times B\f$ rotation frequency [kRad/s]'''
        return self._get_data("omgeb")

    @property
    def er(self):
        r'''! Radial electric field \f$E_r\f$ [kV/m]'''
        return self._get_data("er")

    @property
    def kpol(self):
        r'''! \f$K_{\mathrm{POL}} = V_{\mathrm{POL}}/B_p\f$ [km/s/T]'''
        return self._get_data("kpol")

    @property
    def omghb(self):
        r'''! Hahm-Burrell \f$E\times B\f$ shearing rate'''
        return self._get_data("omghb")

    @property
    def ion_species(self):
        r'''! Ion species dict with `'N'`, `'Z'`, `'A'` arrays, or `None`'''
        return self._raw.get("N Z A")

    # --- Construction helpers ---

    @classmethod
    def new(cls):
        r'''! Create an empty PFile (no profiles loaded from disk)

        @result New empty `PFile` instance
        '''
        obj = object.__new__(cls)
        obj._raw = OrderedDict()
        return obj

    def set_profile(self, key, psinorm, data, derivative=None, units=None):
        r'''! Add or replace a profile

        @param key Profile name (e.g. `"ne"`, `"te"`)
        @param psinorm Normalised poloidal flux grid
        @param data Profile values on *psinorm*
        @param derivative Derivative d(data)/d(psinorm); if `None`, computed
                          via `np.gradient`
        @param units Unit label; if `None`, looked up from `PFILE_UNITS`
        '''
        psinorm = np.asarray(psinorm, dtype=float)
        data = np.asarray(data, dtype=float)
        if derivative is None:
            derivative = np.gradient(data, psinorm)
        else:
            derivative = np.asarray(derivative, dtype=float)
        if units is None:
            units = PFILE_UNITS.get(key, "")
        self._raw[key] = {
            "psinorm": psinorm,
            "data": data,
            "derivative": derivative,
            "units": units,
            "deriv_label": f"d{key}/dpsiN",
        }

    def set_ion_species(self, N, Z, A):
        r'''! Set the ion species block

        @param N Atomic number for each species
        @param Z Charge state for each species
        @param A Mass number for each species
        '''
        self._raw["N Z A"] = {
            "N": np.asarray(N, dtype=float),
            "Z": np.asarray(Z, dtype=float),
            "A": np.asarray(A, dtype=float),
        }

    def compute_derivatives(self):
        r'''! Recompute d(data)/d(psinorm) for all profiles in place'''
        for key, val in self._raw.items():
            if key == "N Z A":
                continue
            val["derivative"] = np.gradient(val["data"], val["psinorm"])

    # --- Physics computations ---

    def compute_pressure(self):
        r'''! Compute total pressure from density and temperature profiles

        Uses \f$p_{\mathrm{tot}} = c_{NT}\,(n_e T_e + (n_i + n_{z1}) T_i) + p_b\f$
        with \f$c_{NT}\f$ converting from (10\f$^{20}\f$/m\f$^3\f$ \f$\cdot\f$ keV)
        to kPa.  Requires `ne`, `te`, `ni`, `ti` on the same psinorm grid.
        `nz1` and `pb` default to zero if absent.  Stores the result as `ptot`.
        '''
        psinorm = self._raw["ne"]["psinorm"]
        ne = self._raw["ne"]["data"]
        te = self._raw["te"]["data"]
        ni = self._raw["ni"]["data"]
        ti = self._raw["ti"]["data"]

        nz1 = self._get_data("nz1")
        if nz1 is None:
            nz1 = np.zeros_like(psinorm)
        pb = self._get_data("pb")
        if pb is None:
            pb = np.zeros_like(psinorm)

        ptot = _NT_TO_KPA * (ne * te + (ni + nz1) * ti) + pb
        self.set_profile("ptot", psinorm, ptot)

    def compute_quasineutrality(self):
        r'''! Compute impurity density \f$n_{z1}\f$ from quasi-neutrality

        \f$n_{z1} = (n_e - n_i - n_b) / Z_{\mathrm{imp}}\f$.
        Requires `ne`, `ni` on the same grid and an `"N Z A"` block with at
        least one impurity species.  `nb` defaults to zero if absent.  Stores
        the result as `nz1`.
        '''
        nza = self._raw.get("N Z A")
        if nza is None:
            raise ValueError("Ion species (N Z A) block required")
        Z_imp = nza["Z"][0]

        psinorm = self._raw["ne"]["psinorm"]
        ne = self._raw["ne"]["data"]
        ni = self._raw["ni"]["data"]
        nb = self._get_data("nb")
        if nb is None:
            nb = np.zeros_like(psinorm)

        nz1 = (ne - ni - nb) / Z_imp
        n_neg = np.count_nonzero(nz1 < 0)
        if n_neg:
            warnings.warn(
                f"Quasi-neutrality produced negative nz1 at {n_neg}/{len(nz1)} "
                f"grid points (min = {nz1.min():.4g})."
            )
        self.set_profile("nz1", psinorm, nz1)

    def compute_zeff(self):
        r'''! Compute the effective charge profile \f$Z_{\mathrm{eff}}\f$

        \f$Z_{\mathrm{eff}} = (n_i Z_{\mathrm{main}}^2 + n_{z1} Z_{\mathrm{imp}}^2
        + n_b Z_{\mathrm{beam}}^2)/n_e\f$.  Charge states are read from the
        `"N Z A"` block using the OMFIT convention: impurities first, then
        main ion, beam ion last.  Species-count-specific handling:
          - 1 species: treated as the impurity (main = beam = main ion
                       hydrogen, \f$Z=1\f$).
          - 2 species: impurity + main ion (no separate beam; beam \f$n_b\f$
                       defaults to zero so its charge state is irrelevant).
          - 3+ species: impurity(ies) first, main ion at [-2], beam at [-1].

        @result Tuple `(psinorm, zeff)` of 1-D ndarrays
        '''
        nza = self._raw.get("N Z A")
        if nza is None:
            raise ValueError("Ion species (N Z A) block required")
        Zs = np.asarray(nza["Z"], dtype=float)
        if Zs.size < 1:
            raise ValueError("N Z A block must contain at least one species")
        Z_imp = Zs[0]
        if Zs.size == 1:
            # Single species taken as the impurity; assume hydrogenic main/beam
            Z_main = 1.0
            Z_beam = 1.0
        elif Zs.size == 2:
            # Impurity + main ion (no separate beam species)
            Z_main = Zs[-1]
            Z_beam = 1.0   # beam density will default to zero, so Z_beam is unused
        else:
            Z_main = Zs[-2]
            Z_beam = Zs[-1]

        psinorm = self._raw["ne"]["psinorm"]
        ne = self._raw["ne"]["data"]
        ni = self._raw["ni"]["data"]
        nz1 = self._get_data("nz1")
        if nz1 is None:
            nz1 = np.zeros_like(psinorm)
        nb = self._get_data("nb")
        if nb is None:
            nb = np.zeros_like(psinorm)

        with np.errstate(divide="ignore", invalid="ignore"):
            zeff = (ni * Z_main**2 + nz1 * Z_imp**2 + nb * Z_beam**2) / ne
        np.nan_to_num(zeff, copy=False, nan=1.0, posinf=1.0, neginf=1.0)
        return psinorm, zeff

    def compute_diamagnetic_rotations(self, psi, nI=None, TI=None):
        r'''! Compute diamagnetic rotation frequencies from kinetic profiles

        For species \f$s\f$:
        \f$\omega_{\mathrm{dia},s} = \frac{1}{n_s Z_s e}\frac{d(n_s T_s)}{d\psi}\f$.
        In p-file units (n in 10\f$^{20}\f$/m\f$^3\f$, T in keV, \f$\psi\f$ in Wb)
        this gives results in kRad/s.  Sets `omgpp`, `ommpp`, `omepp` profiles.

        @param psi Poloidal flux in SI (Weber), same length as the profile grids
        @param nI Impurity density [10\f$^{20}\f$/m\f$^3\f$]; if `None`, uses
                  `nz1` from the p-file when present, otherwise defaults to zero
        @param TI Impurity temperature [keV]; if `None`, uses `ti`
        '''
        psi = np.asarray(psi, dtype=float)
        dpsi = np.gradient(psi)
        psinorm = self._raw["ne"]["psinorm"]

        ne = self._raw["ne"]["data"]
        te = self._raw["te"]["data"]
        ni = self._raw["ni"]["data"]
        ti = self._raw["ti"]["data"]

        if nI is None:
            nz1_entry = self._raw.get("nz1")
            if nz1_entry is not None:
                nI = nz1_entry["data"]
            else:
                nI = np.zeros_like(ni, dtype=float)
        else:
            nI = np.asarray(nI, dtype=float)
        if TI is None:
            TI = ti
        else:
            TI = np.asarray(TI, dtype=float)

        nza = self._raw.get("N Z A")
        if nza is None:
            raise ValueError("Ion species (N Z A) block required")
        Z_imp = nza["Z"][0]

        dpsi_floor = max(1e-4 * np.max(np.abs(dpsi)), 1e-30)
        dpsi_safe = np.where(np.abs(dpsi) > dpsi_floor, dpsi,
                             np.sign(dpsi) * dpsi_floor)
        n_floor = 1e-4 * np.max(np.abs(ne))
        nI_safe = np.maximum(np.abs(nI), n_floor)
        ni_safe = np.maximum(np.abs(ni), n_floor)
        ne_safe = np.maximum(np.abs(ne), n_floor)

        with np.errstate(divide="ignore", invalid="ignore"):
            omgpp = -np.abs(np.gradient(nI * TI) / dpsi_safe / (nI_safe * Z_imp))
            ommpp = -np.abs(np.gradient(ni * ti) / dpsi_safe / (ni_safe * 1.0))
            omepp = np.abs(np.gradient(ne * te) / dpsi_safe / (ne_safe * 1.0))
        np.nan_to_num(omgpp, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        np.nan_to_num(ommpp, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        np.nan_to_num(omepp, copy=False, nan=0.0, posinf=0.0, neginf=0.0)

        self.set_profile("omgpp", psinorm, omgpp)
        self.set_profile("ommpp", psinorm, ommpp)
        self.set_profile("omepp", psinorm, omepp)

    def compute_rotation_decomposition(self, R=None, Bp=None, Bt=None,
                                       psi=None):
        r'''! Compute \f$E\times B\f$ and \f$V\times B\f$ rotation frequencies and derived quantities

        From the diamagnetic terms (`omgpp`, `ommpp`, `omepp`) and the impurity
        \f$V\times B\f$ rotation `omgvb`, computes:
        `omgeb` = `omgvb` + `omgpp`, `ommvb` = `omgeb` - `ommpp`,
        `omevb` = `omgeb` - `omepp`.  If equilibrium data are provided, also
        computes `er = omgeb * R * Bp` and the Hahm-Burrell shearing rate
        `omghb`.

        @param R Midplane major radius [m] on the profile psinorm grid
        @param Bp Poloidal field [T] on the profile psinorm grid
        @param Bt Toroidal field [T] on the profile psinorm grid
        @param psi Poloidal flux in SI (Weber)
        '''
        psinorm = self._raw["omgpp"]["psinorm"]

        omgvb = self._get_data("omgvb")
        if omgvb is None:
            omgvb = np.zeros_like(psinorm)

        omgpp = self._raw["omgpp"]["data"]
        ommpp = self._raw["ommpp"]["data"]
        omepp = self._raw["omepp"]["data"]

        omgeb = omgvb + omgpp
        ommvb = omgeb - ommpp
        omevb = omgeb - omepp

        self.set_profile("omgeb", psinorm, omgeb)
        self.set_profile("ommvb", psinorm, ommvb)
        self.set_profile("omevb", psinorm, omevb)

        if R is not None and Bp is not None:
            R = np.asarray(R, dtype=float)
            Bp = np.asarray(Bp, dtype=float)
            er = omgeb * R * Bp
            self.set_profile("er", psinorm, er)

            if Bt is not None and psi is not None:
                Bt = np.asarray(Bt, dtype=float)
                psi = np.asarray(psi, dtype=float)
                Bt_safe = np.where(np.abs(Bt) > 1e-6, Bt,
                                   np.sign(Bt) * 1e-6)
                from scipy.signal import savgol_filter
                n_pts = len(omgeb)
                win = max(7, int(np.round(n_pts * 0.03)) | 1)
                if win >= n_pts:
                    win = n_pts - (1 - n_pts % 2)
                delta_psi = np.abs(np.median(np.diff(psi)))
                if delta_psi == 0:
                    delta_psi = 1.0
                domgeb_dpsi = savgol_filter(
                    omgeb, win, min(3, win - 1),
                    deriv=1, delta=delta_psi,
                )
                omghb = (R * Bp) ** 2 / Bt_safe * domgeb_dpsi
                np.nan_to_num(omghb, copy=False, nan=0.0,
                              posinf=0.0, neginf=0.0)
                self.set_profile("omghb", psinorm, omghb)

    # --- Remap ---

    def remap(self, psinorm=None, key="ne"):
        r'''! Return a new `PFile` with all profiles on a common grid

        @param psinorm Target grid. If `None`, use the grid from *key*.
                       If `int`, use `np.linspace(0, 1, psinorm)`.
                       If array-like, use as given.
        @param key Profile whose grid to use when *psinorm* is `None`
        @result New `PFile` instance with interpolated profiles on the common grid
        '''
        if psinorm is None:
            if key not in self._raw:
                raise KeyError(f"Profile {key!r} not found for grid reference")
            target = self._raw[key]["psinorm"]
        elif isinstance(psinorm, (int, np.integer)):
            target = np.linspace(0, 1, int(psinorm))
        else:
            target = np.asarray(psinorm, dtype=float)

        new_raw = OrderedDict()
        for k, val in self._raw.items():
            if k == "N Z A":
                new_raw[k] = val.copy()
                continue

            f_data = interpolate.interp1d(
                val["psinorm"], val["data"],
                kind="linear", bounds_error=False,
                fill_value=(val["data"][0], val["data"][-1]),
            )
            f_deriv = interpolate.interp1d(
                val["psinorm"], val["derivative"],
                kind="linear", bounds_error=False,
                fill_value=(val["derivative"][0], val["derivative"][-1]),
            )
            new_raw[k] = {
                "psinorm": target.copy(),
                "data": f_data(target),
                "derivative": f_deriv(target),
                "units": val.get("units", ""),
                "deriv_label": val.get("deriv_label", f"d{k}/dpsiN"),
            }

        obj = object.__new__(PFile)
        obj._raw = new_raw
        return obj

    def __repr__(self):
        profile_keys = [k for k in self._raw if k != "N Z A"]
        return (
            f"PFile({len(profile_keys)} profiles: "
            f"{', '.join(profile_keys[:6])}"
            f"{'...' if len(profile_keys) > 6 else ''})"
        )

    # --- Introspection / discovery ---

    def describe(self):
        r'''! Print a categorised summary of all profiles and derived quantities

        Groups profiles by physics category (kinetic, rotation, impurity, etc.),
        showing grid length, data range, units, and whether a descriptive label
        from `PFILE_DESCRIPTIONS` is known.  The ion species block is also
        printed if present.  Access any profile as `pf.ne`, `pf.te`, etc. (for
        common named profiles), or generically via `pf['key']` /
        `pf._get_data('key')`.
        '''
        # Group profile keys by physics category
        categories = [
            ("Kinetic profiles",        ["ne", "te", "ni", "ti", "nb", "pb", "ptot"]),
            ("Rotation / Er",           ["omeg", "omegp", "omgvb", "omgpp", "omgeb",
                                          "er", "ommvb", "ommpp", "omevb", "omepp",
                                          "kpol", "omghb"]),
            ("Impurity / beam",         ["nz1", "vtor1", "vpol1",
                                          "nz2", "vtor2", "vpol2"]),
        ]
        present = set(self._raw.keys())

        print(f"PFile ({len(self._raw)} entries)")
        print("=" * 78)

        seen = set()
        for section, keys in categories:
            section_keys = [k for k in keys if k in present]
            if not section_keys:
                continue
            print(f"\n{section}")
            print("-" * len(section))
            for key in section_keys:
                seen.add(key)
                entry = self._raw[key]
                data = entry["data"]
                units = entry.get("units", "")
                desc = PFILE_DESCRIPTIONS.get(key, "")
                unit_str = f" [{units}]" if units else ""
                print(f"  {key:10s} {desc:48s}{unit_str}")
                print(f"  {'':10s}   -> len={len(data)}, "
                      f"range=[{np.min(data):.4g}, {np.max(data):.4g}]")

        # Any extra / unknown keys
        extras = [k for k in self._raw if k not in seen and k != "N Z A"]
        if extras:
            print("\nOther profiles")
            print("-" * len("Other profiles"))
            for key in extras:
                entry = self._raw[key]
                data = entry["data"]
                units = entry.get("units", "")
                unit_str = f" [{units}]" if units else ""
                print(f"  {key:10s} (no registered description){unit_str}")
                print(f"  {'':10s}   -> len={len(data)}, "
                      f"range=[{np.min(data):.4g}, {np.max(data):.4g}]")

        # Ion species block
        if "N Z A" in self._raw:
            print("\nIon species (N, Z, A)")
            print("-" * len("Ion species (N, Z, A)"))
            nza = self._raw["N Z A"]
            for i, (N, Z, A) in enumerate(zip(nza["N"], nza["Z"], nza["A"])):
                print(f"  [{i}]  N={N:g}  Z={Z:g}  A={A:g}")

        # Hint about computed quantities
        print("\nComputable (not stored as profiles until you call these):")
        print("-" * len("Computable (not stored as profiles until you call these):"))
        print("  compute_pressure()           -> fills 'ptot'")
        print("  compute_quasineutrality()    -> fills 'nz1' from (ne - ni - nb)/Z_imp")
        print("  compute_zeff()               -> returns (psinorm, Zeff)")
        print("  compute_diamagnetic_rotations(psi)")
        print("                               -> fills omgpp, ommpp, omepp")
        print("  compute_rotation_decomposition(R, Bp, Bt, psi)")
        print("                               -> fills omgeb, ommvb, omevb, er, omghb")


# ===========================================================================
# Convenience functions
# ===========================================================================

def read_geqdsk(filename, cocos=1, nlevels=None, resample="theta",
                extrapolate_edge=True):
    r'''! Read a GEQDSK file and return a `GEQDSKEquilibrium` object

    @param filename Path to the g-file
    @param cocos COCOS convention index (default 1)
    @param nlevels Number of \f$\hat{\psi}\f$ levels for flux-surface analysis
    @param resample Contour resampling method (`"theta"` or `"arc_length"`)
    @param extrapolate_edge If `True` (default), extrapolate \f$p'\f$ and
                            \f$FF'\f$ at the separatrix when the g-file has
                            them forced to zero
    @result `GEQDSKEquilibrium` instance
    '''
    return GEQDSKEquilibrium(
        filename, cocos=cocos, nlevels=nlevels, resample=resample,
        extrapolate_edge=extrapolate_edge,
    )


def read_pfile(filename):
    r'''! Read an Osborne p-file and return a `PFile` object

    @param filename Path to the p-file
    @result `PFile` instance
    '''
    return PFile(filename)
