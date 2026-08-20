'''! Functions and helpers to run varyped equilibrium scan on TokaMaker instance 

Runs equilibrium scan over range of scaling values using profiles from 
input an g-File and p-File. Pedestal of pressure and edge spike of 
current density profiles are scaled self-consistently, and the results 
are saved to new g-files and p-files generated with TokaMaker.

@authors Kevin Clavijo
@date June 2026
'''

from contextlib import contextmanager
import copy
import datetime as dt
from itertools import product
import os
import pprint

import numpy as np

from scipy.constants import constants
from scipy.integrate import trapezoid
from scipy.optimize import least_squares

from matplotlib import cm, pyplot as plt
from matplotlib.colors import Normalize

from OpenFUSIONToolkit.TokaMaker import eqdsk
from OpenFUSIONToolkit.TokaMaker.util import get_jphi_from_GS

MU0 = 4.0e-7 * np.pi
KEV_TO_KPA = constants.eV * 1.0e20
DENOM_TOL = 1.0e-7
PSI_PAD = 1.0e-3
LCFS_PAD = 1.0e-2
FALLBACK_ERRORS = (AttributeError, IndexError, KeyError, TypeError, ValueError, RuntimeError)
_GET_Q_RAVG_INDEX = {"<R>": 0, "<1/R>": 1, "<1/R^2>": 2, "dV/dPsi": 3}


class _NonConvergenceError(ValueError):
    """Solver exhaustion with the tolerance used by the final attempt."""

    def __init__(self, message, nl_tol):
        super().__init__(message)
        self.nl_tol = float(nl_tol)


#=============================================================================
#                       PRESSURE SCALING
#=============================================================================

def _get_pedestal_parameters(psi_n, pprime):
    if psi_n.ndim != 1 or pprime.shape != psi_n.shape or psi_n.size < 3:
        raise ValueError("psi_n and pprime must be matching one-dimensional profiles")
    if not np.all(np.isfinite(psi_n)) or not np.all(np.isfinite(pprime)):
        raise ValueError("psi_n and pprime must contain only finite values")
    if np.any(np.diff(psi_n) <= 0.0):
        raise ValueError("psi_n must be strictly increasing")

    hr_psi = np.linspace(0.0, 1.0, 1000)
    hr_pprime = np.interp(hr_psi, psi_n, pprime)

    prof_edge = 0.75
    hr_psi_slice = hr_psi[np.where(hr_psi >= 0.75)[0]]
    pprime_slice = hr_pprime[hr_psi >= 0.75]
    pprime_peak = pprime_slice[np.argmax(np.abs(pprime_slice))]
    if np.isclose(pprime_peak, 0.0):
        return None
    pprime_sgn = np.sign(pprime_peak)
    pprime_slice = pprime_sgn*pprime_slice

    ppp = np.gradient(pprime_slice, hr_psi_slice)
    roots_ppp_idx = np.where(np.diff(np.sign(ppp)))[0]
    if roots_ppp_idx.size == 0:
        return None
    elif roots_ppp_idx.size == 1:
        psimax_idx = roots_ppp_idx[-1]
    else:
        pppp = np.gradient(ppp, hr_psi_slice)
        roots_ppp_idx = np.compress(pppp[roots_ppp_idx] < 0.0, roots_ppp_idx )
        if roots_ppp_idx.size == 0:
            return None
        psimax_idx = roots_ppp_idx[np.argmax(pprime_slice[roots_ppp_idx])]
    if pprime_slice[-1] > pprime_slice[psimax_idx]:
        psimax_idx = len(pprime_slice) - 1
    psimax_target = hr_psi_slice[psimax_idx]
    psimax_idx = int(np.argmin(np.abs(psi_n - psimax_target)))
    psimax = psi_n[psimax_idx]
    pprime_max = pprime_sgn * pprime[psimax_idx]

    w = np.array([])
    for s,l in [ (0.7862,4.0), (0.5,2.27) ]:
        z = np.where(np.diff(np.sign((pprime_slice-pprime_max*s))))[0]
        if not z.size:
            continue
        z0 = psimax - hr_psi_slice[z]
        ppp1 = ppp[z]
        z00 = z0.copy()
        z0 = np.compress(ppp1*z00 > 0, z0)
        ppp1 = np.compress(ppp1*z00 > 0, ppp1)
        if not z0.size:
            continue
        ppp1min = np.min(abs(0.75*ppp1))
        z0 = np.compress(abs(ppp1)>ppp1min, z0)
        ppp1 = np.compress(abs(ppp1)>ppp1min, ppp1)
        if not z0.size:
            continue
        w = np.concatenate( ( w, abs(z0)*l ) )
    if not w.size:
        return None
    wid = np.mean(w)
    dwid = 0.1 * wid if w.size == 1 else np.std(w)
    return psimax, pprime_max * pprime_sgn, wid, dwid

def _pressure_scaling_function(psi_n, pedestal_factor, norm_factor, pedestal_start):
    psi_n = np.asarray(psi_n, dtype=float)
    u = np.where(
        psi_n < pedestal_start,
        (pedestal_start - psi_n) / pedestal_start,
        0.0,
    )
    norm_value = float(np.asarray(norm_factor, dtype=float).reshape(-1)[0])
    if pedestal_factor <= 1.0:
        return pedestal_factor * np.cosh(norm_value * u)**0.5
    return pedestal_factor / np.cosh(norm_value * u)**0.5

def _pressure_volume_energy(t_object, psi_n, pressure, geometry=None):
    if geometry is None:
        geometry = fsa_current_geometry(
            t_object, psi_n, want_pprime=False,
        )
    pressure = np.asarray(pressure, dtype=float)
    volume_weight = (
        geometry['dV_dpsi'] * geometry['dpsi_dpsiN']
    )
    energy_factor = float(geometry.get('_energy_factor', 1.0))
    if not np.isfinite(energy_factor) or energy_factor <= 0.0:
        energy_factor = 1.0
    return energy_factor * 1.5 * _trapz(
        volume_weight * pressure, geometry['psi_N'],
    )

def _pressure_objective_function(norm_factor, 
                                 t_object, 
                                 psi_n, 
                                 pressure, 
                                 pedestal_factor, 
                                 pedestal_start, 
                                 geometry, 
                                 energy_target=None):
    #may switch back to TokaMaker's compute_flux_integral method for energy calcs
    w_0 = _pressure_volume_energy(t_object, psi_n, pressure, geometry)
    psf = _pressure_scaling_function(psi_n, pedestal_factor, norm_factor, pedestal_start)
    scaled_pressure = pressure * psf
    w_1 = _pressure_volume_energy(t_object, psi_n, scaled_pressure, geometry)
    target = w_0 if energy_target is None else float(energy_target)
    return np.array([w_1 - target])

def scale_pressure(t_object, psi_n, pressure, pedestal_factor, energy_target=None, geometry=None):
    r'''! Scale a pressure profile while preserving a stored-energy target.

    @param t_object TokaMaker equilibrium object
    @param psi_n Normalized flux coordinate
    @param pressure Pressure profile in Pa
    @param pedestal_factor Requested pedestal scale
    @param energy_target Optional stored-energy target in J
    @param geometry Optional FSA geometry for the energy integral
    @result Tuple containing scaled pressure, scaling function, and energy residual
    '''
    psi_n = np.asarray(psi_n, dtype=float)
    pressure = np.asarray(pressure, dtype=float)
    if psi_n.ndim != 1 or pressure.shape != psi_n.shape:
        raise ValueError("psi_n and pressure must be matching one-dimensional profiles")
    if not np.all(np.isfinite(psi_n)) or not np.all(np.isfinite(pressure)):
        raise ValueError("psi_n and pressure must contain only finite values")
    if not np.isfinite(pedestal_factor):
        raise ValueError("pedestal_factor must be finite")
    if energy_target is not None and not np.isfinite(energy_target):
        raise ValueError("energy_target must be finite when provided")
    if np.any(np.diff(psi_n) <= 0.0):
        raise ValueError("psi_n must be strictly increasing")

    if geometry is None:
        geometry = fsa_current_geometry(
            t_object, psi_n, want_pprime=False,
        )

    if np.isclose(pedestal_factor, 1.0):
        psf = np.ones_like(pressure)
        w_0 = _pressure_volume_energy(t_object, psi_n, pressure, geometry)
        target = w_0 if energy_target is None else float(energy_target)
        return pressure.copy(), psf, w_0 - target

    pprime = np.gradient(pressure, psi_n)
    w_0 = _pressure_volume_energy(t_object, psi_n, pressure, geometry)
    pedestal_parameters = _get_pedestal_parameters(psi_n, pprime)
    if pedestal_parameters is None:
        edge_mask = psi_n >= max(psi_n[0], 0.75)
        if np.any(edge_mask):
            edge_indices = np.flatnonzero(edge_mask)
            psimax = float(psi_n[edge_indices[np.argmax(
                np.abs(pprime[edge_mask])
            )]])
        else:
            psimax = float(psi_n[-1])
        grid_width = float(np.median(np.diff(psi_n)))
        pedestal_parameters = (
            psimax, 0.0,
            max(0.1 * float(psi_n[-1] - psi_n[0]), 2.0 * grid_width),
            0.0,
        )
    psimax, _, wid, _ = pedestal_parameters

    pedestal_start = float(psimax - wid)
    pedestal_start = max(
        float(psi_n[0]) + float(np.median(np.diff(psi_n))),
        pedestal_start,
    )
    if pedestal_start >= psi_n[-1]:
        pedestal_start = float(psi_n[-2])

    def energy_residual(norm_factor):
        return float(_pressure_objective_function(
            [norm_factor], t_object, psi_n, pressure,
            pedestal_factor, pedestal_start, geometry, energy_target,
        )[0])

    norm_candidates = np.concatenate(([0.0], np.logspace(-2.0, 2.0, 100)))
    candidate_residuals = np.asarray([
        energy_residual(norm_factor)
        for norm_factor in norm_candidates
    ])
    finite_candidates = np.isfinite(candidate_residuals)
    if not np.any(finite_candidates):
        raise RuntimeError("Unable to evaluate a finite pressure-energy residual")

    finite_indices = np.flatnonzero(finite_candidates)
    best_index = int(finite_indices[np.argmin(
        np.abs(candidate_residuals[finite_candidates])
    )])
    norm_factor = [float(norm_candidates[best_index])]

    fit = least_squares(
        energy_residual,
        norm_factor,
        bounds=(0.0, 100.0),
        ftol=1.0e-9,
        xtol=1.0e-9,
        gtol=1.0e-9,
        max_nfev=100,
    )
    fit_residual = energy_residual(float(fit.x[0]))
    if np.isfinite(fit_residual) and abs(fit_residual) < abs(
            candidate_residuals[best_index]):
        norm_factor = [float(fit.x[0])]

    psf = _pressure_scaling_function(
        psi_n, pedestal_factor, norm_factor, pedestal_start,
    )
    scaled_pressure = pressure * psf
    w_1 = _pressure_volume_energy(t_object, psi_n, scaled_pressure, geometry)
    target = w_0 if energy_target is None else float(energy_target)
    w_diff = w_1 - target
    return scaled_pressure, psf, w_diff

#=============================================================================
#                       CURRENT SCALING
#=============================================================================

def q_ravg(ravgs, which):
    """Extract a flux-surface average from ``get_q`` output."""
    index = _GET_Q_RAVG_INDEX[which]
    if isinstance(ravgs, dict):
        values = ravgs[which]
    else:
        values = np.asarray(ravgs)
        if values.ndim != 2:
            raise ValueError("legacy ravgs must be a two-dimensional array")
        if values.shape[0] == len(_GET_Q_RAVG_INDEX):
            values = values[index, :]
        elif values.shape[1] == len(_GET_Q_RAVG_INDEX):
            values = values[:, index]
        else:
            raise IndexError(f"legacy ravgs has no axis for {which!r}")
    return np.asarray(values, dtype=float)

def _trapz(y, x):
    return float(trapezoid(np.asarray(y, dtype=float),
                           np.asarray(x, dtype=float)))


def fsa_current_geometry(eq, psi_N, psi_pad=PSI_PAD, want_pprime=True):
    r"""Per-surface FSA geometry for :func:`Ip_fsa_integral`, from ``get_q``.

    Returns a dict with ``psi_N``, ``psi_q`` (the clipped sampling grid),
    ``R_avg``, ``inv_R``, ``inv_R2``, ``dV_dpsi`` (magnitude),
    ``dpsi_dpsiN``, ``pprime`` (``None`` if not requested) and
    ``dA_dpsiN = (V'/2pi) <1/R> |dpsi/dpsi_N|``.
    """
    psi_N = np.asarray(psi_N, dtype=float)
    psi_q = np.ascontiguousarray(
        np.clip(psi_N, float(psi_pad), 1.0 - float(psi_pad)), dtype=float)

    ravgs = eq.get_q(psi=psi_q)[2]
    R_avg = q_ravg(ravgs, "<R>")
    inv_R = q_ravg(ravgs, "<1/R>")
    dV_dpsi = np.abs(q_ravg(ravgs, "dV/dPsi"))
    inv_R2 = None
    if isinstance(ravgs, dict) and "<1/R^2>" in ravgs:
        inv_R2 = q_ravg(ravgs, "<1/R^2>")
    elif not isinstance(ravgs, dict):
        inv_R2 = q_ravg(ravgs, "<1/R^2>")

    if R_avg.size > 1 and np.ptp(R_avg) <= 1.0e-9 * float(np.mean(R_avg)):
        raise RuntimeError(
            "fsa_current_geometry: get_q returned a CONSTANT <R> across all "
            f"{R_avg.size} surfaces ({float(R_avg[0]):.6f} m) -- the surface "
            "tracer collapsed onto the magnetic axis.  This is the silent "
            "failure mode triggered by sampling psi_N = 0 exactly; raise "
            "psi_pad.")
    if not (
        np.all(np.isfinite(R_avg))
        and np.all(np.isfinite(inv_R))
        and (inv_R2 is None or np.all(np.isfinite(inv_R2)))
        and np.all(np.isfinite(dV_dpsi))
    ):
        raise RuntimeError("fsa_current_geometry: get_q returned non-finite "
                           "flux-surface averages")

    bounds = np.asarray(eq.psi_bounds, dtype=float)
    dpsi_dpsiN = abs(float(bounds[1]) - float(bounds[0]))
    if not np.isfinite(dpsi_dpsiN) or dpsi_dpsiN <= 0.0:
        raise RuntimeError(f"fsa_current_geometry: bad psi_bounds {bounds}")

    pprime = None
    if want_pprime:
        # .copy() -- get_profiles writes through its psi argument on some builds
        prof = eq.get_profiles(psi=psi_q.copy())
        pprime = np.asarray(prof[4], dtype=float)

    return {
        "psi_N": psi_N,
        "psi_q": psi_q,
        "R_avg": R_avg,
        "inv_R": inv_R,
        "inv_R2": inv_R2,
        "dV_dpsi": dV_dpsi,
        "dpsi_dpsiN": dpsi_dpsiN,
        "pprime": pprime,
        "dA_dpsiN": dV_dpsi / (2.0 * np.pi) * inv_R * dpsi_dpsiN,
    }


def Ip_fsa_weights(geom, convention="jphi-linterp", pprime_sign=1.0):
 
    g = geom["dV_dpsi"] / (2.0 * np.pi) * geom["dpsi_dpsiN"]

    if convention == "fsa":
        return g * geom["inv_R"], 0.0
    if convention != "jphi-linterp":
        raise ValueError(f"unknown convention {convention!r} "
                         "(want 'jphi-linterp' or 'fsa')")

    inv_R2 = geom["inv_R2"]
    if inv_R2 is None:
        raise ValueError(
            "the 'jphi-linterp' measure needs <1/R^2>, which this OFT build's "
            "get_q does not return (legacy positional ravgs).  Use "
            "convention='fsa' on that build.")
    pprime = geom["pprime"]
    if pprime is None:
        raise ValueError("the 'jphi-linterp' measure needs P'; rebuild the "
                         "geometry with want_pprime=True")

    inv_R = np.asarray(geom["inv_R"], dtype=float)
    inv_R_safe = np.where(
        np.abs(inv_R) > DENOM_TOL,
        inv_R,
        np.where(inv_R >= 0.0, DENOM_TOL, -DENOM_TOL),
    )
    ratio = inv_R2 / inv_R_safe
    w = g * ratio
    c = _trapz(g * float(pprime_sign) * pprime
               * (1.0 - geom["R_avg"] * ratio), geom["psi_N"])
    return w, c


def Ip_fsa_integral(eq, psi_N, j_profile, convention="jphi-linterp",
                    psi_pad=PSI_PAD, pprime_sign=1.0, geom=None):
    r"""Plasma current [A] carried by a toroidal current profile *j_profile*."""
    if geom is None:
        geom = fsa_current_geometry(
            eq, psi_N, psi_pad=psi_pad,
            want_pprime=(convention == "jphi-linterp"),
        )
    w, c = Ip_fsa_weights(geom, convention=convention,
                           pprime_sign=pprime_sign)
    return _trapz(w * np.asarray(j_profile, dtype=float), geom["psi_N"]) + c


def eq_jphi_profile(geom, convention="jphi-linterp", eq=None,
                    pprime_sign=1.0):
    r"""Return the equilibrium's own current profile in *convention*."""
    FFp = geom.get("FFp")
    if FFp is None:
        if eq is None:
            raise ValueError("eq_jphi_profile needs eq (or geom['FFp'])")
        prof = eq.get_profiles(psi=geom["psi_q"].copy())
        FFp = np.asarray(prof[1], dtype=float) * np.asarray(prof[2], dtype=float)
    FFp = float(pprime_sign) * np.asarray(FFp, dtype=float)
    pprime = float(pprime_sign) * np.asarray(geom["pprime"], dtype=float)

    if convention == "jphi-linterp":
        return get_jphi_from_GS(
            FFp, pprime, geom["R_avg"], geom["inv_R"],
        )
    if convention == "fsa":
        if geom["inv_R2"] is None:
            raise ValueError("the 'fsa' equilibrium profile needs <1/R^2>")
        return (pprime + geom["inv_R2"] * FFp / MU0) / geom["inv_R"]
    raise ValueError(f"unknown convention {convention!r}")

def _get_current_width(psi_n, jtor, min_width=0.03):

    if psi_n.ndim != 1 or jtor.shape != psi_n.shape or psi_n.size < 2:
        raise ValueError("psi_n and jtor must be matching one-dimensional profiles")
    if not np.all(np.isfinite(psi_n)) or not np.all(np.isfinite(jtor)):
        raise ValueError("psi_n and jtor must contain only finite values")
    if np.any(np.diff(psi_n) <= 0.0):
        raise ValueError("psi_n must be strictly increasing")

    idx = int(np.argmin(np.abs(psi_n - 0.8)))
    j_edge = jtor[idx:]

    peak_rel = int(np.argmax(np.abs(j_edge)))
    peak_idx = idx + peak_rel
    peak = psi_n[peak_idx]

    wid = max(1.0 - peak, min_width)
    return peak, wid

def _edge_basis(psi_n, s, w):

    if w <= 0.0:
        raise ValueError("Edge basis width must be greater than zero")
    w_eff = 1.2*w

    u = np.where(psi_n > s, (psi_n - s) / w_eff, 0.0)
    b = np.where(psi_n > s, 1.0 - 1.0 / np.cosh(np.clip(u, 0.0, 8))**2, 0.0)

    peak = np.max(b)
    if np.isclose(peak, 0.0):
        raise ValueError("Edge basis has no support on psi_n")

    return b / peak

def _core_basis(psi_n, s, w):

    if w <= 0.0:
        raise ValueError("Core basis width must be greater than zero")
    w_eff = 6.*w

    u = np.where(psi_n < s, (s - psi_n) / w_eff, 0.0)
    b = np.where(psi_n < s, 1.0 - 1.0 / np.cosh(np.clip(u, 0.0, 8))**2, 0.0)

    peak = np.max(b)
    if np.isclose(peak, 0.0):
        raise ValueError("Core basis has no support on psi_n")

    return b / peak

def _scale_current_density(t_object, psi_n, jtor, scale_j,
                           denom_tol=DENOM_TOL, pprime=None,
                           target_current=None):
    r'''! Scale the edge current profile while preserving total plasma current.

    @param t_object TokaMaker equilibrium object
    @param psi_n Normalized flux coordinate
    @param jtor Toroidal current profile
    @param scale_j Requested current scale
    @param denom_tol Minimum normalization denominator
    @param pprime Pressure derivative used by the current measure
    @result Tuple containing scaled current and current integrals
    '''
    psi_n = np.asarray(psi_n, dtype=float)
    jtor = np.asarray(jtor, dtype=float)
    if psi_n.ndim != 1 or jtor.shape != psi_n.shape:
        raise ValueError("psi_n and jtor must be matching one-dimensional profiles")
    if not np.all(np.isfinite(psi_n)) or not np.all(np.isfinite(jtor)):
        raise ValueError("psi_n and jtor must contain only finite values")
    if np.any(np.diff(psi_n) <= 0.0):
        raise ValueError("psi_n must be strictly increasing")
    if not np.isfinite(scale_j):
        raise ValueError("scale_j must be finite")
    if target_current is not None and not np.isfinite(target_current):
        raise ValueError("target_current must be finite")

    if pprime is None:
        fsa_geometry = fsa_current_geometry(t_object, psi_n)
    else:
        fsa_geometry = fsa_current_geometry(t_object, psi_n, want_pprime=False)
        pprime = np.asarray(pprime, dtype=float)
        if pprime.shape != np.asarray(psi_n).shape:
            raise ValueError("pprime must match the psi_n profile shape")
        if not np.all(np.isfinite(pprime)):
            raise ValueError("pprime must contain only finite values")
        fsa_geometry['pprime'] = pprime
    try:
        current_weight, current_offset = Ip_fsa_weights(fsa_geometry)
    except ValueError:
        if fsa_geometry.get('inv_R2') is not None:
            raise
        current_weight, current_offset = Ip_fsa_weights(
            fsa_geometry, convention='fsa',
        )
    current_coordinate = fsa_geometry['psi_N']
    i_base = _trapz(current_weight * jtor, current_coordinate) + current_offset

    psimax, wid1 = _get_current_width(psi_n, jtor)
    spike_start = psimax - wid1
    edge_basis = _edge_basis(psi_n, spike_start, wid1)
    core_basis = _core_basis(psi_n, spike_start, wid1)

    scale_delta = float(scale_j) - 1.0
    edge_integral = _trapz(
        current_weight * jtor * edge_basis, current_coordinate,
    )
    core_integral = _trapz(
        current_weight * jtor * core_basis, current_coordinate,
    )
    if abs(i_base) < denom_tol or abs(core_integral) < denom_tol:
        raise ValueError("Current scaling normalization integral is too small")

    core_coefficient = -scale_delta * edge_integral / core_integral
    current_multiplier = (
        1.0 + scale_delta * edge_basis + core_coefficient * core_basis
    )
    shaped_jtor = jtor * current_multiplier
    weighted_base = i_base - current_offset
    if target_current is not None:
        weighted_base = float(target_current) - current_offset
    weighted_shaped = _trapz(
        current_weight * shaped_jtor, current_coordinate,
    )
    if abs(weighted_shaped) < denom_tol:
        raise ValueError("Current scaling result integral is too small")

    scaled_jtor = shaped_jtor * (weighted_base / weighted_shaped)
    i_final = _trapz(
        current_weight * scaled_jtor, current_coordinate,
    ) + current_offset

    return scaled_jtor, i_base, i_final    


def scale_current_density(t_object, psi_n, jtor, scale_j,
                          denom_tol=DENOM_TOL, pprime=None):
    r'''! Scale current using the input profile's own current as target.'''
    return _scale_current_density(
        t_object, psi_n, jtor, scale_j,
        denom_tol=denom_tol, pprime=pprime,
    )

#=============================================================================
#                       GFILE
#=============================================================================

def _gfile_profile(gfile, key):
    for container in (gfile.averages, gfile.geometry):
        if key in container:
            return np.asarray(container[key])
    raise AttributeError(f"gfile object does not have attribute '{key}'")

#resample the needed gfile profiles onto the normalized psi coordinate being used in equilbrium scan
def resample_gfile(gfile, psi_n, eq_snapshot=None):
    r'''! Resample required g-file profiles onto a common flux coordinate.'''

    resampled_profiles = {'psi_n': np.asarray(psi_n, dtype=float)}
    resampled_profiles['psi_n_gfile'] = np.asarray(gfile.psi_N, dtype=float)
    resampled_profiles['Ip'] = float(gfile.Ip)
    resampled_profiles['j_tor_averaged_direct'], resampled_profiles['jtor_flag'] = (
        _resampled_jtor(gfile, resampled_profiles, eq_snapshot)
    )
    resampled_profiles.update(_resampled_geom(gfile, resampled_profiles, eq_snapshot))
    resampled_profiles.update(_resampled_psi(gfile, resampled_profiles))
    resampled_profiles.update(_resampled_Bs(gfile, resampled_profiles, eq_snapshot))
    return resampled_profiles

def _resampled_psi(gfile, resampled_profiles):
    psi_gfile = np.linspace(gfile.psi_axis, gfile.psi_boundary, len(gfile.psi_N))
    psi = np.interp(resampled_profiles['psi_n'], gfile.psi_N, psi_gfile)
    return {'psi': psi, 'dpsi': np.gradient(psi, resampled_profiles['psi_n'])}


def _resampled_jtor(gfile, resampled_profiles, eq_snapshot):
    if not np.isfinite(gfile.Ip) or gfile.Ip <= 0.0:
        raise ValueError(
            "TokaMaker requires a positive Ip target; input g-file "
            f"contains Ip={gfile.Ip}"
        )
    jtor_flag = False
    try:
        if gfile.j_tor_averaged_direct[0] < 0:
            jtor_tmp = -np.interp(
                resampled_profiles['psi_n'], gfile.psi_N,
                gfile.j_tor_averaged_direct,
            )
        else:
            jtor_tmp = np.interp(
                resampled_profiles['psi_n'], gfile.psi_N,
                gfile.j_tor_averaged_direct,
            )
    
    except FALLBACK_ERRORS as exc:
        if eq_snapshot is None:
            raise RuntimeError("resampled_jtor fallback requires eq_snapshot, but eq_snapshot is None") from exc   
        _,f,fp,_,pprime = eq_snapshot.get_profiles(psi=resampled_profiles['psi_n'])
        _, _, ravgs, _, _, _ = eq_snapshot.get_q(psi=resampled_profiles['psi_n']) # get flux averaged R from equilibrium solution
        R_avg = q_ravg(ravgs, '<R>')
        one_over_R_avg = q_ravg(ravgs, '<1/R>')
        jtor_tmp = get_jphi_from_GS(
            f * fp, pprime, R_avg, one_over_R_avg,
        )
        jtor_flag = True
    return jtor_tmp, jtor_flag

def _resampled_geom(gfile, resampled_profiles, eq_snapshot):
    psi_n = resampled_profiles['psi_n']
    try:
        geom = {
            key: np.interp(psi_n, gfile.psi_N, _gfile_profile(gfile, key))
            for key in ('R', 'Z', 'a', 'kappa', 'delta')
        }
    except FALLBACK_ERRORS as exc:
        if eq_snapshot is None:
            raise RuntimeError("resampled_geom fallback requires eq_snapshot, but eq_snapshot is None") from exc

        geom = {
            'R': _fallback_R(resampled_profiles, eq_snapshot),
            'Z': _fallback_Z(resampled_profiles, eq_snapshot),
        }

        geom_tmp = _fallback_geom_contour(eq_snapshot)
        if geom_tmp is None:
            _, _, _, _, rbounds, zbounds = eq_snapshot.get_q(psi=psi_n, compute_geo=True)
            geom_tmp = _fallback_geom_bounds(rbounds, zbounds)

        geom.update({
            key: np.full_like(psi_n, geom_tmp[key])
            for key in ('a', 'kappa', 'delta')
        })
        geom_flag = True
    else:
        geom_flag = False

    geom['geom_flag'] = geom_flag
    return geom


def _fallback_R(resampled_profiles, eq_snapshot):
    try:
        _,_,ravgs,_,_,_ = eq_snapshot.get_q(psi=resampled_profiles['psi_n'])
        return q_ravg(ravgs, '<R>')
    except FALLBACK_ERRORS:
        _,_,_,_,rbounds,_ = eq_snapshot.get_q(psi=resampled_profiles['psi_n'], compute_geo=True)
        lcfs_r = 0.5 * (rbounds[0, 0] + rbounds[1, 0])
        return np.linspace(eq_snapshot.o_point[0], lcfs_r, len(resampled_profiles['psi_n']))
    
def _fallback_Z(resampled_profiles, eq_snapshot):
    
    if eq_snapshot is None:
        raise ValueError("fallback_Z requires eq_snapshot")

    Z = np.full(len(resampled_profiles['psi_n']), eq_snapshot.o_point[1], dtype=float)
    got_any = False

    for i, psi_val in enumerate(resampled_profiles['psi_n']):
        try:
            contour = eq_snapshot.trace_surf(float(psi_val))
        except FALLBACK_ERRORS:
            contour = None
        if contour is not None and len(contour) > 0:
            Z[i] = np.mean(np.asarray(contour)[:, 1])
            got_any = True

    if not got_any:
        _, _, _, _, _, zbounds = eq_snapshot.get_q(psi=resampled_profiles['psi_n'], compute_geo=True)
        lcfs_z = (zbounds[0, 1] + zbounds[1, 1]) / 2.0
        Z = np.linspace(eq_snapshot.o_point[1], lcfs_z, len(resampled_profiles['psi_n']))

    return Z

def _fallback_geom_contour(eq_snapshot):
    if eq_snapshot is None:
        raise ValueError("fallback_geom_contour requires eq_snapshot")

    psi_candidates = [1.0, 0.999, 0.995, 0.99]

    for psi_val in psi_candidates:
        try:
            contour = eq_snapshot.trace_surf(float(psi_val))
        except FALLBACK_ERRORS:
            contour = None

        if contour is None or len(contour) < 4:
            continue

        contour = np.asarray(contour, dtype=float)
        if contour.ndim != 2 or contour.shape[1] < 2:
            continue

        R = contour[:, 0]
        Z = contour[:, 1]

        if not (np.isclose(R[0], R[-1]) and np.isclose(Z[0], Z[-1])):
            R = np.r_[R, R[0]]
            Z = np.r_[Z, Z[0]]

        try:
            geo = eqdsk._flux_geometry(R, Z)
        except FALLBACK_ERRORS:
            continue

        a = float(geo.get('a', 0.0))
        if not np.isfinite(a) or a <= 0.0:
            continue

        return {
            'R0': float(geo.get('R', eq_snapshot.o_point[0])),
            'Z0': float(geo.get('Z', eq_snapshot.o_point[1])),
            'a': a,
            'kappa': float(geo.get('kappa', 1.0)),
            'delta': float(geo.get('delta', 0.0)),
        }

    return None

def _fallback_geom_bounds(rbounds, zbounds):
    rbounds = np.asarray(rbounds, dtype=float)
    zbounds = np.asarray(zbounds, dtype=float)

    r_min = float(rbounds[0, 0])
    r_max = float(rbounds[1, 0])
    z_min = float(zbounds[0, 1])
    z_max = float(zbounds[1, 1])
    r_center = 0.5 * (r_max + r_min)
    z_center = 0.5 * (z_max + z_min)

    dr = r_max - r_min
    eps = 1.0e-12
    if abs(dr) < eps:
        dr = eps if dr >= 0.0 else -eps

    return {
        'R0': r_center,
        'Z0': z_center,
        'a': 0.5 * dr,
        'kappa': (z_max - z_min) / dr,
        'delta': (
            (r_center - 0.5 * (float(zbounds[1, 0]) + float(zbounds[0, 0])))
            * 2.0 / dr
        ),
    }

def get_init_geom(gfile, resampled_profiles, eq_snapshot=None):
    r'''! Return initial solver geometry from a g-file or equilibrium snapshot.

    @param gfile Input g-file object
    @param resampled_profiles Resampled equilibrium profiles
    @param eq_snapshot Optional fallback equilibrium
    @result Dictionary containing initial geometry parameters
    '''
    ginit_flag = False
    try:
        return {
            'R0': gfile.R_mag,
            'Z0': gfile.Z_mag,
            'a': np.asarray(resampled_profiles['a'])[-1],
            'kappa': np.asarray(resampled_profiles['kappa'])[-1],
            'delta': np.asarray(resampled_profiles['delta'])[-1],
            'ginit_flag': ginit_flag
        }
    except FALLBACK_ERRORS:
        pass
    
    ginit_flag = True
    if eq_snapshot is not None:
        geom_tmp = _fallback_geom_contour(eq_snapshot)
        if geom_tmp is not None:
            geom_tmp['ginit_flag'] = ginit_flag
            return geom_tmp
        _, _, _, _, rbounds, zbounds = eq_snapshot.get_q(psi=resampled_profiles['psi_n'],compute_geo=True)
        geom_tmp = _fallback_geom_bounds(rbounds, zbounds)
        geom_tmp['ginit_flag'] = ginit_flag
        return geom_tmp

    raise RuntimeError("Unable to determine init geometry")


def _resampled_Bs(gfile, resampled_profiles, eq_snapshot):
    psi_n = resampled_profiles['psi_n']
    try:
        Bs = {
            key: np.interp(psi_n, gfile.psi_N, _gfile_profile(gfile, key))
            for key in ('Bt', 'Bp')
        }
    except FALLBACK_ERRORS as exc:
        if eq_snapshot is None:
            raise RuntimeError("resampled_Bs fallback requires eq_snapshot, but eq_snapshot is None") from exc
        Bs = _fallback_B(resampled_profiles, eq_snapshot)
        B_flag = True
    else:
        B_flag = False
    Bs['B_flag'] = B_flag
    return Bs

def _fallback_B(resampled_profiles, eq_snapshot):
    profile_size = len(resampled_profiles['psi_n'])
    Bt = np.zeros(profile_size, dtype=float)
    Bp = np.zeros(profile_size, dtype=float)
    contour_ok = np.zeros(len(resampled_profiles['psi_n']), dtype=bool)

    try:
        b_eval = eq_snapshot.get_field_eval("B")  
        for i, psi_val in enumerate(resampled_profiles['psi_n']):
            contour = eq_snapshot.trace_surf(psi_val)        
            if contour is None or len(contour) == 0:
                continue
            b_vals = np.asarray([b_eval.eval(pt) for pt in contour])     
            if b_vals.ndim != 2 or b_vals.shape[1] < 3:
                continue
            Bt_i = np.mean(b_vals[:, 2])
            Bp_i = np.mean(np.sqrt(b_vals[:, 0]**2 + b_vals[:, 1]**2))
            if np.isfinite(Bt_i) and np.isfinite(Bp_i):
                Bt[i] = Bt_i
                Bp[i] = Bp_i
                contour_ok[i] = True
    except FALLBACK_ERRORS:
        contour_ok[:] = False

    if np.all(contour_ok):
        return {"Bt": Bt, "Bp": Bp}
    
    _, qvals, ravgs, _, _, _ = eq_snapshot.get_q(psi=resampled_profiles['psi_n'])
    _, f, _, _, _ = eq_snapshot.get_profiles(psi=resampled_profiles['psi_n'])

    Bt_fb = f * q_ravg(ravgs, '<1/R>')
    R = q_ravg(ravgs, '<R>')
    R0 = float(eq_snapshot.o_point[0])

    eps = DENOM_TOL
    R_safe = np.where(np.abs(R) < eps, np.where(R >= 0.0, eps, -eps), R)
    q_safe = np.where(np.abs(qvals) < eps, np.where(qvals >= 0.0, eps, -eps), qvals)
    r_eff = np.maximum(np.abs(R - R0), eps)

    Bp_fb = np.abs((r_eff * Bt_fb) / (R_safe * q_safe))
    Bp_fb = np.maximum(Bp_fb, eps)

    if np.any(contour_ok):
        Bt = np.where(contour_ok, Bt, Bt_fb)
        Bp = np.where(contour_ok, Bp, Bp_fb)
    else:
        Bt = Bt_fb
        Bp = Bp_fb

    return {"Bt": Bt, "Bp": Bp}

#=============================================================================
#                       PFILE
#=============================================================================

def profile_exists(pfile, key):
    r'''! Check whether a p-file contains a usable profile.

    @param pfile P-file object
    @param key Profile name
    @result True when the requested profile is present and usable
    '''

    if key == "N Z A":
        if "N Z A" not in pfile:
            return False
        nza = pfile["N Z A"]
        if not isinstance(nza, dict):
            return False
        required_keys = ("N", "Z", "A")
        if not all(k in nza for k in required_keys):
            return False
        n_species = len(np.asarray(nza.get("N", []), dtype=float))
        return n_species >= 2

    if key not in pfile:
        return False
    profile = pfile[key]
    if not isinstance(profile, dict) or "data" not in profile:
        return False
    data = np.asarray(profile["data"], dtype=float)
    return data.size > 1 and not np.all(data == 0)

def _update_impurity_species(pfile, resampled_gfile_profiles, psi_n):

    num_impurity_species = 0
    if profile_exists(pfile, 'N Z A') and len(pfile['N Z A']['N']) > 2:
        num_impurity_species = len(pfile['N Z A']['N']) - 2
    else:
        # No impurities: mark n_imps = 0 and return early
        resampled_gfile_profiles['n_imps'] = 0
        resampled_gfile_profiles['ZIs'] = []
        return

    # Update impurity profiles
    for i in range(num_impurity_species):
        for key in [f'vtor{i+1}', f'vpol{i+1}', f'nz{i+1}']:
            if profile_exists(pfile, key):
                old_psi_n = np.asarray(pfile.psinorm_for(key), dtype=float)
                tmp = np.interp(psi_n, old_psi_n, pfile[key]['data'])
            else:
                tmp = np.zeros_like(psi_n)
            pfile.set_profile(key, psi_n, tmp)
    
    nz_fracs = []
    for i in range(num_impurity_species):
        nz_frac = pfile[f'nz{i+1}']['data'] / (np.abs(pfile['ni']['data']) + DENOM_TOL)
        nz_fracs.append(nz_frac)
    
    ZIs = []
    for i in range(num_impurity_species):
        ZIs.append(pfile['N Z A']['Z'][i])
    Zavg = sum(nz_fracs[i]*ZIs[i] for i in range(len(ZIs)))
    
    # Avoid division by zero on profile array
    Zavg = np.where(np.abs(Zavg) < DENOM_TOL, 1.0, Zavg)
    
    nb = (
        np.asarray(pfile['nb']['data'], dtype=float)
        if profile_exists(pfile, 'nb')
        else np.zeros_like(pfile['ni']['data'], dtype=float)
    )
    nz_tot = (pfile['ne']['data'] - pfile['ni']['data'] - nb) / Zavg

    for i in range(num_impurity_species):
        pfile.set_profile(f'nz{i+1}', psi_n, nz_tot * nz_fracs[i])

    resampled_gfile_profiles['n_imps'] = num_impurity_species
    resampled_gfile_profiles['ZIs'] = ZIs

def _update_main_profiles(pfile, psf, psi_n):

    ne = pfile['ne']['data']*(psf**(2./3.))
    pfile.set_profile('ne',psi_n, ne)

    ni = pfile['ni']['data']*(psf**(2./3.))
    pfile.set_profile('ni', psi_n, ni)

    te = pfile['te']['data']*(psf**(1./3.))
    pfile.set_profile('te', psi_n, te)

    ti = pfile['ti']['data']*(psf**(1./3.))
    pfile.set_profile('ti', psi_n, ti)

def _update_fast_ion_profiles(pfile, psf, psi_n):

    if profile_exists(pfile, 'nb'):
        nb = pfile['nb']['data']*(psf**(2./3.))
    else:
        nb = np.zeros_like(psi_n)
    pfile.set_profile('nb', psi_n, nb)
    
    if profile_exists(pfile, 'pb'):
        pb = pfile['pb']['data']*psf
    else:
        pb = np.zeros_like(psi_n)
    pfile.set_profile('pb', psi_n, pb)

def _update_total_pressure(pfile, resampled_gfile_profiles, psi_n):

    num_impurity_species = resampled_gfile_profiles['n_imps']
    
    p_electron = pfile['ne']['data'] * pfile['te']['data'] * KEV_TO_KPA
    nzs = sum(pfile[f'nz{i+1}']['data'] for i in range(num_impurity_species))

    p_ions = (pfile['ni']['data'] + nzs) * pfile['ti']['data'] * KEV_TO_KPA
    p_beam = pfile['pb']['data'] 

    pfile.set_profile('ptot', psi_n, (p_electron + p_ions + p_beam))

def _safe_divide(numerator, denominator, psi_n, eps=1.0e-6):
    num = np.asarray(numerator, dtype=float)
    den = np.asarray(denominator, dtype=float)
    psi = np.asarray(psi_n, dtype=float)

    out = np.full_like(num, np.nan, dtype=float)

    valid = np.isfinite(num) & np.isfinite(den) & (np.abs(den) > eps)
    np.divide(num, den, out=out, where=valid)

    bad = ~np.isfinite(out)
    if np.any(bad):
        good = ~bad
        if np.count_nonzero(good) >= 2:
            out[bad] = np.interp(psi[bad], psi[good], out[good])
        elif np.count_nonzero(good) == 1:
            out[bad] = out[good][0]
        else:
            out[:] = 0.0

    return out


def _update_rotation_profiles(pfile, resampled_gfile_profiles, psi_n):

    Bp = resampled_gfile_profiles['Bp']
    Bt = resampled_gfile_profiles['Bt']
    R = resampled_gfile_profiles['R']
    has_vpol = profile_exists(pfile, 'vpol1')
    has_vtor = profile_exists(pfile, 'vtor1')
    vpol1 = np.asarray(pfile['vpol1']['data'], dtype=float) if has_vpol else None
    vtor1 = np.asarray(pfile['vtor1']['data'], dtype=float) if has_vtor else None

    if has_vpol:
        kpol = _safe_divide(vpol1, Bp, psi_n)
    elif profile_exists(pfile, 'kpol'):
        kpol = np.interp(psi_n, pfile.psinorm_for('kpol'), pfile['kpol']['data'])
    else:
        kpol = np.zeros_like(psi_n, dtype=float)
    pfile.set_profile('kpol', psi_n, kpol)

    if has_vtor:
        omeg = _safe_divide(vtor1, R, psi_n)
    elif profile_exists(pfile, 'omeg'):
        omeg = np.interp(psi_n, pfile.psinorm_for('omeg'), pfile['omeg']['data'])
    else:
        omeg = np.zeros_like(psi_n, dtype=float)
    pfile.set_profile('omeg', psi_n, omeg)

    if has_vpol:
        omegp = _safe_divide(vpol1 * Bt, R * Bp, psi_n)
    elif profile_exists(pfile, 'omegp'):
        omegp = np.interp(psi_n, pfile.psinorm_for('omegp'), pfile['omegp']['data'])
    else:
        omegp = np.zeros_like(psi_n, dtype=float)
    pfile.set_profile('omegp', psi_n, omegp)

def _effective_impurity_density(pfile, num_imps, ZIs):
    
    # nI_eff * Z_ref = sum_i(nz_i * Z_i)
    
    if num_imps <= 0 or len(ZIs) == 0:
        return None

    z_ref = float(ZIs[0])
    if np.isclose(z_ref, 0.0):
        return None

    charge_density = np.zeros_like(pfile['ni']['data'], dtype=float)
    for i in range(num_imps):
        key = f'nz{i+1}'
        if key in pfile:
            charge_density += np.asarray(pfile[key]['data'], dtype=float) * float(ZIs[i])

    return charge_density / z_ref


def _update_diamagnetic_and_rotation_decomposition(pfile, resampled_gfile_profiles):

    psi_n = np.asarray(resampled_gfile_profiles['psi_n'], dtype=float)

    # Seed VxB term in the same way as previous varyped implementation.
    if profile_exists(pfile, 'omeg') and profile_exists(pfile, 'omegp'):
        pfile.set_profile('omgvb', psi_n, pfile['omeg']['data'] + pfile['omegp']['data'])
    else:
        pfile.set_profile('omgvb', psi_n, np.zeros_like(psi_n))

    num_imps = int(resampled_gfile_profiles.get('n_imps', 0))
    ZIs = resampled_gfile_profiles.get('ZIs', [])
    nI_eff = _effective_impurity_density(pfile, num_imps, ZIs)

    psi = np.asarray(resampled_gfile_profiles['psi'], dtype=float)
    ti = np.asarray(pfile['ti']['data'], dtype=float) if profile_exists(pfile, 'ti') else None

    # Robust diamagnetic terms from pfile.py.
    pfile.compute_diamagnetic_rotations(psi=psi, nI=nI_eff, TI=ti)

    # Robust decomposition + er + omghb from pfile.py.
    pfile.compute_rotation_decomposition(
        R=np.asarray(resampled_gfile_profiles['R'], dtype=float),
        Bp=np.asarray(resampled_gfile_profiles['Bp'], dtype=float),
        Bt=np.asarray(resampled_gfile_profiles['Bt'], dtype=float),
        psi=psi,
    )


def make_updated_pfile(pfile, resampled_gfile_profiles, psf, psi_n):
    r'''! Return a deep-copied p-file with all scan profiles updated.

    The returned p-file updates kinetic, fast-ion, impurity, rotation, and
    derived force-balance profiles on the supplied normalized-flux grid.
    '''
    
    pfile_copy = copy.deepcopy(pfile)

    _update_main_profiles(pfile_copy, psf, psi_n)
    _update_fast_ion_profiles(pfile_copy, psf, psi_n)
    _update_impurity_species(pfile_copy, resampled_gfile_profiles, psi_n)
    _update_total_pressure(pfile_copy, resampled_gfile_profiles, psi_n)
    _update_rotation_profiles(pfile_copy, resampled_gfile_profiles, psi_n)
    _update_diamagnetic_and_rotation_decomposition(pfile_copy, resampled_gfile_profiles)

    return pfile_copy

#=============================================================================
#                       EQUILIBRIUM SCAN
#=============================================================================

def _scan_output_name(path, timestamp):
    date_part = timestamp.strftime('%y%m%d')
    time_part = timestamp.second * 100 + timestamp.microsecond // 10000
    for offset in range(10000):
        candidate_time = (time_part + offset) % 10000
        result_name = f'{date_part}_{candidate_time:04d}_OFT_tokamaker'
        out_dir = os.path.join(path, f'varyped{result_name}')
        if not os.path.exists(out_dir):
            return result_name, out_dir
    raise FileExistsError("Unable to allocate a unique VARYPED output directory")

@contextmanager
def _temporary_solver_settings(t_object, urf, nl_tol, maxits):
    original_settings = (t_object.settings.urf, t_object.settings.nl_tol, t_object.settings.maxits)
    try:
        t_object.settings.urf = urf
        t_object.settings.nl_tol = nl_tol
        t_object.settings.maxits = maxits
        t_object.update_settings()
        yield
    finally:
        t_object.settings.urf, t_object.settings.nl_tol, t_object.settings.maxits = original_settings
        t_object.update_settings()


def _solve_baseline_with_retries(t_object, seed_eq):
    last_nl_tol = getattr(t_object.settings, 'nl_tol', 1.0e-6)
    for attempt in range(6):
        nl_tol = min(5.0e-6, (attempt + 1) * 1.0e-6)
        try:
            with _temporary_solver_settings(t_object, 0.3, nl_tol, 125):
                t_object.solve()
            if hasattr(t_object, 'get_stats'):
                stats = t_object.get_stats()
                for key in ('Ip', 'W_MHD'):
                    if not np.isfinite(float(stats[key])):
                        raise RuntimeError(
                            f"TokaMaker returned a non-finite {key}"
                        )
            return
        except FALLBACK_ERRORS + (FloatingPointError, OverflowError) as error:
            last_nl_tol = nl_tol
            t_object.replace_eq(source_eq=seed_eq)
            if attempt == 5:
                raise _NonConvergenceError(str(error), last_nl_tol) from error

    raise _NonConvergenceError(
        "Baseline equilibrium solve exhausted retry attempts", last_nl_tol,
    )

def _prepare_solver_profiles(t_object, ptot, jtor, psi_n, target_scales,
                     Ip_target, energy_target=None, pressure_geometry=None,
                     energy_geometry_factor=1.0, pressure_amplitude=1.0):
    if pressure_geometry is None:
        pressure_geometry = fsa_current_geometry(
            t_object, psi_n, want_pprime=False,
        )
    else:
        pressure_geometry = dict(pressure_geometry)
    if np.isfinite(energy_geometry_factor) and energy_geometry_factor > 0.0:
        pressure_geometry['_energy_factor'] = float(energy_geometry_factor)
    ptot_scaled, psf, _ = scale_pressure(
        t_object, psi_n, ptot, target_scales[0],
        energy_target=energy_target, geometry=pressure_geometry,
    )
    if np.isfinite(pressure_amplitude) and pressure_amplitude > 0.0:
        ptot_scaled = ptot_scaled * float(pressure_amplitude)
        psf = psf * float(pressure_amplitude)
    pprime_scaled = np.gradient(ptot_scaled, psi_n)
    if (np.isclose(target_scales[0], 1.0)
            and np.isclose(target_scales[1], 1.0)):
        jtor_scaled = np.asarray(jtor, dtype=float).copy()
    else:
        jtor_scaled, _, _ = _scale_current_density(
            t_object, psi_n, jtor, target_scales[1],
            pprime=pprime_scaled, target_current=Ip_target,
        )

    pp_prof = {
        'type': 'linterp',
        'x': psi_n,
        'y': pprime_scaled,
    }
    ffp_prof = {
        'type': 'jphi-linterp',
        'x': psi_n,
        'y': jtor_scaled,
    }

    t_object.set_targets(Ip=Ip_target, pax=ptot_scaled[0])
    t_object.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof)
    return psf, target_scales

def _interpolate_scales(target_scales, conv_scales):
    return tuple(
        0.5 * (float(target_value) + float(converged_value))
        for target_value, converged_value
        in zip(target_scales, conv_scales)
    )


def _solve_with_retries(
    t_object,
    ptot,
    jtor,
    psi_n,
    target_scales,
    conv_scales,
    conv_eq,
    Ip_target,
    max_attempts,
    initial_settings=None,
    return_settings=False,
    energy_target=None,
    pressure_geometry=None,
    energy_geometry_factor=1.0,
    pressure_amplitude=1.0,
):
    requested_scales = tuple(target_scales)
    trial_scales = requested_scales
    last_scales = tuple(conv_scales)
    last_eq = conv_eq
    is_recovery = False
    target_attempt = 0
    last_nl_tol = getattr(t_object.settings, 'nl_tol', 1.0e-6)
    last_eq_geometry = None

    while target_attempt < max_attempts:
        solve_settings = initial_settings
        if not is_recovery:
            target_attempt += 1
        if solve_settings is not None:
            _, initial_nl_tol, _ = solve_settings
            if is_recovery:
                solve_settings = (0.3, 1.0e-5, 125)
            elif float(initial_nl_tol) <= 1.0e-6:
                nl_tol = (1.0e-6, 2.0e-6, 5.0e-6)[
                    min(target_attempt - 1, 2)
                ]
                solve_settings = (0.3, nl_tol, 125)
            else:
                solve_settings = (0.3, float(initial_nl_tol), 125)
        if solve_settings is not None:
            last_nl_tol = float(solve_settings[1])
        try:
            use_predicted_geometry = (
                pressure_geometry is not None
                and not is_recovery
                and target_attempt == 1
            )
            if use_predicted_geometry:
                trial_pressure_geometry = pressure_geometry
            else:
                if last_eq_geometry is None:
                    last_eq_geometry = fsa_current_geometry(
                        last_eq, psi_n, want_pprime=False,
                    )
                trial_pressure_geometry = last_eq_geometry
            prepared = _prepare_solver_profiles(
                t_object, ptot, jtor, psi_n, trial_scales, Ip_target,
                energy_target=energy_target,
                pressure_geometry=trial_pressure_geometry,
                energy_geometry_factor=energy_geometry_factor,
                pressure_amplitude=pressure_amplitude,
            )
            if solve_settings is None:
                t_object.solve()
            else:
                with _temporary_solver_settings(t_object, *solve_settings):
                    t_object.solve()
            if hasattr(t_object, 'get_stats'):
                stats = t_object.get_stats()
                for key in ('Ip', 'W_MHD'):
                    if not np.isfinite(float(stats[key])):
                        raise RuntimeError(
                            f"TokaMaker returned a non-finite {key}"
                        )
        except FALLBACK_ERRORS + (FloatingPointError, OverflowError) as error:
            failed_nl_tol = (
                solve_settings[1]
                if solve_settings is not None
                else getattr(t_object.settings, 'nl_tol', 1.0e-6)
            )
            if not is_recovery and target_attempt >= max_attempts:
                t_object.replace_eq(source_eq=last_eq)
                raise _NonConvergenceError(str(error), failed_nl_tol) from error
            t_object.replace_eq(source_eq=last_eq)
            if is_recovery:
                is_recovery = False
                trial_scales = requested_scales
                print(f"RETRYING TARGET: {requested_scales}")
            else:
                recovery_scales = _interpolate_scales(
                    requested_scales, last_scales,
                )
                if (np.allclose(recovery_scales, requested_scales)
                        or np.allclose(recovery_scales, last_scales)):
                    is_recovery = False
                    trial_scales = requested_scales
                    print(f"RETRYING TARGET: {requested_scales}")
                else:
                    is_recovery = True
                    trial_scales = recovery_scales
                    print(f"SETTING MIDPOINT TARGET: {trial_scales}")
            continue

        if np.allclose(trial_scales, requested_scales):
            if return_settings:
                return (*prepared, solve_settings)
            return prepared

        last_scales = tuple(trial_scales)
        last_eq = t_object.copy_eq()
        last_eq_geometry = None
        trial_scales = requested_scales
        is_recovery = False

    raise _NonConvergenceError(
        "Scaled equilibrium solve exhausted retry attempts", last_nl_tol,
    )


def _geometry_extrapolation_factor(anchor_scales, midpoint_scales,
                                   target_scales):
    anchor = np.asarray(anchor_scales, dtype=float)
    midpoint = np.asarray(midpoint_scales, dtype=float)
    target = np.asarray(target_scales, dtype=float)
    midpoint_step = np.linalg.norm(midpoint - anchor)
    target_step = np.linalg.norm(target - anchor)
    if midpoint_step <= DENOM_TOL or target_step <= DENOM_TOL:
        return 1.0
    return float(np.clip(
        target_step / midpoint_step, 1.0, 2,
    ))


def _predicted_pressure_geometry(anchor_eq, midpoint_eq, psi_n,
                                 anchor_scales, midpoint_scales,
                                 target_scales):
    anchor_geometry = fsa_current_geometry(
        anchor_eq, psi_n, want_pprime=False,
    )
    midpoint_geometry = fsa_current_geometry(
        midpoint_eq, psi_n, want_pprime=False,
    )
    factor = _geometry_extrapolation_factor(
        anchor_scales, midpoint_scales, target_scales,
    )
    geometry = {}
    for key, midpoint_value in midpoint_geometry.items():
        anchor_value = anchor_geometry.get(key)
        if anchor_value is None or midpoint_value is None:
            geometry[key] = midpoint_value
            continue
        midpoint_array = np.asarray(midpoint_value, dtype=float)
        anchor_array = np.asarray(anchor_value, dtype=float)
        if midpoint_array.shape != anchor_array.shape:
            geometry[key] = midpoint_value
            continue
        predicted = anchor_array + factor * (midpoint_array - anchor_array)
        if key == 'dV_dpsi':
            predicted = np.maximum(np.abs(predicted), DENOM_TOL)
        elif key == 'dpsi_dpsiN':
            predicted = max(abs(float(predicted)), DENOM_TOL)
        geometry[key] = predicted
    return geometry


def _failed_scan_result(scale_p, scale_j, Ip_target, pax_target,
                         nl_tol=np.nan):
    return {
        'scale_p': scale_p,
        'scale_j': scale_j,
        'gfile_name': '',
        'pfile_name': '',
        'Ip_target': Ip_target,
        'Ip_final': np.nan,
        'pax_target': float(pax_target),
        'W_target': np.nan,
        'W_final': np.nan,
        'nl_tol': float(nl_tol),
        'converged': False,
        'eq_snapshot': None,
        'flags': {'jtor_flag': False, 'geom_flag': False, 'B_flag': False},
    }


def _solve_point(
    t_object, ptot, jtor, psi_n, scale_p, scale_j,
    anchor_scales, anchor_eq, Ip_target,
    energy_target=None, energy_geometry_factor=1.0,
    midpoint_scales=None, last_converged_scales=None):

    target_scales = (scale_p, scale_j)
    continuation_scales = (
        anchor_scales
        if last_converged_scales is None
        else last_converged_scales
    )

    if midpoint_scales is None:
        solve_result = _solve_with_retries(
            t_object, ptot, jtor, psi_n, target_scales,
            continuation_scales, anchor_eq, Ip_target,
            return_settings=True,
            max_attempts=3,
            initial_settings=(0.3, 1.0e-6, 20),
            energy_target=energy_target,
            energy_geometry_factor=energy_geometry_factor,
        )
    else:
        print(f"MID-STEP: (SCALE_P = {midpoint_scales[0]}, "
            f"SCALE_J = {midpoint_scales[1]})")

        _solve_with_retries(
            t_object, ptot, jtor, psi_n, midpoint_scales,
            continuation_scales, anchor_eq, Ip_target,
            max_attempts=2,
            initial_settings=(0.3, 1.0e-5, 20),
            energy_target=energy_target,
            energy_geometry_factor=energy_geometry_factor,
        )

        midpoint_eq = t_object.copy_eq()
        pressure_geometry = _predicted_pressure_geometry(
            anchor_eq, midpoint_eq, psi_n,
            continuation_scales, midpoint_scales, target_scales,
        )
        solve_result = _solve_with_retries(
            t_object, ptot, jtor, psi_n, target_scales,
            midpoint_scales, midpoint_eq, Ip_target,
            return_settings=True,
            max_attempts=3,
            initial_settings=(0.3, 1.0e-6, 20),
            energy_target=energy_target,
            pressure_geometry=pressure_geometry,
            energy_geometry_factor=energy_geometry_factor,
        )

    if energy_target is None or not hasattr(t_object, 'get_stats'):
        return solve_result

    feedback_amplitude = 1.0
    for feedback_attempt in range(2):
        stats = t_object.get_stats()
        actual_energy = float(stats['W_MHD'])
        relative_energy_error = abs(actual_energy - energy_target) / max(
            abs(energy_target), 1.0,
        )
        if not np.isfinite(relative_energy_error) or relative_energy_error <= 1.0e-2:
            break

        feedback_eq = t_object.copy_eq()
        feedback_amplitude *= float(np.clip(
            energy_target / max(actual_energy, 1.0),
            0.5, 2.0,
        ))
        feedback_amplitude = float(np.clip(feedback_amplitude, 0.25, 4.0))
        try:
            feedback_geometry = fsa_current_geometry(
                feedback_eq, psi_n, want_pprime=False,
            )
            solve_result = _solve_with_retries(
                t_object, ptot, jtor, psi_n, target_scales,
                target_scales, feedback_eq, Ip_target,
                return_settings=True,
                max_attempts=3,
                initial_settings=(0.3, 1.0e-6, 100),
                energy_target=energy_target,
                pressure_geometry=feedback_geometry,
                energy_geometry_factor=energy_geometry_factor,
                pressure_amplitude=feedback_amplitude,
            )
        except _NonConvergenceError as error:
            t_object.replace_eq(source_eq=feedback_eq)
            break

    return solve_result


def _save_scan_outputs(
    t_object, pfile_copy, psi_n, psf, eq_snapshot,
    scan_index, yyyymm, out_dir, file_output):
    if not file_output:
        return {'jtor_flag': False, 'geom_flag': False, 'B_flag': False}

    output_gfile_path = os.path.join(out_dir, f'g{yyyymm}.{scan_index:05d}')
    t_object.save_eqdsk(output_gfile_path,
                        nr=len(psi_n),
                        nz=len(psi_n),
                        truncate_eq=False,
                        lcfs_pad=LCFS_PAD,
                        cocos=7)

    output_gfile = eqdsk.read_geqdsk(output_gfile_path, cocos=7)
    output_profiles = resample_gfile(output_gfile, psi_n, eq_snapshot=eq_snapshot)
    
    output_pfile = make_updated_pfile(pfile_copy, output_profiles, psf, psi_n)
    output_pfile_path = os.path.join(out_dir, f'p{yyyymm}.{scan_index:05d}')
    print(f'Saving pFile: {output_pfile_path}\n')
    output_pfile.save(output_pfile_path)

    return {'jtor_flag': output_profiles['jtor_flag'],
            'geom_flag': output_profiles['geom_flag'],
            'B_flag': output_profiles['B_flag']}


def equilibrium_scan(t_object, pfile, gfile, scaling_values, path_to_output=None,
                     file_output=True, result_output=True):
    r'''! Run a VARYPED scan and restore the caller state on exit.'''
    initial_eq = t_object.copy_eq()
    try:
        return _run_equilibrium_scan(
            t_object, pfile, gfile, scaling_values, path_to_output,
            file_output, result_output,
        )
    finally:
        t_object.replace_eq(source_eq=initial_eq)


def _run_equilibrium_scan(t_object, pfile, gfile, scaling_values,
                           path_to_output=None, file_output=True,
                           result_output=True):
    results = []
    result_errors = []

    #ensure scaling_values is array-like
    if isinstance(scaling_values, (int, float)):
        scaling_values = [scaling_values]
    scaling_array = np.asarray(scaling_values, dtype=float)
    if scaling_array.ndim != 1 or scaling_array.size == 0:
        raise ValueError("scaling_values must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(scaling_array)) or np.any(scaling_array <= 0.0):
        raise ValueError("scaling_values must contain finite positive values")
    if np.any(np.diff(scaling_array) < 0.0):
        raise ValueError("scaling_values must be non-decreasing")

    time = dt.datetime.now()
    yyyymm = time.strftime('%Y%m')

    if path_to_output is None:
        path_to_output = os.getcwd()
    else:
        path_to_output = os.path.abspath(path_to_output)
    result_name, out_dir = _scan_output_name(path_to_output, time)
    if file_output or result_output:
        os.makedirs(out_dir, exist_ok=True)

   #Extract psi_n from pfile 
    if not profile_exists(pfile, 'ptot'):
        raise ValueError("pfile must contain 'ptot' (total pressure) profile")
    psi_n_pfile = pfile.psinorm_for('ptot')
    if psi_n_pfile is None:
        raise ValueError("Could not extract psinorm for 'ptot' profile")
    psi_n_pfile = np.asarray(psi_n_pfile, dtype=float)
    
    # Remap all pfile profiles onto ptot's psi_n 
    # ensures all profiles share a common coordinate system
    pfile_copy = pfile.remap(psi_n_pfile)
    
    # Get psi_n from gfile 
    psi_n_gfile = np.asarray(gfile.psi_N, dtype=float)
    
    # Select unified psi_n as one with higher resolution
    # If resolutions are equal, prefer gfile's coordinate system
    if len(psi_n_pfile) > len(psi_n_gfile):
        psi_n = psi_n_pfile
    else:
        psi_n = psi_n_gfile
        pfile_copy = pfile_copy.remap(psi_n)

    # eqdsk.PFile stores ptot in kPa; TokaMaker targets and P' profiles use Pa.
    ptot = np.asarray(pfile_copy['ptot']['data'], dtype=float) * 1.0e3
    
    # Resample all needed gfile profiles onto unified psi_n.
    input_eq = t_object.copy_eq()
    gfile_profiles = resample_gfile(gfile, psi_n, input_eq)
    refit_jtor = gfile_profiles['j_tor_averaged_direct']

    init_geom = get_init_geom(gfile, gfile_profiles, input_eq)
    R0 = init_geom['R0']
    Z0 = init_geom['Z0']
    a = init_geom['a']
    kappa = init_geom['kappa']
    delta = init_geom['delta']

    # Get Ip target
    Ip_target = gfile.Ip
    pax_target = float(ptot[0])

    t_object.init_psi(R0, Z0, a, kappa, delta)

    init_pp_prof = {'type':'linterp',
                    'y': np.gradient(ptot, psi_n),
                    'x': psi_n}
    
    init_ffp_prof = {'type':'jphi-linterp',
                    'y': refit_jtor,
                    'x': psi_n}
    
    t_object.set_targets(Ip=Ip_target, pax=pax_target)

    t_object.set_profiles(ffp_prof=init_ffp_prof, pp_prof=init_pp_prof)
    baseline_seed = t_object.copy_eq()
    
    print("=============== STARTING BASELINE SOLVE ===============")
    _solve_baseline_with_retries(t_object, baseline_seed)
    baseline_solve = t_object.copy_eq()
    W_target = float(baseline_solve.get_stats()['W_MHD'])
    energy_geometry_factor = 1.0
    try:
        baseline_geometry = fsa_current_geometry(
            baseline_solve, psi_n, want_pprime=False,
        )
        baseline_pressure = np.asarray(
            baseline_solve.get_profiles(psi=psi_n.copy())[3],
            dtype=float,
        )
        baseline_fsa_energy = _pressure_volume_energy(
            baseline_solve, psi_n, baseline_pressure, baseline_geometry,
        )
        if (np.isfinite(baseline_fsa_energy)
                and baseline_fsa_energy > 0.0
                and np.isfinite(W_target)
                and W_target > 0.0):
            energy_geometry_factor = W_target / baseline_fsa_energy
    except FALLBACK_ERRORS:
        energy_geometry_factor = 1.0

    results.append({'input_gfile': gfile,
                    'input_pfile': pfile,
                    'baseline_eq': baseline_solve ,
                    'scaling_values': scaling_values,
                    'g_init_flag': init_geom['ginit_flag']
                    })
    
    teqs = 0
    print("\n=============== STARTING VARYPED SCAN ===============")
    print(f'Scaling Values Array: {scaling_values}')
    
    scan_points = product(scaling_values, scaling_values)
    # The baseline equilibrium is the unscaled (1, 1) continuation state.
    last_conv = (1.0, 1.0)
    last_conv_eq = t_object.copy_eq()
    row_anchor_scales = last_conv
    row_anchor_eq = last_conv_eq
    converged_anchors = [(last_conv, last_conv_eq)]

    prev_scale_p = None

    for scale_p, scale_j in scan_points:
        print(f"SOLVE #{teqs}: (SCALE_P = {scale_p}, SCALE_J = {scale_j})")

        new_pressure_row = scale_p != prev_scale_p
        pressure_transition = (
            prev_scale_p is not None and new_pressure_row
        )
        prev_scale_p = scale_p
        target_scales = np.asarray((scale_p, scale_j), dtype=float)
        pressure_midpoint = None
        if teqs > 0:
            if pressure_transition:
                pressure_midpoint = (
                    0.5 * (float(last_conv[0]) + float(scale_p)),
                    float(scale_j),
                )
            else:
                pressure_midpoint = (
                    float(scale_p),
                    0.5 * (float(last_conv[1]) + float(scale_j)),
                )
        anchor_scales, anchor_eq = min(
            converged_anchors,
            key=lambda item: np.linalg.norm(
                np.log(target_scales / np.asarray(item[0], dtype=float))
            ),
        )
        candidate_anchors = [(anchor_scales, anchor_eq)]
        for candidate in (
            (row_anchor_scales, row_anchor_eq),
            (last_conv, last_conv_eq),
            ((1.0, 1.0), results[0]['baseline_eq']),
        ):
            if not any(np.allclose(candidate[0], existing[0])
                       for existing in candidate_anchors):
                candidate_anchors.append(candidate)

        solve_error = None
        solved = None
        for candidate_scales, candidate_eq in candidate_anchors:
            t_object.replace_eq(source_eq=candidate_eq)
            try:
                solved = _solve_point(
                    t_object, ptot, refit_jtor, psi_n, scale_p, scale_j,
                    candidate_scales, candidate_eq, Ip_target,
                    energy_target=W_target,
                    energy_geometry_factor=energy_geometry_factor,
                    midpoint_scales=pressure_midpoint,
                    last_converged_scales=last_conv,
                )
            except _NonConvergenceError as error:
                solve_error = error
                continue
            except (RuntimeError, ValueError) as error:
                solve_error = _NonConvergenceError(
                    str(error),
                    float(getattr(t_object.settings, 'nl_tol', 5.0e-6)),
                )
                continue
            else:
                anchor_scales = candidate_scales
                anchor_eq = candidate_eq
                break

        if solved is None:
            error = solve_error
            t_object.replace_eq(source_eq=last_conv_eq)
            results.append(_failed_scan_result(
                scale_p, scale_j, Ip_target, pax_target,
                nl_tol=error.nl_tol,
            ))
            result_errors.append(error.nl_tol)
            print(
                f"NON-CONVERGED: (SCALE_P = {scale_p}, "
                f"SCALE_J = {scale_j}); nl_tol = {error.nl_tol}"
            )
            teqs += 1
            continue

        psf, solved_scales, solve_settings = solved

        last_conv = solved_scales
        last_conv_eq = t_object.copy_eq()
        converged_anchors.append((tuple(solved_scales), last_conv_eq))
        if new_pressure_row:
            row_anchor_scales = last_conv
            row_anchor_eq = last_conv_eq
        eq_snapshot = t_object.copy_eq()

        flags = _save_scan_outputs(
            t_object, pfile_copy, psi_n, psf, eq_snapshot,
            teqs, yyyymm, out_dir, file_output,
        )
        nl_tol = (
            float(solve_settings[1])
            if solve_settings is not None
            else float(getattr(t_object.settings, 'nl_tol', 1.0e-6))
        )
        result_errors.append(nl_tol)
        eq_stats = eq_snapshot.get_stats()
        results.append({
            'scale_p': scale_p,
            'scale_j': scale_j,
            'gfile_name': f'g{yyyymm}.{teqs:05d}',
            'pfile_name': f'p{yyyymm}.{teqs:05d}',
            'Ip_target': Ip_target,
            'Ip_final': eq_stats['Ip'],
            'pax_target': float(pax_target),
            'W_target': W_target,
            'W_final': float(eq_stats['W_MHD']),
            'nl_tol': nl_tol,
            'converged': True,
            'eq_snapshot': eq_snapshot,
            'flags': flags,
        })

        teqs += 1

    results[0]['out_dir'] = out_dir
    results[0]['result_filename'] = result_name
    print("\nPre-run Info:")
    pprint.pprint(results[0])
    for i in range(1, len(results)):
        print(
            f"\nVaryped Equilibrium Solve #{i - 1} "
            f"({results[i]['scale_p']}, {results[i]['scale_j']}) Results:"
        )
        pprint.pprint(results[i])
    if result_output:
        make_results_file(results, result_errors=result_errors)
    
    return results


def make_results_file(results, result_errors=None):
    r'''! Write the VARYPED text results file.

    @param results Scan result dictionaries, including the baseline entry.
    @param result_errors Optional private values for the ``err`` column.
    '''
    time = dt.datetime.now()
    time_run = time.ctime()

    file_dir = results[0]['out_dir']
    file_name = 'results_'+results[0]['result_filename']
    if result_errors is not None and len(result_errors) != len(results) - 1:
        raise ValueError(
            "result_errors must contain one value for each scan result"
        )

    lines = [f'OFT TokaMaker VARYPED Results; Run on: {time_run}', 
             '  i pres w nustar den temp cur ne_shift te_shift pe_shift err convrg']

    for i, res in enumerate(results[1:]):
        error_value = (
            result_errors[i]
            if result_errors is not None
            else 1.0e-6
        )
        try:
            error_value = float(error_value)
        except (TypeError, ValueError):
            error_value = 1.0e-6
        if not np.isfinite(error_value):
            error_value = 1.0e-6
        converged = int(bool(res.get('converged', False)))
        line = (
            f"  {i} {res['scale_p']} 1.0 1.0 1.0 1.0 "
            f"{res['scale_j']} 0.00 0.00 0.00 "
            f"{error_value:.6g} {converged}"
        )
        lines.append(line)
    
    txt = "\n".join(lines) + "\n"
    with open(os.path.join(file_dir, file_name), 'w') as file:
        file.write(txt)

def summarize_scan(results, n_best=5):
    r'''! Print a compact run summary for a VARYPED scan.

    @param results Result structure returned by :func:`equilibrium_scan`.
    @param n_best Number of best converged points to print.
    '''
    try:
        import pandas as pd
    except ImportError as error:
        raise ImportError(
            "summarize_scan requires pandas to build its quality table"
        ) from error

    if not isinstance(results, (list, tuple)) or not results:
        raise ValueError("results must be the non-empty output of equilibrium_scan")

    scan_points = [point for point in results[1:] if isinstance(point, dict)]

    def _stats(values):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            return np.nan, np.nan, np.nan
        return float(np.mean(values)), float(np.median(values)), float(np.max(values))

    def _relative_error(final_value, target_value):
        if not np.isfinite(final_value) or not np.isfinite(target_value):
            return np.nan
        return abs(final_value - target_value) / max(abs(target_value), 1.0)

    ip_target = next(
        (float(point.get('Ip_target', np.nan)) for point in scan_points
         if np.isfinite(float(point.get('Ip_target', np.nan)))),
        np.nan,
    )
    w_target = next(
        (float(point.get('W_target', np.nan)) for point in scan_points
         if np.isfinite(float(point.get('W_target', np.nan)))),
        np.nan,
    )

    rows = []
    for index, point in enumerate(scan_points):
        ip_final = float(point.get('Ip_final', np.nan))
        w_final = float(point.get('W_final', np.nan))
        rows.append({
            'scan_index': index,
            'scale_p': float(point.get('scale_p', np.nan)),
            'scale_j': float(point.get('scale_j', np.nan)),
            'converged': bool(point.get('converged', False)),
            'nl_tol': float(point.get('nl_tol', np.nan)),
            'Ip_final': ip_final,
            'Ip_rel_err': _relative_error(ip_final, ip_target),
            'W_final': w_final,
            'W_rel_err': _relative_error(w_final, w_target),
        })

    columns = [
        'scan_index', 'scale_p', 'scale_j', 'converged', 'nl_tol',
        'Ip_final', 'Ip_rel_err', 'W_final', 'W_rel_err',
    ]
    quality_table = pd.DataFrame(rows, columns=columns)
    if quality_table.empty:
        best_points = quality_table.copy()
        converged_table = quality_table.copy()
    else:
        converged_table = quality_table.loc[
            quality_table['converged']
        ].copy()
        best_points = converged_table.sort_values(
            ['W_rel_err', 'Ip_rel_err', 'nl_tol'],
            na_position='last',
        ).reset_index(drop=True)

    ip_mean, ip_median, ip_max = _stats(converged_table['Ip_rel_err'])
    w_mean, w_median, w_max = _stats(converged_table['W_rel_err'])
    nl_mean, nl_median, nl_max = _stats(quality_table['nl_tol'])
    n_points = len(quality_table)
    n_converged = int(quality_table['converged'].sum())

    def _format(value):
        return 'n/a' if not np.isfinite(value) else f'{value:.3e}'

    print('=' * 90)
    print('VARYPED RUN QUALITY SUMMARY')
    print('=' * 90)
    print(f"Scan points:                 {n_points}")
    print(f"Converged:                   {n_converged}")
    print(f"Failed:                      {n_points - n_converged}")
    print(
        f"Convergence fraction:        "
        f"{n_converged / n_points if n_points else np.nan:.3f} "
        f"({100.0 * n_converged / n_points if n_points else np.nan:.1f}%)"
    )
    print(f"Ip target:                   {_format(ip_target)} A")
    print(f"W target:                    {_format(w_target)} J")
    print()
    print(
        f"Current error mean/median/max: "
        f"{_format(ip_mean)} / {_format(ip_median)} / {_format(ip_max)}"
    )
    print(
        f"Energy error mean/median/max:  "
        f"{_format(w_mean)} / {_format(w_median)} / {_format(w_max)}"
    )
    print(
        f"nl_tol mean/median/max:         "
        f"{_format(nl_mean)} / {_format(nl_median)} / {_format(nl_max)}"
    )
    print()
    print('Best converged scan points by energy, then current error:')
    print(best_points.head(n_best)[[
        'scale_p', 'scale_j', 'Ip_rel_err', 'W_rel_err',
        'nl_tol', 'Ip_final', 'W_final',
    ]].to_string(index=False))
    print('=' * 90)


#=============================================================================
#                       PLOTTING
#=============================================================================

def _find_results_file(path):
    for entry in sorted(os.scandir(path), key=lambda item: item.name):
        if entry.is_file() and entry.name.startswith('r'):
            return entry.path
    raise FileNotFoundError(f"No VARYPED results file found in '{path}'")


def _read_scan_values(path, usecols):
    values = np.asarray(
        np.genfromtxt(path, skip_header=2, usecols=usecols),
        dtype=float,
    )
    if np.isscalar(usecols):
        return np.atleast_1d(values)
    values = np.atleast_2d(values)
    return tuple(values[:, index] for index in range(values.shape[1]))


def _converged_scan_results(results):
    return [
        result for result in results[1:]
        if result.get('converged', result.get('eq_snapshot') is not None)
    ]


def plot_pressure_results(results, savefig=False):
    r'''! Plot pressure and pressure-derivative scan results.'''
    scaling_values = np.asarray(results[0]["scaling_values"])
    vmin=(scaling_values.min())
    vmax=(scaling_values.max())
    norm = Normalize(vmin=vmin, vmax=vmax)    
    cmap = cm.plasma
    fig, ax = plt.subplots()

    for result in _converged_scan_results(results):
        color = cmap(norm(float(result['scale_p'])))    
        psi_n, _, _, ptot, _ = result['eq_snapshot'].get_profiles(
            npsi=514, psi_pad=PSI_PAD)
        ax.plot(psi_n, ptot, color = color)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_ticks(scaling_values)
    cbar.set_label("Pressure Scale Factor")

    baseline = results[0]['baseline_eq']
    psi_b, _, _, ptot_b, pprime_b = baseline.get_profiles(
        npsi=514, psi_pad=PSI_PAD)

    ax.plot(psi_b, ptot_b, color = 'black', label = 'Input Pressure Profile', linestyle = '--', linewidth = 2.5)

    plt.grid()
    plt.legend()
    plt.xlabel('Normalized ψ')
    plt.ylabel('Total Pressure [Pa]')
    plt.title('Scaled Pressures')
    if savefig:
        fig_out_dir = results[0]['out_dir']
        plt.savefig(os.path.join(fig_out_dir, 'scaled_pressure_scan.png'))
    plt.show()
    plt.close(fig)

    fig, ax = plt.subplots()

    for result in _converged_scan_results(results):
        color = cmap(norm(float(result['scale_p'])))    
        psi_n, _, _, _, pprime = result['eq_snapshot'].get_profiles(
            npsi=514, psi_pad=PSI_PAD)
        ax.plot(psi_n, -pprime, color = color)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_ticks(scaling_values)
    cbar.set_label("Pressure Scale Factor")

    ax.plot(psi_b, -pprime_b, color = 'black', label = 'Initial Pressure Derivative Profile', linestyle = '--', linewidth = 2.5)

    plt.grid()
    plt.legend()
    plt.xlabel('Normalized ψ')
    plt.ylabel("P' [Pa/Wb]")
    plt.title('Scaled Pressure Derivatives')
    if savefig:
        fig_out_dir = results[0]['out_dir']
        plt.savefig(os.path.join(fig_out_dir, 'scaled_pprime_scan.png'))
    plt.show()
    plt.close(fig)
    
def plot_current_results(results, savefig=False):
    r'''! Plot toroidal current-density scan results.'''
    scaling_values = np.asarray(results[0]["scaling_values"])
    vmin=(scaling_values.min())
    vmax=(scaling_values.max())
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.plasma
    fig, ax = plt.subplots()


    for result in _converged_scan_results(results):
        color = cmap(norm(float(result['scale_j'])))    
        
        psi_n, jtor = _jtor_from_eq(result['eq_snapshot'])

        ax.plot(psi_n, jtor, color = color)
    
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_ticks(scaling_values)
    cbar.set_label("Current Density Scale Factor")

    baseline = results[0]['baseline_eq']
    psi_b, baseline_jtor = _jtor_from_eq(baseline)

    ax.plot(psi_b, baseline_jtor, color = 'black', label = 'Initial Current Density Profile', linestyle = '--', linewidth = 2.5)
    
    plt.grid()
    plt.legend()
    plt.xlabel('Normalized ψ')
    plt.ylabel('Current Density [Ampere/m^2]')
    plt.title('Scaled Current Densities')
    if savefig:
        fig_out_dir = results[0]['out_dir']
        plt.savefig(os.path.join(fig_out_dir, 'scaled_current_scan.png'))
    plt.show()
    plt.close(fig)


def plot_varyped_results(results, savefig=False):
    r'''! Plot pressure and current profiles from a VARYPED scan.

    @param results Result structure returned by :func:`equilibrium_scan`
    @param savefig Save generated figures in the scan output directory?
    '''
    plot_pressure_results(results, savefig=savefig)
    plot_current_results(results, savefig=savefig)


def _jtor_from_eq(eq, npsi=514, psi_pad=PSI_PAD):
    psi_n, f, fp, _, pprime = eq.get_profiles(
        npsi=npsi, psi_pad=psi_pad,
    )
    _, _, ravgs, _, _, _ = eq.get_q(
        npsi=npsi, psi_pad=psi_pad,
    )
    R_avg = q_ravg(ravgs, '<R>')
    one_over_R_avg = q_ravg(ravgs, '<1/R>')
    return psi_n, get_jphi_from_GS(
        f * fp, pprime, R_avg, one_over_R_avg,
    )


def plot_pfiles(dir_name, profiles, scaling_range=None, savefig=False):
    r'''! Plot selected profiles from p-files in a scan directory.'''

    if profiles == 'all':
        profiles = ['all']
    elif not isinstance(profiles, list):
        profiles = [profiles]

    cwd = os.getcwd()
    path = os.path.join(cwd, dir_name)
    pfiles = []
    for entry in sorted(os.scandir(path), key=lambda entry: entry.name):
        if entry.is_file() and entry.name.startswith('p'):
            pfiles.append(eqdsk.read_pfile(entry))

    if not pfiles:
        raise ValueError(f"No p-files found in '{path}'")
    results_file = _find_results_file(path)
    scale_values, convergence = _read_scan_values(results_file, (1, 11))
    convergence = convergence.astype(bool)
    scaling_values = scale_values[convergence]
    if len(pfiles) != len(scaling_values):
        raise ValueError(
            f"Found {len(pfiles)} p-files but {len(scaling_values)} "
            "converged scan results"
        )
    
    if scaling_range is not None:
        scale_min, scale_max = sorted(map(float, scaling_range))
        keep = (scaling_values >= scale_min) & (scaling_values <= scale_max)
        pfiles = [pfile for pfile, ok in zip(pfiles, keep) if ok]
        scaling_values = scaling_values[keep]
        if len(pfiles) == 0:
            raise ValueError(f"No pfiles found with scaling values in range [{scale_min}, {scale_max}]")

    
    vmin=(scaling_values.min())
    vmax=(scaling_values.max())
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.plasma
    
    if profiles == ['all']:
        profiles = [key for key in pfiles[0].keys if key != 'N Z A']
    
    print("Plotting: ",list(profiles))
  
    for profile in profiles:
        fig, ax = plt.subplots()
        for i, pfile in enumerate(pfiles):
            color = cmap(norm(float(scaling_values[i]))) 
            ax.plot(pfile.psinorm_for(f'{profile}'), pfile[f'{profile}']['data'], color=color)

        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.02)
        cbar.set_ticks(scaling_values)
        cbar.set_label("Scale Factor")
        plt.grid()
        plt.xlabel('Normalized ψ')
        plt.ylabel(f'{profile}')
        plt.title(f'\'{profile}\' Profiles')
        if savefig:
            plt.savefig(os.path.join(path, f'{profile}_scan.png'))
        plt.show()
        plt.close(fig)


def plot_gfiles(
    dir_name, profiles, container=None, scaling_range=None, savefig=False,
    cocos=None,
):
    r'''! Plot selected profiles from g-files in a scan directory.'''

    if profiles == 'all':
        profiles = ['all']
    elif not isinstance(profiles, list):
        profiles = [profiles]
    container_dicts = ['geometry', 'averages', 'midplane']
    not_plot = ['Ip', 'R_mag', 'Z_mag', 'R_center', 'B_center', 'psi_axis', 'psi_boundary', 'cocos', 'psi_RZ', 'psi_N_RZ', 'limiter_R', 'limiter_Z', 'boundary_R', 'boundary_Z', 'R_grid', 'Z_grid', 'li', 'betas', 'contours']

    cwd = os.getcwd()
    path = os.path.join(cwd, dir_name)
    if cocos is None:
        directory_name = os.path.basename(os.path.normpath(path))
        cocos = 7 if directory_name.startswith('varyped') else 1
    gfiles = []
    for entry in sorted(os.scandir(path), key = lambda e: e.name):
        if entry.is_file() and entry.name.startswith('g'):
            gfiles.append(eqdsk.read_geqdsk(entry, cocos=cocos))

    if not gfiles:
        raise ValueError(f"No g-files found in '{path}'")
    results_file = _find_results_file(path)
    scale_p_all, scale_j_all, convergence = _read_scan_values(
        results_file, (1, 6, 11)
    )
    convergence = convergence.astype(bool)
    scale_p_values = scale_p_all[convergence]
    scale_j_values = scale_j_all[convergence]
    if len(gfiles) != len(scale_p_values):
        raise ValueError(
            f"Found {len(gfiles)} g-files but {len(scale_p_values)} "
            "converged scan results"
        )

    if profiles == ['all'] and container is not None:
        c = getattr(gfiles[0], container)
        profiles = list(c.keys())
    elif profiles == ['all']:
        profiles = [key for key in gfiles[0].keys if key not in container_dicts + not_plot]

    print("Plotting: ", list(profiles))

    for profile in profiles:
        profile_prefix = profile[0].lower()
        uses_current_scale = (
            profile_prefix in ('q', 'j', 'i', 'v', 'd', 'k')
            or (container == 'averages' and profile.startswith('Bp'))
        )
        scaling_values = (
            scale_j_values if uses_current_scale else scale_p_values
        )
        plot_gfiles_subset = gfiles

        if scaling_range is not None:
            scale_min, scale_max = sorted(map(float, scaling_range))
            keep = (scaling_values >= scale_min) & (scaling_values <= scale_max)
            plot_gfiles_subset = [gfile for gfile, ok in zip(gfiles, keep) if ok]
            scaling_values = scaling_values[keep]
            if not plot_gfiles_subset:
                raise ValueError(
                    f"No g-files found with scaling values in range "
                    f"[{scale_min}, {scale_max}]"
                )

 

        vmin = scaling_values.min()
        vmax = scaling_values.max()
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = cm.plasma

        fig, ax = plt.subplots()
        for i, gfile in enumerate(plot_gfiles_subset):
            color = cmap(norm(float(scaling_values[i])))
            if container is not None:
                values = getattr(gfile, container)[profile]
            else:
                values = getattr(gfile, profile)
            if profile_prefix == 'j' and values[0] < 0:
                values = -values
            ax.plot(gfile.psi_N, values, color=color)

        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.02)
        cbar.set_ticks(scaling_values)
        cbar.set_label("Scale Factor")
        plt.grid()
        plt.xlabel('Normalized ψ')
        plt.ylabel(profile)
        plt.title(f'\'{profile}\' Profiles')
        if savefig:
            plt.savefig(os.path.join(path, f'{profile}_scan.png'))
        plt.show()
        plt.close(fig)