#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-only
"""Compute a fixed- or free-boundary Grad-Shafranov equilibrium with TokaMaker
from a JSON case definition.

Usage:
  python run_equilibrium.py case.json [--output-dir DIR] [--nthreads N] [--no-plot]

See `../reference.md` for the full case-file schema and `../examples/` for
complete examples.
"""
import argparse
import json
import os
import sys

import numpy


def _load_tokamaker():
    """Import the OpenFUSIONToolkit python package, honoring OFT_ROOTPATH.

    Returns:
        Dictionary of the TokaMaker symbols used by this script.
    """
    root_path = os.getenv('OFT_ROOTPATH')
    if root_path is not None:
        sys.path.append(os.path.join(root_path, 'python'))
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain, load_gs_mesh
    from OpenFUSIONToolkit.TokaMaker.util import (create_isoflux,
                                                  create_isoflux_xpts,
                                                  create_power_flux_fun,
                                                  create_spline_flux_fun,
                                                  read_eqdsk,
                                                  xpoints_from_moments)
    return {
        'OFT_env': OFT_env,
        'TokaMaker': TokaMaker,
        'gs_Domain': gs_Domain,
        'load_gs_mesh': load_gs_mesh,
        'create_isoflux': create_isoflux,
        'create_isoflux_xpts': create_isoflux_xpts,
        'create_power_flux_fun': create_power_flux_fun,
        'create_spline_flux_fun': create_spline_flux_fun,
        'read_eqdsk': read_eqdsk,
        'xpoints_from_moments': xpoints_from_moments,
    }


def resolve_path(path, base_dir):
    """Resolve a path from the case file relative to the case file's directory.

    Args:
        path: Path as given in the case file.
        base_dir: Directory containing the case file.

    Returns:
        Absolute path.
    """
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(base_dir, path))


def build_lcfs(lcfs_spec, base_dir, api):
    """Build the LCFS contour and X-point list from the case "lcfs" block.

    Args:
        lcfs_spec: "lcfs" block of the case file.
        base_dir: Directory containing the case file.
        api: Dictionary of TokaMaker symbols from _load_tokamaker().

    Returns:
        Tuple of (contour [npts,2], xpoints [nx,2] or None, moments dictionary).
    """
    lcfs_type = lcfs_spec.get('type', 'shape')
    if lcfs_type == 'points':
        contour = numpy.asarray(lcfs_spec['points'], dtype=numpy.float64)
        xpoints = lcfs_spec.get('xpoints')
        if xpoints is not None:
            xpoints = numpy.asarray(xpoints, dtype=numpy.float64)
        return contour, xpoints, moments_from_contour(contour)
    if lcfs_type == 'eqdsk':
        eq_in = api['read_eqdsk'](resolve_path(lcfs_spec['filename'], base_dir))
        contour = numpy.asarray(eq_in['rzout'], dtype=numpy.float64)
        xpoints = lcfs_spec.get('xpoints')
        if xpoints is not None:
            xpoints = numpy.asarray(xpoints, dtype=numpy.float64)
        return contour, xpoints, moments_from_contour(contour)
    if lcfs_type != 'shape':
        raise ValueError('Unknown LCFS type "{0}"'.format(lcfs_type))
    #
    r0 = float(lcfs_spec['R0'])
    z0 = float(lcfs_spec.get('Z0', 0.0))
    a = float(lcfs_spec['a'])
    kappa_u = float(lcfs_spec['kappa_upper'] if 'kappa_upper' in lcfs_spec else lcfs_spec['kappa'])
    kappa_l = float(lcfs_spec['kappa_lower'] if 'kappa_lower' in lcfs_spec else lcfs_spec['kappa'])
    delta_u = float(lcfs_spec['delta_upper'] if 'delta_upper' in lcfs_spec else lcfs_spec['delta'])
    delta_l = float(lcfs_spec['delta_lower'] if 'delta_lower' in lcfs_spec else lcfs_spec['delta'])
    npts = int(lcfs_spec.get('npts', 80))
    moments = {'R0': r0, 'Z0': z0, 'a': a, 'kappa_upper': kappa_u,
               'kappa_lower': kappa_l, 'delta_upper': delta_u, 'delta_lower': delta_l}
    #
    xpt_spec = lcfs_spec.get('xpoints', 'none')
    if isinstance(xpt_spec, (list, tuple)):
        # Explicit X-points: use the smooth Miller contour and the given saddles
        contour = api['create_isoflux'](npts, r0, z0, a, kappa_u, delta_u, kappa_l, delta_l)
        return contour, numpy.asarray(xpt_spec, dtype=numpy.float64), moments
    if xpt_spec in (None, False, 'none'):
        contour = api['create_isoflux'](npts, r0, z0, a, kappa_u, delta_u, kappa_l, delta_l)
        return contour, None, moments
    if xpt_spec not in ('auto', 'both', 'upper', 'lower'):
        raise ValueError('Unknown "xpoints" value "{0}"'.format(xpt_spec))
    # Diverted shape: X-pointed contour from the shaping moments
    squareness = lcfs_spec.get('squareness', {})
    contour = api['create_isoflux_xpts'](
        npts, r0, z0, a, kappa_u, delta_u, kappa_l, delta_l,
        r_inner_mid=lcfs_spec.get('r_inner_mid'),
        zeta_outer_upper=squareness.get('outer_upper', 0.0),
        zeta_outer_lower=squareness.get('outer_lower', 0.0),
        zeta_inner_upper=squareness.get('inner_upper', 0.0),
        zeta_inner_lower=squareness.get('inner_lower', 0.0))
    all_xpts = api['xpoints_from_moments'](r0, z0, a, kappa_u, delta_u, kappa_l, delta_l)
    if xpt_spec == 'upper':
        xpoints = all_xpts[0:1, :]
    elif xpt_spec == 'lower':
        xpoints = all_xpts[1:2, :]
    else:
        xpoints = all_xpts
    return contour, xpoints, moments


def moments_from_contour(contour):
    """Estimate Miller-like shaping moments from an explicit LCFS contour.

    Used only to provide sensible defaults for init_psi().

    Args:
        contour: LCFS point list [npts,2].

    Returns:
        Dictionary of approximate shaping moments.
    """
    r_min, r_max = contour[:, 0].min(), contour[:, 0].max()
    z_min, z_max = contour[:, 1].min(), contour[:, 1].max()
    r0 = 0.5*(r_max + r_min)
    a = 0.5*(r_max - r_min)
    z0 = 0.5*(z_max + z_min)
    r_zmax = contour[contour[:, 1].argmax(), 0]
    r_zmin = contour[contour[:, 1].argmin(), 0]
    return {
        'R0': r0, 'Z0': z0, 'a': a,
        'kappa_upper': (z_max - z0)/a, 'kappa_lower': (z0 - z_min)/a,
        'delta_upper': (r0 - r_zmax)/a, 'delta_lower': (r0 - r_zmin)/a,
    }


def build_profile(prof_spec, api):
    """Build a flux-function profile dictionary from a case-file profile block.

    Args:
        prof_spec: Profile block (None or {'type': 'default'} to leave unset).
        api: Dictionary of TokaMaker symbols from _load_tokamaker().

    Returns:
        Profile dictionary for set_profiles(), or None.
    """
    if prof_spec is None:
        return None
    prof_type = prof_spec.get('type', 'power')
    if prof_type == 'default':
        return None
    if prof_type == 'power':
        return api['create_power_flux_fun'](int(prof_spec.get('npts', 40)),
                                            float(prof_spec['alpha']),
                                            float(prof_spec['gamma']))
    if prof_type == 'linterp':
        return {'type': 'linterp',
                'x': numpy.asarray(prof_spec['x'], dtype=numpy.float64),
                'y': numpy.asarray(prof_spec['y'], dtype=numpy.float64)}
    if prof_type == 'spline':
        return api['create_spline_flux_fun'](int(prof_spec.get('npts', 40)),
                                             numpy.asarray(prof_spec['x'], dtype=numpy.float64),
                                             numpy.asarray(prof_spec['y'], dtype=numpy.float64))
    raise ValueError('Unknown profile type "{0}"'.format(prof_type))


def setup_mesh(mygs, case, base_dir, lcfs_contour, api):
    """Load an existing mesh or generate a fixed-boundary mesh from the LCFS.

    Args:
        mygs: TokaMaker instance.
        case: Full case dictionary.
        base_dir: Directory containing the case file.
        lcfs_contour: LCFS point list [npts,2].
        api: Dictionary of TokaMaker symbols from _load_tokamaker().

    Returns:
        True if the mesh defines coils/conductors (free-boundary capable).
    """
    mesh_file = case.get('mesh_file')
    if mesh_file is not None:
        mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = \
            api['load_gs_mesh'](resolve_path(mesh_file, base_dir))
        mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
        mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
        return True
    if case.get('mode', 'free') == 'free':
        raise ValueError('"mesh_file" is required for free-boundary cases')
    # Fixed boundary: mesh the plasma region bounded by the LCFS
    mesh_dx = case.get('mesh_dx')
    if mesh_dx is None:
        # Default to ~40 cells across the minor diameter
        r_extent = lcfs_contour[:, 0].max() - lcfs_contour[:, 0].min()
        mesh_dx = r_extent/40.0
    gs_mesh = api['gs_Domain']()
    gs_mesh.define_region('plasma', float(mesh_dx), 'plasma')
    gs_mesh.add_polygon(lcfs_contour, 'plasma')
    mesh_pts, mesh_lc, _ = gs_mesh.build_mesh()
    mygs.setup_mesh(mesh_pts, mesh_lc)
    return False


def apply_coil_settings(mygs, coil_spec):
    """Apply VSC definition, current bounds, and current regularization.

    Args:
        mygs: TokaMaker instance.
        coil_spec: "coils" block of the case file.
    """
    vsc = coil_spec.get('vsc')
    if vsc:
        mygs.set_coil_vsc(vsc)
    #
    bounds_spec = coil_spec.get('bounds')
    if bounds_spec:
        default_bound = bounds_spec.get('default')
        coil_bounds = {}
        if default_bound is not None:
            coil_bounds = {name: list(default_bound) for name in mygs.coil_sets}
        for name, bound in bounds_spec.items():
            if name == 'default':
                continue
            if name not in mygs.coil_sets:
                raise KeyError('Unknown coil "{0}" in "bounds"'.format(name))
            coil_bounds[name] = list(bound)
        mygs.set_coil_bounds(coil_bounds)
    #
    reg_spec = coil_spec.get('regularization', {})
    default_weight = float(reg_spec.get('default_weight', 1.0E-2))
    weights = reg_spec.get('weights', {})
    targets = reg_spec.get('targets', {})
    for name in list(weights) + list(targets):
        if name not in mygs.coil_sets:
            raise KeyError('Unknown coil "{0}" in "regularization"'.format(name))
    reg_terms = []
    for name in mygs.coil_sets:
        reg_terms.append(mygs.coil_reg_term({name: 1.0},
                                            target=float(targets.get(name, 0.0)),
                                            weight=float(weights.get(name, default_weight))))
    if vsc:
        # Constrain the virtual VSC coil; a small weight allows vertical adjustment,
        # a large weight effectively disables it in favor of the real VS coil set
        reg_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0,
                                            weight=float(reg_spec.get('vsc_weight', 1.0E-2))))
    if reg_terms:
        mygs.set_coil_reg(reg_terms=reg_terms)


def init_psi_args(case, moments):
    """Determine init_psi() arguments from the case file and LCFS moments.

    Args:
        case: Full case dictionary.
        moments: Shaping moments of the LCFS.

    Returns:
        Argument dictionary for init_psi() (empty for a uniform-current start).
    """
    init_spec = case.get('init', {})
    if init_spec.get('uniform', False):
        return {}
    kappa_avg = 0.5*(moments['kappa_upper'] + moments['kappa_lower'])
    return {
        'r0': float(init_spec.get('R0', moments['R0'])),
        'z0': float(init_spec.get('Z0', moments['Z0'])),
        # A slightly smaller, rounder seed converges more reliably than the target shape
        'a': float(init_spec.get('a', 0.8*moments['a'])),
        'kappa': float(init_spec.get('kappa', 1.0 + 0.6*(kappa_avg - 1.0))),
        'delta': float(init_spec.get('delta', 0.0)),
    }


def write_plot(mygs, filename, has_regions, isoflux_pts, lcfs_contour):
    """Write a flux-surface plot of the converged equilibrium.

    Args:
        mygs: TokaMaker instance.
        filename: Output image path.
        has_regions: Does the mesh define coils/conductors?
        isoflux_pts: Isoflux constraint points to overlay (or None).
        lcfs_contour: Target LCFS contour to overlay when no constraints were used.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(6, 8), constrained_layout=True)
    if has_regions:
        mygs.plot_machine(fig, ax, coil_colormap='seismic', coil_symmap=True,
                          coil_scale=1.0E-6, coil_clabel=r'$I_C$ [MA]')
    mygs.plot_psi(fig, ax, xpoint_color='k', vacuum_nlevels=4)
    if isoflux_pts is not None:
        mygs.plot_constraints(fig, ax, isoflux_color='tab:red', isoflux_marker='.')
    else:
        # Fixed boundary: the target LCFS is the mesh boundary
        closed = numpy.vstack((lcfs_contour, lcfs_contour[0:1, :]))
        ax.plot(closed[:, 0], closed[:, 1], color='tab:red', lw=1.5)
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    fig.savefig(filename, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Compute a TokaMaker G-S equilibrium from a JSON case file')
    parser.add_argument('case_file', help='JSON case definition')
    parser.add_argument('--output-dir', default=None, help='Directory for output files (default: case file directory)')
    parser.add_argument('--nthreads', type=int, default=None, help='OpenMP thread count (default: case file or 2)')
    parser.add_argument('--no-plot', action='store_true', help='Skip writing the flux-surface plot')
    parser.add_argument('--quiet', action='store_true', help='Suppress solver iteration output')
    args = parser.parse_args()

    with open(args.case_file, 'r') as fid:
        case = json.load(fid)
    base_dir = os.path.dirname(os.path.abspath(args.case_file))
    out_spec = case.get('output', {})
    out_dir = args.output_dir or resolve_path(out_spec.get('dir', '.'), base_dir)
    os.makedirs(out_dir, exist_ok=True)
    prefix = out_spec.get('prefix', os.path.splitext(os.path.basename(args.case_file))[0])

    api = _load_tokamaker()
    mode = case.get('mode', 'free')
    if mode not in ('free', 'fixed'):
        raise ValueError('"mode" must be "free" or "fixed"')

    # Boundary definition
    lcfs_contour, xpoints, moments = build_lcfs(case['lcfs'], base_dir, api)

    # Solver setup
    nthreads = args.nthreads or int(case.get('nthreads', 2))
    myOFT = api['OFT_env'](nthreads=nthreads)
    mygs = api['TokaMaker'](myOFT)
    has_regions = setup_mesh(mygs, case, base_dir, lcfs_contour, api)

    for key, value in case.get('settings', {}).items():
        setattr(mygs.settings, key, value)
    # `mode` and the CLI flags win over anything in the "settings" block
    mygs.settings.free_boundary = (mode == 'free')
    if args.quiet:
        mygs.settings.pm = False
    mygs.update_settings()

    if 'F0' in case:
        F0 = float(case['F0'])
    elif ('B0' in case) and ('R0_tor' in case):
        F0 = float(case['B0'])*float(case['R0_tor'])
    else:
        raise ValueError('Either "F0" or both "B0" and "R0_tor" must be given')
    mygs.setup(order=int(case.get('order', 2)), F0=F0)

    # Coils (free boundary only)
    if mode == 'free':
        apply_coil_settings(mygs, case.get('coils', {}))

    # Global targets
    targets = dict(case.get('targets', {}))
    if not targets:
        raise ValueError('At least one global target (eg. "Ip") is required')
    mygs.set_targets(**targets)

    # Shape constraints (free boundary only; for fixed boundary the mesh *is* the shape)
    isoflux_pts = None
    if mode == 'free':
        isoflux_pts = lcfs_contour
        if xpoints is not None:
            isoflux_pts = numpy.vstack((lcfs_contour, xpoints))
        mygs.set_isoflux_constraints(isoflux_pts)
        if xpoints is not None:
            mygs.set_saddle_constraints(xpoints)

    # Flux functions
    profiles = case.get('profiles', {})
    ffp_prof = build_profile(profiles.get('ffp'), api)
    pp_prof = build_profile(profiles.get('pp'), api)
    if (ffp_prof is not None) or (pp_prof is not None):
        mygs.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof)

    # Initialize and solve
    mygs.init_psi(**init_psi_args(case, moments))
    try:
        _, nl_its = mygs.solve(return_its=True)
    except ValueError as solve_error:
        print('\nERROR: G-S solve failed: {0}'.format(solve_error), file=sys.stderr)
        print('       See "Troubleshooting" in the skill reference '
              '(start with "settings": {"urf": 0.3, "maxits": 100} and a softer "init" shape).',
              file=sys.stderr)
        return 1
    if nl_its < 0:
        print('\nERROR: G-S solve did not converge after {0} iterations'.format(-nl_its), file=sys.stderr)
        return 1
    print('\nConverged in {0} nonlinear iterations'.format(nl_its))

    # Outputs
    mygs.print_info()
    stats = mygs.get_stats()
    coil_currents = {}
    if mode == 'free':
        currents, _ = mygs.get_coil_currents()
        coil_currents = {key: float(value) for key, value in currents.items()}
        print('\nCoil currents [MA]:')
        for key, current in coil_currents.items():
            print('  {0:12} {1:10.3F}'.format(key+':', current/1.0E6))

    written = []
    if out_spec.get('stats', True):
        stats_file = os.path.join(out_dir, prefix+'_stats.json')
        payload = {key: (value.tolist() if isinstance(value, numpy.ndarray) else float(value))
                   for key, value in stats.items()}
        payload['nl_iterations'] = int(nl_its)
        if coil_currents:
            payload['coil_currents'] = coil_currents
        with open(stats_file, 'w') as fid:
            json.dump(payload, fid, indent=2)
        written.append(stats_file)
    if out_spec.get('eqdsk', True):
        eqdsk_file = os.path.join(out_dir, 'g'+prefix)
        mygs.save_eqdsk(eqdsk_file, run_info=out_spec.get('run_info', 'TokaMaker')[:40])
        written.append(eqdsk_file)
    if out_spec.get('plot', True) and (not args.no_plot):
        plot_file = os.path.join(out_dir, prefix+'_psi.png')
        write_plot(mygs, plot_file, has_regions, isoflux_pts, lcfs_contour)
        written.append(plot_file)

    print('\nWrote:')
    for path in written:
        print('  '+path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
