#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! General utility and supporting functions for TokaMaker

@authors Chris Hansen
@date April 2024
@ingroup doxy_oft_python
'''
import struct
import numpy
from collections import OrderedDict
from .._interface import *
from ..util import read_fortran_namelist

## @cond
tokamaker_eval_green = ctypes_subroutine(oftpy_lib.tokamaker_eval_green,
    [c_int, ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,1), c_double, c_double, ctypes_numpy_array(float64,1)])
## @endcond


def create_isoflux(npts, r0, z0, a, kappa, delta, kappaL=None, deltaL=None):
    r'''! Create isoflux points using simple analytic form

    @param npts Number of points to sample (evenly spaced in \f$\theta\f$)
    @param r0 Major radial position for magnetic axis
    @param z0 Vertical position for magnetic axis
    @param a Minor radius
    @param kappa Elongation (upper only if kappaL is set)
    @param delta Triangularity (upper only if deltaL is set)
    @param kappaL Lower elongation (default: kappa)
    @param deltaL Lower triangularity (default: delta)
    @result Point list [npts,2]
    '''
    kappaU = kappa
    deltaU = delta
    if kappaL is None:
        kappaL = kappaU
    if deltaL is None:
        deltaL = deltaU
    x0 = numpy.r_[r0, z0]
    isoflux_pts = numpy.zeros((npts,2))
    for i in range(npts):
        theta = i*2.0*numpy.pi/npts
        delta = ((deltaU + deltaL) + (deltaU - deltaL)*numpy.sin(theta))/2
        kappa = ((kappaU + kappaL) + (kappaU - kappaL)*numpy.sin(theta))/2
        isoflux_pts[i,:] = x0 + a*numpy.r_[numpy.cos(theta+numpy.arcsin(delta)*numpy.sin(theta)),kappa*numpy.sin(theta)]
    return isoflux_pts


def xpoints_from_moments(r0, z0, a, kappa_upper, delta_upper,
                         kappa_lower=None, delta_lower=None):
    r'''! Compute X-point locations from Miller shaping moments

    Upper: \f$(r_0 - a\,\delta_U,\; z_0 + a\,\kappa_U)\f$.
    Lower: \f$(r_0 - a\,\delta_L,\; z_0 - a\,\kappa_L)\f$.

    @param r0 Major radial position for magnetic axis
    @param z0 Vertical position for magnetic axis
    @param a Minor radius
    @param kappa_upper Upper elongation
    @param delta_upper Upper triangularity
    @param kappa_lower Lower elongation (default: kappa_upper)
    @param delta_lower Lower triangularity (default: delta_upper)
    @result X-point array, shape (2, 2): [[R_upper, Z_upper], [R_lower, Z_lower]]
    '''
    if kappa_lower is None:
        kappa_lower = kappa_upper
    if delta_lower is None:
        delta_lower = delta_upper
    return numpy.array([
        [r0 - a * delta_upper, z0 + a * kappa_upper],
        [r0 - a * delta_lower, z0 - a * kappa_lower],
    ])


def _intor_arc(r0, z0_arc, a, kappa_eff, theta_x, r_x, zeta, n, endpoint=False, sign=1):
    r'''! Generate one half of the INTOR outer arc

    Arc spans \f$\theta\f$ from 0 (outboard midplane) to
    \f$\mathrm{sign}\cdot\theta_x\f$ (X-point).
    sign=+1 gives the upper half (Z > z0_arc); sign=-1 gives the lower half.
    delta_eff is solved so the arc passes exactly through r_x at theta_x.

    @param r0 Major radial position for magnetic axis
    @param z0_arc Vertical center for this arc half (typically the magnetic axis Z)
    @param a Minor radius
    @param kappa_eff Elongation for this half
    @param theta_x Poloidal angle at the X-point
    @param r_x R coordinate of the X-point
    @param zeta INTOR squareness parameter
    @param n Number of points to generate
    @param endpoint Whether to include the final endpoint (default: False)
    @param sign +1 for upper half, -1 for lower half (default: +1)
    @result r, z arrays for the arc segment
    '''
    arccos_val = numpy.arccos(numpy.clip((r_x - r0) / a, -1.0, 1.0))
    delta_eff = (arccos_val - theta_x + zeta * numpy.sin(2 * theta_x)) / max(numpy.sin(theta_x), 1e-8)
    theta = sign * numpy.linspace(0.0, theta_x, n, endpoint=endpoint)
    r = r0 + a * numpy.cos(theta + delta_eff * numpy.sin(theta) - zeta * numpy.sin(2 * theta))
    z = z0_arc + a * kappa_eff * numpy.sin(theta + zeta * numpy.sin(2 * theta))
    return r, z


def _inner_arc(r_x, z_x, r_inner_mid, z_inner_mid, zeta, n):
    r'''! Generate an inverted-INTOR inner arc from X-point to inner midplane

    Uses \f$R(\alpha) = r_x - a_\mathrm{in}\cos(\alpha - \zeta\sin 2\alpha)\f$,
    \f$Z(\alpha) = z_\mathrm{mid} + a_\mathrm{in}\kappa_\mathrm{in}\sin(\alpha + \zeta\sin 2\alpha)\f$
    with \f$\alpha\f$ from \f$\pi/2\f$ (X-point) to 0 (midplane, excluded).
    The tangent is vertical at the midplane, so joining upper and lower halves is kink-free.

    @param r_x R coordinate of the X-point
    @param z_x Z coordinate of the X-point
    @param r_inner_mid R at the inner midplane
    @param z_inner_mid Z at the inner midplane
    @param zeta INTOR squareness parameter
    @param n Number of points (midplane endpoint excluded)
    @result r, z arrays for the arc segment
    '''
    a_in = max(r_x - r_inner_mid, 1e-10)
    kappa_in = (z_x - z_inner_mid) / a_in
    alpha = numpy.linspace(numpy.pi / 2, 0.0, n, endpoint=False)
    r = r_x - a_in * numpy.cos(alpha - zeta * numpy.sin(2 * alpha))
    z = z_inner_mid + a_in * kappa_in * numpy.sin(alpha + zeta * numpy.sin(2 * alpha))
    return r, z


def create_isoflux_xpts(npts, r0, z0, a, kappa_upper, delta_upper,
                        kappa_lower=None, delta_lower=None,
                        r_inner_mid=None,
                        zeta_outer_upper=0.0, zeta_outer_lower=0.0,
                        zeta_inner_upper=0.0, zeta_inner_lower=0.0):
    r'''! Create isoflux boundary points for single-null configurations with X-points

    Generates a closed boundary contour suitable for use with
    \ref OpenFUSIONToolkit.TokaMaker.TokaMaker.set_isoflux_constraints "set_isoflux_constraints()"
    and \ref OpenFUSIONToolkit.TokaMaker.TokaMaker.set_saddle_constraints "set_saddle_constraints()".
    The outer (low-field-side) arc uses the INTOR analytic formula centered on the magnetic
    axis; the inner arcs use an inverted-INTOR formula that gives a smooth vertical tangent
    at the inner midplane.  Supports non-up-down-symmetric equilibria via independent upper
    and lower shaping moments and squareness values.

    Use \ref xpoints_from_moments "xpoints_from_moments()" with the same kappa/delta
    arguments to obtain the matching saddle-constraint array.

    @param npts Total number of boundary points
    @param r0 Major radial position for magnetic axis
    @param z0 Vertical position for magnetic axis
    @param a Minor radius
    @param kappa_upper Upper elongation (sets X-point height and outer arc shape)
    @param delta_upper Upper triangularity
    @param kappa_lower Lower elongation (default: kappa_upper)
    @param delta_lower Lower triangularity (default: delta_upper)
    @param r_inner_mid R at the inner midplane (default: r0 - a)
    @param zeta_outer_upper INTOR squareness of the upper outer corner (default: 0)
    @param zeta_outer_lower INTOR squareness of the lower outer corner (default: 0)
    @param zeta_inner_upper INTOR squareness of the upper inner corner (default: 0)
    @param zeta_inner_lower INTOR squareness of the lower inner corner (default: 0)
    @result Point list [npts, 2]
    '''
    if npts < 3:
        raise ValueError("\"npts\" must be at least 3")
    kU = float(kappa_upper)
    dU = float(delta_upper)
    kL = float(kappa_lower) if kappa_lower is not None else kU
    dL = float(delta_lower) if delta_lower is not None else dU
    zou = float(zeta_outer_upper)
    zol = float(zeta_outer_lower)
    ziu = float(zeta_inner_upper)
    zil = float(zeta_inner_lower)

    r_xu = r0 - a * dU
    z_xu = z0 + a * kU
    r_xl = r0 - a * dL
    z_xl = z0 - a * kL
    zsurf = (z_xu + z_xl) / 2.0

    if r_inner_mid is None:
        r_inner_mid = r0 - a

    # Outer arc: INTOR centered on z0; each half spans θ = π/2 by the moment definition
    n_seed = max(8 * npts, 800)
    n_su = n_seed // 2
    n_sl = n_seed - n_su

    r_seed_u, z_seed_u = _intor_arc(r0, z0, a, kU, numpy.pi / 2, r_xu, zou, n_su, endpoint=True, sign=+1)
    r_seed_l, z_seed_l = _intor_arc(r0, z0, a, kL, numpy.pi / 2, r_xl, zol, n_sl, endpoint=True, sign=-1)

    r_arc = numpy.concatenate([r_seed_l[::-1], r_seed_u[1:]])
    z_arc = numpy.concatenate([z_seed_l[::-1], z_seed_u[1:]])
    r_arc[0],  z_arc[0]  = r_xl, z_xl
    r_arc[-1], z_arc[-1] = r_xu, z_xu

    outer_segs = numpy.hypot(numpy.diff(r_arc), numpy.diff(z_arc))
    outer_len  = float(outer_segs.sum())

    # Estimate inner arc lengths for proportional point allocation
    n_est = 400
    r_ui_e, z_ui_e = _inner_arc(r_xu, z_xu, r_inner_mid, zsurf, ziu, n_est)
    r_li_e, z_li_e = _inner_arc(r_xl, z_xl, r_inner_mid, zsurf, zil, n_est)
    len_ui = float(numpy.sum(numpy.hypot(numpy.diff(r_ui_e), numpy.diff(z_ui_e))))
    len_li = float(numpy.sum(numpy.hypot(numpy.diff(r_li_e), numpy.diff(z_li_e))))
    inner_len = len_ui + len_li

    n_outer   = max(2, round(npts * outer_len / (outer_len + inner_len)))
    n_inner   = npts - n_outer
    n_inner_u = max(1, round(n_inner * len_ui / inner_len))
    n_inner_l = n_inner - n_inner_u

    # Resample outer arc; upper X-point is the last point
    outer_cum = numpy.concatenate([[0.0], numpy.cumsum(outer_segs)])
    t_sample  = numpy.linspace(0.0, outer_cum[-1], n_outer, endpoint=True)
    r_outer   = numpy.interp(t_sample, outer_cum, r_arc)
    z_outer   = numpy.interp(t_sample, outer_cum, z_arc)
    r_outer[-1], z_outer[-1] = r_xu, z_xu

    # Upper inner arc: X-point → midplane (skip X-point, already last in outer arc)
    r_iu, z_iu = _inner_arc(r_xu, z_xu, r_inner_mid, zsurf, ziu, n_inner_u + 1)
    r_iu, z_iu = r_iu[1:], z_iu[1:]

    # Lower inner arc: midplane → X-point (reversed; drop the final X-point,
    # which is already the first point of the outer arc, to avoid a duplicate)
    r_il_raw, z_il_raw = _inner_arc(r_xl, z_xl, r_inner_mid, zsurf, zil, n_inner_l + 1)
    r_il = r_il_raw[::-1][:-1]
    z_il = z_il_raw[::-1][:-1]

    return numpy.column_stack([numpy.concatenate([r_outer, r_iu, r_il]),
                               numpy.concatenate([z_outer, z_iu, z_il])])


def create_spline_flux_fun(npts,x,y,axis_bc=[1,0.0],edge_bc=[1,0.0],normalize=True):
    r'''! Build cubic spline flux function

    @param npts Number of points for definition
    @param x Location of spline "knots" in normalized flux
    @param y Value of flux function at spline "knots"
    @param axis_bc [SciPy BC specification](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html) on axis (\f$ \hat{\psi} = 0 \f$)
    @param edge_bc [SciPy BC specification](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html) on LCFS (\f$ \hat{\psi} = 1 \f$)
    @result Flux function definition dictionary
    '''
    from scipy.interpolate import CubicSpline
    prof = CubicSpline(x,y,bc_type=[axis_bc,edge_bc])
    x_sample = numpy.linspace(0.0,1.0,npts)
    prof = {
        'type': 'linterp',
        'x': x_sample,
        'y': prof(x_sample)
    }
    if normalize:
        prof['y'] /= prof['y'][0]
    return prof



def create_power_flux_fun(npts,alpha,gamma):
    r'''! Build power law flux function of the form \f$ (1-\hat{\psi}^{\alpha})^{\gamma} \f$

    @param npts Number of points for definition
    @param alpha Inner exponent
    @param gamma Outer exponent
    @result Flux function definition dictionary
    '''
    psi_pts = numpy.linspace(0.0,1.0,npts)
    return {
        'type': 'linterp',
        'x': psi_pts,
        'y': numpy.power(1.0-numpy.power(psi_pts,alpha),gamma)
    }


def read_eqdsk(filename):
    '''! Read gEQDSK file

    @param filename Path to gEQDSK file
    @result Dictionary containing gEQDSK information
    '''
    def read_1d(fid, n):
        j = 0
        output = numpy.zeros((n,))
        for i in range(n):
            if j == 0:
                line = fid.readline()
            output[i] = line[j:j+16]
            j += 16
            if j == 16*5:
                j = 0
        return output

    def read_2d(fid, n, m):
        j = 0
        output = numpy.zeros((n, m))
        for k in range(n):
            for i in range(m):
                if j == 0:
                    line = fid.readline()
                output[k, i] = line[j:j+16]
                j += 16
                if j == 16*5:
                    j = 0
        return output
    # Read-in data
    eqdsk_obj = {}
    with open(filename, 'r') as fid:
        # Get sizes
        line = fid.readline()
        eqdsk_obj['case'] = line[:48]
        split_line = line[48:].split()
        eqdsk_obj['nr'] = int(split_line[-2])
        eqdsk_obj['nz'] = int(split_line[-1])
        # Read header content
        line_keys = [['rdim',  'zdim',  'rcentr',  'rleft',  'zmid'],
                     ['raxis', 'zaxis', 'psimag', 'psibry', 'bcentr'],
                     ['ip',    'skip',  'skip',   'skip',   'skip'],
                     ['skip',  'skip',  'skip',   'skip',   'skip']]
        for i in range(4):
            line = fid.readline()
            for j in range(5):
                if line_keys[i][j] == 'skip':
                    continue
                line_seg = line[j*16:(j+1)*16]
                eqdsk_obj[line_keys[i][j]] = float(line_seg)
        # Read flux profiles
        keys = ['fpol', 'pres', 'ffprim', 'pprime']
        for key in keys:
            eqdsk_obj[key] = read_1d(fid, eqdsk_obj['nr'])
        # Read PSI grid
        eqdsk_obj['psirz'] = read_2d(fid, eqdsk_obj['nz'],
                                        eqdsk_obj['nr'])
        # Read q-profile
        eqdsk_obj['qpsi'] = read_1d(fid, eqdsk_obj['nr'])
        # Read limiter count
        line = fid.readline()
        eqdsk_obj['nbbs'] = int(line.split()[0])
        eqdsk_obj['nlim'] = int(line.split()[1])
        # Read outer flux surface
        eqdsk_obj['rzout'] = read_2d(fid, eqdsk_obj['nbbs'], 2)
        # Read limiting corners
        eqdsk_obj['rzlim'] = read_2d(fid, eqdsk_obj['nlim'], 2)
    return eqdsk_obj


def read_ifile(filename):
    '''! Read i-file inverse equilibrium file

    @param filename Path to file
    @result Dictionary containing i-file information
    '''
    def read_array(content,offset,var_type,count):
        if var_type in ("i", "f"):
            var_size = 4
        elif var_type in ("l", "d"):
            var_size = 8
        else:
            raise ValueError("Invalid variable type")
        array_size = var_size*count
        #
        head_val = struct.unpack_from("i",content,offset=offset)
        if head_val[0] != array_size:
            raise ValueError("Dataframe size does not match array size")
        offset += 4
        body_val = struct.unpack_from("="+var_type*count,content,offset=offset)
        offset += array_size
        tail_val = struct.unpack_from("i",content,offset=offset)
        if head_val[0] != tail_val[0]:
            raise ValueError("Head and tail values disagree")
        offset += 4
        return numpy.array(body_val), offset
    #
    with open(filename, 'rb') as handle:
        content = handle.read()
    out_dict = {}
    offset = 0
    sizes, offset = read_array(content,offset,"i",2)
    out_dict["npsi"] = sizes[0]
    out_dict["ntheta"] = sizes[1]
    var_type = "d"
    try:
        out_dict["psi"], offset = read_array(content,offset,var_type,sizes[0])
    except ValueError:
        try:
            var_type = "f"
            out_dict["psi"], offset = read_array(content,offset,var_type,sizes[0])
        except ValueError:
            raise ValueError("Unable to determine float point datatype")
    out_dict["f"], offset = read_array(content,offset,var_type,sizes[0])
    out_dict["p"], offset = read_array(content,offset,var_type,sizes[0])
    out_dict["q"], offset = read_array(content,offset,var_type,sizes[0])
    R, offset = read_array(content,offset,var_type,sizes[0]*sizes[1])
    Z, offset = read_array(content,offset,var_type,sizes[0]*sizes[1])
    out_dict["R"] = R.reshape(sizes)
    out_dict["Z"] = Z.reshape(sizes)
    return out_dict


def eval_green(x,xc):
        r'''! Evaluate Green's function for a toroidal filament

        @param x Observation point [2]
        @param xc Coil location [:,2]
        @result \f$\psi(x)\f$ due to a coil with unit current [A] at xc
        '''
        n = x.shape[0]
        vals = numpy.zeros((n,),dtype=numpy.float64)
        r = x[:,0].copy()
        z = x[:,1].copy()
        tokamaker_eval_green(c_int(n),r,z,
            c_double(xc[0]),c_double(xc[1]),vals)
        return vals*mu0


def compute_forces_components(tMaker_obj,psi,cell_centered=False):
    r'''! Compute terms needed for evaluating forces in passively conducting regions

    @param tMaker_obj TokaMaker equilibrium object
    @param psi \f$ \psi \f$ corresponding to desired currents
    @param cell_centered Evaluate at cell centers instead of node points?
    @result J_cond, B_cond, mask, R
    '''
    # Get conductor currents at cell centers
    mask, J_cond = tMaker_obj.get_conductor_currents(psi,cell_centered=cell_centered)
    
    # Find points inside conducting regions
    pt_mask = numpy.zeros((tMaker_obj.r.shape[0],), dtype=numpy.int32)
    pt_mask[tMaker_obj.lc[mask,:]] = 1
    
    # Set psi and evaluate B-field in conducting regions
    psi_save = tMaker_obj.get_psi(normalized=False)
    tMaker_obj.set_psi(psi)
    field_eval = tMaker_obj.get_field_eval('B')
    B_cond = numpy.zeros((tMaker_obj.r.shape[0],3))
    for i in range(tMaker_obj.r.shape[0]):
        if pt_mask[i] == 0:
            continue
        B_cond[i,:] = field_eval.eval(tMaker_obj.r[i,:2])
    tMaker_obj.set_psi(psi_save) # Reset psi

    if cell_centered:
        # Convert to cell centered
        Bv_cond = (B_cond[tMaker_obj.lc[:,0],:] + B_cond[tMaker_obj.lc[:,1],:] + B_cond[tMaker_obj.lc[:,2],:])/3.0
        rcc = (tMaker_obj.r[tMaker_obj.lc[:,0],:] + tMaker_obj.r[tMaker_obj.lc[:,1],:] + tMaker_obj.r[tMaker_obj.lc[:,2],:])/3.0
        return J_cond, Bv_cond, mask, rcc
    else:
        return J_cond, B_cond, mask, tMaker_obj.r
    
def read_mhdin(path, e_coil_names=None, f_coil_names=None):
    r'''Read mhdin.dat file.

    @param path Path to file
    @param e_coil_names Names of E coils (hardcoded, generates indexed names if None)
    @param f_coil_names Names of F coils (hardcoded, generates indexed names if None)
    @result machine_dict Dictionary containing coil coordinates and turns, loop names, and probe names and angles.
    @result raw Dictionary containing all other data from mhdin.dat
    '''
    raw = read_fortran_namelist(path)
    machine_dict = OrderedDict()
    
    # Expand later
    keys = ['MPNAM2', 'LPNAME']
    for key in keys:
        names = raw[key].replace("'", " ")
        names = names.split()
        machine_dict[key] = names        
        
    e_coil_dict = OrderedDict()

    for i in range(len(raw['ECID'])):
        idx = int(raw['ECID'][i]) - 1
        e_coil_name = "ECOIL{:03d}".format(idx + 1)
        if e_coil_names:
            e_coil_name = e_coil_names[idx]
        if e_coil_name not in e_coil_dict:
            e_coil_dict[e_coil_name] = []
        e_coil_dict[e_coil_name].append([float(raw['RE'][i]), float(raw['ZE'][i]), float(raw['WE'][i]), float(raw['HE'][i])])
    machine_dict['ECOIL'] = e_coil_dict

    f_coil_dict = OrderedDict()
    for i in range(len(raw['FCID'])):
        f_coil_name = "FCOIL{:03d}".format(i + 1)
        if f_coil_names:
            f_coil_name = f_coil_names[i]
        f_coil_dict[f_coil_name] = [float(raw['RF'][i]), float(raw['ZF'][i]), float(raw['WF'][i]), float(raw['HF'][i]), float(raw['TURNFC'][i])]
    machine_dict['FCOIL'] = f_coil_dict

    probe_angle_dict = OrderedDict()
    i = 0
    probe_angles = raw['AMP2']
    for probe_name in machine_dict['MPNAM2']:
        probe_angle_dict[probe_name] = float(probe_angles[i])
        i = i + 1
    machine_dict['AMP2'] = probe_angle_dict

    return machine_dict, raw

def read_kfile(path, machine_dict, e_coil_names=None, f_coil_names=None):
    r'''Read k-file.

    @param path Path to file
    @param e_coil_names Names of E coils (hardcoded, generates indexed names if None)
    @param f_coil_names Names of F coils (hardcoded, generates indexed names if None)
    @param machine_dict Result from read_mhdin (contents of mhdin.dat file)
    @result probes_dict Dictionary containing probe values and weights (0 if not selected).
    @result loops_dict Dictionary containing loop values and weights (0 if not selected).
    @result e_coil_dict Dictionary containing E copil values and weights (0 if not selected).
    @result f_coil_dict Dictionary containing F coil values and weights (0 if not selected).
    @result raw Dictionary containing all other data from k-file.
    '''
    raw = read_fortran_namelist(path)
    
    probe_names = machine_dict['MPNAM2']
    probe_vals = raw['EXPMP2']
    probe_weights = raw['FWTMP2']
    probes_dict = OrderedDict()
    for i in range(len(probe_names)):
        probes_dict[probe_names[i]] = [probe_vals[i], probe_weights[i]]

    loop_names = machine_dict['LPNAME']
    loop_vals = raw['COILS']
    loop_weights = raw['FWTSI']
    loops_dict = OrderedDict()
    for i in range(len(loop_names)):
        loops_dict[loop_names[i]] = [loop_vals[i], loop_weights[i]]
        
    if e_coil_names is None:
        e_coil_names = sorted(machine_dict['ECOIL'].keys())
    e_coil_vals = raw['ECURRT']
    e_coil_weights = raw['FWTEC']
    e_coil_dict = OrderedDict()
    for i in range(len(e_coil_names)):
        e_coil_dict[e_coil_names[i]] = [e_coil_vals[i], e_coil_weights[i]]
    
    if f_coil_names is None:
        f_coil_names = sorted(machine_dict['FCOIL'].keys())
    f_coil_vals = raw['BRSP']
    f_coil_weights = raw['FWTFC']
    f_coil_dict = OrderedDict()

    for i in range(len(f_coil_names)):
        f_coil_dict[f_coil_names[i]] = [f_coil_vals[i], f_coil_weights[i]]
    return probes_dict, loops_dict, e_coil_dict, f_coil_dict, raw

def get_jphi_from_GS(ffprime, pprime, R_avg, one_over_R_avg):
    r'''! Calculate j_phi profile from Grad-Shafranov equation
    @param ffprime FF'(psi_N) profile
    @param pprime P'(psi_N) profile
    @param R_avg <R>(psi_N) profile
    @param one_over_R_avg <1/R>(psi_N) profile
    Returns:
    j_phi(\psi_N) profile
    '''
    return ffprime * (one_over_R_avg / mu0) + R_avg * pprime