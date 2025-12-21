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