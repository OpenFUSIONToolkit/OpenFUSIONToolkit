'''! Python interface for TokaMaker Grad-Shafranov functionality

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''
import numpy
from ..util import *

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
    try:
        from scipy.interpolate import CubicSpline
    except ImportError:
        print("Spline flux function requires SciPy")
        raise
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
    r'''! Build power law flux function of the form \f$ ((1-\hat{\psi})^{\alpha})^{\gamma} \f$

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