import numpy
import h5py
from ..util import write_native_mesh


def write_ThinCurr_mesh(filename, r, lc, reg, holes=[], closures=[], pmap=None, nfp=None):
    r'''Create a native HDF5 mesh file for ThinCurr from the given mesh information

    @param filename Filename for mesh file
    @param r Points list [np,3]
    @param lc Cell list [nc,3] (1-based)
    @param reg Region list [nc]
    @param holes List of node sets for hole definitions
    @param closures List of closures
    @param pmap Point mapping for periodic meshes (single surface only)
    @param nfp Number of field periods for periodic meshes (single surface only)
    '''
    write_native_mesh(filename, r, lc, reg, nodesets=holes, sidesets=[closures,])
    with h5py.File(filename, 'r+') as h5_file:
        if pmap is not None:
            h5_file.create_dataset('thincurr/periodicity/pmap', data=pmap, dtype='i4')
        if nfp is not None:
            h5_file.create_dataset('thincurr/periodicity/nfp', data=[nfp,], dtype='i4')


def write_periodic_mesh(filename, r, lc, reg, tnodeset, pnodesets, pmap=None, nfp=1, include_closures=True):
    if include_closures:
        closures = numpy.arange(nfp)*int(lc.shape[0]/nfp)+1
    else:
        closures = []
    if nfp == 1:
        write_ThinCurr_mesh(filename, r, lc, reg,
           holes=[tnodeset+1, pnodesets[0]+1], closures=closures)
    else:
        write_ThinCurr_mesh(filename, r, lc, reg,
           holes=[tnodeset+1] + [pnodeset+1 for pnodeset in pnodesets], closures=closures, pmap=pmap, nfp=nfp)


def build_regcoil_grid(filename, field_suffix, ntheta, nphi, full_torus=False):
    r'''Build a uniform grid from a REGCOIL definition

    @param field_suffix Suffix for netCDF fields (eg. "plasma" or "coil")
    @param ntheta Number of points in the \f$ \theta \f$ (poloidal) direction
    @param nphi Number of points in the \f$ \phi \f$ (toroidal) direction
    @param full_torus Construct grid for the full torus (default: one field period)
    @result 
    '''
    try:
        import netCDF4
    except ImportError:
        print('"netCDF4" required to load REGCOIL files')
        raise
    with netCDF4.Dataset(filename) as file:
        rmnc = file['rmnc_{0}'.format(field_suffix)][:]
        zmns = file['zmns_{0}'.format(field_suffix)][:],
        xm = file['xm_{0}'.format(field_suffix)][:]
        xn = file['xn_{0}'.format(field_suffix)][:]
        if full_torus:
            nfp = 1
        else:
            nfp = file['nfp'][0]
    #
    theta_span = numpy.linspace(0.0,2.0*numpy.pi,ntheta,endpoint=False)
    zeta_span = numpy.linspace(0.0,2.0*numpy.pi/nfp,nphi,endpoint=(nfp>1))
    rgrid = numpy.zeros((zeta_span.shape[0],theta_span.shape[0],3))
    for j in range(zeta_span.shape[0]):
        for i in range(theta_span.shape[0]):
            r = numpy.sum(rmnc*numpy.cos(xm*theta_span[i] - xn*zeta_span[j]))
            z = numpy.sum(zmns*numpy.sin(xm*theta_span[i] - xn*zeta_span[j]))
            rgrid[j,i,:] = [r*numpy.cos(zeta_span[j]), r*numpy.sin(zeta_span[j]), z]
    return rgrid, nfp


def build_torus_bnorm_grid(filename,nsample,nphi,resample_type='theta',use_spline=False):
    r'''Build a uniform grid from a toroidal surface B-norm definition file

    @param filename Filename of B-norm poloidal mode definition
    @param nsample Number of points in the \f$ \theta \f$ (poloidal) direction
    @param nphi Number of points in the \f$ \phi \f$ (toroidal) direction
    @param resample_type Construct grid for the full torus (default: one field period)
    @param use_spline 
    @result 
    '''
    print('Loading toroidal plasma mode')
    print('  filename = {0}'.format(filename))
    # Read DCON mode from modified "surf.out" file produced by match
    data_lines = []
    with open(filename,'r') as fid:
        for line in fid:
            if not line.strip().startswith('#'):
                data_lines.append(line)
    (npts, nmode) = [int(val) for val in data_lines[0].split()]
    mode_in = numpy.zeros((npts+2,7))
    r0 = numpy.zeros((2,))
    for i in range(npts):
        vals = [float(val) for val in data_lines[i+1].split()]
        mode_in[i+1,:3]=vals[:3]
        mode_in[i+1,6]=vals[3]
    r0 = mode_in[:,:2].sum(axis=0)
    if (abs(mode_in[npts,:]-mode_in[1,:])).max()<1.E-6:
        npts = npts-1 # First and last point are the same
        r0 = r0-mode_in[npts,:2] # Remove one copy from center sum
    r0 /= npts
    #
    print('  N        = {0}'.format(nmode))
    print('  # of pts = {0}'.format(npts))
    print('  R0       = ({0:.4E}, {1:.4E})'.format(*r0))
    # Sort input grid to consistent ordering
    theta_tmp = numpy.zeros((npts,))
    rhat = numpy.zeros((2,))
    for i in range(1,npts+1):
        theta_tmp[i-1]=numpy.arctan2(mode_in[i,1]-r0[1],mode_in[i,0]-r0[0])
        if theta_tmp[i-1] < 0.0:
            theta_tmp[i-1]=2.0*numpy.pi+theta_tmp[i-1]
        i_n=i+1
        if i == npts:
            i_n=0
        rhat[0]=-(mode_in[i_n,1]-mode_in[i,1])
        rhat[1]=mode_in[i_n,0]-mode_in[i,0]
        rhat /= numpy.linalg.norm(rhat)
        mode_in[i,4:6]=rhat
    if theta_tmp[0] > numpy.pi:
        theta_tmp[0] = 2.0*numpy.pi-theta_tmp[0]
    ind_reorder = theta_tmp.argsort()
    mode_tmp = mode_in.copy()
    for i in range(1,npts+1):
      mode_in[i,:] = mode_tmp[ind_reorder[i-1]+1,:]
      mode_in[i,3] = theta_tmp[ind_reorder[i-1]]
    mode_in[0,:]=mode_in[npts,:]; mode_in[0,3]=mode_in[npts,3]-2.0*numpy.pi
    mode_in[npts+1,:]=mode_in[1,:]; mode_in[npts+1,3]=mode_in[1,3]+2.0*numpy.pi
    #
    mode_resample = numpy.zeros((nsample,6))
    resample_grid = numpy.zeros((nsample,))
    # Setup grid based on sampling type
    if resample_type == "theta":
        for i in range(nsample):
            resample_grid[i]=i*2.0*numpy.pi/nsample
    elif resample_type == "arc_len":
        mode_in[0,3]=0.0
        for i in range(npts+1):
            mode_in[i+1,3]=mode_in[i,3]+numpy.linalg.norm(mode_in[i+1,:2]-mode_in[i,:2])
        for i in range(nsample):
            resample_grid[i]=i*mode_in[npts,3]/nsample
    else:
        print("Invalid resampling type")
        return None
        #raise ValueError("Invalid resampling type")
    if use_spline:
        try:
            from scipy.interpolate import CubicSpline
        except ImportError:
            print("SciPy required for spline interpolation")
            raise
        splines = [
            CubicSpline(mode_in[:npts+2,3], mode_in[:npts+2,0], axis=0, bc_type='not-a-knot'),
            CubicSpline(mode_in[:npts+2,3], mode_in[:npts+2,1], axis=0, bc_type='not-a-knot'),
            CubicSpline(mode_in[:npts+2,3], mode_in[:npts+2,2], axis=0, bc_type='not-a-knot'),
            CubicSpline(mode_in[:npts+2,3], mode_in[:npts+2,4], axis=0, bc_type='not-a-knot'),
            CubicSpline(mode_in[:npts+2,3], mode_in[:npts+2,5], axis=0, bc_type='not-a-knot'),
            CubicSpline(mode_in[:npts+2,3], mode_in[:npts+2,6], axis=0, bc_type='not-a-knot')
        ]
        for i in range(nsample):
            for j, spline in enumerate(splines):
                mode_resample[i,j]=spline(resample_grid[i])
    else:
        for j in range(3):
            mode_resample[:,j]=numpy.interp(resample_grid,mode_in[:npts+2,3],mode_in[:npts+2,j])
            mode_resample[:,3+j]=numpy.interp(resample_grid,mode_in[:npts+2,3],mode_in[:npts+2,4+j])
    print('  Mode pair sums {0:.4E} {1:.4E}'.format(mode_resample[:,2].sum(),mode_resample[:,5].sum()))
    #
    phi_grid = numpy.linspace(0.0,2.0*numpy.pi/nmode,nphi,endpoint=(nmode>1))
    r = numpy.zeros((nphi,nsample,3))
    bnorm = numpy.zeros((2,nphi,nsample))
    for i in range(nsample):
      i_n=i+1
      if i == nsample-1:
          i_n=0
      for j in range(nphi):
          r[j,i,:]=[
              mode_resample[i,0]*numpy.cos(phi_grid[j]),
              mode_resample[i,0]*numpy.sin(phi_grid[j]),
              mode_resample[i,1]
          ]
          bnorm[0,j,i]=mode_resample[i,2]*numpy.sin(nmode*phi_grid[j]) \
              + mode_resample[i,5]*numpy.cos(nmode*phi_grid[j])
          bnorm[1,j,i]=mode_resample[i,2]*numpy.cos(nmode*phi_grid[j]) \
              - mode_resample[i,5]*numpy.sin(nmode*phi_grid[j])
    return r, bnorm, nmode


def build_periodic_mesh(r_grid,nfp):
    r'''Build triangular mesh for the full surface from a uniform grid of a single field period (toroidal)

    @param r_grid Uniform grid [nphi,ntheta,3] (\f$ \phi \f$ and \f$ \theta \f$ vary along the first and second dimension respectively)
    @param nfp Number of field periods
    @result point list [np,3], cell list [nc,3], toroidal nodeset, poloidal nodesets, periodicity map
    '''
    if nfp == 1:
        lc, r, tnodeset, pnodeset = build_triangles_from_grid(r_grid)
        return r, lc, tnodeset, [pnodeset,], None
    nphi = r_grid.shape[0]
    ntheta = r_grid.shape[1]
    r_base = r_grid.copy().reshape((nphi*ntheta,3))
    r_full = r_base[:-ntheta,:].copy()
    r_map = numpy.arange((nphi-1)*ntheta, dtype=numpy.int32)
    rotation = 2.0*numpy.pi/nfp
    for i in range(1,nfp):
        r_rotated = r_base[:-ntheta,:].copy()
        for j in range(r_rotated.shape[0]):
            theta = numpy.arctan2(r_rotated[j,1],r_rotated[j,0])
            if theta < 0.0:
                theta += 2*numpy.pi
            r = numpy.sqrt(numpy.power(r_rotated[j,0],2)+numpy.power(r_rotated[j,1],2))
            r_rotated[j,0] = r*numpy.cos(theta + rotation*i)
            r_rotated[j,1] = r*numpy.sin(theta + rotation*i)
        r_full = numpy.vstack((r_full,r_rotated))
        r_map = numpy.hstack((r_map,numpy.arange((nphi-1)*ntheta, dtype=numpy.int32)))
    lc, r, tnodeset, pnodeset = build_triangles_from_grid(r_full.reshape((int(r_full.shape[0]/(ntheta)),ntheta,3)))
    pnodesets = [pnodeset]
    for i in range(1,nfp):
        pnodesets.append(pnodeset+i*(nphi-1)*ntheta)
    return r, lc, tnodeset, pnodesets, r_map


def build_triangles_from_grid(data_grid,wrap_n=True,wrap_m=True):
    r'''Build triangles from a uniform grid of points

    @param data_grid Uniform point grid [n,m,3]
    @param wrap_n Wrap grid in n-direction?
    @param wrap_m Wrap grid in m-direction?
    @result point list [n*m,3], cell list [nc,3], nodeset in the n-direction, nodeset in the m-direction
    '''
    n = data_grid.shape[0]
    m = data_grid.shape[1]
    triangles = []
    nodeset_n = []
    nodeset_m = []
    #
    if wrap_n:
        nwrap = n
    else:
        nwrap = n-1
    if wrap_m:
        mwrap = m
    else:
        mwrap = m-1
    for i in range(nwrap):
        i1 = (i + 1) % n 
        nodeset_n.append(i*m)
        for j in range(mwrap):
            # Wrap around the indices when at the borders
            j1 = (j + 1) % m  
            if i == 0:
                nodeset_m.append(j)
            # Get the indices of the vertices of the two triangles
            # that divide the current grid square
            tri1 = [i * m + j, i * m + j1, i1 * m + j]
            tri2 = [i1 * m + j1, i1 * m + j, i * m + j1]
    
            # Add the two triangles to the list
            triangles += [tri1, tri2]
    if not wrap_n:
        nodeset_n.append((n-1)*m)
    if not wrap_m:
        nodeset_m.append(m-1)
    nodeset_n = numpy.asarray(nodeset_n)
    nodeset_m = numpy.asarray(nodeset_m)
    return data_grid.reshape((n*m,3)), numpy.asarray(triangles), nodeset_n, nodeset_m