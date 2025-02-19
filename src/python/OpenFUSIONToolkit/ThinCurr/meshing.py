import numpy
import h5py

import numpy.linalg
from ..util import write_native_mesh


def write_ThinCurr_mesh(filename, r, lc, reg, holes=[], closures=[], pmap=None, nfp=None):
    r'''! Create a native HDF5 mesh file for ThinCurr from the given mesh information

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


def build_ThinCurr_dummy(center,size=1.0,nsplit=0):
    '''! Build simple square dummy mesh for ThinCurr (1 active node)
    
    @param center Center of mesh [3]
    @param size Physical size of dummy mesh in X and Y
    @param nsplit Number of refinement iterations to perform on grid (starting mesh is `np=5`, `nc=4`)
    @returns `r` Point list [:,3], lc Cell list [:,3]
    '''
    r = numpy.array([
        [-size/2.0,-size/2.0,0.0],
        [size/2.0,-size/2.0,0.0],
        [size/2.0,size/2.0,0.0],
        [-size/2.0,size/2.0,0.0],
        [0.0,0.0,0.0]
    ])
    lc = numpy.array([
        [0,1,4],
        [1,2,4],
        [2,3,4],
        [3,0,4]
    ])
    for i in range(3):
        r[:,i] += center[i]
    # Refine grid as desired
    for i in range(nsplit):
        lc_new = []
        r_new = [r_old for r_old in r]
        for j in range(lc.shape[0]):
            new_inds = [0,0,0]
            r_candidates = [
                (r[lc[j,0],:]+r[lc[j,1],:])/2.0,
                (r[lc[j,1],:]+r[lc[j,2],:])/2.0,
                (r[lc[j,0],:]+r[lc[j,2],:])/2.0
            ]
            for k in range(3):
                for k2 in range(r.shape[0],len(r_new)):
                    if numpy.linalg.norm(r_new[k2]-r_candidates[k]) < 1.E-10:
                        new_inds[k] = k2
                        break
                else:
                    r_new.append(r_candidates[k])
                    new_inds[k] = len(r_new)-1
            lc_new.append([lc[j,0], new_inds[0], new_inds[2]])
            lc_new.append([new_inds[0], lc[j,1], new_inds[1]])
            lc_new.append([new_inds[1], lc[j,2], new_inds[2]])
            lc_new.append([new_inds[0], new_inds[1], new_inds[2]])
        lc = numpy.array(lc_new)
        r = numpy.array(r_new)
    return r, lc


def build_regcoil_grid(filename, field_suffix, ntheta, nphi, full_torus=False):
    r'''! Build a uniform grid from a REGCOIL definition

    @param field_suffix Suffix for netCDF fields (eg. "plasma" or "coil")
    @param ntheta Number of points in the \f$ \theta \f$ (poloidal) direction
    @param nphi Number of points (per field period) in the \f$ \phi \f$ (toroidal) direction
    @param full_torus Construct grid for the full torus (default: one field period)
    @result `rgrid` Structed phi-theta grid of points [nphi,ntheta,3], `nfp` Number of field periods
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
        nfp = file['nfp'][0]
        if full_torus:
            nphi *= nfp
            nfp = 1
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
    r'''! Build a uniform grid from a toroidal surface B-norm definition file

    @param filename Filename of B-norm poloidal mode definition
    @param nsample Number of points in the \f$ \theta \f$ (poloidal) direction
    @param nphi Number of points in the \f$ \phi \f$ (toroidal) direction
    @param resample_type Construct grid for the full torus (default: one field period)
    @param use_spline Fit a spline to the boundary to produce a smoother representation?
    @result `rgrid` Structed phi-theta grid of points [nphi,ntheta,3], `bnorm` Bn on same grid [nphi,ntheta], `nfp` Number of field periods
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
    sin_sum = mode_resample[:,2].sum(); sin_sum_abs = max(1.E-14,abs(mode_resample[:,2]).sum())
    cos_sum = mode_resample[:,5].sum(); cos_sum_abs = max(1.E-14,abs(mode_resample[:,5]).sum())
    print('  Mode pair sums {0:.4E} {1:.4E}'.format(sin_sum,cos_sum))
    if (abs(sin_sum/sin_sum_abs) > 1.E-4) or (abs(cos_sum/cos_sum_abs) > 1.E-4):
        print('Warning: Large net flux present in one or more modes! This may indicate an error.')
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


class ThinCurr_periodic_toroid:
    '''! Helper class for working with periodic toroid models'''
    def __init__(self,r_grid,nfp,ntheta,nphi):
        r'''! Build triangular mesh for the full surface from a uniform grid of a single field period (toroidal)

        @param r_grid Uniform grid [nphi,ntheta,3] (\f$ \phi \f$ and \f$ \theta \f$ vary along the first and second dimension respectively)
        @param nfp Number of field periods
        @param ntheta Number of node points in the poloidal direction
        @param nphi Number of node points (per field period) in the toroidal direction
        '''
        self.nfp = nfp
        self.nphi = nphi
        self.ntheta = ntheta
        if nfp == 1:
            self.r, self.lc, self.tnodeset, pnodeset = build_triangles_from_grid(r_grid)
            self.pnodesets = [pnodeset,]
            self.r_map = numpy.s_[:]
        else:
            nphi = r_grid.shape[0]
            ntheta = r_grid.shape[1]
            r_base = r_grid.copy().reshape((nphi*ntheta,3))
            r_full = r_base[:-ntheta,:].copy()
            self.r_map = numpy.arange((nphi-1)*ntheta, dtype=numpy.int32)
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
                self.r_map = numpy.hstack((self.r_map,numpy.arange((nphi-1)*ntheta, dtype=numpy.int32)))
            self.r, self.lc, self.tnodeset, pnodeset = build_triangles_from_grid(r_full.reshape((int(r_full.shape[0]/(ntheta)),ntheta,3)))
            self.pnodesets = [pnodeset]
            for i in range(1,nfp):
                self.pnodesets.append(pnodeset+i*(nphi-1)*ntheta)
    
    def plot_mesh(self,fig,equal_aspect=True,surf_alpha=0.1,surf_cmap='viridis'):
        '''! Plot mesh and holes

        @param fig Figure to use for plots (must be empty)
        @param equal_aspect Set plot aspect ratio to be more physically realistic
        @param surf_alpha Transparency of surface in hole plot (left plot)
        @param surf_cmap Colormap for mesh plot (right plot)
        '''
        ax = fig.add_subplot(1,2,1, projection='3d')
        ax.plot(self.r[self.tnodeset,0], self.r[self.tnodeset,1], self.r[self.tnodeset,2], c='tab:red')
        for pnodeset in self.pnodesets:
            pnodeset_tmp = numpy.append(pnodeset, (pnodeset[0],))
            ax.plot(self.r[pnodeset_tmp,0], self.r[pnodeset_tmp,1], self.r[pnodeset_tmp,2], c='tab:blue')
        ax.plot_trisurf(self.r[:,0], self.r[:,1], self.r[:,2], triangles=self.lc, color=[0.0,0.0,0.0,surf_alpha])
        if equal_aspect:
            mesh_range = numpy.ptp(self.r,axis=0)
            ax.set_box_aspect(mesh_range)
        ax = fig.add_subplot(1,2,2, projection='3d')
        ax.plot_trisurf(self.r[:,0], self.r[:,1], self.r[:,2], triangles=self.lc, cmap=surf_cmap, antialiased=False)
        if equal_aspect:
            ax.set_box_aspect(mesh_range)
    
    def write_to_file(self,filename,reg=None,include_closures=True):
        '''! Save mesh to file in ThinCurr format

        @param filename Filename for mesh file
        @param reg Region list [nc,]
        @param include_closures Include closure defitions in file?
        '''
        if include_closures:
            closures = numpy.arange(self.nfp)*int(self.lc.shape[0]/self.nfp)+1
        else:
            closures = []
        if reg is None:
            reg = numpy.ones((self.lc.shape[0],))
        if self.nfp == 1:
            write_ThinCurr_mesh(filename, self.r, self.lc+1, reg,
                holes=[self.tnodeset+1, self.pnodesets[0]+1], closures=closures)
        else:
            write_ThinCurr_mesh(filename, self.r, self.lc+1, reg,
                holes=[self.tnodeset+1] + [pnodeset+1 for pnodeset in self.pnodesets], closures=closures, pmap=self.r_map, nfp=self.nfp)

    def condense_matrix(self,matrix,axis=None):
        '''! Condense matrix to unique DOF only by combining poloidal nodesets

        @param matrix Initial matrix in full ThinCurr representation (includes copies of poloidal nodesets)
        @param axis Axis to act on (both if `None`)
        @returns Condensed matrix
        '''
        # Condense model to single mode period
        if self.nfp > 1:
            nelems_new = matrix.shape[0]-self.nfp+1
            if axis is None:
                matrix_new = numpy.zeros((nelems_new,nelems_new))
                matrix_new[:-1,:-1] = matrix[:-self.nfp,:-self.nfp]
                matrix_new[:-1,-1] = matrix[:-self.nfp,-self.nfp:].sum(axis=1)
                matrix_new[-1,:-1] = matrix[-self.nfp:,:-self.nfp].sum(axis=0)
                matrix_new[-1,-1] = matrix[-self.nfp:,-self.nfp:].sum(axis=None)
            elif axis == 0:
                matrix_new = numpy.zeros((nelems_new,matrix.shape[1]))
                matrix_new[:-1,:] = matrix[:-self.nfp,:]
                matrix_new[-1,:] = matrix[-self.nfp:,:].sum(axis=0)
            elif axis == 1:
                matrix_new = numpy.zeros((matrix.shape[0],nelems_new))
                matrix_new[:,:-1] = matrix[:,:-self.nfp]
                matrix_new[:,-1] = matrix[:,-self.nfp:].sum(axis=1)
            else:
                raise ValueError('Invalid value for "axis", must be (None,0,1)')
            return matrix_new
        else:
            return matrix
    
    def nodes_to_unique(self,vector,tor_val=0.0,pol_val=0.0,remove_closure=True):
        '''! Maps a vector of values on all nodes within a single field period
        to unique DOFs

        @param vector Vector of values on nodes
        @param tor_val Value of toroidal hole element
        @param pol_val Value of poloidal hole element
        @param remove_closure Remove closure element?
        '''
        if self.nfp > 1:
            vector = vector[:-self.pnodesets[0].shape[0]]
        if remove_closure:
            return numpy.r_[vector[1:],tor_val,pol_val]
        else:
            return numpy.r_[vector,tor_val,pol_val]

    def expand_vector(self,vector):
        '''! Expand vector from unique DOF by duplicating poloidal nodesets

        @param vector Vector of values for unique DOF
        @returns Expanded vector
        '''
        if self.nfp == 1:
            return vector
        output = numpy.zeros((vector.shape[0]+(self.nfp-1)))
        output[:-self.nfp+1] = vector
        output[-self.nfp+1:] = output[-self.nfp]
        return output

    def unique_to_nodes_2D(self,data):
        '''! Maps a vector of values for unique DOFs to all nodes within a single field period
        and reshape to a structured 2D grid

        @param vector Vector of values on nodes
        @param tor_val Value of toroidal hole element
        @param pol_val Value of poloidal hole element
        @param remove_closure Remove closure element?
        '''
        data_out = numpy.zeros((self.nphi,self.ntheta))
        if self.nfp > 1:
            data_out[:-1,:] = numpy.r_[0.0,data[:-2]].reshape((self.nphi-1,self.ntheta))
            data_out[-1,:] = data_out[0,:] + data[-1]
        else:
            data_out = numpy.r_[0.0,data[:-2]].reshape((self.nphi,self.ntheta))
        return data_out.transpose()


def build_triangles_from_grid(data_grid,wrap_n=True,wrap_m=True):
    r'''! Build triangles from a uniform grid of points

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