'''! Python interface for TokaMaker Grad-Shafranov functionality

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''
import json
import math
import numpy
from .._interface import *

## @cond
class triangle_struct(c_struct):
    r'''! Triangle library I/O structure
    '''
    _fields_ = [
        ("pointlist", c_double_ptr),
        ("pointattributelist", c_double_ptr),
        ("pointmarkerlist", c_int_ptr),
        ("numberofpoints", c_int),
        ("numberofpointattributes", c_int),
        ("trianglelist", c_int_ptr),
        ("triangleattributelist", c_double_ptr),
        ("trianglearealist", c_double_ptr),
        ("neighborlist", c_int_ptr),
        ("numberoftriangles", c_int),
        ("numberofcorners", c_int),
        ("numberoftriangleattributes", c_int),
        ("segmentlist", c_int_ptr),
        ("segmentmarkerlist", c_int_ptr),
        ("numberofsegments", c_int),
        ("holelist", c_double_ptr),
        ("numberofholes", c_int),
        ("regionlist", c_double_ptr),
        ("numberofregions", c_int),
        ("edgelist", c_int_ptr),
        ("edgemarkerlist", c_int_ptr),
        ("normlist", c_double_ptr),
        ("numberofedges", c_int)
    ]

oft_triangulate = ctypes_subroutine(oft_triangle_lib.triangulate,
    [c_char_p, ctypes.POINTER(triangle_struct), ctypes.POINTER(triangle_struct), ctypes.POINTER(triangle_struct)], c_int)
## @endcond

def run_triangle(alpha):
    # # Old version using Python "triangle" package
    # try:
    #     import triangle as tr
    # except ImportError:
    #     print('Meshing requires "triangle" python library')
    #     return None
    # return tr.triangulate(alpha,'pqaA')

    # New version using internal interface
    in_struct = triangle_struct()
    point_list = numpy.ascontiguousarray(alpha['vertices'], dtype=numpy.float64)
    in_struct.pointlist = point_list.ctypes.data_as(c_double_ptr)
    in_struct.numberofpoints = point_list.shape[0]
    segments = numpy.ascontiguousarray(alpha['segments'], dtype=numpy.int32)
    in_struct.segmentlist = segments.ctypes.data_as(c_int_ptr)
    in_struct.numberofsegments = segments.shape[0]
    regions = numpy.ascontiguousarray(alpha['regions'], dtype=numpy.float64)
    in_struct.regionlist = regions.ctypes.data_as(c_double_ptr)
    in_struct.numberofregions = regions.shape[0]
    out_struct = triangle_struct()
    vor_struct = triangle_struct()
    setting_string = c_char_p(b'pqaAzQ')
    errval = oft_triangulate(setting_string, ctypes.byref(in_struct), ctypes.byref(out_struct), ctypes.byref(vor_struct))
    if errval != 0:
        raise ValueError("Meshing failed!")
    vertices = numpy.ctypeslib.as_array(out_struct.pointlist,shape=(out_struct.numberofpoints,2))
    triangles = numpy.ctypeslib.as_array(out_struct.trianglelist,shape=(out_struct.numberoftriangles,3))
    triangle_attributes = numpy.ctypeslib.as_array(out_struct.triangleattributelist,shape=(out_struct.numberoftriangles,))
    return {'vertices': vertices, 'triangles': triangles, 'triangle_attributes': triangle_attributes}


class gs_Domain:
    '''! Grad-Sharvanov domain definitions for TokaMaker with Triangle library meshing'''
    def __init__(self,rextent=None,zextents=[None,None],rpad=1.2,zpad=[1.2,1.2],json_filename=None):
        '''! Create a new Grad-Shafranov domain

        @param region_list List of @ref oftpy.Region objects that define mesh
        @param merge_thresh Distance threshold for merging nearby points
        '''
        self._r = None
        self._lc = None
        self._reg = None
        if json_filename is not None:
            with open(json_filename, 'r') as fid:
                input_dict = json.load(fid)
            self.rextent = input_dict['rextent']
            self.zextents = input_dict['zextents']
            self.rpad = input_dict['rpad']
            self.zpad = input_dict['zpad']
            self.rmax = input_dict['rmax']
            self.zmin = input_dict['zmin']
            self.zmax = input_dict['zmax']
            self.boundary_reg = input_dict['boundary_reg']
            self.reg_type_counts = input_dict['reg_type_counts']
            self.region_info = input_dict['region_info']
            self.regions = []
            for region in input_dict['regions']:
                region.append(Region(load_dict=region))
        else:
            self.rextent = rextent
            self.zextents = zextents
            self.rpad = None
            if self.rextent is None:
                self.rpad = rpad
            self.zpad = [None,None]
            if self.zextents[0] is None:
                self.zpad[0] = zpad[0]
            if self.zextents[1] is None:
                self.zpad[1] = zpad[1]
            self.rmax = 0.0
            self.zmin = 0.0
            self.zmax = 0.0
            self.boundary_reg = None
            self.regions = []
            self.reg_type_counts = {
                'plasma': 0,
                'vacuum': 0,
                'boundary': 0,
                'conductor': 0,
                'coil': 0,
            }
            self.region_info = {}
            self._extra_reg_defs = []
    
    def define_region(self,name,dx,reg_type,eta=None,noncontinuous=None,nTurns=None,coil_set=None,allow_xpoints=False):
        '''! Define a new region and its properties (geometry is given in a separate call)

        @param name Name of region
        @param dx Target mesh size for region
        @param reg_type Type of region ("plasma", "vacuum", "boundary", "conductor", or "coil")
        @param eta Resistivity for "conductor" regions (raises error if region is other type)
        @param nTurns Number of turns for "coil" regions (raises error if region is other type)
        @param allow_xpoints Allow X-points in this region (for non-plasma regions only)
        '''
        if (dx is None) or (dx < 0.0):
            raise ValueError('"dx" must have a non-negative value')
        name = name.upper()
        if name[0] == '#':
            raise ValueError('Invalid value for "name", cannot start with "#"')
        if (name in self.region_info):
            raise KeyError('Region already exists!')
        next_id = -1
        if reg_type == 'plasma':
            next_id = 1
            allow_xpoints = True
        elif reg_type == 'vacuum':
            pass
        elif reg_type == 'boundary':
            self.boundary_reg = name
        elif reg_type in ('conductor', 'coil'):
            pass
        else:
            raise ValueError("Unknown region type")
        if next_id < 0:
            next_id = len(self.region_info) - self.reg_type_counts['plasma'] + 2
        self.reg_type_counts[reg_type] += 1
        self.region_info[name] = {
            'id': next_id,
            'dx': dx,
            'count': 0,
            'type': reg_type,
            'allow_xpoints': allow_xpoints
        }
        if eta is not None:
            if reg_type != 'conductor':
                raise ValueError('Resistivity specification only valid for "conductor" regions')
            else:
                self.region_info[name]['eta'] = eta
        else:
            if reg_type == 'conductor':
                raise ValueError('Resistivity not specified for "conductor" region')
        if noncontinuous is not None:
            if reg_type != 'conductor':
                raise ValueError('Non-contiguous specification only valid for "conductor" regions')
            else:
                self.region_info[name]['noncontinuous'] = noncontinuous
        if nTurns is not None:
            if reg_type != 'coil':
                raise ValueError('"nTurns" specification only valid for "coil" regions')
            else:
                self.region_info[name]['nturns'] = nTurns
        if coil_set is not None:
            if reg_type != 'coil':
                raise ValueError('"coil_set" specification only valid for "coil" regions')
            else:
                if coil_set[0] == '#':
                    raise ValueError('Invalid value for "coil_set", cannot start with "#"')
                self.region_info[name]['coil_set'] = coil_set
        

    def add_annulus(self,inner_countour,inner_name,outer_contour,annulus_name,parent_name=None,angle_tol=30.0,sliver_tol=120.0,small_thresh=None):
        '''! Add annular geometry defining region boundaries to the mesh

        @param inner_contour Curve defining inner boundary
        @param inner_name Name of region enclosed by the inner boundary
        @param outer_contour Curve defining outer boundary
        @param annulus_name Name of annular region between inner and outer boundaries
        @param parent_name Name of region beyond the outer boundary
        @param angle_tol Corner tolerance used when resampling curve at desired resolution
        @param sliver_tol Tolerance used for "sliver" region warnings
        @param small_thresh Tolerance used for "small" curve warnings (default: dx/2)
        '''
        inner_countour = numpy.array(inner_countour)
        outer_contour = numpy.array(outer_contour)
        inner_name = inner_name.upper()
        annulus_name = annulus_name.upper()
        inner_reg = self.region_info.get(inner_name, None)
        annulus_reg = self.region_info.get(annulus_name, None)
        if inner_reg is None:
            raise KeyError('Region "{0}" not defined'.format(inner_name))
        else:
            inner_dx = inner_reg.get('dx', None)
            inner_dx_curve = inner_dx
            if inner_dx is None:
                raise ValueError('Resolution for region "{0}" not defined'.format(inner_name))
            if inner_countour[:,0].min() < 0.0:
                raise ValueError('Negative radial value detected in "inner_countour"')
        if annulus_reg is None:
            raise KeyError('Region "{0}" not defined'.format(annulus_name))
        else:
            annulus_dx = annulus_reg.get('dx', None)
            outer_dx_curve = annulus_dx
            if annulus_dx is None:
                raise ValueError('Resolution for region "{0}" not defined'.format(annulus_name))
            else:
                inner_dx_curve = min(inner_dx_curve,annulus_dx)
            if outer_contour[:,0].min() < 0.0:
                raise ValueError('Negative radial value detected in "outer_contour"')
        if parent_name is not None:
            parent_name = parent_name.upper()
            parent_reg = self.region_info.get(parent_name, None)
            if parent_reg is None:
                raise KeyError('Region "{0}" not defined'.format(parent_name))
            else:
                parent_dx = parent_reg.get('dx', None)
                if parent_dx is None:
                    raise ValueError('Resolution for region "{0}" not defined'.format(parent_name))
                else:
                    outer_dx_curve = min(outer_dx_curve,parent_dx)
        # Add inner region
        maxes = inner_countour.max(axis=0)
        self.rmax = max(self.rmax,maxes[0])
        self.zmax = max(self.zmax,maxes[1])
        self.zmin = min(self.zmin,inner_countour[:,1].min())
        self.regions.append(Region(inner_countour,inner_dx,inner_dx_curve,angle_tol,sliver_tol,small_thresh,inner_reg["id"]))
        inner_reg["count"] += 1
        # Add outer region
        maxes = outer_contour.max(axis=0)
        self.rmax = max(self.rmax,maxes[0])
        self.zmax = max(self.zmax,maxes[1])
        self.zmin = min(self.zmin,outer_contour[:,1].min())
        self.regions.append(Region(outer_contour,annulus_dx,outer_dx_curve,angle_tol,sliver_tol,small_thresh,annulus_reg["id"]))
        annulus_reg["count"] += 1
    
    def add_polygon(self,contour,name,parent_name=None,angle_tol=30.0,sliver_tol=120.0,small_thresh=None):
        '''! Add polygon geometry defining region boundaries to the mesh

        @param contour Curve defining polygon
        @param name Name of region enclosed by the polygon
        @param parent_name Name of region outside the polygon
        @param angle_tol Corner tolerance used when resampling curve at desired resolution
        @param sliver_tol Tolerance used for "sliver" region warnings
        @param small_thresh Tolerance used for "small" curve warnings (default: dx/2)
        '''
        contour = numpy.array(contour)
        name = name.upper()
        reg = self.region_info.get(name, None)
        if reg is None:
            raise KeyError('Region "{0}" not defined'.format(name))
        else:
            dx = reg.get('dx', None)
            dx_curve = dx
            if dx is None:
                raise ValueError('Resolution for region "{0}" not defined'.format(name))
            if contour[:,0].min() < 0.0:
                raise ValueError('Negative radial value detected in "inner_countour"')
        if parent_name is not None:
            parent_name = parent_name.upper()
            parent_reg = self.region_info.get(parent_name, None)
            if parent_reg is None:
                raise KeyError('Region "{0}" not defined'.format(parent_name))
            else:
                parent_dx = parent_reg.get('dx', None)
                if parent_dx is None:
                    raise ValueError('Resolution for region "{0}" not defined'.format(parent_name))
                else:
                    dx = min(dx,parent_dx)
        # Add region
        maxes = contour.max(axis=0)
        self.rmax = max(self.rmax,maxes[0])
        self.zmax = max(self.zmax,maxes[1])
        self.zmin = min(self.zmin,contour[:,1].min())
        self.regions.append(Region(contour,dx,dx_curve,angle_tol,sliver_tol,small_thresh,reg["id"]))
        reg["count"] += 1

    def add_rectangle(self,rc,zc,w,h,name,parent_name=None, rot=None):
        '''! Add rectangular geometry defining region boundaries to the mesh

        @param rc Radial center of rectangle
        @param zc Vertical center of rectangle
        @param w Width of rectangle (radial direction)
        @param h Height of the rectangle (vertical direction)
        @param name Name of region enclosed by the polygon
        @param parent_name Name of region outside the polygon
        @param rot Rotation of rectangle (degrees)
        '''
        contour = numpy.asarray([
            [-w/2.0, -h/2.0],
            [+w/2.0, -h/2.0],
            [+w/2.0, +h/2.0],
            [-w/2.0, +h/2.0]
        ])

        if rot is not None:
            rot = numpy.deg2rad(rot)
            rotmat = numpy.asarray([numpy.cos(rot), -numpy.sin(rot), numpy.sin(rot), numpy.cos(rot)]).reshape((2,2))

            contour = numpy.dot(contour,rotmat.T)
        
        contour[:,0] += rc
        contour[:,1] += zc

        self.add_polygon(contour,name,parent_name)
    
    def add_enclosed(self,in_point,name):
        name = name.upper()
        reg = self.region_info.get(name, None)
        if reg is None:
            raise KeyError('Region "{0}" not defined'.format(name))
        else:
            dx = reg['dx']
            id = reg['id']
            self._extra_reg_defs.append([in_point[0], in_point[1], id, dx*dx/2.0])
            reg["count"] += 1
    
    def get_coils(self):
        '''! Get dictionary describing coil regions in domain

        @result Dictionary of coil regions and attributes
        '''
        coil_list = {}
        coil_id = 0
        for key in self.region_info:
            if self.region_info[key]['type'] == 'coil':
                coil_list[key] = {
                    'reg_id': self.region_info[key]['id'],
                    'coil_id': coil_id,
                    'nturns': self.region_info[key].get('nturns',1),
                    'coil_set': self.region_info[key].get('coil_set',key),
                    'allow_xpoints': self.region_info[key].get('allow_xpoints',False)
                }
                coil_id += 1
        return coil_list

    def get_conductors(self):
        '''! Get dictionary describing conducting regions in domain

        @result Dictionary of conducting regions and attributes
        '''
        cond_list = {}
        cond_id = 0
        vac_id = 0
        for key in self.region_info:
            if self.region_info[key]['type'] == 'conductor':
                cond_list[key] = {
                    'reg_id': self.region_info[key]['id'],
                    'cond_id': cond_id,
                    'eta': self.region_info[key]['eta'],
                    'noncontinuous': self.region_info[key].get('noncontinuous',False),
                    'allow_xpoints': self.region_info[key].get('allow_xpoints',False)
                }
                cond_id += 1
            elif self.region_info[key]['type'] in ('vacuum','boundary'):
                cond_list[key] = {
                    'reg_id': self.region_info[key]['id'],
                    'vac_id': vac_id,
                    'allow_xpoints': self.region_info[key].get('allow_xpoints',False)
                }
                vac_id += 1
        return cond_list
    
    def build_mesh(self,debug=False,merge_thresh=1.E-4,require_boundary=True,setup_only=False):
        '''! Build mesh for specified domains

        @result Meshed representation (pts[np,2], tris[nc,3], regions[nc])
        '''
        # Check for single plasma region
        if self.reg_type_counts['plasma'] > 1:
            raise ValueError('More than one plasma region specified')
        elif self.reg_type_counts['plasma'] == 0:
            raise ValueError('No plasma region specified')
        else:
            # Make sure a boundary exists if we have regions other than plasma
            if ((self.reg_type_counts['vacuum'] > 0) or (self.reg_type_counts['coil'] > 0) or (self.reg_type_counts['conductor'] > 0)) and require_boundary:
                if self.boundary_reg is None:
                    raise ValueError('No boundary region specified')
                # Check or set extents
                if self.rextent is None:
                    self.rextent = self.rpad*self.rmax
                else:
                    if self.rmax > self.rextent:
                        raise ValueError('User specified "rextent" does not enclose all regions')
                if self.zextents[0] is None:
                    self.zextents[0] = self.zpad[0]*self.zmin
                else:
                    if self.zmin < self.zextents[0]:
                        raise ValueError('User specified "zextents[0]" does not enclose all regions')
                if self.zextents[1] is None:
                    self.zextents[1] = self.zpad[1]*self.zmax
                    if self.zmax > self.zextents[1]:
                        raise ValueError('User specified "zextents[1]" does not enclose all regions')
                # Create boundary region
                vac_dx = self.region_info[self.boundary_reg]['dx']
                if vac_dx is None:
                    raise ValueError('Resolution for region "vacuum" not defined')
                vac_contour = numpy.asarray([
                    [0.0,           self.zextents[0]],
                    [self.rextent, self.zextents[0]],
                    [self.rextent, self.zextents[1]],
                    [0.0,           self.zextents[1]]
                ])
                self.regions.append(Region(vac_contour,vac_dx,vac_dx,id=self.region_info[self.boundary_reg]['id']))
                self.region_info[self.boundary_reg]['count'] += 1
        # Check for undefined regions
        for key in self.region_info:
            if self.region_info[key]['count'] == 0:
                raise KeyError('Region "{0}" defined but never created'.format(key))
        # Re-index regions
        reg_reorder = [-1 for _ in self.region_info]
        reg_reorder[0] = 1
        reg_id = 1
        for reg_type in ('boundary', 'vacuum','conductor','coil'):
            for key in self.region_info:
                if self.region_info[key]['type'] == reg_type:
                    reg_id += 1
                    reg_reorder[self.region_info[key]['id']-1] = reg_id
                    self.region_info[key]['id'] = reg_id
        for region in self.regions:
            region._id = reg_reorder[region._id-1]
        for point_def in self._extra_reg_defs:
            point_def[2] = reg_reorder[point_def[2]-1]
        # Generate mesh
        self.mesh = Mesh(self.regions,debug=debug,extra_reg_defs=self._extra_reg_defs,merge_thresh=merge_thresh)
        if not setup_only:
            self._r, self._lc, self._reg = self.mesh.get_mesh()
        return self._r, self._lc, self._reg
    
    def save_json(self,filename):
        '''! Create a JSON file containing a description of the mesh 

        @param filename Path to create JSON file
        '''
        output_dict = {
            'rextent': self.rextent,
            'zextents': self.zextents,
            'rpad': self.rpad,
            'zpad': self.zpad,
            'rmax': self.rmax,
            'zmin': self.zmin,
            'zmax': self.zmax,
            'boundary_reg': self.boundary_reg,
            'reg_type_counts': self.reg_type_counts,
            'regions': [],
            'region_info': self.region_info

        }
        for region in self.regions:
            output_dict['regions'].append(region.get_dict())
        with open(filename, 'w+') as fid:
            fid.write(json.dumps(output_dict))
    
    def plot_topology(self,fig,ax,linewidth=None):
        '''! Plot mesh topology

        @param fig Figure to add to (unused)
        @param ax Axes to add to (must be scalar)
        @param linewidth Line width for plots
        '''
        for region in self.regions:
            region.plot_segments(fig,ax,linewidth=linewidth)
        # Format plot
        ax.set_aspect('equal','box')
        ax.set_xlabel('R (m)')
        ax.set_ylabel('Z (m)')
    
    def plot_mesh(self,fig,ax,lw=0.5,show_legends=True,col_max=10,split_coil_sets=False):
        '''! Plot machine geometry

        @param fig Figure to add to (unused)
        @param ax Axes to add to (must be scalar, [2], or [2,2])
        @param lw Width of lines in calls to "triplot()"
        @param show_legends Show legends for plots with more than one region?
        @param col_max Maximum number of entries per column in each legend
        '''
        if self._r is None:
            raise ValueError('"plot_mesh()" can only be called after "build_mesh()"')
        # Get format type from shape of axis object
        format_type = -1
        try:
            if (ax.shape[0] == 2):
                format_type = 1
                try:
                    if (ax.shape[1] == 2):
                        format_type = 2
                    else:
                        format_type = -1
                except:
                    pass
            else:
                format_type = -1
        except:
            format_type = 0
        if format_type < 0:
            raise ValueError("Axes for plotting must be scalar, [2], or [2,2]")
        # Set appropriate axes based on format type
        if format_type == 0:
            plasma_axis = ax
            cond_axis = ax
            coil_axis = ax
            vac_axis = ax
            ax_flat = [ax]
        elif format_type == 1:
            plasma_axis = ax[0]
            vac_axis = ax[0]
            cond_axis = ax[1]
            coil_axis = ax[1]
            ax_flat = ax.flatten()
        else:
            plasma_axis = ax[1,0]
            cond_axis = ax[0,1]
            coil_axis = ax[1,1]
            vac_axis = ax[0,0]
            ax_flat = ax.flatten()
        # Get region count
        nregs = self._reg.max()
        reg_mark = numpy.zeros((nregs,))
        reg_mark[0] = 1
        # Plot the plasma region
        plasma_axis.triplot(self._r[:,0],self._r[:,1],self._lc[self._reg==1,:],lw=lw,label='Plasma')
        # Plot conductor regions
        nCond = 0
        for key, cond in self.get_conductors().items():
            if 'vac_id' in cond:
                continue
            nCond += 1
            reg_mark[cond['reg_id']-1] = 1
            cond_axis.triplot(self._r[:,0],self._r[:,1],self._lc[self._reg==cond['reg_id'],:],lw=lw,label=key)
        # Plot coil regions
        coil_colors = {}
        for key, coil in self.get_coils().items():
            reg_mark[coil['reg_id']-1] = 1
            if split_coil_sets:
                leg_key = key
            else:
                leg_key = coil.get('coil_set',key)
            if leg_key not in coil_colors:
                lines, _ = coil_axis.triplot(self._r[:,0],self._r[:,1],self._lc[self._reg==coil['reg_id'],:],lw=lw,label=leg_key)
                coil_colors[leg_key] = lines.get_color()
            else:
                coil_axis.triplot(self._r[:,0],self._r[:,1],self._lc[self._reg==coil['reg_id'],:],lw=lw,color=coil_colors[leg_key])
        nCoil = len(coil_colors)
        # Plot the vacuum regions
        nVac = 0
        for i in range(nregs):
            if reg_mark[i] == 0:
                nVac += 1
                vac_axis.triplot(self._r[:,0],self._r[:,1],self._lc[self._reg==i+1,:],lw=lw,label='Vacuum_{0}'.format(nVac))
        # Format plots
        for ax_tmp in ax_flat:
            ax_tmp.set_aspect('equal','box')
            ax_tmp.set_xlabel('R (m)')
            ax_tmp.set_ylabel('Z (m)')
        if show_legends:
            if format_type == 0:
                ncols = max(1,math.floor((1+nCond+nCoil+nVac)/col_max))
                plasma_axis.legend(bbox_to_anchor=(1.05,0.5), loc='center left', ncols=ncols)
            elif format_type == 1:
                ncols = max(1,math.floor((1+nVac)/col_max))
                plasma_axis.legend(bbox_to_anchor=(1.05,0.5), loc='center left', ncols=ncols)
                ncols = max(1,math.floor((nCond+nCoil)/col_max))
                cond_axis.legend(bbox_to_anchor=(1.05,0.5), loc='center left', ncols=ncols)
            elif format_type == 2:
                ncols = max(1,math.floor((nCond)/col_max))
                cond_axis.legend(bbox_to_anchor=(1.05,0.5), loc='center left', ncols=ncols)
                ncols = max(1,math.floor((nCoil)/col_max))
                coil_axis.legend(bbox_to_anchor=(1.05,0.5), loc='center left', ncols=ncols)


def save_gs_mesh(pts,tris,regions,coil_dict,cond_dict,filename,use_hdf5=True):
        '''! Save G-S mesh to file in HDF5 format

        @param pts[np,2] Vertex list
        @param tris[nc,3] Cell list
        @param regions[nc] Region list
        @param coil_dict Coil region dictionary
        @param cond_dict Conducting region dictionary
        @param filename Path to create HDF5 mesh file
        '''
        if use_hdf5:
            import h5py
            coil_json = json.dumps(coil_dict)
            cond_json = json.dumps(cond_dict)
            with h5py.File(filename, 'w') as h5_file:
                h5_file.create_dataset('mesh/r', data=pts, dtype='f8')
                h5_file.create_dataset('mesh/lc', data=tris, dtype='i4')
                h5_file.create_dataset('mesh/reg', data=regions, dtype='i4')
                string_datatype = h5py.string_dtype('ascii')
                h5_file.create_dataset('mesh/coil_dict', data=coil_json, dtype=string_datatype)
                h5_file.create_dataset('mesh/cond_dict', data=cond_json, dtype=string_datatype)
        else:
            with open(filename, 'w+') as fid:
                fid.write(json.dumps({
                    'mesh': {
                        'r': pts.tolist(),
                        'lc': tris.tolist(),
                        'reg': regions.tolist(),
                        'coil_dict': coil_dict,
                        'cond_dict': cond_dict
                    }
                }))


def load_gs_mesh(filename,use_hdf5=True):
        '''! Load G-S mesh to file in HDF5 format

        @param filename Path to HDF5 mesh file
        @result pts[np,2], tris[nc,3], regions[nc], coil_dict, cond_dict
        '''
        if use_hdf5:
            import h5py
            with h5py.File(filename, 'r') as h5_file:
                pts = numpy.asarray(h5_file['mesh/r'])
                tris = numpy.asarray(h5_file['mesh/lc'])
                regions = numpy.asarray(h5_file['mesh/reg'])
                coil_dict = json.loads(h5_file['mesh/coil_dict'][()])
                cond_dict = json.loads(h5_file['mesh/cond_dict'][()])
        else:
            with open(filename, 'r') as fid:
                input_dict = json.load(fid)
            pts = numpy.asarray(input_dict['mesh']['r'])
            tris = numpy.asarray(input_dict['mesh']['lc'])
            regions = numpy.asarray(input_dict['mesh']['reg'])
            coil_dict = input_dict['mesh']['coil_dict']
            cond_dict = input_dict['mesh']['cond_dict']
        return pts, tris, regions, coil_dict, cond_dict


class Mesh:
    '''! Mesh builder class for triangle library'''
    def __init__(self,region_list,merge_thresh=1.E-4,debug=False,extra_reg_defs=[]):
        '''! Initialize Mesh builder object

        @param region_list List of @ref oftpy.Region objects that define mesh
        @param merge_thresh Distance threshold for merging nearby points
        '''
        self._merge_thresh = merge_thresh
        self._unique_points = []
        self._reg_seg_map = []
        self._segments = []
        print('Assembling regions:')
        # Build list of unique points
        for ireg, region in enumerate(region_list):
            local_seg_map = []
            reg_pt_map = [i+len(self._unique_points) for i in range(region._points.shape[0])]
            ilocal = len(self._unique_points)-1
            for tmp_pts in region._segments:
                tmp_pt_map = [reg_pt_map[i] for i in tmp_pts]
                for ipt, reg_id in enumerate(tmp_pts):
                    reg_pt = region._points[reg_id,:]
                    if reg_pt_map[reg_id] > len(self._unique_points)-1:
                        for jpt, unique_pt in enumerate(self._unique_points):
                            if numpy.linalg.norm(reg_pt-unique_pt) < merge_thresh:
                                if debug:
                                    print('  Merging points:',ireg,reg_pt_map[reg_id],jpt,reg_pt,unique_pt)
                                tmp_pt_map[ipt] = jpt
                                reg_pt_map[reg_id] = jpt
                                break
                        else:
                            ilocal += 1
                            tmp_pt_map[ipt] = ilocal
                            reg_pt_map[reg_id] = ilocal
                            self._unique_points.append(reg_pt)
                for iseg, segment in enumerate(self._segments):
                    nOverlap = 0
                    for test_pt in segment[0]:
                        if test_pt in tmp_pt_map:
                            nOverlap += 1
                    if (nOverlap > 1) and (len(tmp_pts) == len(segment[0])): # Full overlap
                        # Look forward
                        for i,test_pt in enumerate(segment[0]):
                            if tmp_pt_map[i] != test_pt:
                                break
                        else: # Matched segment
                            if debug:
                                print('  Merging curve segments:',ireg,iseg)
                            segment[1] = min(segment[1],region._dx_curve)
                            segment[2] = min(segment[2],region._small_thresh)
                            local_seg_map.append(iseg)
                            break
                        # Look backward
                        for i,test_pt in enumerate(segment[0]):
                            if tmp_pt_map[-i-1] != test_pt:
                                break
                        else: # Matched segment
                            if debug:
                                print('  Merging curve segments:',ireg,iseg)
                            segment[1] = min(segment[1],region._dx_curve)
                            segment[2] = min(segment[2],region._small_thresh)
                            local_seg_map.append(-iseg)
                            break
                    elif  (nOverlap > 1): # Partial match
                        if debug:
                            print('  Merging partially overlapping curve segments:',ireg,iseg)
                        if len(tmp_pts) < len(segment[0]):
                            overlap_start = len(segment[0])
                            overlap_end = -1
                            first_pt = -1
                            for i, test_pt in enumerate(segment[0]):
                                if test_pt in tmp_pt_map:
                                    overlap_start = min(overlap_start,i)
                                    overlap_end = max(overlap_end,i)
                                    if first_pt < 0 :
                                        first_pt = tmp_pt_map.index(test_pt)
                            segment_split = [iseg, -1, -1]
                            if overlap_start > 0:
                                self._segments.append([segment[0][:overlap_start+1], segment[1], segment[2]])
                                segment_split[1] = len(self._segments)-1
                            if overlap_end < len(segment[0])-1:
                                self._segments.append([segment[0][overlap_end:], segment[1], segment[2]])
                                segment_split[2] = len(self._segments)-1
                            segment[0] = segment[0][overlap_start:overlap_end+1]
                            segment[1] = min(segment[1],region._dx_curve)
                            segment[2] = min(segment[2],region._small_thresh)
                            #
                            for reg_seg_map in self._reg_seg_map:
                                try:
                                    ifound = reg_seg_map.index(iseg)
                                    if segment_split[1] >= 0:
                                        reg_seg_map.insert(ifound,segment_split[1])
                                        ifound += 1
                                    if segment_split[2] >= 0:
                                        reg_seg_map.insert(ifound+1,segment_split[2])
                                except ValueError:
                                    ifound = -1
                                    pass
                                try:
                                    ifound = reg_seg_map.index(-iseg)
                                    if ifound >= 0:
                                        if segment_split[2] >= 0:
                                            reg_seg_map.insert(ifound,-segment_split[2])
                                            ifound += 1
                                        if segment_split[1] >= 0:
                                            reg_seg_map.insert(ifound+1,-segment_split[1])
                                except ValueError:
                                    pass
                            #
                            if first_pt == overlap_start:
                                local_seg_map.append(iseg)
                                break
                            else:
                                local_seg_map.append(-iseg)
                                break
                        else:
                            overlap_start = len(tmp_pts)
                            overlap_end = -1
                            first_pt = -1
                            for i, test_pt in enumerate(tmp_pt_map):
                                if test_pt in segment[0]:
                                    overlap_start = min(overlap_start,i)
                                    overlap_end = max(overlap_end,i)
                                    if first_pt < 0 :
                                        first_pt = segment[0].index(test_pt)
                            if overlap_start > 0:
                                self._segments.append([tmp_pt_map[:overlap_start+1], region._dx_curve, region._small_thresh])
                            if overlap_end < len(segment[0])-1:
                                self._segments.append([tmp_pt_map[overlap_end:], region._dx_curve, region._small_thresh])
                            segment[1] = min(segment[1],region._dx_curve)
                            segment[2] = min(segment[2],region._small_thresh)
                            #
                            if first_pt == overlap_start:
                                local_seg_map.append(iseg)
                                break
                            else:
                                local_seg_map.append(-iseg)
                                break
                else:
                    self._segments.append([tmp_pt_map, region._dx_curve, region._small_thresh])
                    local_seg_map.append(len(self._segments)-1)
            self._reg_seg_map.append(local_seg_map)
        # Resample segments
        pts_out = self._unique_points
        self._unique_points = numpy.array(self._unique_points)
        self._resampled_segments = []
        for segment in self._segments:
            dx = segment[1]
            pts_tmp = self._unique_points[segment[0],:]
            seg_tmp = [segment[0][0],]
            dl = numpy.zeros(len(segment[0]))
            for i in range(len(segment[0])-1):
                dl[i+1] = dl[i] + numpy.linalg.norm(pts_tmp[i+1,:]-pts_tmp[i,:])
            if int(dl[-1]/dx) >= 2:
                for dl_samp in numpy.linspace(0.0,dl[-1],int(dl[-1]/dx)+1)[1:-1]:
                    pts_out.append([
                        numpy.interp(dl_samp,dl,pts_tmp[:,0]),
                        numpy.interp(dl_samp,dl,pts_tmp[:,1])
                    ])
                    seg_tmp.append(len(pts_out)-1)
            elif dl[-1] < segment[2]:
                print("  Warning: small feature (dl={0:.2E}) detected at point {1} ({2}, {3})".format(dl[-1], segment[0][0], *pts_tmp[0,:]))
            seg_tmp.append(segment[0][-1])
            self._resampled_segments.append(seg_tmp)
        # Reindex points
        reindex = -numpy.ones((len(pts_out),),dtype=numpy.int32)
        pt_count = 0
        for segment in self._resampled_segments:
            for i in range(len(segment)-1):
                if reindex[segment[i]] < 0:
                    pt_count += 1
                    reindex[segment[i]] = pt_count-1
                if reindex[segment[i+1]] < 0:
                    pt_count += 1
                    reindex[segment[i+1]] = pt_count-1
                segment[i]=reindex[segment[i]]
            segment[-1]=reindex[segment[-1]]
        self._resampled_points = numpy.zeros((pt_count,2))
        for i in range(len(pts_out)):
            if reindex[i] >= 0:
                self._resampled_points[reindex[i],:] = pts_out[i]
        #
        print('  # of unique points    = {0}'.format(self._resampled_points.shape[0]))
        print('  # of unique segments  = {0}'.format(len(self._segments)))
        # Build list of region definitions
        self._reg_defs = extra_reg_defs
        for ireg, region in enumerate(region_list):
            imin = -1
            dmin = 1.E99
            dmax = -1.0
            reg_points = []
            for segment in self._reg_seg_map[ireg]:
                if segment < 0:
                    tmp_pt_map = [entry for entry in reversed(self._resampled_segments[-segment])]
                else:
                    tmp_pt_map = self._resampled_segments[segment]
                for ipt, unique_id in enumerate(tmp_pt_map):
                    if ipt == len(tmp_pt_map)-1:
                        continue
                    reg_pt = self._resampled_points[unique_id]
                    reg_points.append(reg_pt)
                    dmin_loc = 1.E99
                    for jpt, unique_pt in enumerate(self._resampled_points):
                        if unique_id == jpt:
                            continue
                        dmin_loc = min(dmin_loc,numpy.linalg.norm(reg_pt-unique_pt))
                    dmin = min(dmin,dmin_loc)
                    if (dmin_loc > dmax):
                        dmax = dmin_loc
                        imin = len(reg_points)-1
            region._resampled_points = numpy.asarray(reg_points)
            in_pt = region.get_in_point(imin,dmin)
            self._reg_defs.append([in_pt[0], in_pt[1], region._id, region._dx_vol*region._dx_vol/2.0])
    
    def get_mesh(self):
        '''! Generate mesh using triangle

        @result pts[np,2], tris[nc,3], regions[nc] 
        '''
        print('Generating mesh:')
        resampled_flat = []
        for segment in self._resampled_segments:
            resampled_flat += [[segment[i], segment[i+1]] for i in range(len(segment)-1)]
        alpha = dict(vertices=self._resampled_points, segments=resampled_flat, regions=self._reg_defs)
        beta = run_triangle(alpha)
        regions = beta['triangle_attributes'].astype('int32').ravel()
        print('  # of points  = {0}'.format(len(beta['vertices'])))
        print('  # of cells   = {0}'.format(len(beta['triangles'])))
        print('  # of regions = {0}'.format(regions.max()))
        return beta['vertices'], beta['triangles'], regions


class Region:
    '''! Region class for @ref oftpy.Mesh class'''
    def __init__(self,points,dx=None,dx_curve=None,angle_tol=30.0,sliver_tol=120.0,small_thresh=None,id=0,load_dict=None):
        '''! Create Region object from a closed bounding curve

        @param points List of points forming region boundary [:,2]
        @param dx Target mesh resolution inside region
        @param angle_tol Corner tolerance used when resampling curve at desired resolution
        @param sliver_tol Tolerance used for "sliver" region warnings
        @param small_thresh Tolerance used for "small" curve warnings (default: dx/2)
        @param id Region id number (default: 0)
        '''
        if load_dict is not None:
            self._points = load_dict['points']
            self._segments = load_dict['segments']
            self._id = load_dict['id']
            self._dx_curve = load_dict['dx_curve']
            self._dx_vol = load_dict['dx_vol']
        else:
            if numpy.any(points[:,0] < 0):
                raise ValueError("Point with negative radial position detected!")
            self._points = points
            self._id = id
            if dx is None:
                raise ValueError('No target mesh size set')
            else:
                self._dx_vol = dx
                if dx_curve is not None:
                    self._dx_curve = dx_curve
                else:
                    self._dx_curve = dx
            if numpy.linalg.norm(self._points[0,:]-self._points[-1,:]) < dx/1.E3:
                self._points = self._points[:-1,:]
            nv = self._points.shape[0]
            # Detect corner points and split contour
            keep_tol = numpy.cos(numpy.pi*angle_tol/180)
            sliver_tol = numpy.cos(numpy.pi*sliver_tol/180)
            keep_points = [0]
            for i in range(nv):
                if (i == nv-1):
                    tang = self._points[0,:] - self._points[i,:]
                else:
                    tang = self._points[i+1,:] - self._points[i,:]
                tang_norm = numpy.linalg.norm(tang)
                if tang_norm < self._dx_curve/1.E3:
                    print("  Warning: repeated points detected at point {0} ({1}, {2})".format(i, *self._points[i,:]))
                    continue
                tang /= tang_norm
                if i > 0:
                    angle = numpy.dot(tang,tangp)
                    if angle < keep_tol:
                        keep_points.append(i)
                        if angle < sliver_tol:
                            print("  Warning: sliver (angle={0:.1F}) detected at point {1} ({2}, {3})".format(180.-numpy.arccos(angle)*180./numpy.pi, i, *self._points[i,:]))
                tangp = tang
            keep_points.append(nv+1)
            # Index segments
            k=1
            self._segments = []
            seg_tmp = []
            for i in range(nv):
                if i >= keep_points[k]:
                    self._segments.append(seg_tmp + [i,])
                    seg_tmp = [i,]
                    k += 1
                else:
                    seg_tmp.append(i)
            self._segments.append(seg_tmp + [0,])
        # Get small feature threshold
        if small_thresh is None:
            small_thresh = self._dx_curve/2.0
        self._small_thresh = small_thresh
        self._resampled_points = None
    
    def get_resampled_points(self):
        '''! Get resampled points for bounding curve

        @result Point list [np,2]
        '''
        # Resample curves to approximate target dx
        if self._resampled_points is None:
            pts_out = []
            for segment in self._segments:
                pts_tmp = self._points[segment,:]
                dl = numpy.zeros(len(segment))
                for i in range(len(segment)-1):
                    dl[i+1] = dl[i] + numpy.linalg.norm(pts_tmp[i+1,:]-pts_tmp[i,:])
                pts_out.append(pts_tmp[0,:])
                if int(dl[-1]/self._dx_curve) >= 2:
                    for dl_samp in numpy.linspace(0.0,dl[-1],int(dl[-1]/self._dx_curve)+1)[1:-1]:
                        pts_out.append([
                            numpy.interp(dl_samp,dl,pts_tmp[:,0]),
                            numpy.interp(dl_samp,dl,pts_tmp[:,1])
                        ])
                elif dl[-1] < self._small_thresh:
                    print("  Warning: small feature (dl={0:.2E}) detected at point {1} ({2}, {3})".format(dl[-1], segment[0], *pts_tmp[0,:]))
            self._resampled_points = numpy.asarray(pts_out)
        return self._resampled_points

    def check_in_poly(self,pt):
        '''! Check if point is inside region (polygon-based approach)

        @param pt Point to check [2]
        @result Is pt inside the region?
        '''
        ncuts = 0
        for j in range(self._resampled_points.shape[0]-1):
            if (pt[1]-self._resampled_points[j,1])*(pt[1]-self._resampled_points[j+1,1]) <= 0.0:
                if (pt[1]-self._resampled_points[j,1]) == 0.0:
                    continue
                xInter = (self._resampled_points[j+1,0]-self._resampled_points[j,0])* \
                    (pt[1]-self._resampled_points[j,1])/(self._resampled_points[j+1,1]-self._resampled_points[j,1]) \
                    + self._resampled_points[j,0]
                if pt[0] <= xInter:
                    ncuts +=1
        if (pt[1]-self._resampled_points[-1,1])*(pt[1]-self._resampled_points[0,1]) <= 0.0:
            if (pt[1]-self._resampled_points[-1,1]) == 0.0:
                return (ncuts % 2 == 1)
            xInter = (self._resampled_points[0,0]-self._resampled_points[-1,0])* \
                (pt[1]-self._resampled_points[-1,1])/(self._resampled_points[0,1]-self._resampled_points[-1,1]) \
                + self._resampled_points[-1,0]
            if pt[0] <= xInter:
                ncuts +=1
        return (ncuts % 2 == 1)

    def get_in_point(self,i,dx):
        '''! Get a suitable point for defining the "inside" of the region for triangle 

        @param i Index of desired nearest boundary point
        @param dx Offset distance from bounding curve
        @result Point used to define region [2]
        '''
        dx = min(dx,self._dx_curve)/4.0
        if i==self._resampled_points.shape[0]-1:
            that2 = self._resampled_points[0,:] - self._resampled_points[i,:]
        else:
            that2 = self._resampled_points[i+1,:] - self._resampled_points[i,:]
        that1 = self._resampled_points[i,:] - self._resampled_points[i-1,:]
        that1 /= numpy.linalg.norm(that1)
        that2 /= numpy.linalg.norm(that2)
        # nhat = dx*numpy.r_[-that1[1], that1[0]]
        # dx2 = dx-(nhat[1]*that2[0]-nhat[0]*that2[1])
        # nhat += dx2*numpy.r_[-that2[1],that2[0]]
        nhat = dx*(numpy.r_[-that1[1], that1[0]] + numpy.r_[-that2[1], that2[0]])
        pt_out = self._resampled_points[i,:] + nhat
        # Check if inside contour
        if not self.check_in_poly(pt_out):
            pt_out = self._resampled_points[i,:] - nhat
        return pt_out
    
    def get_segments(self):
        segments = []
        for i in range(len(self._segments)):
            segments.append(self._points[self._segments[i],:])
        return segments
    
    def plot_segments(self,fig,ax,linewidth=None):
        '''! Plot boundary curve
        
        @param fig Figure to add curves to
        @param ax Axis to add curves to
        @param linewidth Line width for plots
        '''
        for i in range(len(self._segments)):
            ax.plot(self._points[self._segments[i],0],self._points[self._segments[i],1],linewidth=linewidth)
    
    def get_json(self):
        return {
            'points': self._points,
            'segments': self._segments,
            'id': self._id,
            'dx_curve': self._dx_curve,
            'dx_vol': self._dx_vol
        }