#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Core definitions for ThinCurr thin-wall E-M functionality

@authors Chris Hansen
@date March 2024
@ingroup doxy_oft_python
'''
import ctypes
import numpy
import h5py
from ._interface import *
from ..io import build_XDMF


class ThinCurr():
    '''! ThinCurr thin-wall E-M model class'''
    def __init__(self,OFT_env):
        '''! Initialize ThinCurr object

        @param OFT_env OFT runtime environment object (See @ref OpenFUSIONToolkit._core.OFT_env "OFT_env")
        '''
        # Create OFT execution environment
        self._oft_env = OFT_env
        self._oft_env.oft_in_groups['mesh_options'] = {'cad_type': "0"}
        self._oft_env.oft_in_groups['thincurr_hodlr_options'] = {
            'target_size': '1200',
            'aca_min_its': '20',
            'L_svd_tol': '1.E-8',
            'L_aca_rel_tol': '0.05',
            'B_svd_tol': '1.E-3',
            'B_aca_rel_tol': '0.05',
        }
        self._oft_env.update_oft_in()
        ## Thin-wall model object
        self.tw_obj = c_void_p()
        ## Number of regions in mesh
        self.nregs = -1
        ## Number of points in mesh
        self.np = -1
        ## Number of edges in mesh
        self.ne = -1
        ## Number of cells in mesh
        self.nc = -1
        ## Number of vertices that are active in the thin-wall model
        self.np_active = -1
        ## Number of hole elements in the thin-wall model
        self.nholes = -1
        ## Number of V-coil elements in the thin-wall model
        self.n_vcoils = -1
        ## Total number of active DOFs in the thin-wall model
        self.nelems = -1
        ## Number of V-coil elements in the thin-wall model
        self.n_icoils = -1
        ## Model inductance matrix (dense)
        self.Lmat = None
        ## Model inductance matrix (HODLR)
        self.Lmat_hodlr = c_void_p()
        ## Model resistance matrix
        self.Rmat = None
        ## Mesh vertices [np,3] (last column should be all zeros)
        self.r = None
        ## Mesh triangles [nc,3]
        self.lc = None
        ## Mesh regions [nc]
        self.reg = None
        ## Pointer to XML element in Fortran (FoX)
        self._xml_ptr = c_void_p()
        ## I/O basepath for plotting/XDMF output
        self._io_basepath = "."

    def setup_model(self,r=None,lc=None,reg=None,mesh_file=None,pmap=None,xml_filename=None,jumper_start=0):
        '''! Setup ThinCurr model

        @param r Point list `(np,3)`
        @param lc Cell list `(nc,3)`
        @param reg Region tag `(nc,)`
        @param mesh_file File containing model in native mesh format
        @param pmap Point map for periodic grids
        @param xml_filename Path to XML file for model
        @param jumper_start Index of first jumper nodeset in meshfile (positive values Fortran style, negative values Python style)
        '''
        if self.nregs != -1:
            raise ValueError('Mesh already setup, delete or create new instance for new model')
        if xml_filename is not None:
            cfilename = self._oft_env.path2c(xml_filename)
            oftpy_load_xml(cfilename,ctypes.byref(self._xml_ptr))
        if mesh_file is not None:
            if (r is not None) or (lc is not None) or (reg is not None):
                raise ValueError('Specification of "mesh_file" is incompatible with specification of "r", "lc", and "reg"')
            idummy = c_int(-1)
            rfake = numpy.ones((1,1),dtype=numpy.float64)
            lcfake = numpy.ones((1,1),dtype=numpy.int32)
            regfake = numpy.ones((1,),dtype=numpy.int32)
            if pmap is None:
                pmap = -numpy.ones((1,),dtype=numpy.int32)
            else:
                pmap = numpy.ascontiguousarray(pmap, dtype=numpy.int32)
            sizes = numpy.zeros((9,),dtype=numpy.int32)
            cfilename = self._oft_env.path2c(mesh_file)
            error_string = self._oft_env.get_c_errorbuff()
            thincurr_setup(cfilename,idummy,rfake,idummy,lcfake,regfake,pmap,c_int(jumper_start),ctypes.byref(self.tw_obj),sizes,error_string,self._xml_ptr)
            if error_string.value != b'':
                raise Exception(error_string.value.decode())
        elif r is not None:
            if lc is None:
                raise ValueError('"r" and "lc" must be both be specified')
            if jumper_start != 0:
                raise ValueError('"jumper_start" not supported with manual mesh specification')
            np = c_int(r.shape[0])
            nc = c_int(lc.shape[0])
            r = numpy.ascontiguousarray(r, dtype=numpy.float64)
            lc = numpy.ascontiguousarray(lc, dtype=numpy.int32)
            if reg is None:
                reg = numpy.ones((nc.value,),dtype=numpy.int32)
            else:
                reg = numpy.ascontiguousarray(reg, dtype=numpy.int32)
            if pmap is None:
                pmap = -numpy.ones((1,),dtype=numpy.int32)
            else:
                # pmap = numpy.ascontiguousarray(pmap, dtype=numpy.int32)
                raise ValueError('"pmap" not supported with manual mesh specification')
            sizes = numpy.zeros((9,),dtype=numpy.int32)
            cfilename = self._oft_env.path2c('')
            error_string = self._oft_env.get_c_errorbuff()
            thincurr_setup(cfilename,np,r,nc,lc+1,reg,pmap,c_int(jumper_start),ctypes.byref(self.tw_obj),sizes,error_string,self._xml_ptr)
            if error_string.value != b'':
                raise Exception(error_string.value.decode())
        else:
            raise ValueError('Mesh filename (native format) or mesh values (r, lc) required')
        self.np = sizes[0]
        self.ne = sizes[1]
        self.nc = sizes[2]
        self.nregs = sizes[3]
        self.np_active = sizes[4]
        self.nholes = sizes[5]
        self.n_vcoils = sizes[6]
        self.nelems = sizes[7]
        self.n_icoils = sizes[8]
    
    def setup_io(self,basepath=None,save_debug=False,legacy_hdf5=False):
        '''! Setup XDMF+HDF5 I/O for 3D visualization

        @param basepath Path to root directory to use for I/O
        @param save_debug Save model debug information?
        @param legacy_hdf5 Use legacy HDF5 format (required for VisIt)
        '''
        # tw_ptr,basepath,save_debug,error_str
        if basepath is None:
            basepath_c = self._oft_env.path2c('')
            self._io_basepath = "."
        else:
            if basepath[-1] != '/':
                basepath += '/'
            self._io_basepath = basepath[:-1]
            basepath_c = self._oft_env.path2c(basepath)
        error_string = self._oft_env.get_c_errorbuff()
        thincurr_setup_io(self.tw_obj,basepath_c,c_bool(save_debug),c_bool(legacy_hdf5),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
    
    def reconstruct_current(self,potential,centering='cell'):
        '''! Reconstruct current field on mesh

        @param potential Current potential
        @param centering Desired field centering ('cell', 'vertex')
        @result Current field on selected centering `(:,3)`
        '''
        if centering == 'cell':
            curr = numpy.zeros((self.nc,3), dtype=numpy.float64)
            cent_key = 1
        elif centering == 'vertex':
            curr = numpy.zeros((self.np,3), dtype=numpy.float64)
            cent_key = 2
        else:
            raise ValueError('Unknown centering, should be "cell" or "vertex"')
        potential = numpy.ascontiguousarray(potential, dtype=numpy.float64)
        thincurr_recon_curr(self.tw_obj,potential,curr,cent_key)
        return curr/mu0
    
    def reconstruct_Bfield(self,potential,coil_currs=None):
        '''! Reconstruct magnetic field on original grid

        @param weights Reduced model basis weights
        @param coil_currs Coil currents [A]
        @result Magnetic field on original grid `(:,3)`
        '''
        potential = numpy.ascontiguousarray(potential, dtype=numpy.float64)
        coil_currs = numpy.ascontiguousarray(coil_currs, dtype=numpy.float64)
        field = numpy.zeros((self.np,3), dtype=numpy.float64)
        if self.Lmat_hodlr:
            thincurr_recon_field(self.tw_obj,potential,coil_currs,field,self.Lmat_hodlr)
        else:
            thincurr_recon_field(self.tw_obj,potential,coil_currs,field,c_void_p())
        return field/mu0
    
    def save_current(self,potential,tag):
        '''! Save current field from ThinCurr to plot files

        @param potential Current potential to save
        @param tag Name of field in plot files
        '''
        if potential.shape[0] != self.nelems:
            raise IndexError('Incorrect shape of "potential", should be [nelems]')
        potential = numpy.ascontiguousarray(potential, dtype=numpy.float64)
        ctag = self._oft_env.path2c(tag)
        thincurr_save_field(self.tw_obj,potential,ctag)
    
    def save_scalar(self,field,tag):
        '''! Save scalar field to plot files

        @param field Pointwise data to save
        @param tag Name of field in plot files
        '''
        if field.shape[0] != self.np:
            raise IndexError('Incorrect shape of "field", should be [np]')
        field = numpy.ascontiguousarray(field, dtype=numpy.float64)
        ctag = self._oft_env.path2c(tag)
        thincurr_save_scalar(self.tw_obj,field,ctag)
    
    def build_XDMF(self,repeat_static=False,pretty=False):
        '''! Build XDMF plot metadata files for model

        @param repeat_static Repeat static fields (0-th timestep) in all timesteps?
        @param pretty Use pretty printing (indentation) in XDMF files?
        '''
        return build_XDMF(path=self._io_basepath,repeat_static=repeat_static,pretty=pretty)
    
    def scale_va(self,data,div_flag=False):
        '''! Scale a vertex array by vertex areas (eg. B_n -> flux)

        @param data Data to scale
        @param div_flag Divide by vertex areas instead?
        @result Scaled data `(:)`
        '''
        data_in = numpy.ascontiguousarray(data.copy(), dtype=numpy.float64)
        thincurr_scale_va(self.tw_obj,data_in,div_flag)
        return data_in
    
    def compute_Lmat(self,cache_file=None,use_hodlr=False):
        '''! Compute the self-inductance matrix for this model

        @param cache_file Path to cache file to store/load matrix
        @param use_hodlr Use HODLR compression for the matrix
        @result Self-inductance matrix (`(:,:)` if `use_hodlr=False` else reference to HODLR object)
        '''
        if cache_file is None:
            cache_string = self._oft_env.path2c("")
        else:
            cache_string = self._oft_env.path2c(cache_file)
        Lmat_loc = c_void_p()
        error_string = self._oft_env.get_c_errorbuff()
        thincurr_Lmat(self.tw_obj,use_hodlr,ctypes.byref(Lmat_loc),cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        if use_hodlr:
            self.Lmat_hodlr = Lmat_loc
        else:
            self.Lmat = numpy.ctypeslib.as_array(ctypes.cast(Lmat_loc, c_double_ptr),shape=(self.nelems,self.nelems))
    
    def compute_Bmat(self,cache_file=None):
        '''! Compute magnetic field reconstruction operators for this model

        @param cache_file Path to cache file to store/load matrix
        @result Element B-field reconstruction matrix (`(:,:)` if HODLR is available else `None`)
        @result Icoil B-field reconstruction matrix `(:,:)`
        '''
        if cache_file is None:
            cache_string = self._oft_env.path2c("")
        else:
            cache_string = self._oft_env.path2c(cache_file)
        Bmat_loc = c_void_p()
        Bdr_ptr = c_void_p()
        error_string = self._oft_env.get_c_errorbuff()
        if self.Lmat_hodlr:
            thincurr_Bmat(self.tw_obj,self.Lmat_hodlr,ctypes.byref(Bmat_loc),ctypes.byref(Bdr_ptr),cache_string,error_string)
        else:
            thincurr_Bmat(self.tw_obj,c_void_p(),ctypes.byref(Bmat_loc),ctypes.byref(Bdr_ptr),cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        if self.Lmat_hodlr:
            return None, numpy.ctypeslib.as_array(ctypes.cast(Bdr_ptr, c_double_ptr),shape=(3,self.n_icoils,self.np))
        else:
            return numpy.ctypeslib.as_array(ctypes.cast(Bmat_loc, c_double_ptr),shape=(3,self.nelems,self.np)), \
                numpy.ctypeslib.as_array(ctypes.cast(Bdr_ptr, c_double_ptr),shape=(3,self.n_icoils,self.np))
    
    def compute_Mcoil(self,cache_file=None):
        '''! Compute the mutual inductance between passive (mesh+vcoils) and active elements (icoils)

        @param cache_file Path to cache file to store/load matrix
        @result Mutual inductance matrix `(:,:)`
        '''
        if cache_file is None:
            cache_string = self._oft_env.path2c("")
        else:
            cache_string = self._oft_env.path2c(cache_file)
        Mc_loc = c_void_p()
        error_string = self._oft_env.get_c_errorbuff()
        thincurr_Mcoil(self.tw_obj,ctypes.byref(Mc_loc),cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return numpy.ctypeslib.as_array(ctypes.cast(Mc_loc, c_double_ptr),shape=(self.n_icoils,self.nelems))

    def compute_Msensor(self,sensor_file=None,cache_file=None):
        '''! Compute the mutual inductance between model and sensors

        @param sensor_file Path to file contatining flux loop definitions
        @param cache_file Path to cache file to store/load matrix
        @result Mutual inductance matrix between model and sensors `(:,:)`
        @result Mutual inductance matrix between icoils and sensors `(:,:)`
        @result Internal ThinCurr sensor object
        '''
        if cache_file is None:
            cache_string = self._oft_env.path2c("")
        else:
            cache_string = self._oft_env.path2c(cache_file)
        if sensor_file is None:
            sensor_string = self._oft_env.path2c("none")
        else:
            sensor_string = self._oft_env.path2c(sensor_file)
        Ms_loc = c_void_p()
        Msc_loc = c_void_p()
        nsensors = c_int()
        njumpers = c_int()
        sensor_loc = c_void_p()
        error_string = self._oft_env.get_c_errorbuff()
        thincurr_Msensor(self.tw_obj,sensor_string,ctypes.byref(Ms_loc),ctypes.byref(Msc_loc), 
                         ctypes.byref(nsensors),ctypes.byref(njumpers),ctypes.byref(sensor_loc),cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        sensor_names = []
        for i in range(nsensors.value):
            sensor_name = create_string_buffer(b"",40)
            error_string = create_string_buffer(b"",200)
            thincurr_get_sensor_name(sensor_loc,c_int(i+1),sensor_name,error_string)
            if error_string.value != b'':
                raise Exception(error_string.value.decode())
            sensor_names.append(sensor_name.value.decode().strip())
        return numpy.ctypeslib.as_array(ctypes.cast(Ms_loc, c_double_ptr),shape=(self.nelems,nsensors.value)), \
               numpy.ctypeslib.as_array(ctypes.cast(Msc_loc, c_double_ptr),shape=(self.n_icoils,nsensors.value)), \
               {'names': sensor_names, 'ptr': sensor_loc}
    
    def get_eta_values(self):
        '''! Get resistivity values for model

        @returns `eta_values` Resistivity values for model [nregs]
        '''
        eta_values = numpy.zeros((self.nregs,), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        thincurr_get_eta(self.tw_obj,eta_values,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return eta_values
    
    def set_eta_values(self,eta_values=None):
        '''! Set resistivity values for model (overrides those in XML if specified)

        @param eta_values New resistivity values for model [nregs]
        '''
        if eta_values.shape[0] != self.nregs:
            raise IndexError('Incorrect shape of "eta_values", should be [nregs]')
        eta_values = numpy.ascontiguousarray(eta_values, dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        thincurr_set_eta(self.tw_obj,eta_values,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())

    def compute_Rmat(self,copy_out=False):
        '''! Compute the resistance matrix for this model

        @param copy_out Copy matrix to python and store in `self.Rmat`?
        '''
        if copy_out:
            self.Rmat = numpy.zeros((self.nelems,self.nelems), dtype=numpy.float64)
            Rmat_tmp = self.Rmat
        else:
            Rmat_tmp = numpy.zeros((1,1), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        thincurr_curr_Rmat(self.tw_obj,c_bool(copy_out),Rmat_tmp,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
    
    def cross_coupling(self,model2,cache_file=None):
        '''! Compute the mutual inductance between this and another ThinCurr model

        @param model2 The second model for mutual calculation
        @param cache_file Path to cache file to store/load matrix
        @result Mutual inductance matrix `(:,:)`
        '''
        Mmat = numpy.zeros((self.nelems,model2.nelems), dtype=numpy.float64)
        if cache_file is None:
            cache_string = self._oft_env.path2c("")
        else:
            cache_string = self._oft_env.path2c(cache_file)
        error_string = self._oft_env.get_c_errorbuff()
        thincurr_cross_coupling(self.tw_obj,model2.tw_obj,Mmat,cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return Mmat
    
    def cross_eval(self,model2,field):
        '''! Compute the voltage/flux induced on another ThinCurr model from a current structure on this model

        @param model2 The second model for mutual calculation
        @param field One or more current fields
        @result Flux on `model2` from `field` on `self` `(field.shape[0],:)`
        '''
        nrhs = field.shape[0]
        if field.shape[1] != self.nelems:
            raise IndexError('Incorrect shape of "field", should be [nelems]')
        vec_out = numpy.zeros((nrhs,model2.nelems), dtype=numpy.float64)
        vec_in = numpy.ascontiguousarray(field.copy(), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        thincurr_cross_eval(self.tw_obj,model2.tw_obj,c_int(nrhs),vec_in,vec_out,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return vec_out
    
    def get_regmat(self):
        '''! Compute the current regularization matrix for this model

        @result Regularization matrix `(:,:)`
        '''
        Rmat = numpy.zeros((self.nelems,3*self.nc), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        thincurr_curr_regmat(self.tw_obj,Rmat,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return Rmat
    
    def get_eigs(self,neigs,direct=False):
        '''! Compute eigenmodes of this model

        @param neigs Number of eigenvalues to compute
        @param direct Use direct solver?
        @result Eigenvalues `(neigs)`
        @result Eigenvectors `(neigs,:)`
        '''
        eig_vals = numpy.zeros((neigs,), dtype=numpy.float64)
        eig_vecs = numpy.zeros((neigs,self.nelems), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        if self.Lmat_hodlr:
            thincurr_eigenvalues(self.tw_obj,c_bool(direct),c_int(neigs),eig_vals,eig_vecs,self.Lmat_hodlr,error_string)
        else:
            thincurr_eigenvalues(self.tw_obj,c_bool(direct),c_int(neigs),eig_vals,eig_vecs,c_void_p(),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return eig_vals, eig_vecs

    def compute_freq_response(self,fdriver=None,vdriver=None,freq=0.0,fr_limit=0,direct=False):
        '''! Compute frequence response of the current model to given driver voltages

        @param fdriver Real/Imaginary driver flux pair (M*I)
        @param vdriver Real/Imaginary driver voltage pair
        @param freq Frequency for response calculation [Hz] (unused if `fr_limit!=0`)
        @param fr_limit Frequency limit for calculation (0: none, 1: inductive, 2: resistive)
        @param direct Use direct solver?
        @result Real/Imaginary eddy current response
        '''
        if fdriver is not None:
            if vdriver is not None:
                raise ValueError('"fdriver" and "vdriver" cannot be specified simultaneously')
            # Assume fixed current V = -i*omega*M*I
            if fr_limit == 0:
                omega = 2.0*numpy.pi*freq
            else:
                omega = 1.0
            vdriver = fdriver.copy()
            vdriver[0,:] = omega*fdriver[1,:]
            vdriver[1,:] = -omega*fdriver[0,:]
        result = numpy.ascontiguousarray(vdriver.copy(), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        if self.Lmat_hodlr:
            thincurr_freq_response(self.tw_obj,c_bool(direct),c_int(fr_limit),c_double(freq),result,self.Lmat_hodlr,error_string)
        else:
            thincurr_freq_response(self.tw_obj,c_bool(direct),c_int(fr_limit),c_double(freq),result,c_void_p(),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return result
    
    def run_td(self,dt,nsteps,coil_currs=None,coil_volts=None,full_volts=None,direct=False,
               status_freq=10,plot_freq=10,sensor_obj=None,sensor_values=None,lin_tol=1.E-6,timestep_cn=True):
        '''! Perform a time-domain simulation

        @param dt Time step for simulation
        @param nsteps Number of steps to take
        @param coil_currs Current vs time array for Icoils `(:,n_icoils+1)` (first column is time)
        @param coil_volts Voltage vs time array for Vcoils `(:,n_vcoils+1)` (first column is time)
        @param full_volts Voltage vs time array for Vcoils `(:,nelems+1)` (first column is time)
        @param direct Use direct solver?
        @param status_freq Frequency to print status information
        @param plot_freq Frequency to save plot files
        @param sensor_obj Sensor object to use
        @param lin_tol Tolerance for linear solver when `direct=False`
        @param timestep_cn Use Crank-Nicolson timestep?
        '''
        vec_ic = numpy.zeros((self.nelems,), dtype=numpy.float64)
        if coil_currs is None:
            ncurr = c_int(0)
            coil_currs = numpy.zeros((1,1), dtype=numpy.float64)
        else:
            if coil_currs.shape[1]-1 != self.n_icoils:
                raise ValueError("# of currents in waveform does not match # of icoils")
            ncurr = c_int(coil_currs.shape[0])
            coil_currs = numpy.ascontiguousarray(coil_currs.transpose(), dtype=numpy.float64)
        volts_full = c_bool(False)
        sensor_ptr = c_void_p()
        if sensor_obj is not None:
            sensor_ptr = sensor_obj['ptr']
        sensor_values_ptr = c_void_p()
        if coil_volts is None:
            if full_volts is not None:
                if full_volts.shape[1]-1 != self.nelems:
                    raise ValueError("# of voltages in waveform does not match # of elements")
                nvolt = c_int(full_volts.shape[0])
                coil_volts = numpy.ascontiguousarray(full_volts.transpose(), dtype=numpy.float64)
                volts_full = c_bool(True)
                if sensor_values is not None:
                    if sensor_values.shape[0] != full_volts.shape[0]:
                        raise ValueError("# of timepoints in sensor waveform must match voltage waveform")
                    if sensor_values.shape[1]-1 != len(sensor_obj['names']):
                        raise ValueError("# of values in sensor waveform does not match number of sensors")
                    sensor_values = numpy.ascontiguousarray(sensor_values.transpose(), dtype=numpy.float64)
                    sensor_values_ptr = sensor_values.ctypes.data_as(c_double_ptr)
            else:
                nvolt = c_int(0)
                coil_volts = numpy.zeros((1,1), dtype=numpy.float64)
        else:
            if full_volts is not None:
                raise ValueError('"full_volts" and "coil_volts" cannot be used simultaneously, include coil voltages in "full_volts"')
            if coil_volts.shape[1]-1 != self.n_vcoils:
                raise ValueError("# of voltages in waveform does not match # of vcoils")
            nvolt = c_int(coil_volts.shape[0])
            coil_volts = numpy.ascontiguousarray(coil_volts.transpose(), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        if self.Lmat_hodlr:
            thincurr_time_domain(self.tw_obj,c_bool(direct),c_double(dt),c_int(nsteps),c_double(lin_tol),c_bool(timestep_cn),
                                 c_int(status_freq),c_int(plot_freq),vec_ic,sensor_ptr,ncurr,coil_currs,nvolt,coil_volts,volts_full,
                                 sensor_values_ptr,self.Lmat_hodlr,error_string)
        else:
            thincurr_time_domain(self.tw_obj,c_bool(direct),c_double(dt),c_int(nsteps),c_double(lin_tol),c_bool(timestep_cn),
                                 c_int(status_freq),c_int(plot_freq),vec_ic,sensor_ptr,ncurr,coil_currs,nvolt,coil_volts,volts_full,
                                 sensor_values_ptr,c_void_p(),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
    
    def plot_td(self,nsteps,compute_B=False,rebuild_sensors=False,plot_freq=10,sensor_obj=None,sensor_values=None):
        '''! Perform a time-domain simulation

        @param nsteps Number of steps to take
        @param compute_B Compute B-field on grid vertices
        @param rebuild_sensors Recompute sensor signals (overwriting if present)
        @param plot_freq Frequency to load plot files
        @param sensor_obj Sensor object to use
        '''
        sensor_ptr = c_void_p()
        if sensor_obj is not None:
            sensor_ptr = sensor_obj['ptr']
            if sensor_values is not None:
                if sensor_values.shape[1]-1 != len(sensor_obj['names']):
                    raise ValueError("# of values in sensor waveform does not match number of sensors")
                nsensor = c_int(sensor_values.shape[0])
                sensor_values = numpy.ascontiguousarray(sensor_values.transpose(), dtype=numpy.float64)
        if sensor_values is None:
            nsensor = c_int(0)
            sensor_values = numpy.zeros((1,1), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        if self.Lmat_hodlr:
            thincurr_time_domain_plot(self.tw_obj,c_bool(compute_B),c_bool(rebuild_sensors),c_int(nsteps),c_int(plot_freq),sensor_ptr,
                                      sensor_values,nsensor,self.Lmat_hodlr,error_string)
        else:
            thincurr_time_domain_plot(self.tw_obj,c_bool(compute_B),c_bool(rebuild_sensors),c_int(nsteps),c_int(plot_freq),sensor_ptr,
                                      sensor_values,nsensor,c_void_p(),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())

    def build_reduced_model(self,basis_set,filename='tCurr_reduced.h5',compute_B=False,sensor_obj=None):
        r'''! Build reduced model by projecting full model onto defined basis set of currents

        @param basis_set Basis set for projection `(nBasis,:)`
        @param filename Filename for saving reduced model
        @param compute_B Compute B-field reconstruction operators for reduced model?
        @param sensor_obj Sensor object to use
        @result Reduced model (see \ref ThinCurr_reduced)
        '''
        basis_set = numpy.ascontiguousarray(basis_set, dtype=numpy.float64)
        nbasis = c_int(basis_set.shape[0])
        sensor_ptr = c_void_p()
        if sensor_obj is not None:
            sensor_ptr = sensor_obj['ptr']
        cfilename = self._oft_env.path2c(filename)
        error_string = self._oft_env.get_c_errorbuff()
        if self.Lmat_hodlr:
            thincurr_reduce_model(self.tw_obj,cfilename,nbasis,basis_set,c_bool(compute_B),sensor_ptr,self.Lmat_hodlr,error_string)
        else:
            thincurr_reduce_model(self.tw_obj,cfilename,nbasis,basis_set,c_bool(compute_B),sensor_ptr,c_void_p(),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return ThinCurr_reduced(filename)


class ThinCurr_reduced:
    '''! Reduced ThinCurr thin-wall eddy current model class'''
    def __init__(self, filename):
        '''! Initialize Reduced ThinCurr object

        @param filename File containing reduced model
        '''
        with h5py.File(filename,'r') as file:
            mu0_scale = mu0 # Prior to addition of version flag magnetic units were saved for coils
            if 'ThinCurr_Version' in file:
                mu0_scale = 1.0 
            ## Current potential basis set
            self.Basis = numpy.asarray(file['Basis'])
            ## Self-inductance matrix for reduced model
            self.L = numpy.asarray(file['L'])
            ## Resistance matrix for reduced model
            self.R = numpy.asarray(file['R'])
            ## B-field reconstruction operator
            self.B = None
            if 'Bx' in file:
                self.B = [numpy.asarray(file['Bx']), numpy.asarray(file['By']), numpy.asarray(file['Bz'])]
            ## Model-sensor mutual inductance matrix 
            self.Ms = None
            if 'Ms' in file:
                self.Ms = numpy.asarray(file['Ms'])
            ## Model-coil mutual inductance matrix
            self.Mc = None
            ## B-field coil reconstruction operator
            self.Bc = None
            if 'Mc' in file:
                self.Mc = mu0_scale*numpy.asarray(file['Mc'])
                if 'Bx_c' in file:
                    self.Bc = mu0_scale*numpy.asarray([file['Bx_c'], file['By_c'], file['Bz_c']])
            ## Coil-sensor mutual inductance matrix
            self.Msc = None
            if 'Msc' in file:
                self.Msc = mu0_scale*numpy.asarray(file['Msc'])
    
    def reconstruct_potential(self,weights):
        r'''! Reconstruct full current potential on original grid

        @note To reconstruct the full current pass the output to @ref ThinCurr.reconstruct_current "reconstruct_current()"
        with the original model

        @param weights Reduced model basis weights
        @result Full current potential on original grid `(:)`
        '''
        return numpy.dot(weights,self.Basis)
    
    def reconstruct_Bfield(self,weights,coil_currs=None):
        '''! Reconstruct magnetic field on original grid

        @param weights Reduced model basis weights
        @param coil_currs Coil currents [A]
        @result Magnetic field on original grid `(:,3)`
        '''
        if (self.B is None):
            raise ValueError('Magnetic field reconstruction operator not part of this model')
        if coil_currs is not None:
            if (self.Bc is None):
                raise ValueError('Coil magnetic field reconstruction operator not part of this model')
            if coil_currs.shape[0] != self.Bc[0].shape[0]:
                raise ValueError('Size of "coil_currs" does not match number of coils in reduced model')
        Bfield = numpy.zeros((self.B[0].shape[1],3))
        for i in range(3):
            Bfield[:,i] = numpy.dot(weights,self.B[i])
            if self.Bc is not None:
                Bfield[:,i] += numpy.dot(coil_currs,self.Bc[i])
        return Bfield

    def get_eigs(self):
        '''! Compute eigenmodes for reduced model

        @result Eigenvalues `(:)`
        @result Eigenvectors `(:,:)`
        '''
        eig_vals, eig_vecs = numpy.linalg.eig(numpy.dot(numpy.linalg.inv(self.R),self.L))
        sort_inds = (-eig_vals).argsort()
        return eig_vals[sort_inds], eig_vecs[sort_inds,:]

    def run_td(self,dt,nsteps,coil_currs,status_freq=10,plot_freq=10):
        r'''! Perform a time-domain simulation

        @param dt Time step for simulation
        @param nsteps Number of steps to take
        @param coil_currs Current vs time array for Icoils `(:,n_icoils+1)` (first index is time)
        @param status_freq Frequency to print status information
        @param plot_freq Frequency to save plot files
        @result Sensor signals dictionary
          - `time` Timebase [s] `(:)`, 
          - `sensors` Sensor signals `(:,:)`
        @result Current history dictionary
          - `time` Timebase [s] `(:)`
          - `curr` Reduced model basis weights `(:,:)`,
          - `coil_curr` Coil currents [A / \f$ \mu_0 \f$ ] `(:,:)`
        '''
        Lforward = self.L - (dt/2.0)*self.R
        Lbackward = numpy.linalg.inv(self.L + (dt/2.0)*self.R)
        #
        vec_time = []
        vec_hist = []
        coil_hist = []
        sen_time = []
        sen_hist = []
        pot_tmp = numpy.zeros((self.L.shape[0],))
        t = 0.0
        print('Timestep {0} {1:12.4E} {2:12.4E}'.format(0,t,numpy.linalg.norm(pot_tmp)))
        vec_time.append(t)
        curr_tmp = numpy.zeros((coil_currs.shape[1]-1,))
        for j in range(coil_currs.shape[1]-1):
            curr_tmp[j] = numpy.interp(t,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1],right=coil_currs[-1,j+1])
        coil_hist.append(curr_tmp)
        vec_hist.append(pot_tmp)
        if self.Ms is not None:
            sen_tmp = numpy.dot(pot_tmp,self.Ms)
            if self.Msc is not None:
                sen_tmp += numpy.dot(curr_tmp,self.Msc)
            sen_time.append(t)
            sen_hist.append(sen_tmp)
        #
        for i in range(nsteps):
            rhs = numpy.dot(pot_tmp,Lforward)
            if self.Mc is not None:
                dcurr_tmp = numpy.zeros((coil_currs.shape[1]-1,))
                for j in range(coil_currs.shape[1]-1):
                    dcurr_tmp[j] = numpy.interp(t+dt/4.0,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1],right=coil_currs[-1,j+1])
                    dcurr_tmp[j] -= numpy.interp(t-dt/4.0,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1],right=coil_currs[-1,j+1])
                    dcurr_tmp[j] += numpy.interp(t+dt*5.0/4.0,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1],right=coil_currs[-1,j+1])
                    dcurr_tmp[j] -= numpy.interp(t+dt*3.0/4.0,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1],right=coil_currs[-1,j+1])
                rhs -= numpy.dot(dcurr_tmp,self.Mc)
            pot_tmp = numpy.dot(rhs,Lbackward)
            t += dt
            if ((i+1) % status_freq) == 0:
                print('Timestep {0} {1:12.4E} {2:12.4E}'.format(i+1,t,numpy.linalg.norm(pot_tmp)))
            curr_tmp = numpy.zeros((coil_currs.shape[1]-1,))
            for j in range(coil_currs.shape[1]-1):
                curr_tmp[j] = numpy.interp(t,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1],right=coil_currs[-1,j+1])
            if ((i+1) % plot_freq) == 0:
                vec_time.append(t)
                coil_hist.append(curr_tmp)
                vec_hist.append(pot_tmp)
            if self.Ms is not None:
                sen_tmp = numpy.dot(pot_tmp,self.Ms)
                if self.Msc is not None:
                    sen_tmp += numpy.dot(curr_tmp,self.Msc)
                sen_time.append(t)
                sen_hist.append(sen_tmp)
        #
        sensor_obj = {
            'time': numpy.array(sen_time),
            'sensors': numpy.array(sen_hist)
        }
        curr_obj = {
            'time': numpy.array(vec_time),
            'curr': numpy.array(vec_hist),
            'coil_curr': numpy.array(coil_hist)
        }
        return sensor_obj, curr_obj