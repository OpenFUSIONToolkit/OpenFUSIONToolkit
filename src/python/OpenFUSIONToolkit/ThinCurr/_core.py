'''! Python interface for ThinCurr thin-wall eddy current functionality

@authors Chris Hansen
@date March 2024
@ingroup doxy_oft_python
'''

#
# Python interface for TokaMaker Grad-Shafranov functionality
import ctypes
import numpy
from ..util import *
from ._interface import *

oft_in_template = """&runtime_options
 debug={DEBUG_LEVEL}
/

&mesh_options
 meshname='none'
 {MESH_TYPE}
/

{MESH_DEF}

&thincurr_hodlr_options
 target_size=1200
 aca_min_its=20
 L_svd_tol=1.E-8
 L_aca_rel_tol=0.05
 B_svd_tol=1.E-3
 B_aca_rel_tol=0.05
/
"""

class ThinCurr():
    '''! ThinCurr thin-wall eddy current model class'''
    def __init__(self,debug_level=0,nthreads=2):
        '''! Initialize ThinCurr object

        @param debug_level Level of debug printing (0-3)
        '''
        ## Input file settings
        self._psin_dict = {'DEBUG_LEVEL': debug_level, 'MESH_TYPE': '', 'MESH_DEF': ''}
        self._update_psin()
        oft_init(c_int(nthreads))
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

    def _update_psin(self):
        '''! Update input file (`oftpyin`) with current settings'''
        with open('oftpyin','w+') as fid:
            fid.write(oft_in_template.format(**self._psin_dict))

    def setup_model(self,r=None,lc=None,reg=None,mesh_file=None,pmap=None,xml_filename=None):
        '''! Setup ThinCurr model

        @param r Point list [np,3]
        @param lc Cell list [nc,3]
        @param reg Region tag [nc,]
        @param mesh_file File containing model in native mesh format
        @param pmap Point map for periodic grids
        @param xml_filename Path to XML file for model
        '''
        if self.nregs != -1:
            raise ValueError('Mesh already setup, delete or create new instance for new model')
        if xml_filename is not None:
            filename = c_char_p(xml_filename.encode())
            oftpy_load_xml(filename,ctypes.byref(self._xml_ptr))
        nregs = c_int()
        if mesh_file is not None:
            idummy = c_int(-1)
            rfake = numpy.ones((1,1),dtype=numpy.float64)
            lcfake = numpy.ones((1,1),dtype=numpy.int32)
            regfake = numpy.ones((1,),dtype=numpy.int32)
            if pmap is None:
                pmap = -numpy.ones((1,),dtype=numpy.int32)
            else:
                pmap = numpy.ascontiguousarray(pmap, dtype=numpy.int32)
            sizes = numpy.zeros((8,),dtype=numpy.int32)
            filename = c_char_p(mesh_file.encode())
            cstring = c_char_p(b""*200)
            thincurr_setup(filename,idummy,rfake,idummy,lcfake,regfake,pmap,ctypes.byref(self.tw_obj),sizes,cstring,self._xml_ptr)
            if cstring.value != b'':
                raise Exception(cstring.value.decode())
            self.np = sizes[0]
            self.ne = sizes[1]
            self.nc = sizes[2]
            self.np_active = sizes[3]
            self.nholes = sizes[4]
            self.n_vcoils = sizes[5]
            self.nelems = sizes[6]
            self.n_icoils = sizes[7]
        elif r is not None:
            raise ValueError('Specifying mesh values not yet supported')
            # r = numpy.ascontiguousarray(r, dtype=numpy.float64)
            # lc = numpy.ascontiguousarray(lc, dtype=numpy.int32)
            # np = c_int(r.shape[0])
            # nc = c_int(lc.shape[0])
            # if reg is None:
            #     reg = numpy.ones((nc.value,),dtype=numpy.int32)
            # else:
            #     reg = numpy.ascontiguousarray(reg, dtype=numpy.int32)
            # if pmap is None:
            #     pmap = -numpy.ones((1,),dtype=numpy.int32)
            # else:
            #     pmap = numpy.ascontiguousarray(pmap, dtype=numpy.float64)
            # sizes = numpy.zeros((8,),dtype=numpy.int32)
            # filename = c_char_p(b"")
            # cstring = c_char_p(b""*200)
            # thincurr_setup(filename,np,r,nc,lc,reg,ctypes.byref(self.tw_obj),sizes,cstring)
            # if cstring.value != b'':
            #     raise Exception(cstring.value)
            # self.np = sizes[0]
            # self.ne = sizes[1]
            # self.nc = sizes[2]
            # self.np_active = sizes[3]
            # self.nholes = sizes[4]
            # self.n_vcoils = sizes[5]
            # self.nelems = sizes[6]
            # self.n_icoils = sizes[7]
        else:
            raise ValueError('Mesh filename (native format) or mesh values required')
        self.nregs = nregs.value
    
    def setup_io(self,basepath=None,save_debug=False):
        '''! Setup XDMF+HDF5 I/O for 3D visualization

        @param basepath Path to root directory to use for I/O
        @param save_debug Save model debug information?
        '''
        # tw_ptr,basepath,save_debug,error_str
        if basepath is None:
            basepath_c = c_char_p(b'')
            self._io_basepath = "."
        else:
            if basepath[-1] != '/':
                basepath += '/'
            self._io_basepath = basepath[:-1]
            basepath_c = c_char_p(basepath.encode())
        cstring = c_char_p(b""*200)
        thincurr_setup_io(self.tw_obj,basepath_c,c_bool(save_debug),cstring)
        if cstring.value != b'':
            raise Exception(cstring.value.decode())
    
    def save_current(self,potential,tag):
        '''! Save current field from ThinCurr to plot files

        @param field Pointwise data to save
        @param tag Name of field in plot files
        '''
        # if potential.shape[0] != self.np:
        #     raise ValueError('Incorrect shape of "psi", should be [np]')
        potential = numpy.ascontiguousarray(potential, dtype=numpy.float64)
        cstring = c_char_p(tag.encode())
        thincurr_save_field(self.tw_obj,potential,cstring)
    
    def save_scalar(self,field,tag):
        '''! Save scalar field to plot files

        @param field Pointwise data to save
        @param tag Name of field in plot files
        '''
        # if field.shape[0] >= self.np:
        #     raise ValueError('Incorrect shape of "psi", should be [np]')
        field = numpy.ascontiguousarray(field, dtype=numpy.float64)
        ctag = c_char_p(tag.encode())
        thincurr_save_scalar(self.tw_obj,field,ctag)
    
    def build_XDMF(self,repeat_static=False,pretty=False):
        '''! Build XDMF plot metadata files for model

        @param repeat_static Repeat static fields (0-th timestep) in all timesteps?
        @param pretty Use pretty printing (indentation) in XDMF files?
        '''
        build_XDMF(path=self._io_basepath,repeat_static=repeat_static,pretty=pretty)
    
    def scale_va(self,data,div_flag=False):
        '''! Scale a vertex array by vertex areas (eg. B_n -> flux)

        @param data Data to scale
        @param div_flag Divide by vertex areas instead?
        '''
        data_in = numpy.ascontiguousarray(data.copy(), dtype=numpy.float64)
        thincurr_scale_va(self.tw_obj,data_in,div_flag)
        return data_in
    
    def compute_Lmat(self,cache_file=None,use_hodlr=False):
        '''! Compute the self-inductance matrix for this model

        @param cache_file Path to cache file to store/load matrix
        @param use_hodlr Use HODLR compression for the matrix
        '''
        if cache_file is None:
            cache_string = c_char_p(b"")
        else:
            cache_string = c_char_p(cache_file.encode())
        Lmat_loc = c_void_p()
        error_string = c_char_p(b""*200)
        thincurr_curr_Lmat(self.tw_obj,use_hodlr,ctypes.byref(Lmat_loc),cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        if use_hodlr:
            self.Lmat_hodlr = Lmat_loc
        else:
            self.Lmat = numpy.ctypeslib.as_array(ctypes.cast(Lmat_loc, c_double_ptr),shape=(self.nelems,self.nelems))
    
    def compute_Mcoil(self,cache_file=None):
        '''! Compute the mutual inductance between passive (mesh+vcoils) and active elements (icoils)

        @param cache_file Path to cache file to store/load matrix
        @result Mutual inductance matrix
        '''
        if cache_file is None:
            cache_string = c_char_p(b"")
        else:
            cache_string = c_char_p(cache_file.encode())
        Mc_loc = c_void_p()
        error_string = c_char_p(b""*200)
        thincurr_Mcoil(self.tw_obj,ctypes.byref(Mc_loc),cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return numpy.ctypeslib.as_array(ctypes.cast(Mc_loc, c_double_ptr),shape=(self.n_icoils,self.nelems))

    def compute_Msensor(self,sensor_file,cache_file=None):
        '''! Compute the mutual inductance between model and sensors

        @param sensor_file Path to file contatining flux loop definitions
        @param cache_file Path to cache file to store/load matrix
        @result Mutual inductance matrix between model and sensors, Mutual inductance matrix between icoils and sensors
        '''
        if cache_file is None:
            cache_string = c_char_p(b"")
        else:
            cache_string = c_char_p(cache_file.encode())
        sensor_string = c_char_p(sensor_file.encode())
        Ms_loc = c_void_p()
        Msc_loc = c_void_p()
        nsensors = c_int()
        sensor_loc = c_void_p()
        error_string = c_char_p(b""*200)
        thincurr_Msensor(self.tw_obj,sensor_string,ctypes.byref(Ms_loc),ctypes.byref(Msc_loc),ctypes.byref(nsensors),ctypes.byref(sensor_loc),cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return numpy.ctypeslib.as_array(ctypes.cast(Ms_loc, c_double_ptr),shape=(self.nelems,nsensors.value)), \
               numpy.ctypeslib.as_array(ctypes.cast(Msc_loc, c_double_ptr),shape=(self.n_icoils,nsensors.value)), sensor_loc
    
    def compute_Rmat(self,copy_out=False):
        '''! Compute the resistance matrix for this model

        @param copy_out Copy matrix to python and store in `self.Rmat`?
        '''
        if copy_out:
            self.Rmat = numpy.zeros((self.nelems,self.nelems), dtype=numpy.float64)
            Rmat_tmp = self.Rmat
        else:
            Rmat_tmp = numpy.zeros((1,1), dtype=numpy.float64)
        error_string = c_char_p(b""*200)
        thincurr_curr_Rmat(self.tw_obj,c_bool(copy_out),Rmat_tmp,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
    
    def cross_coupling(self,model2,cache_file=None):
        '''! Compute the mutual inductance between this and another ThinCurr model

        @param model2 The second model for mutual calculation
        @param cache_file Path to cache file to store/load matrix
        '''
        Mmat = numpy.zeros((self.nelems,model2.nelems), dtype=numpy.float64)
        if cache_file is None:
            cache_string = c_char_p(b"")
        else:
            cache_string = c_char_p(cache_file.encode())
        error_string = c_char_p(b""*200)
        thincurr_cross_coupling(self.tw_obj,model2.tw_obj,Mmat,cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return Mmat
    
    def cross_eval(self,model2,field):
        '''! Compute the voltage/current induced on another ThinCurr model from a current structure on this model

        @param model2 The second model for mutual calculation
        @param field One or more current fields
        '''
        nrhs = field.shape[0]
        vec_out = numpy.zeros((nrhs,model2.nelems), dtype=numpy.float64)
        vec_in = numpy.ascontiguousarray(field.copy(), dtype=numpy.float64)
        error_string = c_char_p(b""*200)
        thincurr_cross_eval(self.tw_obj,model2.tw_obj,c_int(nrhs),vec_in,vec_out,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return vec_out
    
    def get_regmat(self):
        '''! Compute the current regularization matrix for this model
        '''
        Rmat = numpy.zeros((self.nelems,3*self.nc), dtype=numpy.float64)
        error_string = c_char_p(b""*200)
        thincurr_curr_regmat(self.tw_obj,Rmat,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return Rmat
    
    def get_eigs(self,neigs,direct=False):
        '''! Compute eigenmodes of this model

        @param neigs Number of eigenvalues to compute
        @param direct Use direct solver?
        @result Eigenvalues, Eigenvectors
        '''
        eig_vals = numpy.zeros((neigs,), dtype=numpy.float64)
        eig_vecs = numpy.zeros((neigs,self.nelems), dtype=numpy.float64)
        error_string = c_char_p(b""*200)
        if self.Lmat_hodlr:
            thincurr_eigenvalues(self.tw_obj,c_bool(direct),c_int(neigs),eig_vals,eig_vecs,self.Lmat_hodlr,error_string)
        else:
            thincurr_eigenvalues(self.tw_obj,c_bool(direct),c_int(neigs),eig_vals,eig_vecs,c_void_p(),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return eig_vals, eig_vecs

    def compute_freq_response(self,driver,freq=0.0,fr_limit=0,direct=False):
        '''! Compute frequence response of the current model to given driver voltages

        @param driver Real/Imaginary driver voltage pair
        @param freq Frequency for response calculation [Hz] (unused if `fr_limit!=0`)
        @param fr_limit Frequency limit for calculation (0: none, 1: inductive, 2: resistive)
        @param direct Use direct solver?
        @result Real/Imaginary eddy current response
        '''
        result = numpy.ascontiguousarray(driver.copy(), dtype=numpy.float64)
        error_string = c_char_p(b""*200)
        if self.Lmat_hodlr:
            thincurr_freq_response(self.tw_obj,c_bool(direct),c_int(fr_limit),c_double(freq),result,self.Lmat_hodlr,error_string)
        else:
            thincurr_freq_response(self.tw_obj,c_bool(direct),c_int(fr_limit),c_double(freq),result,c_void_p(),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
        return result
    
    def run_td(self,dt,nsteps,coil_currs=None,coil_volts=None,direct=False,status_freq=10,plot_freq=10,sensor_obj=c_void_p()):
        '''! Perform a time-domain simulation

        @param dt Time step for simulation
        @param nsteps Number of steps to take
        @param coil_currs Current vs time array for Icoils [:,n_icoils+1] (first index is time)
        @param coil_volts Voltage vs time array for Vcoils [:,n_vcoils+1] (first index is time)
        @param direct Use direct solver?
        @param status_freq Frequency to print status information
        @param plot_freq Frequency to save plot files
        @param sensor_obj Sensor object to use
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
        if coil_volts is None:
            nvolt = c_int(0)
            coil_volts = numpy.zeros((1,1), dtype=numpy.float64)
        else:
            if coil_volts.shape[1]-1 != self.n_vcoils:
                raise ValueError("# of voltages in waveform does not match # of vcoils")
            nvolt = c_int(coil_volts.shape[0])
            coil_volts = numpy.ascontiguousarray(coil_volts.transpose(), dtype=numpy.float64)
        error_string = c_char_p(b""*200)
        if self.Lmat_hodlr:
            thincurr_time_domain(self.tw_obj,c_bool(direct),c_double(dt),c_int(nsteps),c_int(status_freq),c_int(plot_freq),vec_ic,
                                sensor_obj,ncurr,coil_currs,nvolt,coil_volts,self.Lmat_hodlr,error_string)
        else:
            thincurr_time_domain(self.tw_obj,c_bool(direct),c_double(dt),c_int(nsteps),c_int(status_freq),c_int(plot_freq),vec_ic,
                                sensor_obj,ncurr,coil_currs,nvolt,coil_volts,c_void_p(),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())

    def build_reduced_model(self,basis_set,filename='tCurr_reduced.h5',sensor_obj=c_void_p()):
        r'''! Build reduced model by projecting full model onto defined basis set of currents

        @param basis_set Basis set for projection [nBasis,:]
        @param filename Filename for reduction model
        @param sensor_obj Sensor object to use
        @result Reduced model (see \ref ThinCurr_reduced)
        '''
        basis_set = numpy.ascontiguousarray(basis_set, dtype=numpy.float64)
        nbasis = c_int(basis_set.shape[0])
        error_string = c_char_p(b""*200)
        if self.Lmat_hodlr:
            thincurr_reduce_model(self.tw_obj,c_char_p(filename.encode()),nbasis,basis_set,sensor_obj,self.Lmat_hodlr,error_string)
        else:
            thincurr_reduce_model(self.tw_obj,c_char_p(filename.encode()),nbasis,basis_set,sensor_obj,c_void_p(),error_string)
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
            self.Basis = numpy.asarray(file['Basis'])
            self.L = numpy.asarray(file['L'])
            self.R = numpy.asarray(file['R'])
            if 'Ms' in file:
                self.Ms = numpy.asarray(file['Ms'])
            else:
                self.Ms = None
            if 'Mc' in file:
                self.Mc = numpy.asarray(file['Mc'])
            else:
                self.Mc = None
            if 'Msc' in file:
                self.Msc = numpy.asarray(file['Msc'])
            else:
                self.Msc = None

    def get_eigs(self):
        '''! Compute eigenmodes for reduced model

        @result Eigenvalues, Eigenvectors
        '''
        eig_vals, eig_vecs = numpy.linalg.eig(numpy.dot(numpy.linalg.inv(self.R),self.L))
        sort_inds = (-eig_vals).argsort()
        return eig_vals[sort_inds], numpy.dot(eig_vecs[sort_inds,:],self.Basis)

    def run_td(self,dt,nsteps,coil_currs,status_freq=10,plot_freq=10):
        '''! Perform a time-domain simulation

        @param dt Time step for simulation
        @param nsteps Number of steps to take
        @param coil_currs Current vs time array for Icoils [:,n_icoils+1] (first index is time)
        @param status_freq Frequency to print status information
        @param plot_freq Frequency to save plot files
        '''
        Lforward = self.L - (dt/2.0)*self.R
        Lbackward = numpy.linalg.inv(self.L + (dt/2.0)*self.R)
        #
        vec_time = []
        vec_hist = []
        sen_time = []
        sen_hist = []
        pot_tmp = numpy.zeros((self.L.shape[0],))
        t = 0.0
        print('Timestep {0} {1:12.4E} {2:12.4E}'.format(0,t,numpy.linalg.norm(pot_tmp)))
        vec_time.append(t)
        vec_hist.append(numpy.dot(pot_tmp,self.Basis))
        if self.Ms is not None:
            sen_tmp = numpy.dot(pot_tmp,self.Ms)
            if self.Msc is not None:
                curr_tmp = numpy.zeros((coil_currs.shape[1]-1,))
                for j in range(coil_currs.shape[1]-1):
                    curr_tmp[j] = numpy.interp(t,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1])
                curr_tmp *= mu0
                sen_tmp += numpy.dot(curr_tmp,self.Msc)
            sen_time.append(t)
            sen_hist.append(sen_tmp)
        #
        for i in range(nsteps):
            rhs = numpy.dot(pot_tmp,Lforward)
            if self.Mc is not None:
                dcurr_tmp = numpy.zeros((coil_currs.shape[1]-1,))
                for j in range(coil_currs.shape[1]-1):
                    dcurr_tmp[j] = numpy.interp(t+dt/4.0,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1])
                    dcurr_tmp[j] -= numpy.interp(t-dt/4.0,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1])
                    dcurr_tmp[j] += numpy.interp(t+dt*5.0/4.0,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1])
                    dcurr_tmp[j] -= numpy.interp(t+dt*3.0/4.0,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1])
                dcurr_tmp *= mu0
                rhs -= numpy.dot(dcurr_tmp,self.Mc)
            pot_tmp = numpy.dot(rhs,Lbackward)
            t += dt
            if ((i+1) % status_freq) == 0:
                print('Timestep {0} {1:12.4E} {2:12.4E}'.format(i+1,t,numpy.linalg.norm(pot_tmp)))
            if ((i+1) % plot_freq) == 0:
                vec_time.append(t)
                vec_hist.append(numpy.dot(pot_tmp,self.Basis))
            if self.Ms is not None:
                sen_tmp = numpy.dot(pot_tmp,self.Ms)
                if self.Msc is not None:
                    curr_tmp = numpy.zeros((coil_currs.shape[1]-1,))
                    for j in range(coil_currs.shape[1]-1):
                        curr_tmp[j] = numpy.interp(t,coil_currs[:,0],coil_currs[:,j+1],left=coil_currs[0,j+1])
                    curr_tmp *= mu0
                    sen_tmp += numpy.dot(curr_tmp,self.Msc)
                sen_time.append(t)
                sen_hist.append(sen_tmp)
        #
        return numpy.array(sen_time), numpy.array(sen_hist), numpy.array(vec_time), numpy.array(vec_hist)