#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Core definitions for ThickCurr thick-wall electromagnetic functionality

@authors Chris Hansen
@date November 2026
@ingroup doxy_oft_python
'''
import os
import numpy
from ._interface import *
from ..io import build_XDMF, histfile
from ..ThinCurr.sensor import flux_loop, Mirnov


def _gauss_legendre(n):
    '''! Gauss-Legendre points and weights on [0,1]'''
    x, w = numpy.polynomial.legendre.leggauss(n)
    return (x+1.0)/2.0, w/2.0


def _read_sensor_file(filename):
    '''! Read sensors from a ThinCurr sensor file (eg. `floops.loc`)

    @param filename Sensor file
    @result List of sensor objects
    '''
    sensors = []
    with open(filename, 'r') as fid:
        lines = [line.strip() for line in fid if line.strip() and not line.strip().startswith('#')]
    nsensors = int(lines[0])
    i = 1
    for _ in range(nsensors):
        np_str, scale_str, name = lines[i].split()[:3]
        npts = int(np_str)
        if npts < 0:
            pt = numpy.array([float(val) for val in lines[i+1].split()[:3]])
            norm = numpy.array([float(val) for val in lines[i+2].split()[:3]])
            sensors.append(Mirnov(pt, norm, name, dx=None))
            i += 3
        else:
            pts = numpy.array([[float(val) for val in line.split()[:3]] for line in lines[i+1:i+1+npts]])
            sensors.append(flux_loop(pts, name, scale=float(scale_str)))
            i += 1+npts
    return sensors


class ThickCurr():
    r'''! ThickCurr thick-wall electromagnetic simulation class

    Evolves the magnetic vector potential \f$ A \f$ in conducting structures and the surrounding (pseudo-)vacuum
    \f[ \frac{\mu_0}{\eta} \frac{\partial A}{\partial t} + \nabla \times \nabla \times A = 0, \f]
    within a domain bounded by a perfectly conducting wall, using H(Curl) + Grad(H^1) finite elements (see
    @ref thick_wall "thick_wall" for details).

    @note As the underlying Fortran implementation stores state at the module level only one ThickCurr object
    may be active (set up) at a time.
    '''
    def __init__(self,OFT_env):
        '''! Initialize ThickCurr object

        @param OFT_env OFT runtime environment object (See @ref OpenFUSIONToolkit._core.OFT_env "OFT_env")
        '''
        ## OFT execution environment
        self._oft_env = OFT_env
        ## Internal ThickCurr object
        self._tc_ptr = c_void_p()
        ## Internal mesh object
        self._mesh_ptr = c_void_p()
        ## Internal XML object
        self._xml_ptr = c_void_p()
        ## Number of regions in mesh
        self.nregs = -1
        ## Sensor objects created by this instance (destroyed with object)
        self._sensor_objs = []
        ## Time step [s] (also the default shift for eigenvalue solves)
        self.dt = None

    def __del__(self):
        '''! Destroy ThickCurr object'''
        if self._oft_env.mpi_finalized:
            return
        for sensor_obj in self._sensor_objs:
            if sensor_obj['ptr']:
                thickcurr_f_destroy_sensors(sensor_obj['ptr'])
                sensor_obj['ptr'] = c_void_p()
        self._sensor_objs = []
        if (not self._tc_ptr) or self._oft_env.mpi_finalized:
            return
        error_string = self._oft_env.get_c_errorbuff()
        thickcurr_f_destroy(self._tc_ptr,error_string)
        self._tc_ptr = c_void_p()

    def print_ascii_logo(self, italic=True):
        '''! Print ThickCurr ASCII logo

        @param italic Print italicized logo?'''
        if italic:
            print(r'''
  ________    _      __   ______
 /_  __/ /_  (_)____/ /__/ ____/_  ___________
  / / / __ \/ / ___/ //_/ /   / / / / ___/ ___/
 / / / / / / / /__/ ,< / /___/ /_/ / /  / /
/_/ /_/ /_/_/\___/_/|_|\____/\__,_/_/  /_/
''')
        else:
            print(r'''
 _____ _     _      _     ____
|_   _| |__ (_) ___| | __/ ___|   _ _ __ _ __
  | | | '_ \| |/ __| |/ / |  | | | | '__| '__|
  | | | | | | | (__|   <| |__| |_| | |  | |
  |_| |_| |_|_|\___|_|\_\\____\__,_|_|  |_|
''')

    def _check_setup(self):
        '''! Check that model has been setup'''
        if not self._tc_ptr:
            raise ValueError('Model has not been setup, call "setup_model()" first')

    @staticmethod
    def _check_error(error_string):
        '''! Raise exception if Fortran call reported an error'''
        if error_string.value != b'':
            raise Exception(error_string.value.decode())

    def setup_model(self,r=None,lc=None,reg=None,mesh_file=None,grid_order=1,order=2,eta_values=None,xml_filename=None,
                    pmg=False,minlev=-1,dt=1.E-3,theta=1.0,eta_vac=-1.0,lin_tol=1.E-9,nbase=1,nlevels=None):
        r'''! Setup ThickCurr model (mesh, finite element spaces, and operators)

        A mesh should be specified by passing "r", "lc", and optionally "reg" or using a "mesh_file".

        The resistivity of each region is specified using `eta_values` or the `region` elements in the
        `thickcurr` element of the XML file (`xml_filename`). Regions with a resistivity \f$ \le 0 \f$ are treated
        as a pseudo-vacuum with resistivity `eta_vac`. The XML file may also contain a `pre` element defining the
        preconditioner for the linear solver (see @ref thick_wall::thickcurr_read_settings "thickcurr_read_settings").

        When running with MPI the mesh is decomposed across processes, which requires at least one level
        after the base level(s) (`nlevels > nbase`). By default `nlevels=nbase` in serial and `nlevels=nbase+1`
        with MPI, so no refinement is performed.

        @param r Mesh point list [np,3]
        @param lc Mesh cell list [nc,4] (base zero)
        @param reg Mesh region list [nc] (base one)
        @param mesh_file Filename containing mesh to load (native format only)
        @param grid_order Order of mesh representation (1: linear, 2: quadratic)
        @param order Order of finite element basis functions
        @param eta_values Resistivity for each region [Ohm-m] [nregs] (<= 0 for vacuum, overrides XML)
        @param xml_filename Path to XML file with `thickcurr` element (resistivity and/or preconditioner)
        @param pmg Use p-multigrid, with linear elements on the finest mesh as the coarse level?
        @param minlev Lowest FE level for multigrid (-1: no multigrid, ignored if `pmg=True`)
        @param dt Time step [s] (can be changed later in @ref run_td)
        @param theta Time centering (1: backward Euler, 0.5: Crank-Nicolson)
        @param eta_vac Resistivity of pseudo-vacuum regions [Ohm-m] (negative: 1E4*max(eta))
        @param lin_tol Absolute tolerance for linear solver
        @param nbase Number of mesh levels before domain decomposition
        @param nlevels Total number of mesh levels (default: see above)
        '''
        if self._tc_ptr:
            raise ValueError('Model already setup')
        if nlevels is None:
            nlevels = nbase+1 if self._oft_env.nprocs > 1 else nbase
        if (self._oft_env.nprocs > 1) and (nlevels <= nbase):
            raise ValueError('"nlevels" must be > "nbase" when running with MPI')
        #---Mesh and ThickCurr settings
        self._oft_env.oft_in_groups['mesh_options'].update({
            'cad_type': '0',
            'nlevels': '{0:d}'.format(nlevels),
            'nbase': '{0:d}'.format(nbase),
            'grid_order': '{0:d}'.format(grid_order)
        })
        self._oft_env.oft_in_groups['thickcurr_options'] = {
            'dt': '{0:.16E}'.format(dt),
            'theta': '{0:.16E}'.format(theta),
            'eta_vac': '{0:.16E}'.format(eta_vac),
            'lin_tol': '{0:.16E}'.format(lin_tol)
        }
        #---Load mesh
        nregs = c_int()
        if mesh_file is not None:
            if not os.path.isfile(mesh_file):
                raise ValueError('Mesh file "{0}" does not exist'.format(mesh_file))
            self._oft_env.oft_in_groups['native_mesh_options'] = {'filename': '"{0}"'.format(mesh_file)}
            self._oft_env.update_oft_in()
            np = c_int(-1)
            rfake = numpy.ones((1,1),dtype=numpy.float64)
            lcfake = numpy.ones((1,1),dtype=numpy.int32)
            regfake = numpy.ones((1,),dtype=numpy.int32)
            oft_setup_vmesh(np,rfake,np,np,lcfake,regfake,ctypes.byref(nregs),ctypes.byref(self._mesh_ptr))
        elif r is not None:
            self._oft_env.update_oft_in()
            r = numpy.ascontiguousarray(r, dtype=numpy.float64)
            lc = numpy.ascontiguousarray(lc, dtype=numpy.int32)
            if (r.ndim != 2) or (r.shape[1] != 3):
                raise ValueError('Point list must have 3 spatial dimensions `[np,3]`')
            if (lc.ndim != 2) or (lc.shape[1] != 4):
                raise ValueError('Cell list must be tetrahedra `[nc,4]`')
            if reg is None:
                reg = numpy.ones((lc.shape[0],),dtype=numpy.int32)
            else:
                reg = numpy.ascontiguousarray(reg, dtype=numpy.int32)
            oft_setup_vmesh(c_int(r.shape[0]),r,c_int(4),c_int(lc.shape[0]),lc+1,reg,ctypes.byref(nregs),ctypes.byref(self._mesh_ptr))
        else:
            raise ValueError('Mesh filename (native format) or mesh values required')
        self.nregs = nregs.value
        #---Load XML
        if xml_filename is not None:
            if not os.path.isfile(xml_filename):
                raise ValueError('XML file "{0}" does not exist'.format(xml_filename))
            oftpy_load_xml(self._oft_env.path2c(xml_filename),ctypes.byref(self._xml_ptr))
            if not self._xml_ptr:
                raise ValueError('Error parsing XML file "{0}"'.format(xml_filename))
        #---Resistivities (from XML if not specified)
        if eta_values is not None:
            eta_values = numpy.ascontiguousarray(eta_values, dtype=numpy.float64)
            if eta_values.shape != (self.nregs,):
                raise ValueError('Size of "eta_values" ({0}) does not match number of regions ({1})'.format(eta_values.size,self.nregs))
            neta = self.nregs
        elif xml_filename is None:
            raise ValueError('Resistivities must be specified using "eta_values" or "xml_filename"')
        else:
            eta_values = numpy.zeros((1,),dtype=numpy.float64)
            neta = 0
        if minlev == 0:
            raise ValueError('Invalid "minlev", use "pmg=True" for p-multigrid')
        error_string = self._oft_env.get_c_errorbuff()
        thickcurr_f_setup(ctypes.byref(self._tc_ptr),self._mesh_ptr,c_int(order),c_int(0 if pmg else minlev),c_int(neta),
                          eta_values,self._xml_ptr,error_string)
        self._check_error(error_string)
        self.dt = dt

    def set_options(self,dt=None,theta=None,lin_tol=None):
        '''! Update run options after setup

        @note Changing `dt` or `theta` requires rebuilding the linear system and preconditioner

        @param dt Time step [s]
        @param theta Time centering (1: backward Euler, 0.5: Crank-Nicolson)
        @param lin_tol Absolute tolerance for linear solver
        '''
        self._check_setup()
        error_string = self._oft_env.get_c_errorbuff()
        thickcurr_f_set_options(self._tc_ptr,c_double(-1.0 if dt is None else dt),c_double(-1.0 if theta is None else theta),
                                c_double(-1.0 if lin_tol is None else lin_tol),c_int(-1),c_int(-1),error_string)
        self._check_error(error_string)
        if dt is not None:
            self.dt = dt

    def get_quadrature_points(self):
        '''! Get locations used to sample fields for @ref set_field

        @note With MPI only points in the local domain are returned

        @result Point locations [npts,3]
        '''
        self._check_setup()
        npts = c_int()
        error_string = self._oft_env.get_c_errorbuff()
        thickcurr_f_get_qpts(self._tc_ptr,ctypes.byref(npts),c_void_p(),error_string)
        self._check_error(error_string)
        pts = numpy.zeros((npts.value,3), dtype=numpy.float64)
        thickcurr_f_get_qpts(self._tc_ptr,ctypes.byref(npts),pts.ctypes.data_as(c_void_p),error_string)
        self._check_error(error_string)
        return pts

    def set_field(self,A_func=None,zero_boundary=True):
        r'''! Set the vector potential (initial condition for @ref run_td)

        The field is projected onto the finite element basis by sampling `A_func` at the quadrature points
        of every cell (see @ref get_quadrature_points), which is called once with all points.

        As tangential components of \f$ A \f$ on the boundary are held fixed by the perfectly conducting wall,
        by default these components are set to zero (`zero_boundary=True`). If `zero_boundary=False` the projected
        boundary values are retained, which corresponds to a static field applied from outside the domain.

        @param A_func Function returning the vector potential [T-m] at a set of points `A_func(pts[n,3]) -> [n,3]`
        (`None` sets a zero field)
        @param zero_boundary Zero tangential components of the vector potential on the boundary?
        @result Magnetic energy of the projected field [J]
        '''
        self._check_setup()
        energy = c_double()
        error_string = self._oft_env.get_c_errorbuff()
        if A_func is None:
            thickcurr_f_set_field(self._tc_ptr,c_void_p(),c_bool(zero_boundary),ctypes.byref(energy),error_string)
        else:
            pts = self.get_quadrature_points()
            vals = numpy.ascontiguousarray(A_func(pts), dtype=numpy.float64)
            if vals.shape != pts.shape:
                raise ValueError('"A_func" must return an array with shape [n,3]')
            thickcurr_f_set_field(self._tc_ptr,vals.ctypes.data_as(c_void_p),c_bool(zero_boundary),ctypes.byref(energy),error_string)
        self._check_error(error_string)
        return energy.value

    def get_energy(self):
        r'''! Get the magnetic energy of the current vector potential

        @result Magnetic energy \f$ \int B^2 / 2 \mu_0 dV \f$ [J]
        '''
        self._check_setup()
        energy = c_double()
        error_string = self._oft_env.get_c_errorbuff()
        thickcurr_f_get_energy(self._tc_ptr,ctypes.byref(energy),error_string)
        self._check_error(error_string)
        return energy.value

    def get_eigs(self,neigs,shift=None,save_plot=False):
        r'''! Compute the slowest decaying eigenmodes (L/R time constants) of the model

        The eigenvalue problem is solved using a shift, which by default is the (effective) time step
        \f$ \theta \Delta t \f$. Convergence is typically much faster with a shift comparable to the slowest decay
        time (see @ref thick_wall::thickcurr_eigs "thickcurr_eigs").

        @param neigs Number of eigenvalues to compute
        @param shift Shift for eigenvalue solve [s] (default: time step)
        @param save_plot Save eigenmodes to plot files (`eigs` subdirectory, mode index as "time")
        @result Decay times (slowest first) [s] `(neigs)`
        '''
        self._check_setup()
        tau = numpy.zeros((neigs,), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        thickcurr_f_eigs(self._tc_ptr,c_int(neigs),tau,c_double(-1.0 if shift is None else shift),c_bool(save_plot),error_string)
        self._check_error(error_string)
        return tau

    def setup_sensors(self,sensors,filename='floops.hist',nquad=3):
        r'''! Setup sensors for time-dependent simulations

        Sensors are defined using the same classes as ThinCurr (@ref OpenFUSIONToolkit.ThinCurr.sensor), where:
          - Flux loops (@ref OpenFUSIONToolkit.ThinCurr.sensor.flux_loop "flux_loop") measure the flux through the
            loop formed by the polyline through the points as \f$ \oint A \cdot dl \f$, which does not require
            the mesh to conform to a surface bounded by the loop. Signals are multiplied by the sensor's scale factor.
          - Point Mirnov sensors (@ref OpenFUSIONToolkit.ThinCurr.sensor.Mirnov "Mirnov" with `dx=None`) measure
            \f$ B \cdot \hat{n} \f$ at a point (finite size Mirnov sensors are treated as flux loops)

        @param sensors List of sensor objects or path to ThinCurr sensor file (eg. `floops.loc`)
        @param filename Filename for sensor history file
        @param nquad Number of Gauss-Legendre points used on each segment of flux loops
        @result Sensor object (dictionary with keys `ptr`, `names`, and `filename`)
        '''
        self._check_setup()
        if isinstance(sensors, str):
            sensors = _read_sensor_file(sensors)
        qx, qw = _gauss_legendre(nquad)
        names = []
        ptsa, wa, idxa = [], [], []
        ptsb, wb, idxb = [], [], []
        for i, sensor in enumerate(sensors):
            if not isinstance(sensor, flux_loop):
                raise ValueError('Invalid sensor object')
            if len(sensor._name) > 39:
                raise ValueError('Sensor name "{0}" too long (max 39 characters)'.format(sensor._name))
            names.append(sensor._name)
            if getattr(sensor, '_point', False):
                ptsb.append(sensor._pts[0,:])
                wb.append(sensor._pts[1,:]*sensor._scale)
                idxb.append(i+1)
                continue
            pts = numpy.asarray(sensor._pts, dtype=numpy.float64)
            for j in range(pts.shape[0]-1):
                dl = pts[j+1,:]-pts[j,:]
                for x, w in zip(qx, qw):
                    ptsa.append(pts[j,:]+x*dl)
                    wa.append(sensor._scale*w*dl)
                    idxa.append(i+1)
        def _pack(pts, w, idx):
            if len(idx) == 0:
                return numpy.zeros((1,3)), numpy.zeros((1,3)), numpy.zeros((1,),dtype=numpy.int32), 0
            return (numpy.ascontiguousarray(pts, dtype=numpy.float64), numpy.ascontiguousarray(w, dtype=numpy.float64),
                    numpy.ascontiguousarray(idx, dtype=numpy.int32), len(idx))
        ptsa, wa, idxa, npa = _pack(ptsa, wa, idxa)
        ptsb, wb, idxb, npb = _pack(ptsb, wb, idxb)
        names_c = c_char_p(b''.join([name.encode().ljust(40, b'\0') for name in names]))
        sensor_ptr = c_void_p()
        error_string = self._oft_env.get_c_errorbuff()
        thickcurr_f_setup_sensors(self._tc_ptr,c_int(len(names)),names_c,c_int(npa),ptsa,wa,idxa,c_int(npb),ptsb,wb,idxb,
                                  self._oft_env.path2c(filename),ctypes.byref(sensor_ptr),error_string)
        self._check_error(error_string)
        sensor_obj = {'ptr': sensor_ptr, 'names': names, 'filename': filename}
        self._sensor_objs.append(sensor_obj)
        return sensor_obj

    def eval_sensors(self,sensor_obj):
        '''! Evaluate sensor signals for the current vector potential

        @param sensor_obj Sensor object (see @ref setup_sensors)
        @result Sensor signals `(nsensors)`
        '''
        self._check_setup()
        vals = numpy.zeros((len(sensor_obj['names']),), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        thickcurr_f_eval_sensors(self._tc_ptr,sensor_obj['ptr'],vals,error_string)
        self._check_error(error_string)
        return vals

    def run_td(self,dt=None,nsteps=None,theta=None,plot_freq=10,sensor_obj=None,lin_tol=None):
        r'''! Perform a time-domain simulation starting from the current vector potential (see @ref set_field)

        On completion the vector potential is updated to the final state, so additional calls continue
        the simulation (restarting the time at zero).

        @note Output files are written to the current working directory on the head (rank 0) process:
          - `thickcurr.hist`: Step index, linear iterations, time, magnetic energy, Ohmic power, and walltime
          - `thickcurr_XXXXX.rst`: Restart files with the vector potential and its time derivative, every `plot_freq` steps
          - Sensor history file if `sensor_obj` is passed (see @ref setup_sensors)

        @param dt Time step [s] (default: current value)
        @param nsteps Number of steps to take
        @param theta Time centering (1: backward Euler, 0.5: Crank-Nicolson) (default: current value)
        @param plot_freq Frequency to save restart (plot) files
        @param sensor_obj Sensor object to record signals (see @ref setup_sensors)
        @param lin_tol Absolute tolerance for linear solver (default: current value)
        @result Time at the end of the simulation [s]
        '''
        self._check_setup()
        if nsteps is None:
            raise ValueError('"nsteps" must be specified')
        error_string = self._oft_env.get_c_errorbuff()
        thickcurr_f_set_options(self._tc_ptr,c_double(-1.0 if dt is None else dt),c_double(-1.0 if theta is None else theta),
                                c_double(-1.0 if lin_tol is None else lin_tol),c_int(nsteps),c_int(plot_freq),error_string)
        self._check_error(error_string)
        if dt is not None:
            self.dt = dt
        sensor_ptr = c_void_p() if sensor_obj is None else sensor_obj['ptr']
        t_final = c_double()
        thickcurr_f_run(self._tc_ptr,sensor_ptr,ctypes.byref(t_final),error_string)
        self._check_error(error_string)
        return t_final.value

    def plot_td(self,nsteps=None,plot_freq=10,sensor_obj=None):
        r'''! Generate plot files from the restart files of a time domain simulation

        Fields are saved to the plot group `thickcurr`, including the vector potential (`A`), magnetic field (`B`),
        and when available the electric field (`E`) and current density (`J`).

        @param nsteps Last step to plot (default: all available)
        @param plot_freq Frequency of restart files (must match `plot_freq` used in @ref run_td)
        @param sensor_obj Sensor object to recompute signals from restart files (see @ref setup_sensors)
        '''
        self._check_setup()
        error_string = self._oft_env.get_c_errorbuff()
        thickcurr_f_set_options(self._tc_ptr,c_double(-1.0),c_double(-1.0),c_double(-1.0),c_int(-1),c_int(plot_freq),error_string)
        self._check_error(error_string)
        sensor_ptr = c_void_p() if sensor_obj is None else sensor_obj['ptr']
        thickcurr_f_plot(self._tc_ptr,sensor_ptr,c_int(0),c_int(10**8 if nsteps is None else nsteps),error_string)
        self._check_error(error_string)

    def build_XDMF(self,repeat_static=False,pretty=False):
        '''! Build XDMF plot metadata files for time-dependent (current directory) and eigenmode (`eigs`) plot files

        @note Only performed on the head (rank 0) process

        @param repeat_static Repeat static fields (0-th timestep) in all timesteps?
        @param pretty Use pretty printing (indentation) in XDMF files?
        @result Plot file objects for time-dependent and eigenmode results (`None` if not present)
        '''
        if not self._oft_env.head_proc:
            return None, None
        td_file = None
        eig_file = None
        if os.path.isfile('oft_xdmf.0001.h5'):
            td_file = build_XDMF(path='.',repeat_static=repeat_static,pretty=pretty)
        if os.path.isfile(os.path.join('eigs','oft_xdmf.0001.h5')):
            eig_file = build_XDMF(path='eigs',repeat_static=repeat_static,pretty=pretty)
        return td_file, eig_file

    @staticmethod
    def load_history(filename='thickcurr.hist'):
        '''! Load history file from a time-dependent simulation

        @param filename History file (`thickcurr.hist` or a sensor history file)
        @result History file object (see @ref OpenFUSIONToolkit.io.histfile "histfile")
        '''
        return histfile(filename)
