#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Functionality for performing equilibrium reconstructions using TokaMaker

@authors Chris Hansen
@date April 2024
@ingroup doxy_oft_python
'''
from warnings import warn
from .._interface import *
from ..util import oft_warning

## @cond
class tokamaker_recon_settings_cstruct(c_struct):
    _fields_ = [("fitF", c_bool),
                ("fitP", c_bool),
                ("fit_Pscale", c_bool),
                ("fit_FFPscale", c_bool),
                ("fitR0", c_bool),
                ("fitZ0", c_bool),
                ("fitCoils", c_bool),
                ("fitF0", c_bool),
                ("fixedCentering", c_bool),
                ("pm", c_bool),
                ("infile", ctypes.c_char_p),
                ("outfile", ctypes.c_char_p)]

# tokamaker_recon_run(tMaker_ptr,vacuum,settings,error_flag)
tokamaker_recon_run = ctypes_subroutine(oftpy_lib.tokamaker_recon_run,
    [c_void_p, c_bool, tokamaker_recon_settings_cstruct, c_int_ptr])

# tokamaker_recon_err(tMaker_ptr,vacuum,settings,error_mat,error_flag)
tokamaker_recon_err = ctypes_subroutine(oftpy_lib.tokamaker_recon_err,
    [c_void_p, c_bool, tokamaker_recon_settings_cstruct, c_double_ptr, c_int_ptr])

# tokamaker_recon_setup(tMaker_ptr,settings,ncons,error_flag)
tokamaker_recon_setup = ctypes_subroutine(oftpy_lib.tokamaker_recon_setup,
    [c_void_p, tokamaker_recon_settings_cstruct, c_int_ptr, c_int_ptr])

# tokamaker_recon_destroy(tMaker_ptr,error_flag)
tokamaker_recon_destroy = ctypes_subroutine(oftpy_lib.tokamaker_recon_destroy,
    [c_void_p, c_int_ptr])
## @endcond

Mirnov_con_id = 1
Ip_con_id = 2
fluxLoop_con_id = 7
dFlux_con_id = 8
Press_con_id = 9
q_con_id = 10
saddle_con_id = 11
coil_current_con_id = 12


class Mirnov_con:
    '''! TokaMaker equilibrium reconstruction Mirnov sensor constraint'''
    def __init__(self, loc=None, phi=0., norm=None, val=None, err=None):
        r'''! Create Mirnov sensor constraint
        
        @param loc Location of Mirnov in R-Z plane [2]
        @param phi Toroidal location [rad] (only meaningful with 3D fields)
        @param norm Unit normal \f$ (R,\phi,Z) \f$ [3]
        @param val Value of constraint
        @param err Error in constraint
        '''
        self.loc = loc
        self.phi = phi
        self.norm = norm
        self.val = val
        self.err = err

    def read(self, file):
        '''! Read Mirnov constraint from file
        
        @param file Open file object containing constraint, must be positioned at start of constraint
        '''
        values = file.readline().split()
        self.loc = (float(values[0]), float(values[1]))
        self.phi = float(values[2])
        values = file.readline().split()
        self.norm = (float(values[0]), float(values[1]), float(values[2]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        '''! Write Mirnov constraint to file
        
        @param file Open file object for saving constraints
        '''
        file.write('{0}\n'.format(Mirnov_con_id))
        file.write(' {0:E} {1:E} {2:E}\n'.format(self.loc[0], self.loc[1], self.phi))
        file.write(' {0:E} {1:E} {2:E}\n'.format(self.norm[0], self.norm[1], self.norm[2]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class Ip_con:
    '''! TokaMaker equilibrium reconstruction plasma current constraint'''
    def __init__(self, val=None, err=None):
        '''! Create plasma current constraint
        
        @param val Value of constraint
        @param err Error in constraint
        '''
        if val <= 0.0:
            raise ValueError("Plasma current constraint must be positive")
        self.val = val
        self.err = err

    def read(self, file):
        '''! Read plasma current constraint from file
        
        @param file Open file object containing constraint, must be positioned at start of constraint
        '''
        values = file.readline().split()
        Ip = float(values[0])
        err = float(values[1])
        if Ip <= 0.0:
            raise ValueError("Invalid value in file: Plasma current constraint must be positive")
        self.val = Ip
        self.err = 1./err

    def write(self, file):
        '''! Write plasma current constraint to file
        
        @param file Open file object for saving constraints
        '''
        file.write('{0}\n'.format(Ip_con_id))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class fluxLoop_con:
    '''! TokaMaker equilibrium reconstruction full flux loop constraint'''
    def __init__(self, loc=None, val=None, err=None):
        '''! Create full flux loop sensor constraint
        
        @param loc Location of Mirnov in R-Z plane [2]
        @param val Value of constraint
        @param err Error in constraint
        '''
        self.loc = loc
        self.val = val
        self.err = err

    def read(self, file):
        '''! Read full flux loop constraint from file
        
        @param file Open file object containing constraint, must be positioned at start of constraint
        '''
        values = file.readline().split()
        self.loc = (float(values[0]), float(values[1]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        '''! Write full flux loop constraint to file
        
        @param file Open file object for saving constraints
        '''
        file.write('{0}\n'.format(fluxLoop_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.loc[0], self.loc[1]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class dFlux_con:
    '''! TokaMaker equilibrium reconstruction diamagnetic flux constraint'''
    def __init__(self, val=None, err=None):
        '''! Create diamagnetic flux constraint
        
        @param val Value of constraint
        @param err Error in constraint
        '''
        self.val = val
        self.err = err

    def read(self, file):
        '''! Read diamagnetic flux constraint from file
        
        @param file Open file object containing constraint, must be positioned at start of constraint
        '''
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        '''! Write diamagnetic flux constraint to file
        
        @param file Open file object for saving constraints
        '''
        file.write('{0}\n'.format(dFlux_con_id))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class Press_con:
    '''! TokaMaker equilibrium reconstruction plasma pressure constraint'''
    def __init__(self, loc=None, val=None, err=None):
        '''! Create plasma pressure constraint
        
        @param loc Location of measurement in R-Z plane [2]
        @param val Value of pressure constraint
        @param err Error in constraint
        '''
        if val <= 0.0:
            raise ValueError("Plasma pressure constraints must be positive")
        self.loc = loc
        self.val = val
        self.err = err

    def read(self, file):
        '''! Read plasma pressure constraint from file
        
        @param file Open file object containing constraint, must be positioned at start of constraint
        '''
        values = file.readline().split()
        self.loc = (float(values[0]), float(values[1]))
        values = file.readline().split()
        pressure = float(values[0])
        if pressure <= 0.0:
            raise ValueError("Invalid value in file: Plasma pressure constraints must be positive")
        self.val = pressure
        self.err = 1./float(values[1])

    def write(self, file):
        '''! Write plasma pressure constraint to file
        
        @param file Open file object for saving constraints
        '''
        file.write('{0}\n'.format(Press_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.loc[0], self.loc[1]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class q_con:
    '''! TokaMaker equilibrium reconstruction local safety factor constraint'''
    def __init__(self, type=None, loc=0., val=None, err=None):
        r'''! Create local safety factor constraint
        
        @param type Type of constraint
        @param loc Location of constraint in \f$ \hat{\psi} \f$
        @param val Value of constraint
        @param err Error in constraint
        '''
        self.type = type
        self.val = val
        self.err = err
        self.loc = loc

    def read(self, file):
        '''! Read local safety factor constraint from file
        
        @param file Open file object containing constraint, must be positioned at start of constraint
        '''
        values = file.readline().split()
        self.type = int(values[0])
        self.loc = float(values[1])
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        '''! Write local safety factor constraint to file
        
        @param file Open file object for saving constraints
        '''
        file.write('{0}\n'.format(q_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.type, self.loc))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class saddle_con:
    '''! TokaMaker equilibrium reconstruction saddle flux loop constraint '''
    def __init__(self, pt1=None, pt2=None, width=None, val=None, err=None):
        '''! Create saddle flux loop sensor constraint
        
        @param pt1 Location of first toroidal saddle leg in R-Z plane [2]
        @param pt2 Location of second toroidal saddle leg in R-Z plane [2]
        @param width Toroidal extent in radians
        @param val Value of saddle loop constraint
        @param err Error in constraint
        '''
        self.pt1 = pt1
        self.pt2 = pt2
        self.width = width
        self.val = val
        self.err = err

    def read(self, file):
        '''! Read saddle flux loop constraint from file
        
        @param file Open file object containing constraint, must be positioned at start of constraint
        '''
        values = file.readline().split()
        self.pt1 = (float(values[0]), float(values[1]))
        values = file.readline().split()
        self.pt2 = (float(values[0]), float(values[1]))
        value = file.readline()
        self.width = float(value)
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        '''! Write saddle flux loop constraint to file
        
        @param file Open file object for saving constraints
        '''
        file.write('{0}\n'.format(saddle_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.pt1[0], self.pt1[1]))
        file.write(' {0:E} {1:E}\n'.format(self.pt2[0], self.pt2[1]))
        file.write(' {0:E}\n'.format(self.width))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class coil_current_con:
    '''! TokaMaker equilibrium reconstruction coil current constraint'''
    def __init__(self, ind=None, val=None, err=None):
        '''! Create coil current constraint
        
        @param ind Index of coil
        @param val Value of constraint [A]
        @param err Error in constraint [A]
        '''
        self.ind = ind
        self.val = val
        self.err = err

    def read(self, file):
        '''! Read coil current constraint from file
        
        @param file Open file object containing constraint, must be positioned at start of constraint
        '''
        values = file.readline()
        self.ind = int(values) - 1
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        '''! Write coil current constraint to file
        
        @param file Open file object for saving constraints
        '''
        file.write('{0}\n'.format(coil_current_con_id))
        file.write(' {0}\n'.format(self.ind+1))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./(self.err)))


con_map = {
    Mirnov_con_id: Mirnov_con,
    Ip_con_id: Ip_con,
    fluxLoop_con_id: fluxLoop_con,
    dFlux_con_id: dFlux_con,
    Press_con_id: Press_con,
    q_con_id: q_con,
    saddle_con_id: saddle_con,
    coil_current_con_id: coil_current_con
}


class tokamaker_recon_settings:
    r'''! TokaMaker reconstruction settings class'''
    def __init__(self):
        ## Adjust \f$ F*F' \f$ parameterization coefficients?
        self.fitF = False
        ## Adjust \f$ P' \f$ parameterization coefficients?
        self.fitP = False
        ## Adjust \f$ P' \f$ scale factor?
        self.fit_Pscale = False
        ## Adjust \f$ F*F' \f$ scale factor?
        self.fit_FFPscale = False
        ## Utilize and adjust \f$ R_0 \f$ constraint?
        self.fitR0 = False
        ## Utilize and adjust \f$ Z_0 \f$ constraint?
        self.fitZ0 = False
        ## Allow adjustment of PF coil currents?
        self.fitCoils = False
        ## Allow adjustment of TF coil current?
        self.fitF0 = False
        ## Do not update centering (initial guess for NL solve) as solve progresses
        self.fixedCentering = False
        ## Relative step size for finite difference Jacobian calculation
        self.dx = 1.E-2
        ## Minimum absolute step size for finite difference Jacobian calculation
        self.dx_min = 1.E-8
        ## Show detailed progress output?
        self.pm = False
        ## File containing constraint definitions
        self.infile = 'fit.in'
        ## File to write output information to
        self.outfile = 'fit.out'
        # Must be added last
        self._initialized = True

    def __setattr__(self, name, value):
        # Check if the attribute is being created for the first time after initialization
        if (name not in self.__dict__) and hasattr(self, '_initialized'):
            raise AttributeError(f"Cannot add new attribute '{name}' to tokamaker_recon_settings")
        super().__setattr__(name, value)
    
    @property
    def fitV0(self):
        warn(
            "`fitV0` is deprecated, use `fitZ0` instead. This function will be removed in a future version.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.fitZ0
    
    @fitV0.setter
    def fitV0(self, value):
        warn(
            "`fitV0` is deprecated, use `fitZ0` instead. This function will be removed in a future version.",
            DeprecationWarning,
            stacklevel=2
        )
        self.fitZ0 = value
    
    def get_c_struct(self,oft_env):
        r'''! Get C struct representation of settings for passing to TokaMaker Fortran API'''
        c_struct_instance = tokamaker_recon_settings_cstruct()
        c_struct_instance.fitF = self.fitF
        c_struct_instance.fitP = self.fitP
        c_struct_instance.fit_Pscale = self.fit_Pscale
        c_struct_instance.fit_FFPscale = self.fit_FFPscale
        c_struct_instance.fitR0 = self.fitR0
        c_struct_instance.fitZ0 = self.fitZ0
        c_struct_instance.fitCoils = self.fitCoils
        c_struct_instance.fitF0 = self.fitF0
        c_struct_instance.fixedCentering = self.fixedCentering
        c_struct_instance.pm = self.pm
        c_struct_instance.infile = oft_env.path2c(self.infile)
        c_struct_instance.outfile = oft_env.path2c(self.outfile)
        return c_struct_instance


class reconstruction():
    '''! TokaMaker equilibrium reconstruction class'''
    def __init__(self,tMaker_obj,in_filename='fit.in',out_filename='fit.out'):
        '''! Create equilibrium reconstruction object
        
        @param tMaker_obj TokaMaker object used for computing G-S equilibria
        @param in_filename Filename to use for reconstruction input
        @param out_filename Filename to use for reconstruction outputs
        '''
        ## Grad-Shafranov object for reconstruction
        self._tMaker_obj = tMaker_obj
        ## Reconstruction specific settings object
        self.settings = tokamaker_recon_settings()
        ## Plasma current constraint
        self._Ip_con = None
        ## Diamagnetic flux constraint 
        self._Dflux_con = None
        ## Flux loop constraints
        self._flux_loops = []
        ## Mirnov sensor constraints
        self._mirnovs = []
        ## Saddle loop constraints
        self._saddles = []
        ## Plasma pressure constraints
        self._pressure_cons = []
        ## Coil current constraints
        self._coil_current_cons = []
        ## Name of constraint file (input for reconstruction)
        self.con_file = in_filename
        ## Name of reconstruction output file
        self.out_file = out_filename
        ## Number of constraints in underlying Fortran object
        self._ncons = 0
        # Update settings
        self.settings.infile = self.con_file
        self.settings.outfile = self.out_file
        # Fit-specific input file settings
        self._tMaker_obj._oft_env.oft_in_groups['gs_fit_options'] = {
            'ftol': '1.E-3',
            'xtol': '1.E-3',
            'gtol': '1.E-3',
            'maxfev': '100',
            'epsfcn': '1.E-3',
            'factor': '1.0',
            'comp_var': 'F',
            'linearized_fit': 'F'
        }
        self._tMaker_obj._oft_env.update_oft_in()
    
    def __del__(self):
        '''! Destroy reconstruction object'''
        try:
            self.destroy_constraints()
        except Exception:
            pass
        self._tMaker_obj = None
        self.settings = None
        self._Ip_con = None
        self._Dflux_con = None
        self._flux_loops = []
        self._mirnovs = []
        self._saddles = []
        self._pressure_cons = []
        self._coil_current_cons = []
        self.con_file = None
        self.out_file = None
        self._ncons = 0
    
    def set_Ip(self,Ip,err):
        '''! Set plasma current constraint
        
        @param Ip Value of plasma current constraint
        @param err Error in constraint
        '''
        self._Ip_con = Ip_con(val=Ip, err=err)

    def set_DFlux(self,DFlux,err):
        '''! Set diamagnetic flux constraint
        
        @param DFlux Value of diamagnetic flux constraint
        @param err Error in constraint
        '''
        self._Dflux_con = dFlux_con(val=DFlux, err=err)
    
    def set_coil_currents(self,targets,errs):
        '''! Set coil current constraints
        
        @param targets Dictionary of coil current targets
        @param errs Dictionary of coil current constraint errors
        '''
        self._coil_current_cons = []
        for name, coil_set in self._tMaker_obj.coil_sets.items():
            if name not in targets:
                raise KeyError(f"Missing coil current target for coil set '{name}'")
            if name not in errs:
                raise KeyError(f"Missing coil current error for coil set '{name}'")
            self._coil_current_cons.append(coil_current_con(ind=coil_set['id'], val=targets[name], err=errs[name]))

    def add_flux_loop(self,loc,val,err):
        '''! Add full poloidal flux loop constraint
        
        @param loc Location of flux loop in R-Z plane [2]
        @param val Value of flux loop constraint
        @param err Error in constraint
        '''
        self._flux_loops.append(fluxLoop_con(loc=loc, val=val, err=err))

    def add_Mirnov(self,loc,norm,val,err):
        '''! Add Mirnov sensor constraint
        
        @param loc Location of Mirnov in R-Z plane [2]
        @param norm Unit normal in R-Z plane [2]
        @param val Value of Mirnov constraint
        @param err Error in constraint
        '''
        self._mirnovs.append(Mirnov_con(loc=loc, norm=norm, val=val, err=err))
    
    def add_saddle(self,p1,p2,width,val,err):
        '''! Add saddle loop constraint
        
        @param p1 Location of first toroidal saddle leg in R-Z plane [2]
        @param p2 Location of second toroidal saddle leg in R-Z plane [2]
        @param width Toroidal extent in radians
        @param val Value of saddle loop constraint
        @param err Error in constraint
        '''
        self._saddles.append(saddle_con(p1=p1, p2=p2, width=width, val=val, err=err))
    
    def add_pressure(self,loc,val,err):
        '''! Add plasma pressure constraint
        
        @param loc Location of measurement in R-Z plane [2]
        @param val Value of pressure constraint
        @param err Error in constraint
        '''
        self._pressure_cons.append(Press_con(loc=loc,val=val,err=err))
    
    def reset_constraints(self):
        '''! Remove all current constraints'''
        self._Ip_con = None
        self._Dflux_con = None
        self._flux_loops = []
        self._mirnovs = []
        self._saddles = []
        self._pressure_cons = []
    
    def write_fit_in(self):
        '''! Create reconstruction input file for specified constraints'''
        constraints = self._flux_loops + self._mirnovs + self._saddles + self._pressure_cons + self._coil_current_cons
        if self._Ip_con is not None:
            constraints.append(self._Ip_con)
        if self._Dflux_con is not None:
            constraints.append(self._Dflux_con)
        ncons = len(constraints)
        with open(self.con_file, 'w+') as fid:
            fid.write('{0:d}\n\n'.format(ncons))
            for con in constraints:
                con.write(fid)
    
    def read_fit_in(self):
        '''! Read constraints from existing input file'''
        self.reset_constraints()
        with open(self.con_file, 'r') as fid:
            ncons = int(fid.readline())
            for _ in range(ncons):
                fid.readline()
                con_type = int(fid.readline())
                new_con_class = con_map[con_type]
                new_con = new_con_class()
                new_con.read(fid)
                if con_type == Ip_con_id:
                    self._Ip_con = new_con
                elif con_type == dFlux_con_id:
                    self._Dflux_con = new_con
                elif con_type == fluxLoop_con_id:
                    self._flux_loops.append(new_con)
                elif con_type == Mirnov_con_id:
                    self._mirnovs.append(new_con)
                elif con_type == saddle_con_id:
                    self._saddles.append(new_con)
                elif con_type == Press_con_id:
                    self._pressure_cons.append(new_con)
                elif con_type == coil_current_con_id:
                    self._coil_current_cons.append(new_con)
                else:
                    raise ValueError("Unknown constraint type")
                
    def plot_constraints(self, fig, ax, plot_flux_loops=True, plot_mirnovs=True, plot_saddles=True, plot_pressure=True, mirnov_scale=1.0, base_zorder=0):
        '''! Plot constraint locations in R-Z plane for current constraints
        
        @param fig Matplotlib figure to plot on
        @param ax Matplotlib axis to plot on
        @param plot_flux_loops Whether to plot flux loop constraints
        @param plot_mirnovs Whether to plot Mirnov sensor constraints
        @param plot_saddles Whether to plot saddle loop constraints
        @param plot_pressure Whether to plot plasma pressure constraints
        @param mirnov_scale Scaling factor for Mirnov sensor quiver plot (default 1.0, smaller values increase arrow size)
        @param base_zorder Base z-order for plotting constraints (default 0.0, higher values will plot on top of lower values)
        '''
        if plot_flux_loops:
            for flux_loop in self._flux_loops:
                ax.plot(flux_loop.loc[0], flux_loop.loc[1], 'ro', zorder=base_zorder+1)
        mirnov_scale = 1.0/mirnov_scale/(0.005*(self._tMaker_obj.r.max(axis=0) - self._tMaker_obj.r.min(axis=0)).max())
        if plot_mirnovs:
            for mirnov in self._mirnovs:
                ax.quiver(mirnov.loc[0], mirnov.loc[1], mirnov.norm[0], mirnov.norm[2], color='g', scale=mirnov_scale, zorder=base_zorder+2)
        if plot_saddles:
            for saddle in self._saddles:
                ax.plot([saddle.pt1[0],saddle.pt2[0]], [saddle.pt1[1],saddle.pt2[1]], 'b+', zorder=base_zorder+1)
        if plot_pressure:
            for press_con in self._pressure_cons:
                ax.plot(press_con.loc[0], press_con.loc[1], 'm.', zorder=base_zorder+1)
    
    def plot_error(self, ax, error_mat=None, error_file=None, coil_ax=None, chi_ax=None):
        '''! Plot constraint values and reconstructed signals for current equilibrium from reconstruction output file
        
        @param ax Matplotlib axis to plot on
        @param error_mat Error matrix, if `None` `eval_error` will be called
        @param error_file File containing error data
        @param coil_ax Matplotlib axis to plot coil current constraints on
        @param chi_ax Matplotlib axis to plot chi-squared values on
        '''
        if len(ax) != 4:
            raise ValueError("Must provide list of 4 axes for plotting flux loop, Mirnov, saddle, and pressure constraint errors")
        if error_file is not None:
            error_mat = numpy.loadtxt(error_file)[:,1:]
        if error_mat is None:
            error_mat = self.eval_error()
        else:
            if error_mat.shape[1] != 4:
                raise ValueError("Error matrix must have 4 columns (weighted error, reconstructed signal, constraint value, 3D correction)")
            if error_mat.shape[0] != self._ncons:
                raise ValueError("Error matrix must have same number of rows as number of constraints in reconstruction")
        # Plot signals for "main" diagnostics
        i = 0
        ax[0].set_title('Flux Loop Constraints')
        ax[0].set_xlabel('Constraint Index')
        ax[0].set_ylabel('Signal [Wb]')
        for j, flux_loop in enumerate(self._flux_loops):
            ax[0].errorbar(j,error_mat[i,2], yerr=flux_loop.err, color='r', capsize=2)
            ax[0].plot(j,error_mat[i,1], 'bx')
            i += 1
        ax[1].set_title('Mirnov Constraints')
        ax[1].set_xlabel('Constraint Index')
        ax[1].set_ylabel('Signal [T]')
        for j, mirnov in enumerate(self._mirnovs):
            ax[1].errorbar(j,error_mat[i,2], yerr=mirnov.err, color='r', capsize=2)
            ax[1].plot(j,error_mat[i,1], 'bx')
            i += 1
        ax[2].set_title('Saddle Constraints')
        ax[2].set_xlabel('Constraint Index')
        ax[2].set_ylabel('Signal [Wb]')
        for j, saddle in enumerate(self._saddles):
            ax[2].errorbar(j,error_mat[i,2], yerr=saddle.err, color='r', capsize=2)
            ax[2].plot(j,error_mat[i,1], 'bx')
            i += 1
        ax[3].set_title('Pressure Constraints')
        ax[3].set_xlabel('Constraint Index')
        ax[3].set_ylabel('Signal [Pa]')
        for j, press_con in enumerate(self._pressure_cons):
            ax[3].errorbar(j,error_mat[i,2], yerr=press_con.err, color='r', capsize=2)
            ax[3].plot(j,error_mat[i,1], 'bx')
            i += 1
        # Plot coil signals
        if (coil_ax is not None) and (len(self._coil_current_cons) > 0):
            ind_to_name = {self._tMaker_obj.coil_sets[coil_set]['id']: coil_set for coil_set in self._tMaker_obj.coil_sets}
            coil_ax.set_title('Coil Current Constraints')
            coil_ax.set_xlabel('Coil Set Index')
            coil_ax.set_ylabel('Current [A]')
            for coil_con in self._coil_current_cons:
                coil_ax.errorbar(coil_con.ind, coil_con.val, yerr=coil_con.err, color='r', capsize=2)
                coil_ax.plot(coil_con.ind, error_mat[i,1], 'bx')
                i += 1
            coil_ax.set_xticks(range(len(ind_to_name)), labels=[ind_to_name[j] for j in sorted(ind_to_name.keys())],
                rotation=45, ha="right", rotation_mode="anchor")
        # Plot error contributions
        if chi_ax is not None:
            err_ind = numpy.cumsum([0,
                                    len(self._flux_loops),
                                    len(self._mirnovs),
                                    len(self._saddles),
                                    len(self._pressure_cons),
                                    len(self._coil_current_cons),
                                    1 if self._Ip_con is not None else 0,
                                    1 if self._Dflux_con is not None else 0])
            names = ['Flux loops', 'Mirnovs', 'Saddle loops', 'Pressure', r'$I_C$', r'$I_p$', r'$\Delta \Phi$']
            chi_ax.set_title(r'Signal $\chi^2$ contributions')
            chi_ax.set_xlabel('Signal Index')
            chi_ax.set_ylabel(r'$\chi^2_i$')
            for i in range(err_ind.shape[0]-1):
                if err_ind[i+1] == err_ind[i]:
                    continue
                chi_ax.plot(numpy.arange(err_ind[i],err_ind[i+1])+1,numpy.power(error_mat[err_ind[i]:err_ind[i+1],0],2),label=names[i])
            chi_ax.legend()
    
    def setup_constraints(self):
        '''! Set up constraints in TokaMaker for current equilibrium without performing reconstruction
        
        @result Error flag
        '''
        if self._ncons > 0:
            self.destroy_constraints()
        # Modify input file
        self.write_fit_in()
        # Run setup
        ncons = c_int()
        error_flag = c_int()
        tokamaker_recon_setup(self._tMaker_obj._tMaker_ptr,self.settings.get_c_struct(self._tMaker_obj._oft_env),ctypes.byref(ncons),ctypes.byref(error_flag))
        if error_flag.value != 0:
            raise ValueError("Constraint setup failed with error code {0:d}".format(error_flag.value))
        self._ncons = ncons.value
    
    def destroy_constraints(self):
        '''! Destroy constraints in TokaMaker'''
        if self._ncons == 0:
            return
        error_flag = c_int()
        tokamaker_recon_destroy(self._tMaker_obj._tMaker_ptr,ctypes.byref(error_flag))
        if error_flag.value != 0:
            raise ValueError("Constraint destruction failed with error code {0:d}".format(error_flag.value))
        self._ncons = 0

    def eval_error(self, vacuum=False, save_to_file=False):
        '''! Evaluate error in current equilibrium for specified constraints without performing reconstruction
        
        @param vacuum Perform vacuum reconstruction
        @param save_to_file Save error matrix to file specified in `settings.outfile`
        @result Error matrix if `save_to_file=False`, otherwise None (error matrix will be saved to file)
        '''
        if self._ncons == 0:
            raise ValueError("Constraints are not set up; call `setup_constraints()` before `eval_error()`. ")
        mat_ptr = None
        if not save_to_file:
            error_mat = numpy.zeros((self._ncons, 4), dtype=numpy.float64)
            mat_ptr = error_mat.ctypes.data_as(c_double_ptr)
        error_flag = c_int()
        tokamaker_recon_err(self._tMaker_obj._tMaker_ptr,c_bool(vacuum),self.settings.get_c_struct(self._tMaker_obj._oft_env),mat_ptr,ctypes.byref(error_flag))
        if error_flag.value != 0:
            raise ValueError("Error evaluation failed with error code {0:d}".format(error_flag.value))
        return error_mat if not save_to_file else None
    
    def setup_get_opt(self,fail_factor=1.E2,pp_target_weight=1.E5,opoint_target_weight=1.E2):
        r''' Setup object for use with general optimization routines, such as `scipy.optimize.minimize`,
        which will call `opt_error` and, optionally, `opt_error_jacobian` with appropriate arguments.
        This function should be called after setting up constraints with `setup_constraints` and before calling optimization routines.

        @param fail_factor Factor to multiply initial error by for failed solves
        @param pp_target_weight Weight for global P' targets when treated as soft constraints and fit_Pscale is True
        @param opoint_target_weight Weight for R0 and Z0 targets when treated as soft constraints and `fitR0` or `fitZ0` is True
        @returns Initial DoF vector for optimization, error vector for initial DoFs
        '''
        x0_nl = self.opt_get_dofs()
        self._err_min = 1.E99
        self._EQ_center = self._tMaker_obj.copy_eq()
        error0  = self.opt_error(x0_nl,self,False)
        self._fail_error = error0*fail_factor

        # Make sure settings are up to date
        soft_targets_active = any([self._tMaker_obj.settings.ffp_target_weight > 0.0, self._tMaker_obj.settings.pp_target_weight > 0.0, self._tMaker_obj.settings.opoint_target_weight > 0.0])
        if (self.settings.fitR0 or self.settings.fitZ0) and (soft_targets_active):
            self._tMaker_obj.settings.opoint_target_weight=opoint_target_weight
        if self.settings.fit_Pscale and soft_targets_active:
            self._tMaker_obj.settings.pp_target_weight=pp_target_weight
        self._tMaker_obj.update_settings()

        return x0_nl, error0
    
    def opt_get_dofs(self):
        r''' Get DoF vector for general optimization routines using `opt_error` and `opt_error_jacobian`.

        @returns DoF vector
        '''
        dofs = []
        if self.settings.fitR0:
            dofs.append(self._tMaker_obj.o_point[0])
        if self.settings.fitZ0:
            dofs.append(self._tMaker_obj.o_point[1])
        if self.settings.fit_Pscale:
            dofs.append(self._tMaker_obj.p_scale)
        if self.settings.fitF:
            f_dofs = self._tMaker_obj.get_profile_dofs('ffp')
            if f_dofs is not None:
                dofs.extend(f_dofs.tolist())
        if self.settings.fitP:
            p_dofs = self._tMaker_obj.get_profile_dofs('pp')
            if p_dofs is not None:
                dofs.extend(p_dofs.tolist())
        return numpy.array(dofs)

    @staticmethod
    def opt_error(cofs,recon_obj,in_jac):
        r''' Compute weighted constraint error for given DoF vector `cofs` and reconstruction object `recon_obj`.
        This function is designed to be called by general optimization routines, such as `scipy.optimize.minimize`.

        @param cofs DoF vector
        @param recon_obj `TokaMaker.reconstruction.reconstruction` object to be used
        @param in_jac Flag indicating whether this function is being called from a Jacobian calculation to make print statements during finite differencing
        @returns Weighted error vector for current DoFs
        '''
        offset = 0
        if recon_obj.settings.fitR0:
            recon_obj._tMaker_obj.set_targets(R0=cofs[offset],retain_previous=True)
            offset += 1
        if recon_obj.settings.fitZ0:
            recon_obj._tMaker_obj.set_targets(Z0=cofs[offset],retain_previous=True)
            offset += 1
        if recon_obj.settings.fit_Pscale:
            recon_obj._tMaker_obj.p_scale=cofs[offset]
            offset += 1
        if recon_obj.settings.fitF:
            f_dofs = recon_obj._tMaker_obj.get_profile_dofs('ffp')
            if f_dofs is not None:
                f_dofs = cofs[offset:offset+f_dofs.shape[0]]
                recon_obj._tMaker_obj.set_profile_dofs('ffp', f_dofs)
                offset += f_dofs.shape[0]
        if recon_obj.settings.fitP:
            p_dofs = recon_obj._tMaker_obj.get_profile_dofs('pp')
            if p_dofs is not None:
                p_dofs = cofs[offset:offset+p_dofs.shape[0]]
                recon_obj._tMaker_obj.set_profile_dofs('pp', p_dofs)
                offset += p_dofs.shape[0]
        
        # Re-solve
        try:
            recon_obj._tMaker_obj.solve()
            if recon_obj._tMaker_obj.alam == 0.0:
                recon_obj._tMaker_obj.replace_eq(recon_obj._EQ_center)
                return recon_obj._fail_error
        except:
            recon_obj._tMaker_obj.replace_eq(recon_obj._EQ_center)
            return recon_obj._fail_error
        
        # Compute error
        err_out = recon_obj.eval_error()
        total_err = numpy.power(numpy.linalg.norm(err_out[:,0]),2)
        if not in_jac:
            print('chi_rms = {0:.4E}; chi_max = {1:.4E}'.format(numpy.sqrt(total_err/err_out.shape[0]), abs(err_out[:,0]).max()))
            if total_err < recon_obj._err_min:
                if not recon_obj.settings.fixedCentering:
                    recon_obj._EQ_center = recon_obj._tMaker_obj.copy_eq()
                recon_obj._err_min = total_err
        return err_out[:,0]


    @staticmethod
    def opt_error_jacobian(cofs,recon_obj,in_jac):
        r''' Compute Jacobian matrix corresponding by differencing `opt_error` for given DoF vector `cofs` and reconstruction object `recon_obj`.

        @param cofs DoF vector
        @param recon_obj `TokaMaker.reconstruction.reconstruction` object to be used
        @param in_jac Ignored flag just to match signature with `opt_error`
        @returns Jacobian matrix for current DoFs
        '''
        recon_obj._tMaker_obj.replace_eq(recon_obj._EQ_center)
        _ = reconstruction.opt_error(cofs,recon_obj,True)
        jac = numpy.zeros((len(recon_obj._fail_error),len(cofs)))
        center_EQ = recon_obj._tMaker_obj.copy_eq()
        # Resolve for center point to capture "one more step" converged result as achieved in small diffs below
        center_err = reconstruction.opt_error(cofs,recon_obj,True)
        for i, cof_val in enumerate(cofs):
            cof_tmp = cofs.copy()
            cof_diff = max(recon_obj.settings.dx_min,recon_obj.settings.dx*abs(cof_val))
            cof_tmp[i] += cof_diff
            pert_err = reconstruction.opt_error(cof_tmp,recon_obj,True)
            jac[:,i] = (pert_err-center_err)/cof_diff
            recon_obj._tMaker_obj.replace_eq(center_EQ)
        return jac

    def reconstruct(self, vacuum=False, linearized_fit=False, maxits=100, eps=1.E-3, ftol=1.E-3, xtol=1.E-3, gtol=1.E-3):
        '''! Reconstruct G-S equation with specified fitting constraints, profiles, etc.
        
        @param vacuum Perform vacuum reconstruction
        @param linearized_fit Use linearized solve for suitable terms
        @param maxits Maximum number of iterations
        @param eps Epsilong factor for finite difference derivative calculations
        @param ftol Stopping condition: termination occurs when both the actual and predicted relative reductions in the sum of squares are at most `ftol`
        @param xtol Stopping condition: termination occurs when the relative error between two consecutive iterates is at most `xtol`
        @param gtol Stopping condition: termination occurs when the cosine of the angle between fvec and any column of the jacobian is at most `gtol` in absolute value
        @result Error flag
        '''
        # Check for possibly conflicting constraints
        if self._tMaker_obj._tMaker_equil.Isoflux_constraints is not None:
            oft_warning('Removing conflicting isoflux constraints from equilibrium object via `.set_isoflux_constraints(None)`')
            self._tMaker_obj.set_isoflux_constraints(None)
        if self._tMaker_obj._tMaker_equil.Psi_constraints[0] is not None:
            oft_warning('Removing conflicting Psi constraints from equilibrium object via `.set_psi_constraints(None,None)`')
            self._tMaker_obj.set_psi_constraints(None,None)
        if self._tMaker_obj._tMaker_equil.Saddle_constraints is not None:
            oft_warning('Removing conflicting saddle targets from equilibrium object via `.set_saddle_constraints(None)`')
            self._tMaker_obj.set_saddle_constraints(None)
        weight_chk = [self._tMaker_obj.settings.ffp_target_weight > 0.0, self._tMaker_obj.settings.pp_target_weight > 0.0, self._tMaker_obj.settings.opoint_target_weight > 0.0]
        if any(weight_chk):
            oft_warning('Removing conflicting soft targets from equilibrium solve via `.settings.*_target_weight = -1.0; .update_settings()`')
            self._tMaker_obj.settings.ffp_target_weight = -1.0
            self._tMaker_obj.settings.pp_target_weight = -1.0
            self._tMaker_obj.settings.opoint_target_weight = -1.0
            self._tMaker_obj.update_settings()
        # Modify input file
        self.write_fit_in()
        self._tMaker_obj._oft_env.oft_in_groups['gs_fit_options']['linearized_fit'] = 'T' if linearized_fit else 'F'
        self._tMaker_obj._oft_env.oft_in_groups['gs_fit_options']['maxfev'] = '{0:d}'.format(maxits)
        self._tMaker_obj._oft_env.oft_in_groups['gs_fit_options']['epsfcn'] = '{0:.5E}'.format(eps)
        self._tMaker_obj._oft_env.oft_in_groups['gs_fit_options']['ftol'] = '{0:.5E}'.format(ftol)
        self._tMaker_obj._oft_env.oft_in_groups['gs_fit_options']['xtol'] = '{0:.5E}'.format(xtol)
        self._tMaker_obj._oft_env.oft_in_groups['gs_fit_options']['gtol'] = '{0:.5E}'.format(gtol)
        self._tMaker_obj._oft_env.update_oft_in()
        # Run reconstruction
        error_flag = c_int()
        tokamaker_recon_run(self._tMaker_obj._tMaker_ptr,c_bool(vacuum),self.settings.get_c_struct(self._tMaker_obj._oft_env),ctypes.byref(error_flag))
        return error_flag.value