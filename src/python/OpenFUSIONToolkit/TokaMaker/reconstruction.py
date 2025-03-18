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
from .._interface import *

class tokamaker_recon_settings_struct(c_struct):
    r'''! TokaMaker reconstruction settings structure

     - `fitI` Adjust \f$ F*F' \f$ parameterization coefficients?
     - `fitP` Adjust \f$ P' \f$ parameterization coefficients?
     - `fitPnorm` Adjust \f$ P' \f$ scale factor?
     - `fitAlam` Adjust \f$ F*F' \f$ scale factor?
     - `fitR0` Utilize and adjust \f$ R_0 \f$ constraint?
     - `fitV0` Utilize and adjust \f$ Z_0 \f$ constraint?
     - `fitCoils` Allow adjustment of PF coil currents?
     - `fitF0` Allow adjustment of TF coil current?
     - `fixedCentering` Do not update centering (initial guess for NL solve) as solve progresses
     - `pm` Show detailed progress output?
     - `infile` File containing constraint definitions
     - `outfile` File to write output information to
    '''
    _fields_ = [("fitI", c_bool),
                ("fitP", c_bool),
                ("fitPnorm", c_bool),
                ("fitAlam", c_bool),
                ("fitR0", c_bool),
                ("fitV0", c_bool),
                ("fitCoils", c_bool),
                ("fitF0", c_bool),
                ("fixedCentering", c_bool),
                ("pm", c_bool),
                ("infile", ctypes.c_char_p),
                ("outfile", ctypes.c_char_p)]


def tokamaker_recon_default_settings(oft_env):
    '''! Initialize reconstruction settings object with default values

    @param oft_env OFT runtime environment 
    @result tokamaker_recon_settings_struct object
    '''
    settings = tokamaker_recon_settings_struct()
    settings.fitI = False
    settings.fitP = False
    settings.fitPnorm = True
    settings.fitAlam = True
    settings.fitR0 = False
    settings.fitV0 = False
    settings.fitCoils = False
    settings.fitF0 = False
    settings.fixedCentering = False
    settings.pm = False
    settings.infile = oft_env.path2c('fit.in')
    settings.outfile = oft_env.path2c('fit.out')
    return settings

## @cond
# tokamaker_recon_run(tMaker_ptr,vacuum,settings,error_flag)
tokamaker_recon_run = ctypes_subroutine(oftpy_lib.tokamaker_recon_run,
    [c_void_p, c_bool, ctypes.POINTER(tokamaker_recon_settings_struct), c_int_ptr])
## @endcond

Mirnov_con_id = 1
Ip_con_id = 2
fluxLoop_con_id = 7
dFlux_con_id = 8
Press_con_id = 9
q_con_id = 10
saddle_con_id = 11


class Mirnov_con:
    '''! TokaMaker equilibrium reconstruction Mirnov sensor constraint'''
    def __init__(self, loc=None, phi=0., norm=None, val=None, err=None):
        '''! Create Mirnov sensor constraint
        
        @param loc Location of Mirnov in R-Z plane [2]
        @param phi Toroidal location [rad] (only meaningful with 3D fields)
        @param norm Unit normal in R-Z plane [2]
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
        self.val = val
        self.err = err

    def read(self, file):
        '''! Read plasma current constraint from file
        
        @param file Open file object containing constraint, must be positioned at start of constraint
        '''
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

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
        self.val = float(values[0])
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


con_map = {
    Mirnov_con_id: Mirnov_con,
    Ip_con_id: Ip_con,
    fluxLoop_con_id: fluxLoop_con,
    dFlux_con_id: dFlux_con,
    Press_con_id: Press_con,
    q_con_id: q_con,
    saddle_con_id: saddle_con
}


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
        self.settings = tokamaker_recon_default_settings(self._tMaker_obj._oft_env)
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
        ## Name of constraint file (input for reconstruction)
        self.con_file = in_filename
        ## Name of reconstruction output file
        self.out_file = out_filename
        # Update settings
        self.settings.infile = self._tMaker_obj._oft_env.path2c(self.con_file)
        self.settings.outfile = self._tMaker_obj._oft_env.path2c(self.out_file)
    
    def __del__(self):
        '''! Destroy reconstruction object'''
        self._tMaker_obj = None
        self.settings = None
        self._Ip_con = None
        self._Dflux_con = None
        self._flux_loops = []
        self._mirnovs = []
        self._saddles = []
        self._pressure_cons = []
    
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
        constraints = self._flux_loops + self._mirnovs + self._saddles + self._pressure_cons
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
                else:
                    raise ValueError("Unknown constraint type")

    def reconstruct(self, vacuum=False):
        '''! Reconstruct G-S equation with specified fitting constraints, profiles, etc.'''
        self.write_fit_in()
        error_flag = c_int()
        tokamaker_recon_run(self._tMaker_obj._tMaker_ptr,c_bool(vacuum),self.settings,ctypes.byref(error_flag))
        return error_flag.value