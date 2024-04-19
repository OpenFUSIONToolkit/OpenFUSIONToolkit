'''! Python interface for TokaMaker Grad-Shafranov functionality

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''
import numpy as np
from ..util import *

class tokamaker_recon_settings_struct(c_struct):
    r'''! TokaMaker reconstruction settings structure

     - `pm` Print 'performance' information (eg. iteration count) during run?
     - `free_boundary` Perform free-boundary calculation?
     - `has_plasma` Include plasma effects in calculation, vacuum otherwise?
     - `limited_only` Do not search for X-points when determining LCFS?
     - `maxits` Maximum NL iteration count for G-S solver
     - `mode` Parallel current source formulation used (0 -> define \f$F'\f$, 1 -> define \f$F*F'\f$)
     - `urf` Under-relaxation factor for NL fixed-point iteration
     - `nl_tol` Convergence tolerance for NL solver
     - `rmin` Minimum magnetic axis major radius, used to catch 'lost' equilibria
     - `lim_zmax` Maximum vertical range for limiter points, can be used to exclude complex diverter regions
     - `limiter_file` File containing additional limiter points not included in mesh (default: 'none')
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
                ("pm", c_bool)]


def tokamaker_recon_default_settings():
    '''! Initialize reconstruction settings object with default values

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
    return settings

## @cond
tokamaker_recon_run = ctypes_subroutine(oftpy_lib.tokamaker_recon_run,
    [c_bool, ctypes.POINTER(tokamaker_recon_settings_struct), c_int_ptr])
## @endcond

Mirnov_con_id = 1
Ip_con_id = 2
fluxLoop_con_id = 7
dFlux_con_id = 8
Press_con_id = 9
q_con_id = 10
saddle_con_id = 11


class Mirnov_con:
    def __init__(self, pt=None, phi=0., norm=None, val=None, err=None):
        self.pt = pt
        self.phi = phi
        self.norm = norm
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt = (float(values[0]), float(values[1]))
        self.phi = float(values[2])
        values = file.readline().split()
        self.norm = (float(values[0]), float(values[1]), float(values[2]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(Mirnov_con_id))
        file.write(' {0:E} {1:E} {2:E}\n'.format(self.pt[0], self.pt[1], self.phi))
        file.write(' {0:E} {1:E} {2:E}\n'.format(self.norm[0], self.norm[1], self.norm[2]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class Ip_con:
    def __init__(self, val=None, err=None):
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(Ip_con_id))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class fluxLoop_con:
    def __init__(self, pt=None, val=None, err=None):
        self.pt = pt
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt = (float(values[0]), float(values[1]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(fluxLoop_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.pt[0], self.pt[1]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class dFlux_con:
    def __init__(self, val=None, err=None):
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(dFlux_con_id))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class Press_con:
    def __init__(self, pt=None, val=None, err=None):
        self.pt = pt
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt = (float(values[0]), float(values[1]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(Press_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.pt[0], self.pt[1]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class q_con:
    def __init__(self, type=None, val=None, err=None, loc=0.):
        self.type = type
        self.val = val
        self.err = err
        self.loc = loc

    def read(self, file):
        values = file.readline().split()
        self.type = int(values[0])
        self.loc = float(values[1])
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(q_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.type, self.loc))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class saddle_con:
    def __init__(self, pt1=None, pt2=None, width=None, val=None, err=None):
        self.pt1 = pt1
        self.pt2 = pt2
        self.width = width
        self.val = val
        self.err = err

    def read(self, file):
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
    def __init__(self,gs_obj,filename=None):
        ## Grad-Shafranov object for reconstruction
        self._gs_obj = gs_obj
        ## Reconstruction specific settings object
        self.settings = tokamaker_recon_default_settings()
        ##
        self._Ip_con = None
        ##
        self._Dflux_con = None
        ##
        self._flux_loops = []
        ##
        self._mirnovs = []
        ##
        self._saddles = []
        ##
        self._pressure_cons = []
        #
        if filename is not None:
            self.read_fit_in(filename)
    
    def __del__(self):
        self._gs_obj = None
        self.settings = None
        self._Ip_con = None
        self._Dflux_con = None
        self._flux_loops = []
        self._mirnovs = []
        self._saddles = []
        self._pressure_cons = []
    
    def set_Ip(self,Ip,err):
        self._Ip_con = Ip_con(val=Ip, err=err)

    def set_DFlux(self,DFlux,err):
        self._Dflux_con = dFlux_con(val=DFlux, err=err)

    def add_flux_loop(self,loc,val,err):
        self._flux_loops.append(fluxLoop_con(pt=loc, val=val, err=err))

    def add_Mirnov(self,loc,norm,val,err):
        self._mirnovs.append(Mirnov_con(pt=loc, norm=norm, val=val, err=err))
    
    def add_saddle(self,p1,p2,width,val,err):
        self._saddles.append(saddle_con(p1=p1, p2=p2, width=width, val=val, err=err))
    
    def add_pressure(self,loc,val,err):
        self._pressure_cons.append(Press_con(pt=loc,val=val,err=err))
    
    def reset_constraints(self):
        self._Ip_con = None
        self._Dflux_con = None
        self._flux_loops = []
        self._mirnovs = []
        self._saddles = []
        self._pressure_cons = []
    
    def write_fit_in(self,filename='fit.in'):
        constraints = self._flux_loops + self._mirnovs + self._pressure_cons
        if self._Ip_con is not None:
            constraints.append(self._Ip_con)
        if self._Dflux_con is not None:
            constraints.append(self._Dflux_con)
        ncons = len(constraints)
        with open(filename, 'w+') as fid:
            fid.write('{0:d}\n\n'.format(ncons))
            for con in constraints:
                con.write(fid)
    
    def read_fit_in(self,filename='fit.in'):
        self.reset_constraints()
        with open(filename, 'r') as fid:
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
        tokamaker_recon_run(c_bool(vacuum),self.settings,ctypes.byref(error_flag))
        return error_flag.value