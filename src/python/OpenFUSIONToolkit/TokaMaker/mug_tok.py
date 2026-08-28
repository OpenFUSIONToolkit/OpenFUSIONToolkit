#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Python interface for the coupled MUG/TokaMaker (mugtok_td) time-dependent solver

An extension of @ref OpenFUSIONToolkit.TokaMaker._core.TokaMaker "TokaMaker" that couples
the MUG 2D MHD solver to an existing TokaMaker equilibrium (see @ref MUGToksim).

@authors Sophia Guizzo
@date July 2026
@ingroup doxy_oft_python
'''
import ctypes
import numpy
from .._interface import *

## @cond

# mugtok_alloc(mugtok_ptr,equil_ptr,error_str)
mugtok_alloc = ctypes_subroutine(oftpy_lib.mugtok_alloc,
    [c_void_ptr_ptr, c_void_p, c_char_p])

# mugtok_setup(mugtok_ptr,dt,lin_tol,nl_tol,mhd_flag,dens,visc,nreg,incomp,toroidal_flow,error_str)
mugtok_setup = ctypes_subroutine(oftpy_lib.mugtok_setup,
    [c_void_p, c_double, c_double, c_double,
     ctypes_numpy_array(int32,1), ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,1),
     c_int, c_bool, c_bool, c_char_p])

# mugtok_step(mugtok_ptr,curr,ncoils,time,dt,nl_its,lin_its,nretry,error_str)
mugtok_step = ctypes_subroutine(oftpy_lib.mugtok_step,
    [c_void_p, ctypes_numpy_array(float64,1), c_int,
     c_double_ptr, c_double_ptr, c_int_ptr, c_int_ptr, c_int_ptr, c_char_p])

# mugtok_get_field(mugtok_ptr,field_id,vals,error_str)
mugtok_get_field = ctypes_subroutine(oftpy_lib.mugtok_get_field,
    [c_void_p, c_int, ctypes_numpy_array(float64,1), c_char_p])

# mugtok_get_pmesh(mugtok_ptr,np,r_loc,nc,lc_loc,reg_loc,error_str)
mugtok_get_pmesh = ctypes_subroutine(oftpy_lib.mugtok_get_pmesh,
    [c_void_p, c_int_ptr, c_double_ptr_ptr, c_int_ptr, c_int_ptr_ptr, c_int_ptr_ptr, c_char_p])

# mugtok_get_gradf(mugtok_ptr,gr,gz,error_str)
mugtok_get_gradf = ctypes_subroutine(oftpy_lib.mugtok_get_gradf,
    [c_void_p, ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,1), c_char_p])

# mugtok_destroy(mugtok_ptr,error_str)
mugtok_destroy = ctypes_subroutine(oftpy_lib.mugtok_destroy,
    [c_void_p, c_char_p])

## @endcond


class MUGToksim():
    '''! Coupled MUG/TokaMaker (mugtok_td) time-dependent simulation.

    Built on top of an existing @ref OpenFUSIONToolkit.TokaMaker._core.TokaMaker "TokaMaker"
    instance, whose mesh, finite-element representation, and active equilibrium are reused.
    The poloidal flux (`psi`) and coil currents live in the TokaMaker
    equilibrium; the MHD fields (velocity, pressure, F) live in this object.
    '''

    def __init__(self, tokamaker):
        '''! Create a combined MUG/TokaMaker simulation from an existing TokaMaker instance

        @param tokamaker A TokaMaker object with an active equilibrium
        '''
        if tokamaker._tMaker_equil is None:
            raise ValueError('TokaMaker object has no active equilibrium')
        ## Parent TokaMaker object (mesh/equilibrium/coils borrowed from it)
        self._tokamaker = tokamaker
        ## OpenFUSIONToolkit environment (See @ref OpenFUSIONToolkit._core.OFT_env "OFT_env")
        self._oft_env = tokamaker._oft_env
        ## C pointer to the Fortran-side combined simulation object
        self._ptr = ctypes.c_void_p()
        ## Cached order-1 (pressure) plotting mesh (loaded on first access)
        self._pr = None
        self._plc = None
        self._preg = None
        ## Per-region MHD flag [nreg] (1 = MHD region); set by `setup()`
        self._mhd_flag = None
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_alloc(ctypes.byref(self._ptr), tokamaker._tMaker_equil.c_ptr, error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def __del__(self):
        '''! Free the Fortran-side object (the borrowed equilibrium/mesh are left intact)'''
        if getattr(self, '_ptr', None) is not None and self._ptr:
            error_string = self._oft_env.get_c_errorbuff()
            mugtok_destroy(self._ptr, error_string)
            self._ptr = ctypes.c_void_p()

    @property
    def r(self):
        '''! Node coordinates [np,3] for the main order fields.

        These are the same tessellated node points as the parent `TokaMaker.r`, so the
        order-2 fields (and `psi`) can be plotted with `self.r`/`self.lc`.
        '''
        return self._tokamaker.r

    @property
    def lc(self):
        '''! Triangle list [nc,3] for the main order fields (same as `TokaMaker.lc`)'''
        return self._tokamaker.lc

    @property
    def reg(self):
        '''! Per-cell region IDs [nc] for the main order fields (same as `TokaMaker.reg`).
        '''
        return self._tokamaker.reg

    def _load_pmesh(self):
        '''! Load and cache the order-1 (pressure) plotting mesh'''
        if self._pr is not None:
            return
        np_loc = ctypes.c_int()
        nc_loc = ctypes.c_int()
        r_loc = c_double_ptr()
        lc_loc = c_int_ptr()
        reg_loc = c_int_ptr()
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_get_pmesh(self._ptr, ctypes.byref(np_loc), ctypes.byref(r_loc),
                         ctypes.byref(nc_loc), ctypes.byref(lc_loc),
                         ctypes.byref(reg_loc), error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        self._pr = numpy.ctypeslib.as_array(r_loc, shape=(np_loc.value, 3))
        self._plc = numpy.ctypeslib.as_array(lc_loc, shape=(nc_loc.value, 3))
        self._preg = numpy.ctypeslib.as_array(reg_loc, shape=(nc_loc.value,))

    @property
    def pressure_r(self):
        '''! Node coordinates [np_p,3] for the (order-1) pressure field'''
        self._load_pmesh()
        return self._pr

    @property
    def pressure_lc(self):
        '''! Triangle list [nc_p,3] for the (order-1) pressure field'''
        self._load_pmesh()
        return self._plc

    @property
    def pressure_reg(self):
        '''! Per-cell region IDs [nc_p] for the (order-1) pressure mesh.'''
        self._load_pmesh()
        return self._preg

    def setup(self, dt, lin_tol, nl_tol, mhd_dict,
              incomp=True, allow_toroidal_flow=False):
        '''! Set up the coupled solver

        @param dt Timestep [s]
        @param lin_tol Linear solver tolerance
        @param nl_tol Non-linear solver tolerance
        @param mhd_dict Dictionary of regions where the MHD (MUG) solve is active,
                        keyed by region name. Each entry must provide:
                          'density'   mass density [kg/m^3]
                          'viscosity' dynamic viscosity [Pa-s]
                        e.g. {'CHANNEL': {'density': 9.4E3, 'viscosity': 1.5E-3}}
                        Every named region must also be a TokaMaker conducting
                        region, since that is where its resistivity comes from.
                        May be empty, in which case no region is evolved by MUG.
        @param incomp Use incompressible flow (default: True). Compressible flow
                      (`incomp=False`) is not currently supported and will throw an error.
        @param allow_toroidal_flow Allow toroidal (phi) flow in the MHD regions (default:False)
            
        '''
        nreg = self._tokamaker.nregs
        if nreg <= 0:
            raise ValueError('TokaMaker regions are not set up (nregs <= 0)')
        mhd_flag = numpy.zeros((nreg,), dtype=numpy.int32)
        dens_reg = -numpy.ones((nreg,), dtype=numpy.float64)
        visc_reg = -numpy.ones((nreg,), dtype=numpy.float64)
        for name, props in mhd_dict.items():
            # MHD regions are conducting regions, so they live in the TokaMaker cond_dict
            if name not in self._tokamaker._cond_dict:
                raise KeyError('MHD region "{0}" is not a TokaMaker conducting region'.format(name))
            for key in ('density', 'viscosity'):
                if key not in props:
                    raise KeyError('MHD region "{0}" is missing "{1}"'.format(name, key))
            dens = props['density']
            visc = props['viscosity']
            if dens <= 0.0:
                raise ValueError('MHD region "{0}" has non-positive density '
                                 '({1})'.format(name, dens))
            i = self._tokamaker._cond_dict[name]['reg_id'] - 1
            mhd_flag[i] = 1
            dens_reg[i] = dens
            visc_reg[i] = visc/dens  # dynamic [Pa-s] -> kinematic [m^2/s]
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_setup(self._ptr, ctypes.c_double(dt), ctypes.c_double(lin_tol), ctypes.c_double(nl_tol),
                     numpy.ascontiguousarray(mhd_flag, dtype=numpy.int32),
                     numpy.ascontiguousarray(dens_reg, dtype=numpy.float64),
                     numpy.ascontiguousarray(visc_reg, dtype=numpy.float64),
                     ctypes.c_int(nreg), ctypes.c_bool(incomp),
                     ctypes.c_bool(allow_toroidal_flow), error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        # Remember which regions are MHD so `get_field` can build per-cell field masks
        self._mhd_flag = mhd_flag
        ## Per-region MHD properties as supplied to `setup()`
        self._mhd_dict = mhd_dict

    def step(self, time, dt, coil_currents=None):
        '''! Advance the coupled solution by one timestep

        @param time Time at the start of the step [s]
        @param dt Timestep size [s]
        @param coil_currents Coil currents as a dict of name->[A]
                             (if `None`, the equilibrium's present currents are used)
        @result new time, new dt, # nonlinear iterations, # linear iterations, # retries
        '''
        c_time = ctypes.c_double(time)
        c_dt = ctypes.c_double(dt)
        nl_its = ctypes.c_int()
        lin_its = ctypes.c_int()
        nretry = ctypes.c_int()
        if coil_currents is None:
            coil_currents, _ = self._tokamaker.get_coil_currents()
        currents = numpy.ascontiguousarray(self._tokamaker.coil_dict2vec(coil_currents),
                                           dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_step(self._ptr, currents, ctypes.c_int(self._tokamaker.ncoils),
                    ctypes.byref(c_time), ctypes.byref(c_dt),
                    ctypes.byref(nl_its), ctypes.byref(lin_its), ctypes.byref(nretry), error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return c_time.value, c_dt.value, nl_its.value, lin_its.value, nretry.value

    def _get_mug_field(self, field_id):
        '''! Fetch a MUG-owned field (by augmented-vector block index) at node points'''
        # Size the receiving array from the mesh the field lives on: pressure is on the
        # (order-1) pressure mesh, every other MUG field is on the main-order mesh. 
        if field_id == 5:  # pressure (order-1 pressure mesh)
            n = self.pressure_r.shape[0]
        else:
            n = self.r.shape[0]
        vals = numpy.zeros((n,), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_get_field(self._ptr, ctypes.c_int(field_id), vals, error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return vals

    def _get_gradf(self):
        '''! Fetch the projected poloidal gradient of F, (dF/dR, dF/dZ), at node points'''
        # gradF lives at the main-order mesh nodes (like F itself)
        n = self.r.shape[0]
        gr = numpy.zeros((n,), dtype=numpy.float64)
        gz = numpy.zeros((n,), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_get_gradf(self._ptr, gr, gz, error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return gr, gz

    def get_field(self, name, cell_centered=False):
        '''! Get the present values of a field at node points

        @param name One of 'velocity', 'pressure', 'F' (MUG-owned),
                    'psi', or 'current density'
        @param cell_centered If True, return one value per mesh cell (averaged over the cell's
                    vertices) instead of one per node
        @result A tuple `(mask, field)` where `mask` is a per-cell boolean mask to apply at
                plot time

        Returned units by field:
          - 'velocity':        [m/s]
          - 'pressure':        [Pa]   (defined only up to a constant per MHD region;
                                       each region shifted so its own min = 0)
          - 'F':               [T*m]  
          - 'psi':             [Wb]   
          - 'current density': [A/m^2]

        The per-cell masks let each field be plotted only where it is physically defined:
          - 'velocity', 'pressure': `True` only in the MHD (MUG) regions.
          - 'F', 'psi': `True` everywhere (defined over the whole mesh).
          - 'current': the conductor mask from the TokaMaker equilibrium.

        '''
        if name in ('velocity', 'pressure'):
            # These are defined only where the MUG solve runs
            if self._mhd_flag is None:
                raise RuntimeError('MHD regions are not defined; call `setup()` first')
            mhd_ids = numpy.flatnonzero(self._mhd_flag) + 1  # region IDs are 1-based

        if name == 'psi':
            field = self._tokamaker.get_psi(normalized=False)
            mask = numpy.ones((self.lc.shape[0],), dtype=bool)
            lc = self.lc
        elif name == 'velocity':
            field = numpy.stack([self._get_mug_field(2),   # R
                                 self._get_mug_field(3),   # phi
                                 self._get_mug_field(4)],  # Z
                                axis=-1)   # (ndof, 3)
            mask = numpy.isin(self.reg, mhd_ids)
            lc = self.lc
        elif name == 'pressure':
            # Pressure is block 5. Shift each region by its own minimum so every
            # region starts at 0, because pressure is defined only up to a constant.
            p = self._get_mug_field(5)
            mask = numpy.isin(self.pressure_reg, mhd_ids)
            field = p.copy()
            for reg_id in mhd_ids:
                cells = (self.pressure_reg == reg_id)
                if not cells.any():
                    continue
                nodes = numpy.unique(self.pressure_lc[cells])
                field[nodes] = p[nodes] - p[nodes].min()
            lc = self.pressure_lc 
        elif name == 'F':
            field = self._get_mug_field(7)
            mask = numpy.ones((self.lc.shape[0],), dtype=bool)
            lc = self.lc
        elif name == 'current density':
            mask, field = self._get_current()
            lc = self.lc
        else:
            raise ValueError('Unknown field "{0}"; expected one of '
                             'velocity, pressure, F, psi, current density'.format(name))

        if cell_centered:
            field = numpy.mean(field[lc], axis=1)
        return mask, field

    def plot_field(self, fig, ax, name, colormap=None, scale=1.0, clabel=None,
                   colorbar=True, shared_scale=False, **kwargs):
        '''! Plot a field on the mesh

        Wraps @ref get_field and draws it with `tripcolor`, applying the field's mask so
        only the region where it is defined is shaded.

        Vector fields are drawn as three separate colour maps, one per component, 
        so `ax` must hold three axes for those.

        @param fig Figure the axes belong to (needed for the colorbars)
        @param ax Axes to draw on. Scalar fields ('pressure', 'F', 'psi') take a single
                  Axes; vector fields ('velocity', 'current density') take a sequence of three
        @param name Field name, as accepted by `get_field`
        @param colormap Colormap to use. By default 'velocity' and 'current density' use
                  'seismic' with limits symmetric about zero, since their sign carries
                  direction; every other field uses 'viridis' over its own range.
        @param scale Multiplier applied to the values before plotting.
                     Set `clabel` to match: e.g. `scale=1.0E-6, clabel=r'$J$ [MA/m$^2$]'`.
        @param clabel Label for the colorbars. `None` builds one per component from the
                  field name and its base units. 
        @param colorbar Draw colorbars (default True)
        @param shared_scale For vector fields, put all three components on a common
                  colour scale so the panels are directly comparable, and draw a single
                  colorbar. (default False)
        '''
        labels = {'velocity': (r'$v_R$', r'$v_\phi$', r'$v_Z$'),
                  'current density': (r'$J_R$', r'$J_\phi$', r'$J_Z$'),
                  'pressure': ('p',), 'psi': (r'$\psi$',), 'F': ('F',)}
        units = {'velocity': 'm/s', 'pressure': 'Pa', 'F': r'T$\cdot$m',
                 'psi': 'Wb', 'current density': r'A/m$^2$'}
        # Signed fields, where direction matters and zero is meaningful
        signed = ('velocity', 'current density')
        if name not in labels:
            raise ValueError('Unknown field "{0}"; expected one of {1}'.format(
                name, ', '.join(sorted(labels))))

        mask, field = self.get_field(name)
        # 'pressure' is a lower-order field on its own node/cell lists
        if name == 'pressure':
            r, lc = self.pressure_r, self.pressure_lc
        else:
            r, lc = self.r, self.lc

        if mask.sum() == 0:
            print('Warning: field "{0}" is not defined in any region'.format(name))
            return None

        if field.ndim == 2:                      # vector: one panel per component
            axes = list(numpy.atleast_1d(numpy.asarray(ax, dtype=object)).ravel())
            if len(axes) < 3:
                raise ValueError('"{0}" is a vector field; `ax` must hold three axes '
                                 '(got {1}) so the R, phi and Z components can be drawn '
                                 'separately'.format(name, len(axes)))
            columns = [field[:, k]*scale for k in range(3)]
        else:                                    # scalar: a single panel
            axes = [ax]
            columns = [field*scale]

        # Colorbar labels: auto from the field name and its base units, or as given
        if clabel is None:
            panel_labels = ['{0} [{1}]'.format(l, units[name]) for l in labels[name]]
        elif isinstance(clabel, str):
            panel_labels = [clabel]*len(columns)
        else:
            panel_labels = list(clabel)

        user_vmin = kwargs.pop('vmin', None)
        user_vmax = kwargs.pop('vmax', None)
        cmap = colormap if colormap is not None else ('seismic' if name in signed else 'viridis')
        drawn_nodes = numpy.unique(lc[mask])

        def _limits(vals):
            d = vals[drawn_nodes]
            if name in signed:
                span = float(numpy.max(numpy.abs(d)))
                return -span, span
            return float(numpy.min(d)), float(numpy.max(d))

        if shared_scale and len(columns) > 1:
            lims = [_limits(c) for c in columns]
            limits = [(min(l[0] for l in lims), max(l[1] for l in lims))]*len(columns)
        else:
            limits = [_limits(c) for c in columns]

        one_bar = shared_scale and len(columns) > 1
        handles = []
        for i, (a, values, label) in enumerate(zip(axes, columns, labels[name])):
            lo, hi = limits[i]
            clf = a.tripcolor(r[:, 0], r[:, 1], lc[mask, :], values, cmap=cmap,
                              vmin=lo if user_vmin is None else user_vmin,
                              vmax=hi if user_vmax is None else user_vmax,
                              **kwargs)
            if colorbar and not one_bar:
                fig.colorbar(clf, ax=a, label=panel_labels[i])
            if one_bar:
                # A single shared colorbar cannot name the components, so title them
                a.set_title(label, fontsize=11)
            a.set_aspect('equal', 'box')
            handles.append(clf)
        if colorbar and one_bar:
            # One scale for the row, so one colorbar
            fig.colorbar(handles[0], ax=axes, label=panel_labels[0])
        return handles

    def save_to_dict(self, cell_centered=False):
        '''! Collect every field at the current time into a single dictionary

        @param cell_centered If True, average each field over its cell vertices to give
                     one value per cell.
        @result Dictionary keyed by field name ('velocity', 'pressure', 'F', 'psi',
                'current density'). Each entry is a sub-dictionary:
                  'field' values, 'mask' per-cell boolean for where the field is defined,
                  'units' string for the field's base units.

        '''
        units = {'velocity': 'm/s', 'pressure': 'Pa', 'F': 'T*m',
                 'psi': 'Wb', 'current density': 'A/m^2'}
        out = {}
        for name in units:
            mask, field = self.get_field(name, cell_centered=cell_centered)
            out[name] = {'field': field, 'mask': mask, 'units': units[name]}
        return out

    def _get_current(self):
        '''! Current density in the conducting structures, J [A/m^2].

        Poloidal components come from the (projected) grad(F): J_pol = (grad F x phi_hat)/(mu0 R);
        the toroidal component and the conductor mask come from the TokaMaker equilibrium via
        `calc_conductor_currents(psi)`.

        @result Tuple `(mask, J)` where `mask` is the per-cell conductor mask and `J` is a
                (ndof, 3) array with columns ordered (R, phi, Z).
        '''
        # Poloidal current from projected grad(F), with 1/R guarded at the axis
        gr, gz = self._get_gradf()          # dF/dR, dF/dZ
        R = self.r[:, 0]
        JR = numpy.zeros_like(gr)
        JZ = numpy.zeros_like(gr)
        m = R > 0.0
        JR[m] = -gz[m] / (mu0 * R[m])       # -(1/(mu0 R)) dF/dZ
        JZ[m] =  gr[m] / (mu0 * R[m])       #  (1/(mu0 R)) dF/dR
        # Toroidal current + conductor mask from the TokaMaker equilibrium
        psi = self._tokamaker.get_psi(normalized=False)
        mask, Jphi = self._tokamaker._tMaker_equil.calc_conductor_currents(psi)
        J = numpy.stack([JR, Jphi, JZ], axis=-1)   # (ndof, 3): columns (R, phi, Z)
        return mask, J
