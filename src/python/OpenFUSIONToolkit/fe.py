#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Access to finite element functionality for Open FUSION Toolkit (OFT) Python interfaces

@authors Chris Hansen
@date September 2026
@ingroup doxy_oft_python
'''
import ctypes
import numpy
from ._interface import *


class Lagrange_surface_field_interpolator():
    '''! Interpolation class for surface Lagrange FE fields'''
    def __init__(self,oft_env,fe_obj,values,int_type,fbary_tol=1.E-8):
        '''! Initialize interpolation object

        @param fe_obj Address of Fortran FE object
        @param values FE weight vector
        @param int_type Interpolation type (should be 1 or 2)
        @param fbary_tol Tolerance for physical to logical mapping
        '''
        if int_type == 1:
            self.dim = 1
            self._int_type = 211
        elif int_type == 2:
            self.dim = 3
            self._int_type = 212
        else:
            raise ValueError('Invalid interpolation type, must be 1 or 2')
        self._dim_eval = self.dim
        self._oft_env = oft_env
        self._fe_obj = fe_obj
        self.fbary_tol = c_double(fbary_tol)
        self.cell = c_int(-1)
        #
        values = numpy.ascontiguousarray(values, dtype=numpy.float64)
        self._int_obj = c_void_p()
        error_string = self._oft_env.get_c_errorbuff()
        oft_get_field_eval(fe_obj,values,self._int_type,ctypes.byref(self._int_obj),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def __del__(self):
        '''Destroy underlying interpolation object'''
        pts_eval = numpy.zeros((1,1), dtype=numpy.float64)
        vals_out = numpy.zeros((1,self._dim_eval), dtype=numpy.float64)
        oft_apply_field_eval(self._fe_obj,self._int_obj,-self._int_type,pts_eval,1,self.fbary_tol,self._dim_eval,vals_out)

    def eval(self,pts):
        '''! Evaluate field at a given location

        @param pts Location for evaluation [3] or [n,3]
        @result Field at evaluation point [self.dim]
        '''
        if pts.ndim == 1:
            in_1d = True
            npts = 1
            pts = pts.reshape((1,pts.shape[0]))
        else:
            in_1d = False
            npts = pts.shape[0]
        if pts.shape[1] != 3:
            raise ValueError('Input points must be 3D coordinates')
        pts = numpy.ascontiguousarray(pts, dtype=numpy.float64)
        vals_out = numpy.zeros((npts,self._dim_eval), dtype=numpy.float64)
        oft_apply_field_eval(self._fe_obj,self._int_obj,self._int_type,pts,npts,self.fbary_tol,self._dim_eval,vals_out)
        if in_1d:
            return vals_out[0,:self.dim].copy()
        else:
            return vals_out[:,:self.dim].copy()


class Lagrange_2D_field_interpolator(Lagrange_surface_field_interpolator):
    '''! Interpolation class for 2D Lagrange FE fields'''
    def __init__(self,oft_env,fe_obj,values,int_type,fbary_tol=1.E-8):
        '''! Initialize interpolation object

        @param fe_obj Address of Fortran FE object
        @param values FE weight vector
        @param int_type Interpolation type (should be 21 or 22)
        @param fbary_tol Tolerance for physical to logical mapping
        '''
        super().__init__(oft_env,fe_obj,values,int_type,fbary_tol)
        if self.dim == 3: # Drop z-component for 2D interpolation
            self._dim_eval = 3
            self.dim = 2

    def eval(self,pts):
        '''! Evaluate field at a given location

        @param pts Location for evaluation [2] or [n,2]
        @result Field at evaluation point [self.dim_return]
        '''
        pts = numpy.ascontiguousarray(pts, dtype=numpy.float64)
        if pts.ndim == 1:
            in_1d = True
            npts = 1
            pts = pts.reshape((1,pts.shape[0]))
        else:
            in_1d = False
            npts = pts.shape[0]
        if pts.shape[1] != 2:
            raise ValueError('Input points must be 2D coordinates')
        pts_eval = numpy.zeros((npts,3), dtype=numpy.float64)
        pts_eval[:,:2] = pts
        vals_out = numpy.zeros((npts,self._dim_eval), dtype=numpy.float64)
        oft_apply_field_eval(self._fe_obj,self._int_obj,self._int_type,pts_eval,npts,self.fbary_tol,self._dim_eval,vals_out)
        if in_1d:
            return vals_out[0,:self.dim].copy()
        else:
            return vals_out[:,:self.dim].copy()