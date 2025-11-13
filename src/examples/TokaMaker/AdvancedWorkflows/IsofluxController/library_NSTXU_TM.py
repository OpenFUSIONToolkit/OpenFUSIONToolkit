#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 14:34:26 2025

@author: mparsons
"""

import sys
OFT_ROOT = '/Applications/OpenFUSIONToolkit/python'
if(not OFT_ROOT in sys.path):
    sys.path.insert(0,OFT_ROOT)

from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk
from omfit_shape_generator import boundaryShape
import numpy as np
import matplotlib.tri as tri





    
##################################################################
def create_mygs_NSTXU(mygs):

    # Solver settings
    mygs.settings.pm=False
    
    # Set up NSTX-U mesh
    mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('NSTXU_mesh.h5')
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    mygs.setup(order=2,F0=1.0*0.85)
    
    # Create vertical stability coil
    mygs.set_coil_vsc({'PF3U': 1.0, 'PF3L': -1.0})

    # Regularization
    coil_reg_mat = np.eye(mygs.ncoils+1, dtype=np.float64)
    coil_reg_weights = np.ones((mygs.ncoils+1,), dtype=np.float64) * 1.0
    mygs.set_coil_reg(coil_reg_mat, reg_weights=coil_reg_weights)
    
    return mygs
    
##################################################################






    
##################################################################
def setup_plasma(mygs, inputs):
    
    
    # Coil current limits
    coil_bounds = {key: [-24.00E+03,24.00E+03] for key in mygs.coil_sets}
    coil_bounds['CS'][0] = inputs['IPOH_IOH'][0][1]
    coil_bounds['CS'][1] = inputs['IPOH_IOH'][0][1] + 1
    coil_bounds['PF1AU'][0] = 0 # Unipolar, positive
    coil_bounds['PF1AL'][0] = 0 # Unipolar, positive
    coil_bounds['PF2U'][0] = 0 # Unipolar, Positive
    coil_bounds['PF2L'][0] = 0 # Unipolar, Positive
    coil_bounds['PF4'][1] = 0 # Unipolar, negative
    coil_bounds['PF5'][1] = 0 # Unipolar, negative
    mygs.set_coil_bounds(coil_bounds)
    
    
    
    # Set up plasma shape
    if(inputs['flag_shape'] == 1):
        
        eqdsk_shape = read_eqdsk(inputs['gfile_shape'])
        isoflux_pts = eqdsk_shape['rzout'].copy()
        mygs.set_isoflux(isoflux=isoflux_pts)
        
    else:
        
        # Target plasma shape
        TPS = PSP(inputs['R0_target'], inputs['Z0_target'], inputs['a_target'], 
                  inputs['kappa_upper_target'], inputs['kappa_lower_target'],
                 inputs['delta_upper_target'], inputs['delta_lower_target'], 
                 inputs['zeta_upper_outer_target'], inputs['zeta_lower_outer_target'], 
                 inputs['zeta_upper_inner_target'], inputs['zeta_lower_inner_target'], 
                 inputs['null_type'])
        
        isoflux_pts = TPS.trace(inputs['N_points'])
        mygs.set_isoflux(isoflux=isoflux_pts)
        x_points = TPS.x_points()
        mygs.set_saddles(x_points)


        
        
    
    # Set up plasma profiles
    if(inputs['flag_profiles'] == 1):
        
        eqdsk_profiles = read_eqdsk(inputs['gfile_profiles'])
        
        ffprim = eqdsk_profiles['ffprim']
        pprime = eqdsk_profiles['pprime']
        psi_eqdsk = np.linspace(0.0,1.0,np.size(ffprim))
        psi_sample = np.linspace(0.0,1.0,50)
        
        psi_prof = np.copy(psi_sample)
        ffp_prof = np.transpose(np.vstack((psi_prof,np.interp(psi_sample,psi_eqdsk,ffprim)))).copy()
        pp_prof = np.transpose(np.vstack((psi_prof,np.interp(psi_sample,psi_eqdsk,pprime)))).copy()
        
        ffp_prof = {'type': 'linterp', 'y': ffp_prof[:,1], 'x': psi_sample}
        pp_prof={'type': 'linterp', 'y': pp_prof[:,1], 'x': psi_sample}
        
    else:
        
        # Set psi_n grid
        x = np.linspace(0,1,100)
        
        # Set FF' profile
        ffp_prof = {'type': 'linterp', 'x': x, 'y': inputs['ffp_2']*x*x + inputs['ffp_1']*x + inputs['ffp_0']}
        
        # Set p' profile
        pp_prof = {'type': 'linterp', 'x': x, 'y': inputs['pp_2']*x*x + inputs['pp_1']*x + inputs['pp_0']}
        
    ffp_prof['y'] /= max(ffp_prof['y'], key=abs) # Normalize profile (not required but convienient)
    pp_prof['y'] /= max(pp_prof['y'], key=abs) # Normalize profile (not required but convienient)
    mygs.set_profiles(ffp_prof=ffp_prof,pp_prof=pp_prof)
    
    # Set plasma parameters
    mygs.init_psi(0.90, 0, 0.55, 1.8, 0.6)
    mygs.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof)
    if(inputs['flag_pres'] == 1):
        mygs.set_targets(Ip = inputs['IPOH_IP'][0,1], pax=inputs['PRES'][0,1])
    else:
        mygs.set_targets(Ip = inputs['IPOH_IP'][0,1], Ip_ratio = 1/inputs['BETAP'][0,1] - 1)
    
    # Solve GS for initial equilibrium and set up time-dependent solver
    try:
        mygs.solve()
        mygs.setup_td(1.0e-9, inputs['lin_tol'], inputs['nl_tol'], pre_plasma=False)
        err_flag = 0
    except ValueError:
        err_flag = -1
    
    
    
    return mygs, err_flag
    
##################################################################










##################################################################

def interp_fluxes(psi, triObj, RZ):
    
    # Interpolate fluxes at each flux loop location
    fz = tri.LinearTriInterpolator(triObj, psi)
    fluxes = fz(RZ[:,0], RZ[:,1]).data

    return fluxes

##################################################################










##################################################################

def update_coilcurrents(mygs, dt, dI_ISO, dI_DVC, IOH_target):
    
    # Get currents
    coil_currents, _ = mygs.get_coil_currents()
    
    # OH Ramp
    coil_currents['CS'] = np.clip(IOH_target, -24.0e3, 24.0e3)
    
    # Limit rate of change of coil currents
    dI_PF = dI_ISO
    dI_PF[4] = dI_PF[4] + dI_DVC
    dI_PF[7] = dI_PF[7] - dI_DVC
    dI_PF = np.clip(dI_PF, -5.0e5*dt, 5.0e5*dt)
    
    # Limit coil currents to power supply range
    coil_currents['PF1AU'] = np.clip(coil_currents['PF1AU'] + dI_PF[0], 0, 24.0e3)
    coil_currents['PF1BU'] = np.clip(coil_currents['PF1BU'] + dI_PF[1], -24.0e3, 24.0e3)
    coil_currents['PF1CU'] = np.clip(coil_currents['PF1CU'] + dI_PF[2], -24.0e3, 24.0e3)
    coil_currents['PF2U'] = np.clip(coil_currents['PF2U'] + dI_PF[3], 0, 24.0e3)
    coil_currents['PF3U'] = np.clip(coil_currents['PF3U'] + dI_PF[4], -24.0e3, 24.0e3)
    coil_currents['PF4'] = np.clip(coil_currents['PF4'] + dI_PF[5], -24.0e3, 0)
    coil_currents['PF5'] = np.clip(coil_currents['PF5'] + dI_PF[6], -24.0e3, 0)
    coil_currents['PF3L'] = np.clip(coil_currents['PF3L'] + dI_PF[7], -24.0e3, 24.0e3)
    coil_currents['PF2L'] = np.clip(coil_currents['PF2L'] + dI_PF[8], 0, 24.0e3)
    coil_currents['PF1CL'] = np.clip(coil_currents['PF1CL'] + dI_PF[9], -24.0e3, 24.0e3)
    coil_currents['PF1BL'] = np.clip(coil_currents['PF1BL'] + dI_PF[10], -24.0e3, 24.0e3)
    coil_currents['PF1AL'] = np.clip(coil_currents['PF1AL'] + dI_PF[11], 0, 24.0e3)
    
    # Set new coil currents
    mygs.set_coil_currents(coil_currents)

    return mygs

##################################################################






    
##################################################################
'''
# Plasma Shape Parameters (Class)
'''
class PSP:
    
    # Initialize with all shape parameters and null type
    def __init__(self, R0, Z0, a, 
                 kappa_upper, kappa_lower,
                 delta_upper, delta_lower,
                 zeta_upper_outer, zeta_lower_outer, 
                 zeta_upper_inner, zeta_lower_inner, null_type):
    
        self.R0 = R0
        self.Z0 = Z0
        self.a = a
        self.kappa_upper = kappa_upper
        self.kappa_lower = kappa_lower
        self.delta_upper = delta_upper
        self.delta_lower = delta_lower
        self.zeta_upper_outer = zeta_upper_outer
        self.zeta_lower_outer = zeta_lower_outer
        self.zeta_upper_inner = zeta_upper_inner
        self.zeta_lower_inner = zeta_lower_inner
        self.null_type = null_type
        
    # Calculate trace of boundary with N_points per quadrant
    def trace(self, N_points):
        
        # Select nulls
        if(self.null_type == 'DN'):
            upnull_flag = True
            lonull_flag = True
        elif(self.null_type == 'LSN'):
            upnull_flag = False
            lonull_flag = True
        elif(self.null_type == 'USN'):
            upnull_flag = True
            lonull_flag = False
        else:
            upnull_flag = False
            lonull_flag = False
            

        # Shape
        R_ISO, Z_ISO, Z_ref = boundaryShape(
                a=self.a,
                eps=self.a/self.R0,
                kapu=self.kappa_upper,
                kapl=self.kappa_lower,
                delu=self.delta_upper,
                dell=self.delta_lower,
                zetaou=self.zeta_upper_outer,
                zetaiu=self.zeta_upper_inner,
                zetail=self.zeta_lower_inner,
                zetaol=self.zeta_lower_outer,
                zoffset=self.Z0,
                doPlot=False,
                npts=N_points,
                upnull=upnull_flag,
                lonull=lonull_flag
            )
        
        shape_trace = np.transpose(np.array([R_ISO, Z_ISO]))
        
        return shape_trace
    
    # Find x_points
    def x_points(self):
        
        # Trace boundary
        shape_trace = self.trace(5)
        
        # Highest point
        Z_top = shape_trace[:,1].max()
        R_top = shape_trace[np.argmax(shape_trace[:,1]), 0]
        
        # Lowest point
        Z_bot = shape_trace[:,1].min()
        R_bot = shape_trace[np.argmin(shape_trace[:,1]), 0]
        
        # Select nulls
        if(self.null_type == 'DN'):
            x_points = np.array([(R_bot, Z_bot), (R_top, Z_top)])
        elif(self.null_type == 'LSN'):
            x_points = np.array([(R_bot, Z_bot)])
        elif(self.null_type == 'USN'):
            x_points = np.array([(R_top, Z_top)])
        else:
            x_points = np.array()
            
        return x_points
    
    
    
##################################################################