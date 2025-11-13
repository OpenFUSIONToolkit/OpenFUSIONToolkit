#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 15:00:23 2025

@author: mparsons
"""

import numpy as np
    
class PCS:
    
    # Initialize all parameters from input dictionary
    def __init__(self, inputs):
        
        
        # Time Parameters
        self.dt = inputs['dt_PCS']
        self.t_max = inputs['t_max']
        self.N_its = int(np.ceil(inputs['t_max'] / inputs['dt_PCS']))
        self.N_iso = int(inputs['dt_ISO']/inputs['dt_PCS'])
        self.t = np.linspace(0, self.t_max, self.N_its+1)
        
        
        # OH Ramp
        self.IP = np.interp(self.t, inputs['IPOH_IP'][:,0], inputs['IPOH_IP'][:,1])
        self.IOH = np.interp(self.t, inputs['IPOH_IOH'][:,0], inputs['IPOH_IOH'][:,1])
        
        
        # Derivative Voltage Control Algorithm
        self.DVC_Gain_P = np.interp(self.t, inputs['DVC_Gain_P'][:,0], inputs['DVC_Gain_P'][:,1])
        self.DVC_Gain_I = np.interp(self.t, inputs['DVC_Gain_I'][:,0], inputs['DVC_Gain_I'][:,1])
        self.DVC_Gain_D = np.interp(self.t, inputs['DVC_Gain_D'][:,0], inputs['DVC_Gain_D'][:,1])
        self.DVC_err_P = np.zeros((self.N_its+1, ))
        self.DVC_err_I = np.zeros((self.N_its+1, ))
        self.DVC_err_D = np.zeros((self.N_its+1, ))
        
        
        # ISOFLUX
        self.ISO_MATRIX_M = np.zeros((12,8))
        self.ISO_MATRIX_M[:,0] = inputs['ISO_SEG01_MM']
        self.ISO_MATRIX_M[:,1] = inputs['ISO_SEG02_MM']
        self.ISO_MATRIX_M[:,2] = inputs['ISO_SEG03_MM']
        self.ISO_MATRIX_M[:,3] = inputs['ISO_SEG04_MM']
        self.ISO_MATRIX_M[:,4] = inputs['ISO_GRID1R_MM']
        self.ISO_MATRIX_M[:,5] = inputs['ISO_GRID1Z_MM']
        self.ISO_MATRIX_M[:,6] = inputs['ISO_GRID2R_MM']
        self.ISO_MATRIX_M[:,7] = inputs['ISO_GRID2Z_MM']
        
        
        self.ISO_CP_R = np.zeros((self.N_its+1,8))
        self.ISO_CP_R[:,0] = np.interp(self.t, inputs['ISO_SEG01_CP'][:,0], inputs['ISO_SEG01_CP'][:,1])
        self.ISO_CP_R[:,1] = np.interp(self.t, inputs['ISO_SEG02_CP'][:,0], inputs['ISO_SEG02_CP'][:,1])
        self.ISO_CP_R[:,2] = np.interp(self.t, inputs['ISO_SEG03_CP'][:,0], inputs['ISO_SEG03_CP'][:,1])
        self.ISO_CP_R[:,3] = np.interp(self.t, inputs['ISO_SEG04_CP'][:,0], inputs['ISO_SEG04_CP'][:,1])
        self.ISO_CP_R[:,4] = np.interp(self.t, inputs['ISO_GRID1_CP'][:,0], inputs['ISO_GRID1_CP'][:,1])
        self.ISO_CP_R[:,5] = np.interp(self.t, inputs['ISO_GRID1_CP'][:,0], inputs['ISO_GRID1_CP'][:,1])
        self.ISO_CP_R[:,6] = np.interp(self.t, inputs['ISO_GRID2_CP'][:,0], inputs['ISO_GRID2_CP'][:,1])
        self.ISO_CP_R[:,7] = np.interp(self.t, inputs['ISO_GRID2_CP'][:,0], inputs['ISO_GRID2_CP'][:,1])
        
        
        self.ISO_CP_Z = np.zeros((self.N_its+1,8))
        self.ISO_CP_Z[:,0] = np.interp(self.t, inputs['ISO_SEG01_CP'][:,0], inputs['ISO_SEG01_CP'][:,2])
        self.ISO_CP_Z[:,1] = np.interp(self.t, inputs['ISO_SEG02_CP'][:,0], inputs['ISO_SEG02_CP'][:,2])
        self.ISO_CP_Z[:,2] = np.interp(self.t, inputs['ISO_SEG03_CP'][:,0], inputs['ISO_SEG03_CP'][:,2])
        self.ISO_CP_Z[:,3] = np.interp(self.t, inputs['ISO_SEG04_CP'][:,0], inputs['ISO_SEG04_CP'][:,2])
        self.ISO_CP_Z[:,4] = np.interp(self.t, inputs['ISO_GRID1_CP'][:,0], inputs['ISO_GRID1_CP'][:,2])
        self.ISO_CP_Z[:,5] = np.interp(self.t, inputs['ISO_GRID1_CP'][:,0], inputs['ISO_GRID1_CP'][:,2])
        self.ISO_CP_Z[:,6] = np.interp(self.t, inputs['ISO_GRID2_CP'][:,0], inputs['ISO_GRID2_CP'][:,2])
        self.ISO_CP_Z[:,7] = np.interp(self.t, inputs['ISO_GRID2_CP'][:,0], inputs['ISO_GRID2_CP'][:,2])
        
        
        self.ISO_Gain_P = np.zeros((self.N_its+1,8))
        self.ISO_Gain_P[:,0] = np.interp(self.t, inputs['ISO_SEG01_GP'][:,0], inputs['ISO_SEG01_GP'][:,1])
        self.ISO_Gain_P[:,1] = np.interp(self.t, inputs['ISO_SEG02_GP'][:,0], inputs['ISO_SEG02_GP'][:,1])
        self.ISO_Gain_P[:,2] = np.interp(self.t, inputs['ISO_SEG03_GP'][:,0], inputs['ISO_SEG03_GP'][:,1])
        self.ISO_Gain_P[:,3] = np.interp(self.t, inputs['ISO_SEG04_GP'][:,0], inputs['ISO_SEG04_GP'][:,1])
        self.ISO_Gain_P[:,4] = np.interp(self.t, inputs['ISO_GRID1R_GP'][:,0], inputs['ISO_GRID1R_GP'][:,1])
        self.ISO_Gain_P[:,5] = np.interp(self.t, inputs['ISO_GRID1Z_GP'][:,0], inputs['ISO_GRID1Z_GP'][:,1])
        self.ISO_Gain_P[:,6] = np.interp(self.t, inputs['ISO_GRID2R_GP'][:,0], inputs['ISO_GRID2R_GP'][:,1])
        self.ISO_Gain_P[:,7] = np.interp(self.t, inputs['ISO_GRID2Z_GP'][:,0], inputs['ISO_GRID2Z_GP'][:,1])
        
        
        self.ISO_Gain_I = np.zeros((self.N_its+1,8))
        self.ISO_Gain_I[:,0] = np.interp(self.t, inputs['ISO_SEG01_GI'][:,0], inputs['ISO_SEG01_GI'][:,1])
        self.ISO_Gain_I[:,1] = np.interp(self.t, inputs['ISO_SEG02_GI'][:,0], inputs['ISO_SEG02_GI'][:,1])
        self.ISO_Gain_I[:,2] = np.interp(self.t, inputs['ISO_SEG03_GI'][:,0], inputs['ISO_SEG03_GI'][:,1])
        self.ISO_Gain_I[:,3] = np.interp(self.t, inputs['ISO_SEG04_GI'][:,0], inputs['ISO_SEG04_GI'][:,1])
        self.ISO_Gain_I[:,4] = np.interp(self.t, inputs['ISO_GRID1R_GI'][:,0], inputs['ISO_GRID1R_GI'][:,1])
        self.ISO_Gain_I[:,5] = np.interp(self.t, inputs['ISO_GRID1Z_GI'][:,0], inputs['ISO_GRID1Z_GI'][:,1])
        self.ISO_Gain_I[:,6] = np.interp(self.t, inputs['ISO_GRID2R_GI'][:,0], inputs['ISO_GRID2R_GI'][:,1])
        self.ISO_Gain_I[:,7] = np.interp(self.t, inputs['ISO_GRID2Z_GI'][:,0], inputs['ISO_GRID2Z_GI'][:,1])
        
        
        self.ISO_Gain_D = np.zeros((self.N_its+1,8))
        self.ISO_Gain_D[:,0] = np.interp(self.t, inputs['ISO_SEG01_GD'][:,0], inputs['ISO_SEG01_GD'][:,1])
        self.ISO_Gain_D[:,1] = np.interp(self.t, inputs['ISO_SEG02_GD'][:,0], inputs['ISO_SEG02_GD'][:,1])
        self.ISO_Gain_D[:,2] = np.interp(self.t, inputs['ISO_SEG03_GD'][:,0], inputs['ISO_SEG03_GD'][:,1])
        self.ISO_Gain_D[:,3] = np.interp(self.t, inputs['ISO_SEG04_GD'][:,0], inputs['ISO_SEG04_GD'][:,1])
        self.ISO_Gain_D[:,4] = np.interp(self.t, inputs['ISO_GRID1R_GD'][:,0], inputs['ISO_GRID1R_GD'][:,1])
        self.ISO_Gain_D[:,5] = np.interp(self.t, inputs['ISO_GRID1Z_GD'][:,0], inputs['ISO_GRID1Z_GD'][:,1])
        self.ISO_Gain_D[:,6] = np.interp(self.t, inputs['ISO_GRID2R_GD'][:,0], inputs['ISO_GRID2R_GD'][:,1])
        self.ISO_Gain_D[:,7] = np.interp(self.t, inputs['ISO_GRID2Z_GD'][:,0], inputs['ISO_GRID2Z_GD'][:,1])
        
        
        self.ISO_err_P = np.zeros((self.N_its+1,8))
        self.ISO_err_I = np.zeros((self.N_its+1,8))
        self.ISO_err_D = np.zeros((self.N_its+1,8))
        
        self.ISO_GRID1_BOUNDS = inputs['ISO_GRID1_BOUNDS']
        self.ISO_GRID2_BOUNDS = inputs['ISO_GRID2_BOUNDS']
        
    
    # Return isoflux control points for time step i
    def get_control_points_i(self, i):
        
        RZ_cp = np.transpose(np.array([self.ISO_CP_R[i,[0,1,2,3]], self.ISO_CP_Z[i,[0,1,2,3]]]))
        
        return RZ_cp
    
    
    # Return isoflux x-point targets for time step i
    def get_xpt_targets_i(self, i):
        
        RZ_xp = np.transpose(np.array([self.ISO_CP_R[i,[4,6]], self.ISO_CP_Z[i,[4,6]]]))
        
        return RZ_xp
    
    
    # Return isoflux gains for time step i
    def get_iso_gains_i(self, i):
        
        P_gains = self.ISO_Gain_P[i,:]
        
        I_gains = self.ISO_Gain_I[i,:]
        
        D_gains = self.ISO_Gain_D[i,:]
        
        return P_gains, I_gains, D_gains
    
    
    # Return isoflux errors for time step i
    def get_iso_errors_i(self, i):
        
        P_errors = self.ISO_err_P[i,:]
        
        I_errors = self.ISO_err_I[i,:]
        
        D_errors = self.ISO_err_D[i,:]
        
        return P_errors, I_errors, D_errors
    
    
    # Extrapolate isoflux errors for time step i
    def extrap_iso_errors_i(self, i):
        
        P_errors_prev, I_errors_prev, D_errors_prev = self.get_iso_errors_i(i-1)
        
        errors_P = P_errors_prev
        errors_I = I_errors_prev
        errors_D = D_errors_prev
        
        return errors_P, errors_I, errors_D
        
    
    # Calculate isoflux errors for time step i
    def calc_iso_errors_i(self, P_errors_i, i):
        
        P_errors_prev, I_errors_prev, _ = self.get_iso_errors_i(i-1)
        
        errors_P = P_errors_i
        errors_I = I_errors_prev + P_errors_i * self.dt
        errors_D = (P_errors_i - P_errors_prev) / self.dt
        
        return errors_P, errors_I, errors_D
        
        

    # Update PID errors for time step i
    def update_iso_errors_i(self, errors_P, errors_I, errors_D, i):
        
        self.ISO_err_P[i,:] = errors_P
        self.ISO_err_I[i,:] = errors_I
        self.ISO_err_D[i,:] = errors_D
        
    
    # Calculate the isolux P vector for time step i
    def calc_pvector(self, i):
        
        # Get PID errors
        errors_P, errors_I, errors_D = self.get_iso_errors_i(i)
        
        # Get PID gains
        P_gains, I_gains, D_gains = self.get_iso_gains_i(i)
        
        # Calculate P vector
        pvector = errors_P * P_gains + errors_I * I_gains + errors_D * D_gains
        
        return pvector
    
    
    # Calculate the isoflux commands for time step i
    def calc_iso_commands(self, i):
        
        # Calculate P vector
        pvector = self.calc_pvector(i)
        
        # Map P vector to coil commands
        MP = np.matmul(self.ISO_MATRIX_M, pvector)
        
        return MP
        