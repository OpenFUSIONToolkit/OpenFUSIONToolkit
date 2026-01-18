import struct
import sys
import os
import numpy as np
import pyvista
from matplotlib.path import Path
import matplotlib.patches as patches
import json

    
def solve_filaments(self, timesteps, torus, PsiFull, total_current, error_matrix, num_coils, reg_factor_fil = 1.E-8, reg_factor_wall= 1.E-2):
            
    '''
    Function solves filament and wall currents over time 
    
    Parameters:
    ----------
    timesteps : int
        number of timesteps corresponding to the indexes!
    
    torus : reduced thincurr object
        should be a reduced thincurr object 
        
    PsiFull: float array
        magnetic sensor signals where each row is a different magnetic sensor and each column is a measurement in time. 
        PsiFull = [ [sensor1 time0, time1, time2, time3...], [sensor2 time0, time1, time2, time3]... ] should be in this shape! 
        The number of columns should be equal to the number of timesteps!
    
    total_current: float array
        total plasma current measurement should just be in the shape [measurement1, measurement2, measurement3] where each column
        is at a different time step, make sure that the number of columns also matches the number of timesteps. 
    
    reg_factor_fil : float
        regularization strength for the filaments
        
    reg_factor_wall : float
        regularization strength for the wall currents
        
    error_matrix: float array
        gives the sensor error for each point in time for each of the sensors. shape should be [[mag1 error at t0, mag2 error at t0, ...,   current error at t0], [mag1 error at t1, mag2 error at t1, .... ,current error at t1]]. so number of columns matches the number of magnetic sensors + 1 (for the total current), and the number of rows correspond to the number of timesteps 
    
    num_coils : int
        number of coils that are not filaments, i.e they do not contribute to the total plasma current but are within the Msc matrix
    
        
    Returns:
    ----------
    an array with the wall currents and the filaments currents over the entire time range. 
    '''
                
    # Ms corresponds to the wall mode weighting, Msc corresponds to the coil mode weighting 
    Ms = np.append(torus.Ms, np.zeros((torus.Ms.shape[0],1)),axis=1) #adding an extra row for the current constraint, set to 0
    new_col = np.ones((torus.Msc.shape[0], 1)) #adding an extra row for the current constraint 
    new_col[-num_coils:] = 0  # set the nonfilament coils rows to zero (like vcoils, oh coils, f coils etc.) 
    Msc = np.append(torus.Msc, new_col, axis=1) #appending the extra row for the current constraint, i.e filament currents add to total plasma current
    combined_matrix = np.vstack((Ms, Msc))
    combined_matrix = combined_matrix.T #stacking both matrices, wall and coils, and tranposing for shape matching 
    num_Ms = torus.Ms.shape[0]
    num_Msc = torus.Msc.shape[0]
    
    solutions_for_each_time = []
    for t in range(timesteps): #iterating through by every time step
        error_timestep = np.tile(error_matrix[t,:], (np.shape(combined_matrix)[1], 1)) # getting the errors of each 
        # sensor at that specific time step and tiling it such that it matches the shape of combined_matrix (before transpose) 
        A = combined_matrix*error_timestep.T #now multiplying the error to the combined_matrix (transposed for shape matching) 
        
        B = np.append(PsiFull[:,t], total_current[t]) #stacking the sensor signals and the total current at time t
        B = B*error_matrix[t,:].T #multiplying by the error matrix as well. 
        
        reg_identity_fil = reg_factor_fil*np.eye(A.shape[1]) # incorporating the regularization for filaments and wall
        reg_identity_wall = reg_factor_wall*np.eye(A.shape[1])
        A = np.vstack([A, reg_identity_wall[:num_Ms, :], reg_identity_fil[:num_Msc, :]]) 
        B = np.concatenate([B, np.zeros(A.shape[1])])
        AtA = A.T @ A
        AtB = A.T @ B
        solution = np.linalg.solve(AtA, AtB)    # least squares solver
        solutions_for_each_time.append(solution)
        
    solutions_for_each_time = np.array(solutions_for_each_time)       
    return solutions_for_each_time #return wall and filament solutions 