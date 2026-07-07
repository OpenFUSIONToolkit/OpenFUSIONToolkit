#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Core definitions for JAMfit - filament reconstruction
@authors Jamie Xia
@date July 2026
'''

## IMPORTING EXTERNAL LIBRARIES ##
import os
import numpy
import matplotlib.pyplot as plt

import math
from matplotlib.path import Path
from ._core import ThinCurr, ThinCurr_reduced
from .._core import OFT_env
from .sensor import save_sensors
from ..io import histfile
from ..TokaMaker._core import TokaMaker
from ..TokaMaker.meshing import load_gs_mesh

# ===============================
# JAMfit Base Class
# ===============================

class JAMfit():
    '''! Main class for JAMfit filament reconstruction.
    Manages setup, synthetic data generation, reduced model creation,
    and reconstruction of plasma filament currents using the
    both ThinCurr and TokaMaker.
    '''

    def __init__(self, xml_file, thincurr_meshfile, nthreads=None, oft_env=None):
        '''! Initialize the JAMfit object.

        @param xml_file str, path to the XML configuration file
        @param thincurr_meshfile str, path to the ThinCurr mesh file
        @param nthreads int, number of threads to use (required if oft_env is None)
        @param oft_env OFT_env, an existing OFT environment instance (optional)
        '''
        if oft_env is not None:
            self.myOFT = oft_env
        else:
            if nthreads is None:
                raise ValueError("Either oft_env or nthreads must be provided")
            self.myOFT = OFT_env(nthreads=nthreads)

        self.xml_file = xml_file
        self.thincurr_meshfile = thincurr_meshfile
        self.reduced_created_flag = False

    # =====================================
    # JAMfit Creation Relevant Classes
    # =====================================
    def set_xml(self, xml_file):
        '''! Update the XML configuration file path.

        @param xml_file str, path to the new XML configuration file
        '''
        self.xml_file = xml_file

    def set_sensors(self, sensor_array, floops_path='floops.loc'):
        '''! Save a sensor array to file and return the file path.

        @param sensor_array list, array of ThinCurr sensor objects (e.g. Mirnov, flux loops)
        @param floops_path str, output file path for sensor locations (default: 'floops.loc')
        @result str, path to the saved sensor file
        '''
        save_sensors(sensor_array, filename=floops_path)
        return floops_path
    
    def setup_JAMfit(self, floops_path, plot_files = None, use_legacy_io=False, hodlr_path = 'full_HOLDR_M.save'):
        '''! Set up the ThinCurr model, I/O, and sensor/coil mutual inductance matrices.
        
        @param floops_path str, path to the sensor locations file
        '''
        self.torus = ThinCurr(self.myOFT)
        self.torus.setup_model(mesh_file=self.thincurr_meshfile, xml_filename=self.xml_file)
        if plot_files is not None:
            self.torus.setup_io(basepath=plot_files, legacy_hdf5=use_legacy_io)
        else:
            self.torus.setup_io() 
        self.Msensor, self.Msc, self.sensor_obj = self.torus.compute_Msensor(floops_path)
        self.torus.compute_Mcoil(cache_file=hodlr_path)
        self.torus.compute_Lmat(use_hodlr=True, cache_file=hodlr_path)
        self.torus.compute_Rmat(copy_out=True)
        print('JAMfit setup complete.')

    def setup_fil_timeseries(self, time_array, totalip, coil_currs, r_list, z_list, sigma_r, sigma_z, rgrid, zgrid):
        '''! Build a combined coil + plasma current array for a time-dependent run.
        Generates a Gaussian plasma current distribution at each time step and
        appends it to the coil currents for use with run_td.

        @param time_array numpy.ndarray, array of time points
        @param totalip numpy.ndarray, total plasma current at each time point
        @param coil_currs numpy.ndarray, coil currents with shape (ntimes, ncoils)
        @param r_list list, R position of plasma centroid at each time step
        @param z_list list, Z position of plasma centroid at each time step
        @param sigma_r float, Gaussian spread in the R direction
        @param sigma_z float, Gaussian spread in the Z direction
        @param rgrid numpy.ndarray, R coordinates of the filament grid
        @param zgrid numpy.ndarray, Z coordinates of the filament grid
        @result numpy.ndarray, combined time + coil + plasma current array for run_td
        '''
        coil_curr_wtime = numpy.hstack([time_array, coil_currs])
        plasma_curr_wtime = setup_synthetic_current(time_array, totalip, sigma_r, sigma_z, r_list, z_list, rgrid, zgrid)
        plasma_curr = plasma_curr_wtime[:, 1:]
        final_coil_currs = numpy.hstack((coil_curr_wtime, plasma_curr))
        return final_coil_currs

    def gen_synthetic_data(self, coil_currs, dt, nsteps, verbose = False, hodlr_path = 'full_HOLDR_L.save',s_freq = 10, p_freq=10):
        '''! Run a synthetic time-dependent simulation and compute sensor signals.
        Computes the inductance and resistance matrices, runs the time-dependent
        simulation, plots results, and saves the B matrix and sensor history.

        @param coil_currs numpy.ndarray, combined coil + plasma current array (with time column)
        @param dt float, time step size in seconds
        @param nsteps int, number of time steps
        @param verbose bool, if True, returns additional plot data (default: False)
        @param hodlr_path str, path to save/load the HODLR matrices (default: 'full_HOLDR_L.save')
        @param s_freq int, frequency of sensor data output (default: 10)
        @param p_freq int, frequency of plotting during the run for plot_data (build_XDMF) (default: 10)
        @result returns the hist_file and optionally plot_data if verbose is True
        '''
        self.torus.run_td(dt, nsteps, coil_currs=coil_currs, sensor_obj=self.sensor_obj, status_freq= s_freq, plot_freq=p_freq)
        self.torus.plot_td(nsteps, sensor_obj=self.sensor_obj)
        hist_file = histfile(os.path.join(self.torus._io_basepath, 'floops.hist'))      
        if verbose: 
            plot_data = self.torus.build_XDMF()
            return hist_file, plot_data
        else: 
            return hist_file 

    def create_from_runTD_top_modes(self, num_modes, reduced_filename, coil_currs, dt, nsteps, initial_num_eigs=50, verbose=False, s_freq = 10, p_freq = 10):
        '''! Build a reduced model using the dominant eigenmodes from a full run.
        Computes eigenvalues, runs a preliminary reduced model to identify the
        most active modes by current amplitude, then builds a final reduced model
        using only those dominant modes.

        @param num_modes int, number of dominant modes to retain in the reduced model
        @param reduced_filename str, output filename for the reduced model (HDF5)
        @param coil_currs numpy.ndarray, combined coil + plasma current array (with time column)
        @param nsteps int, number of time steps for the preliminary run
        @param dt float, time step size in seconds
        @param num_eigs int, number of eigenmodes to compute initially (default: 50)
        @param verbose bool, if True, plots mode amplitudes and prints mode info (default: False)
        @param s_freq int, frequency of sensor data output (default: 10)
        @param p_freq int, frequency of plotting during the run for plot_data (build_XDMF) (default: 10)
        @result ThinCurr_reduced, the constructed reduced model object
        @result numpy.ndarray, sensor measurements from the time-dependent run
        @result dict, currents from the time-dependent run
        @result numpy.ndarray, eigenvectors corresponding to the selected dominant modes (if verbose is True)
        @result list, indices of the selected dominant modes (if verbose is True)
        @result list, maximum weight amplitudes of the selected dominant modes (if verbose is True)
        '''
        self.eig_vals, self.eig_vecs = self.torus.get_eigs(initial_num_eigs, False)
        torus_first_reduced = self.torus.build_reduced_model(
            self.eig_vecs, filename='first_reduced_model_temp.h5', sensor_obj=self.sensor_obj
        )
        sensors_measurement, currents = torus_first_reduced.run_td(dt, nsteps, coil_currs, status_freq=s_freq, plot_freq=p_freq)

        temp_curr = currents['curr']
        temp_curr = temp_curr[:, 0:initial_num_eigs]  # Only consider the modes we computed
        max_weights = [abs(temp_curr[:, i]).sum() for i in range(temp_curr.shape[1])] # total sum over time 

        top_modenum_indices = sorted(range(len(max_weights)), key=lambda i: max_weights[i], reverse=True)[:num_modes]

        eig_inds = []
        weight_amplitude = []

        if verbose:
            fig, ax = plt.subplots(1, 1)
            self.torus.build_XDMF()

        for i in range(temp_curr.shape[1]):
            if i in top_modenum_indices:
                eig_inds.append(i)
                weight_amplitude.append(max_weights[i])
                if verbose:
                    ax.semilogy(currents['time'], abs(currents['curr'][:, i]), label=f'Mode {i}')
                    print(f'Saved mode {i} has max weight {max_weights[i]}')
            else:
                if verbose:
                    ax.semilogy(currents['time'], abs(currents['curr'][:, i]), color='gray', alpha=0.3)
                    print("eig_inds:", eig_inds)
                    print("eig_vecs shape:", self.eig_vecs.shape)
                    print("eig_vecs[eig_inds] shape:", self.eig_vecs[eig_inds, :].shape)
        self.torus_reduced = self.create_reduced_model(self.eig_vecs[eig_inds, :], reduced_filename, compute_B=False)
        self.reduced_created_flag = True
        if verbose:
             return self.torus_reduced, sensors_measurement, currents, self.eig_vecs[eig_inds, :], eig_inds, weight_amplitude
        else:
            return self.torus_reduced, sensors_measurement, currents
    
    
    def create_reduced_model(self, eig_vecs, reduced_filename, compute_B=False):
        '''! Build a reduced model using specified eigenvectors.
        
        @param eig_vecs numpy.ndarray, array of eigenvectors to use for the reduced model
        @param reduced_filename str, output filename for the reduced model (HDF5)
        @param compute_B bool, if True, computes the B matrix for the reduced model (default: False)
        @result ThinCurr_reduced, the constructed reduced model object
        '''
        self.torus_reduced = self.torus.build_reduced_model(
            eig_vecs, filename=reduced_filename, compute_B=compute_B, sensor_obj=self.sensor_obj
        )
        print(f"Reduced model created with {eig_vecs.shape[0]} modes")
        self.reduced_created_flag = True
        return self.torus_reduced
    
    def add_freq_eigenvalues(self, specific_fil_array): 
        '''! Augment the eigenvector basis with frequency-response vectors for specific filaments.
        For each filament index provided, computes the steady-state frequency
        response driven by that filament's mutual inductance and appends the
        result to the stored eigenvector matrix.

        @param specific_fil_array list, indices of filaments to compute frequency responses for
        @result numpy.ndarray, updated eigenvector matrix including frequency-response vectors, appended to the end of the existing eigenvectors
        '''
        from IPython.display import clear_output
        if self.eig_vecs is None:
            raise ValueError("Eigenvalues and eigenvectors have not been computed yet.")

        eig_vecs_wfreq = numpy.copy(self.eig_vecs)
        for target_fil in specific_fil_array:
            Mcoil = self.torus.compute_Mcoil()
            driver = numpy.zeros((2, self.torus.nelems))
            driver[0, :] = Mcoil[target_fil, :]
            result = self.torus.compute_freq_response(driver, freq=1.E3)
            eig_vecs_wfreq = numpy.concatenate((eig_vecs_wfreq, result[:1, :]), axis=0)
            clear_output(wait=True)

        self.eig_vecs = eig_vecs_wfreq
        return eig_vecs_wfreq

    # ============================================
    # JAMfit Reconstruction Relevant Classes
    # ============================================

    def initialize_reduced_model(self, reduced_filename):
        '''! Load a previously saved reduced model from file.
        
        @param reduced_filename str, path to the HDF5 reduced model file
        @result str, confirmation message on successful load
        '''
        self.torus_reduced = ThinCurr_reduced(reduced_filename)
        self.reduced_created_flag = True
        return "Reduced model initialized from file."

    def run_reconstruction_lstsq(self, Psi_at_time, ip_at_time, num_non_fil_coils, coil_curr_at_time, ip_weight, sigma, reg_factor_fil, reg_factor_wall, num_sensors = None):
        '''! Run the filament current reconstruction with the lstsq method.
        
        @param Psi_at_time numpy.ndarray, sensor flux measurements at time
        @param ip_at_time float, total plasma current measurement at time
        @param num_non_fil_coils int, number of non-filament coils in the system
        @param ip_weight float, weight for the total plasma current in the reconstruction
        @param sigma numpy.ndarray, array of standard deviations for each sensor measurement (for weighting)
        @param reg_factor_fil float, regularization factor for filament currents
        @param reg_factor_wall float, regularization factor for wall currents
        @result tuple, containing the filament solution, wall solution, residual, Ax, and B
        '''
        if not self.reduced_created_flag:
            raise ValueError("Reduced model has not been created yet. Please create or initialize a reduced model before running reconstruction.")
        if num_sensors is None: 
            num_sensors = self.torus_reduced.Ms.shape[1] 
        # intializing the Ms and Msc matrices with the appropriate weighting and scaling based off sigma
        Ms_weighted = self.torus_reduced.Ms[:, :num_sensors]/sigma[:]
        Msc_weighted = self.torus_reduced.Msc[:, :num_sensors]/sigma[:]
        Msc_weighted_fil = Msc_weighted[num_non_fil_coils:, :]

        # This section of code scales the total ip constraint row of the matrix to ensure it has a comparable influence on the 
        # least squares solution as the magnetic measurements, based on the provided ip_weight and the magnitude of ip_at_time
        ip_row_norm = abs(ip_at_time)
        ip_col_ms = numpy.zeros((Ms_weighted.shape[0], 1))
        ip_col_msc = numpy.ones((Msc_weighted_fil.shape[0], 1)) * (ip_weight / ip_row_norm)

        # appending the ip constraint as an additional row to the Ms and Msc matrices, with appropriate scaling
        Ms_final = numpy.append(Ms_weighted, ip_col_ms, axis=1)
        Msc_final = numpy.append(Msc_weighted_fil, ip_col_msc, axis=1)
        combined_matrix_A = numpy.vstack((Ms_final, Msc_final)).T

        # here we prepare the the B vector by subtacting the non plasma filament contribution from the sensor measurements
        # we also scale by sigma here as well (to ensure magnetic sensor signals are normalized to each other - one sensor doesnt dominate)
        # finally we append the ip cosntraint row as well 
        B_weighted = (Psi_at_time/ sigma[:]) - coil_curr_at_time @ self.torus_reduced.Msc[:num_non_fil_coils, :num_sensors] 
        B_weighted = numpy.append(B_weighted, [ip_at_time * ip_weight/ip_row_norm])  # ip_at_time/ip_row_norm reduces to sign(ip_at_time) — magnitude is encoded in A side column

        # here we apply tikonov regularization to both the filament and wall component of the lstq 
        # note that we must use the unmodified shapes of Ms to construct the identity matrices for regularization
        # the first half of the identity matrix corresponds to the wall currents and the second half corresponds to the filament currents, so we slice accordingly when stacking them below
        # we then add the reg rows to B as well to match the extra rows that we stacked to A 
        num_to_solve = Ms_final.shape[0] + Msc_final.shape[0]
        reg_identity = numpy.eye(num_to_solve)
        A = numpy.vstack([combined_matrix_A, reg_factor_wall * reg_identity[0:Ms_weighted.shape[0], :], reg_factor_fil * reg_identity[Ms_weighted.shape[0]:, :]])
        B= numpy.concatenate([B_weighted, numpy.zeros(A.shape[1])])

        # solving the least squares problem
        AtA = A.T @ A
        AtB = A.T @ B
        solution = numpy.linalg.solve(AtA, AtB)
        Ax = numpy.dot(A, solution)
        residual = numpy.sqrt(numpy.sum((B - Ax)**2))
        solution_wall = solution[:Ms_weighted.shape[0]]
        solution_fil = solution[Ms_weighted.shape[0]:]

        return solution_fil, solution_wall, residual, Ax, B


    def prepare_tsvd_laplace(self, sigma, num_non_fil_coils, num_sensors, rgrid, zgrid, nModes, verbose = False):
        '''! Prepare matrices and projections for the svd + laplacian reconstruction method.
        
        @param sigma numpy.ndarray, array of standard deviations for each sensor measurement (for weighting)
        @param num_non_fil_coils int, number of non-filament coils in the system
        @param num_sensors int, number of real sensors in the system
        @param rgrid numpy.ndarray, R coordinates of the filament grid
        @param zgrid numpy.ndarray, Z coordinates of the filament grid
        @param nModes int, number of SVD modes to truncate to
        @param verbose bool, if True, plots singular values (default: False)
        ''' 
        
        # getting laplacina matrix for smoothing purposes
        lap_mat, N = get_laplace_matrix(rgrid, zgrid, verbose)
        num_Ms = self.torus_reduced.Ms.shape[0]

        # intializing the Ms and Msc matrices with the appropriate weighting and scaling based off sigma
        Msc_fil_weighted =  self.torus_reduced.Msc[num_non_fil_coils:, :num_sensors]/sigma[:]
        Msc_coils_weighted = self.torus_reduced.Msc[:num_non_fil_coils, :num_sensors]/sigma[:]
        Ms_weighted = self.torus_reduced.Ms[:, :num_sensors]/sigma[:]

        # break the problem into just the filaments and find svd modes (we do not solve for the shaping coil currents during the reconstruction)
        U, S, Vh = numpy.linalg.svd(Msc_fil_weighted, full_matrices=False)

        # intializing the least squares matrix for the truncated SVD solution 
        ls_mat_fil = Msc_fil_weighted.T @ U[:, :nModes]
        ls_mat_wall = Ms_weighted.T 
        ls_mat = numpy.hstack([ls_mat_wall, ls_mat_fil])

        # projecting the laplacian onto the TSVD space to get a reg term that smooths in the physical filament space
        lap_proj = lap_mat @ U[:, :nModes]

        # getting the ip constraint row in the TSVD space, since U is already normalized to the sensor signals, we only need to normalize the ip row to itself
        ip_row_fil = U[:, :nModes].sum(axis=0, keepdims=True)
        ip_row_fil_norm = numpy.linalg.norm(ip_row_fil)
        ip_row = numpy.hstack([numpy.zeros((1, num_Ms)), ip_row_fil/ip_row_fil_norm])


        if verbose: 
            plt.figure(figsize=(8, 5))
            plt.semilogy(S/S[0], marker='o')
            plt.title('Normalized Singular Values of plasma filament contribution to sensor signals')
            plt.xlabel('Mode Index')
            plt.ylabel('Singular Value (log scale)')
            plt.grid(True)
            plt.show()
        return Msc_coils_weighted, Ms_weighted, U[:, :nModes], ls_mat, ls_mat_fil, lap_proj, ip_row, N
    

    
    def run_reconstruction_tsvd_laplace(self, Psi_at_time, ip_at_time, coil_curr_at_time, Msc_coils, Ms, U_trun, ls_mat, ls_mat_fil, lap_proj, ip_row, N, sigma, lam=None, lap_lam=1e-8, reg_wall=1e-5):
        '''! Run the filament current reconstruction with the svd + laplacian method.
        
        @param Psi_at_time numpy.ndarray, sensor flux measurements at time
        @param ip_at_time numpy.ndarray, plasma current measurements at time
        @param coil_curr_at_time numpy.ndarray, coil currents at time
        @param Msc_coils numpy.ndarray, coil matrix
        @param Ms numpy.ndarray, wall matrix
        @param U_trun numpy.ndarray, truncated singular vectors
        @param ls_mat numpy.ndarray, least squares matrix
        @param ls_mat_fil numpy.ndarray, filament least squares matrix
        @param lap_proj numpy.ndarray, laplacian projection matrix
        @param ip_row numpy.ndarray, ip constraint row, normalized to both sigma and itself
        @param N int, number of filaments
        @param sigma numpy.ndarray, standard deviations for each sensor (for weighting)
        @param lam float, ip weight (usually calculated automatically but can be set manually)
        @param lap_lam float, laplacian regularization parameter
        @param reg_wall float, wall regularization parameter
        @result tuple, containing the filament solution, wall solution, Ax, and diagnostics dictionary
        '''
        num_Ms = Ms.shape[0]

        # taking out coil contributions from the magnetic sensor signals 
        B_weighted = (Psi_at_time)/sigma[:] - coil_curr_at_time @ Msc_coils

        if lam is None: 
            # calculating the weight of the ip constraint row based on the magnitudal difference between magnetic sensor signals and the total plasma current
            # ensures that they are on the same order of magnitude for the least squares solution 
            # Compare typical sensor signal magnitude to IP magnitude
            magnitude_diff = math.floor(math.log10(numpy.mean(numpy.abs(B_weighted)) / (abs(ip_at_time) + 1e-30)))
            lam = 100 * 10**magnitude_diff # note that I multiply by 100 to give the ip constraint slighly more weight as the magnetic sensor signals have more rows over the singular total plasma current row
        
        ip_scale = numpy.linalg.norm(U_trun.sum(axis=0)) # this is for scaling the ip constraint row to the svd space for the totalip on the B side of Ax=B
        
        # constructing the A matrix by stacking the svd matrix with the regularization rows  
        reg_mat = numpy.vstack([
            ls_mat,                                    # measurements
            lam * ip_row,                                                        # Ip constraint 
            numpy.hstack([reg_wall * numpy.eye(num_Ms), numpy.zeros((num_Ms, ls_mat_fil.shape[1]))]),      # wall Tikhonov
            numpy.hstack([numpy.zeros((N, num_Ms)), lap_lam * lap_proj])                   # fil Laplacian
        ])

        # constructing the B vector by stacking the measurement vector with the ip constraint and zeros for the regularization rows
        psi_reg = numpy.concatenate([
            B_weighted,
            numpy.array([lam * ip_at_time / ip_scale]),
            numpy.zeros(num_Ms),
            numpy.zeros(N)
        ])

        # solving using least squares explicitly of Ax=B
        AtA = reg_mat.T @ reg_mat
        AtB = reg_mat.T @ psi_reg
        curr_weights = numpy.linalg.solve(AtA, AtB)

        # projecting the solution back to the physical space 
        curr_expand = U_trun @ curr_weights[num_Ms:]  # only the filament part contributes to the expansion
        wall_expand = curr_weights[:num_Ms]  # wall coefficients


        # calculating diagnostics
        Ax          = ls_mat @ curr_weights  
        ip_reconstructed = curr_expand.sum()
        ip_actual        = ip_at_time
        ip_error_pct     = 100 * numpy.abs(ip_reconstructed - ip_actual) / (abs(ip_actual) + 1e-30)
        fit_residual     = numpy.linalg.norm(Ax - B_weighted)
        fit_residual_nonorm = Ax - B_weighted

        diagnostics = {
            'lam':              lam,
            'lap_lam':          lap_lam,
            'ip_reconstructed': ip_reconstructed,
            'ip_actual':        ip_actual,
            'ip_error_pct':     ip_error_pct,
            'fit_residual':     fit_residual,
            'curr_weights':     curr_weights,
            'fit_residual_nonorm': fit_residual_nonorm,
            'ip_row':           ip_row,
        }

        return curr_expand, wall_expand, Ax, diagnostics
    
    def run_reconstruction_laplace(self, Psi_at_time, ip_at_time, coil_curr_at_time, sigma, 
                                    num_non_fil_coils, num_real_sensors, rgrid, zgrid, 
                                    lam=None, lap_lam=1e-6, reg_wall=1e-5, sigma_r =1e-1, sigma_z =1e-1, gaussian = False, verbose=False):
        '''! Run filament current reconstruction using direct Laplacian regularization.
        Solves directly in the physical filament space (no SVD projection), applying
        Laplacian smoothing to filaments and Tikhonov regularization to wall currents.

        @param Psi_at_time numpy.ndarray, sensor flux measurements at time
        @param ip_at_time float, total plasma current measurement at time
        @param coil_curr_at_time numpy.ndarray, coil currents at time
        @param sigma numpy.ndarray, standard deviations for each sensor (for weighting)
        @param num_non_fil_coils int, number of non-filament coils in the system
        @param rgrid numpy.ndarray, R coordinates of the filament grid
        @param zgrid numpy.ndarray, Z coordinates of the filament grid
        @param lam float, ip constraint weight (auto-calculated if None)
        @param lap_lam float, Laplacian regularization parameter (default: 1e-6)
        @param reg_wall float, Tikhonov regularization for wall currents (default: 1e-5)
        @param sigma_r float, standard deviation in the R direction for the Gaussian weighting (optional)
        @param sigma_z float, standard deviation in the Z direction for the Gaussian weighting (optional)
        @param gaussian bool, if True uses Gaussian weighting (default: False)
        @param verbose bool, if True prints debug info (default: False)
        @result tuple, containing the filament solution, wall solution, Ax, and diagnostics dictionary
        '''

        # intializing the Ms and Msc matrices with the appropriate weighting and scaling based off sigma
        Ms_weighted       = self.torus_reduced.Ms[:, :num_real_sensors] / sigma[:]                          # (n_wall, n_sensors)
        Msc_weighted      = self.torus_reduced.Msc[:, :num_real_sensors] / sigma[:]                         # (n_coils_total, n_sensors)
        Msc_coils_weighted         = Msc_weighted[:num_non_fil_coils, :num_real_sensors]                       # (n_shaping_coils, n_sensors)
        Msc_fil_weighted  = Msc_weighted[num_non_fil_coils:, :num_real_sensors]                       # (n_fil, n_sensors)
        n_fil             = Msc_fil_weighted.shape[0]
        num_Ms = Ms_weighted.shape[0]

        # taking out coil contributions from the magnetic sensor signals 
        B_weighted = (Psi_at_time)/sigma[:] - coil_curr_at_time @ Msc_coils_weighted

        # auto-calculating the weight of the ip constraint row if not provided 
        if lam is None: 
            # calculating the weight of the ip constraint row based on the magnitudal difference between magnetic sensor signals and the total plasma current
            # ensures that they are on the same order of magnitude for the least squares solution 
            # Compare typical sensor signal magnitude to IP magnitude
            magnitude_diff = math.floor(math.log10(numpy.mean(numpy.abs(B_weighted)) / (abs(ip_at_time) + 1e-30)))
            lam = 100 * 10**magnitude_diff # note that I multiply by 100 to give the ip constraint slighly more weight as the magnetic sensor signals have more rows over the singular total plasma current row

        # each filament contributes equally to Ip, so the row is all-ones for filaments, zeros for wall
        ip_row_fil  = numpy.ones((1, n_fil))
        ip_row_norm = numpy.linalg.norm(ip_row_fil)
        ip_row      = numpy.hstack([numpy.zeros((1, num_Ms)), ip_row_fil / ip_row_norm]) # (1, n_wall + n_fil)

        if gaussian: 
            lap_mat, N = get_gaussian_lap_mat(rgrid, zgrid, sigma_r, sigma_z)                # (n_fil, n_fil)
        # getting laplacian matrix for smoothing the filament solution in space 
        lap_mat, N = get_laplace_matrix(rgrid, zgrid, verbose=verbose)                # (n_fil, n_fil)

        # building A matrix's measurements corresponding to sensors: [wall | filament] columns ---
        meas_block = numpy.hstack([Ms_weighted.T, Msc_fil_weighted.T])                # (n_sensors, n_wall + n_fil)

        # Ip constraint row
        ip_block = lam * ip_row                                                        # (1, n_wall + n_fil)

        # Tikhonov on wall
        wall_reg = reg_wall * numpy.hstack([
            numpy.eye(num_Ms),
            numpy.zeros((num_Ms, n_fil))
        ])                                                                             # (n_wall, n_wall + n_fil)

        # Laplacian on filaments
        fil_lap = lap_lam * numpy.hstack([
            numpy.zeros((N, num_Ms)),
            lap_mat
        ])                                                                             # (n_fil, n_wall + n_fil)
        # constructing full A matrix by stacking the measurement block with the regularization blocks
        A = numpy.vstack([meas_block, ip_block, wall_reg, fil_lap])

        # building B matrix by stacking the sensor measurements (with coil contributions removed) with the ip constraint and zeros for regularization rows
        psi_reg = numpy.concatenate([
            B_weighted,
            numpy.array([lam * ip_at_time / ip_row_norm]),
            numpy.zeros(num_Ms),
            numpy.zeros(N)
        ])

        # Solving Ax=B
        AtA = A.T @ A
        AtB = A.T @ psi_reg
        solution = numpy.linalg.solve(AtA, AtB)
        # extracting wall and filament currents from the solution vector
        curr_wall = solution[:num_Ms]
        curr_fil  = solution[num_Ms:]

        # Diagnostics
        Ax = meas_block @ solution
        ip_reconstructed  = curr_fil.sum()
        ip_error_pct      = 100 * numpy.abs(ip_reconstructed - ip_at_time) / (abs(ip_at_time) + 1e-30)
        fit_residual      = numpy.linalg.norm(Ax - B_weighted)
        fit_residual_nonorm = Ax - B_weighted

        diagnostics = {
            'lam':                  lam,
            'lap_lam':              lap_lam,
            'reg_wall':             reg_wall,
            'ip_reconstructed':     ip_reconstructed,
            'ip_actual':            ip_at_time,
            'ip_error_pct':         ip_error_pct,
            'fit_residual':         fit_residual,
            'fit_residual_nonorm':  fit_residual_nonorm,
            'solution':             solution,
        }

        if verbose:
            print(f"lam = {lam:.3e}, lap_lam = {lap_lam:.3e}, reg_wall = {reg_wall:.3e}")
            print(f"Ip reconstructed: {ip_reconstructed:.3f} A | Ip actual: {ip_at_time:.3f} A | Error: {ip_error_pct:.2f}%")
            print(f"Fit residual: {fit_residual:.4e}")

        return curr_fil, curr_wall, Ax, diagnostics

    # =========================================================
    # JAMfit Reconstruction Post Processesing and Visualization
    # =========================================================

    def get_wall_psi_tidx(self, num_sensors, solution_wall_at_time): 
        '''! Compute the wall contribution to the flux at a given time index.
        Uses the wall current solution and the wall mutual inductance matrix to calculate
        
        @param num_sensors int, number of real sensors to filter out 
        @param solution_wall_at_time numpy.ndarray, wall current potential solution at the given time index
        @result numpy.ndarray, wall contribution to the flux at the given time index'''
        wall_psi_probes = self.torus_reduced.Ms[:, num_sensors:]/(2*numpy.pi)
        wall_psi_at_time = solution_wall_at_time @ wall_psi_probes
        return wall_psi_at_time 
    
    def post_process_tidx(self, filaments_at_time, wall_psi, coil_curr_dict, rgrid, zgrid, meshfile_tokamaker, B0, R0, myOFT, verbose = False):
        '''! Post-process the results at a given time index to compute plasma parameters and visualize.
        
        @param filaments_at_time numpy.ndarray, filament currents at the given time index
        @param coil_curr_dict dict, dictionary of coil currents at the given time index
        @param rgrid numpy.ndarray, R coordinates of the filament mesh
        @param zgrid numpy.ndarray, Z coordinates of the filament mesh
        @param meshfile_tokamaker str, path to the Tokamaker mesh file for equilibrium reconstruction
        @param wall_psi numpy.ndarray, precomputed wall contribution to the flux
        @param B0 float, reference magnetic field strength for equilibrium reconstruction
        @param R0 float, reference major radius for equilibrium reconstruction
        @param myOFT OFT_env, the Open FUSION Toolkit environment instance
        @param verbose bool, if True, plots the equilibrium and LCFS (default: False)
        @result tuple, containing LCFS points, limiting points, psi at LCFS, total psi, q values, q95, internal inductance, current centroid, area centroid, and area'''

        #intialize tokamaker 
        mygs = TokaMaker(myOFT)
        mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(meshfile_tokamaker)
        mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
        mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
        mygs.setup(order=2, F0= B0 * R0)# F0 = B0 * R0 
        limiter = mygs.lim_contour

        # calculate psi from plasmas
        fil_points = list(zip(rgrid, zgrid))
        psi_fil = []
        for filcount, (r, z) in enumerate(fil_points):
            mygs.set_coil_currents()
            if filaments_at_time[filcount] > 0:
                mygs.set_targets(Ip=filaments_at_time[filcount])
                mygs.init_psi(r, z, 0.3, 1.0, 0.0) 
                psi_fil.append(mygs.get_psi(False))
        psi_fil = numpy.array(psi_fil)
        psi_total_fil = numpy.sum(psi_fil, axis=0)
        ip = numpy.sum(filaments_at_time)
        
        # Calculate psi from coils 
        mygs.set_coil_currents(coil_curr_dict)
        psi_vf_eq_obj = mygs.vac_solve()
        psi_vf = psi_vf_eq_obj.get_psi(normalized=False) 


        # Summate Psis 
        total_psi = psi_total_fil + psi_vf + wall_psi
        mygs.set_psi(total_psi, update_bounds = True)

        # getting relevant values 
        lcfs_points = mygs.trace_surf(1) # Trace LCFS
        if lcfs_points is not None:
            psiatlcfs = mygs.psinorm_to_absolute(1)

        if lcfs_points is None: 
            lcfs_points = mygs.trace_surf(0.99)
            psiatlcfs = mygs.psinorm_to_absolute(0.99) 

        limiting_pts = [] 
        q_vals = None 
        q95 = numpy.nan 
        internal_inductance = numpy.nan 
        current_cent = calc_current_centroid(rgrid, zgrid, filaments_at_time)
        area_cent = numpy.nan 
        area = numpy.nan 
     
        if lcfs_points is not None:
            lim_R, lim_Z, _ = find_limiting_point(lcfs_points, limiter, touch_tol=0.0005)
            limiting_pts.append([lim_R, lim_Z])
            _, q_vals, _, _, _, _= mygs.get_q() 

            try:
                internal_inductance = mygs.get_stats(beta_Ip = ip)['l_i']
            except ZeroDivisionError:
                print(f"Couldn't compute li, Ip was {ip}, skipping this one")
                internal_inductance = numpy.nan

            _, q95, _ , _, _, _ = mygs.get_q(psi=0.95)
            area, area_cent = calc_lcfs_geo(lcfs_points) 
        else: 
            psiatlcfs = None 
            limiting_pts.append(None)

        if verbose: 
            fig, ax = plt.subplots(1,1) 
            mygs.plot_machine(fig,ax, cond_color='blue') 
            mygs.plot_psi(fig,ax ,plasma_nlevels = 75, normalized=False)
            plt.gca().set_aspect('equal', adjustable='box')
            if lcfs_points is not None:
                ax.plot(lcfs_points[:,0], lcfs_points[:,1], 'r--', label='LCFS')
    
        return lcfs_points, limiting_pts, psiatlcfs, total_psi, q_vals, q95, internal_inductance, current_cent, area_cent, area

# ===============================
# JAMfit Helper Functions
# ===============================

def setup_synthetic_current(timepoints, ip_list, sigma_r, sigma_z, r0, z0, rgrid, zgrid):
    '''! Generate synthetic filament currents using a Gaussian plasma distribution.
    At each time step, spreads the total plasma current across the filament grid
    using a 2D Gaussian centered at (r0, z0) with widths (sigma_r, sigma_z).

    @param timepoints list or numpy.ndarray, array of time values
    @param ip_list list or numpy.ndarray, total plasma current at each time step
    @param sigma_r float, Gaussian width in the R direction
    @param sigma_z float, Gaussian width in the Z direction
    @param r0 list or numpy.ndarray, R position of the plasma centroid at each time step
    @param z0 list or numpy.ndarray, Z position of the plasma centroid at each time step
    @param rgrid numpy.ndarray, R coordinates of the filament mesh
    @param zgrid numpy.ndarray, Z coordinates of the filament mesh
    @result numpy.ndarray, shape (ntimes, 1 + nfilaments), time column followed by filament currents
    '''
    coil_curr = []
    for i in range(len(timepoints)):
        gaussian_raw = numpy.exp(
            -((rgrid - r0[i])**2 / (2 * sigma_r**2) + (zgrid - z0[i])**2 / (2 * sigma_z**2))
        )
        gaussian_values = ip_list[i] * (gaussian_raw / numpy.sum(gaussian_raw))
        coil_curr.append(gaussian_values)

    coil_curr = numpy.array(coil_curr)
    time_column = numpy.array(timepoints).reshape(-1, 1)
    coil_curr = numpy.hstack((time_column, coil_curr))
    return coil_curr

def interpolate_total_current(coil_currs, nsteps, verbose=False):
    '''! Interpolate total plasma current to a higher-resolution time grid.
    Sums the sensor currents at each time step and interpolates the total
    onto a finer time grid using linear interpolation.

    @param coil_currs numpy.ndarray, shape (ntimes, nsensors+1), first column is time
    @param nsteps int, number of high-resolution steps between first and last time
    @param verbose bool, if True, plots the original and interpolated total current (default: False)
    @result tuple of (high_res_time, total_current_high_res) as numpy.ndarrays
    '''
    times = coil_currs[:, 0]
    sensor_currents = coil_currs[:, 1:]

    high_res_time = numpy.linspace(times[0], times[-1], nsteps + 1)

    interpolated_currents = numpy.array([
        numpy.interp(high_res_time, times, sensor)
        for sensor in sensor_currents.T
    ]).T  # shape: (nsteps+1, nsensors)

    total_current_high_res = numpy.sum(interpolated_currents, axis=1)

    if verbose:
        plt.figure(figsize=(8, 5))
        plt.scatter(times, numpy.sum(sensor_currents, axis=1), color='red', label='Original Data')
        plt.plot(high_res_time, total_current_high_res, color='blue', label='Interpolated Total Current')
        plt.xlabel('Time [s]')
        plt.ylabel('Total Current [A]')
        plt.title('Total Current vs Time with Higher Resolution')
        plt.legend()
        plt.grid(True)
        plt.show()

    return high_res_time, total_current_high_res


def get_laplace_matrix(rgrid, zgrid, verbose=False):
    """
    FD Laplacian for D-shaped grids.

    @param rgrid numpy.ndarray, R coordinates of the filament grid
    @param zgrid numpy.ndarray, Z coordinates of the filament grid
    @param verbose bool, if True, prints debug information (default: False)
    rgrid, zgrid: 1D arrays of valid point coordinates inside the limiter
    """
    # Infer the full rectangular grid from unique values
    nr_arr = numpy.sort(numpy.unique(rgrid))
    nz_arr = numpy.sort(numpy.unique(zgrid))
    dr = nr_arr[1] - nr_arr[0] if len(nr_arr) > 1 else 1.0
    dz = nz_arr[1] - nz_arr[0] if len(nz_arr) > 1 else 1.0

    if verbose:
        print(f"dr: {dr}, dz: {dz}, nr: {len(nr_arr)}, nz: {len(nz_arr)}")
        print(f"Valid points: {len(rgrid)} out of {len(nr_arr)*len(nz_arr)} rectangular")

    # Map r,z values to grid indices
    r_to_idx = {round(r, 8): i for i, r in enumerate(nr_arr)}
    z_to_idx = {round(z, 8): j for j, z in enumerate(nz_arr)}

    # Build lookup: (ir, iz) -> index in the valid point list
    N = len(rgrid)
    valid_set = {}
    for k in range(N):
        ir = r_to_idx[round(rgrid[k], 8)]
        iz = z_to_idx[round(zgrid[k], 8)]
        valid_set[(ir, iz)] = k

    if verbose:
        print(f"Built valid_set with {len(valid_set)} points")

    lap_mat = numpy.zeros((N, N))

    for k in range(N):
        ir = r_to_idx[round(rgrid[k], 8)]
        iz = z_to_idx[round(zgrid[k], 8)]

        for dir_r, dir_z, weight in [
            ( 1,  0, 1.0/dr**2),
            (-1,  0, 1.0/dr**2),
            ( 0,  1, 1.0/dz**2),
            ( 0, -1, 1.0/dz**2),
        ]:
            neighbor = (ir + dir_r, iz + dir_z)
            if neighbor in valid_set:
                # Interior neighbor — standard FD
                lap_mat[k, valid_set[neighbor]] = weight
                lap_mat[k, k] -= weight
            # else: boundary point, Neumann BC, skip

    return lap_mat, N


def find_limiting_point(lcfs_points, limiter, touch_tol=0.005):
    
    dists = numpy.min(
        numpy.linalg.norm(lcfs_points[:, None, :] - limiter[None, :, :], axis=2),
        axis=1
    )
    
    idx = numpy.argmin(dists)
    min_dist = dists[idx]
    
    if min_dist > touch_tol:
        return None, None, None
    
    return lcfs_points[idx, 0], lcfs_points[idx, 1], min_dist

    
def calc_lcfs_geo(lcfs_points):
    """
    Calculate area and geometric centroid of the LCFS.

    @param lcfs_points numpy.ndarray, shape (N, 2), array of LCFS points with columns [R, Z]
    @result tuple, area enclosed by LCFS, R_c, and Z_c where R_c and Z_c are the coordinates of the geometric centroid of the LCFS.
    """
    if lcfs_points is None:
        return None, None, None

    R = numpy.array(lcfs_points[:, 0])
    Z = numpy.array(lcfs_points[:, 1])

    # Ensure contour is closed
    if not (numpy.isclose(R[0], R[-1]) and numpy.isclose(Z[0], Z[-1])):
        R = numpy.append(R, R[0])
        Z = numpy.append(Z, Z[0])

    cross = R[:-1] * Z[1:] - R[1:] * Z[:-1]    # Calculates area using Shoelace method, Shoelace terms

    # Area
    area = 0.5 * numpy.abs(numpy.sum(cross))

    # Centroid
    R_c = numpy.abs(numpy.sum((R[:-1] + R[1:]) * cross)) / (6 * area)
    Z_c = numpy.abs(numpy.sum((Z[:-1] + Z[1:]) * cross)) / (6 * area)

    return area, (R_c, Z_c) 

def calc_current_centroid(R_fil, Z_fil, I):
    '''! Calculate the current centroid of a filament grid using signed current weights
        This variant uses signed currents as weights, so opposing currents can partially
        cancel out the centroid position. Useful for computing the net moment of the
        current distribution.
        
        Calculation: R_c = sum(I_i * R_i) / sum(I_i)
                    Z_c = sum(I_i * Z_i) / sum(I_i)
        
        @param R 1D or 2D array of R (radial) coordinates of filament points `[m]`
        @param Z 1D or 2D array of Z (vertical) coordinates of filament points `[m]`
        @param I 1D or 2D array of currents at each filament point `[A]`, must match shape of R and Z
        @result R coordinate of the current centroid `[m]`
        @result Z coordinate of the current centroid `[m]`
        @result Net current (algebraic sum) `[A]`
    '''
    R = numpy.asarray(R_fil)  # Ensure inumpyuts are numpy arrays
    Z = numpy.asarray(Z_fil)
    I = numpy.asarray(I)
    
    # Verify shapes match
    if not (R.shape == Z.shape == I.shape):
        raise ValueError(
            f"R, Z, and I must have the same shape. "
            f"Got R: {R.shape}, Z: {Z.shape}, I: {I.shape}"
        )
    
    R_flat = R.ravel()
    Z_flat = Z.ravel()
    I_flat = I.ravel()
    
    I_net = numpy.sum(I_flat)
    if I_net == 0:
        raise ValueError("Net current is zero - cannot calculate centroid")
    
    R_centroid = numpy.sum(I_flat * R_flat) / I_net
    Z_centroid = numpy.sum(I_flat * Z_flat) / I_net

    return (R_centroid, Z_centroid)



def get_inside_limiter_pts(meshfile_tokamaker, myOFT, stride=1, verbose=False):
    ''' ! Get the points of the filament grid that are inside the limiter contour defined in the Tokamaker mesh file. For computing flux contribution from wall

    @param meshfile_tokamaker str, path to the Tokamaker mesh file
    @param myOFT OFT_env, the Open FUSION Toolkit environment instance
    @param stride int, subsampling stride for the points inside the limiter (default: 1)
    @param verbose bool, if True, plots the grid points and limiter contour (default: False)
    @result tuple, containing the points inside the limiter, a boolean mask of points inside the limiter, and the full grid points
    '''
        
    mygs = TokaMaker(myOFT)
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(meshfile_tokamaker)
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.setup(order=2, F0=1 * 1)

    r_pts_grid = mygs.r[:, 0]
    z_pts_grid = mygs.r[:, 1]
    r_pts_lim = mygs.lim_contour[:, 0]
    z_pts_lim = mygs.lim_contour[:, 1]
    lim_path = Path(numpy.column_stack([r_pts_lim, z_pts_lim]))
    all_pts = numpy.column_stack([r_pts_grid, z_pts_grid])

    inside_mask = lim_path.contains_points(all_pts)  # full inside-limiter mask
    inside_idx = numpy.where(inside_mask)[0]

    # subsample the indices that actually get a flux loop
    sparse_idx = inside_idx[::stride]
    sparse_mask = numpy.zeros_like(inside_mask)
    sparse_mask[sparse_idx] = True

    inside_lim_pts = numpy.column_stack((mygs.r[sparse_idx, 0], mygs.r[sparse_idx, 1]))

    if verbose:
        fig, ax = plt.subplots()
        ax.scatter(r_pts_grid, z_pts_grid, c=sparse_mask, cmap='RdYlGn', s=10, alpha=0.7)
        ax.plot(r_pts_lim, z_pts_lim, 'k-', label='Limiter Contour')
        ax.set_xlabel('R')
        ax.set_ylabel('Z')
        ax.set_title(f'Grid Points With Flux Loops (stride={stride})')
        ax.legend()
        plt.colorbar(ax.collections[0], ax=ax, label='Has flux loop')
        plt.tight_layout()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    return inside_lim_pts, sparse_mask, inside_mask, mygs.r


def get_gaussian_lap_mat(rgrid, zgrid, sigma_r=8e-2, sigma_z=8e-2):
    '''! Construct a Gaussian-weighted Laplacian regularization matrix for the filament grid (Experimental).
    Each filament point is connected to its neighbors with weights based on a Gaussian function of the distance
    @param rgrid numpy.ndarray, R coordinates of the filament grid
    @param zgrid numpy.ndarray, Z coordinates of the filament grid
    @param sigma_r float, standard deviation in the R direction for the Gaussian weighting
    @param sigma_z float, standard deviation in the Z direction for the Gaussian weighting
    @result lap_mat numpy.ndarray, the constructed Gaussian-weighted Laplacian matrix 
    '''
    N = len(rgrid)
    dr = rgrid[:, None] - rgrid[None, :]  # (N, N)
    dz = zgrid[:, None] - zgrid[None, :]  # (N, N)
    
    W = numpy.exp(-0.5 * ((dr / sigma_r)**2 + (dz / sigma_z)**2))
    numpy.fill_diagonal(W, 0)  # no self-coupling
    
    D = numpy.diag(W.sum(axis=1))
    lap_mat = D - W
    
    return lap_mat, N
