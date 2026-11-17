#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Core definitions for JAMfit - plasma current reconstruction with filaments and a ThinCurr wall model

JAMfit reconstructs the toroidal plasma current as a set of axisymmetric current filaments, together with the
induced currents in the vacuum vessel, from magnetic sensor signals. The wall is described by a ThinCurr
reduced model (a small set of current basis vectors, usually the dominant L/R eigenmodes). The fitted filament
and wall currents are then passed to TokaMaker to trace the LCFS and compute q, q95, li and centroids.

Examples (in `examples/ThinCurr/filament_reconstruct/`): `one_filament_ex.ipynb` runs both steps below on a
small torus; `vde_diiid_create_ex.ipynb` and `vde_diiid_reconstruct_ex.ipynb` run them for a synthetic DIII-D VDE.

Typical workflow:
  1. Create:
     @ref JAMfit.__init__ "JAMfit()" -> @ref JAMfit.set_sensors "set_sensors()" -> @ref JAMfit.setup_JAMfit "setup_JAMfit()"
     -> @ref JAMfit.gen_synthetic_data "gen_synthetic_data()" (optional synthetic signals)
     -> @ref JAMfit.create_from_runTD_top_modes "create_from_runTD_top_modes()" (reduced model)
  2. Reconstruct:
     @ref JAMfit.initialize_reduced_model "initialize_reduced_model()" -> one of the solvers
     (@ref JAMfit.run_reconstruction_lstsq "run_reconstruction_lstsq()",
     @ref JAMfit.run_reconstruction_tsvd_laplace "run_reconstruction_tsvd_laplace()",
     @ref JAMfit.run_reconstruction_laplace "run_reconstruction_laplace()")
     -> @ref JAMfit.post_process_tidx "post_process_tidx()"

Conventions:
  - SI units throughout: currents `[A]`, lengths `[m]`, times `[s]`, flux `[Wb]`, TokaMaker psi `[Wb/rad]`.
  - Icoil ordering in the ThinCurr XML: all non-filament coils (shaping/PF coils) first, then one Icoil
    per plasma filament. `num_non_fil_coils` is the number of coils before the first filament.
  - Sensor ordering: the real magnetic sensors first, then (optionally) in-vessel flux probes used only to
    evaluate the wall flux for TokaMaker. `num_sensors` / `num_real_sensors` is the number of real sensors.
  - Filament grids are given as flat arrays `rgrid`, `zgrid` of shape `(n_fil,)`, in the same order as the
    filament Icoils.

@authors Jamie Xia
@date Nov 2026
'''

## IMPORTING EXTERNAL LIBRARIES ##
import math
import os

import matplotlib.pyplot as plt
import numpy
from matplotlib.path import Path

## IMPORTING OFT MODULES ##
from .._core import OFT_env
from ..io import histfile
from ..TokaMaker._core import TokaMaker
from ..TokaMaker.meshing import load_gs_mesh
from ._core import ThinCurr, ThinCurr_reduced
from .sensor import save_sensors

# ==========================================================
# JAMfit Base Class
# ==========================================================

class JAMfit():
    '''! Main class for JAMfit filament reconstruction.
    Manages setup, synthetic data generation, reduced model creation,
    and reconstruction of plasma filament currents using both
    ThinCurr and TokaMaker.

    Matrix names used throughout (all from the reduced model, `self.torus_reduced`):
      - `Ms`  `(n_wall, n_sensors_total)`: signal on each sensor per unit amplitude of each wall basis vector.
      - `Msc` `(n_coils_total, n_sensors_total)`: signal on each sensor per amp in each Icoil (`[Wb/A]` for
        flux loops); rows `[:num_non_fil_coils]` are the non-filament coils, rows `[num_non_fil_coils:]` the
        filaments.
      - `sigma` `(n_sensors,)`: expected noise on each sensor, in the same units as the sensor signal.
        "Weighted" quantities are divided by `sigma` column-wise, so every sensor enters the fit in units of
        its own noise.

    Main attributes (and the method that sets them):
      - `myOFT` (OFT_env), `mygs` (TokaMaker): @ref __init__
      - `torus` (ThinCurr), `Msensor`, `Msc`, `sensor_obj`, `Mcoil`, L and R matrices: @ref setup_JAMfit
      - `eig_vals`, `eig_vecs`, `freq_vecs`, `freq_labels`: @ref create_from_runTD_top_modes
      - `torus_reduced` (ThinCurr_reduced), `reduced_created_flag`: @ref create_reduced_model,
        @ref create_from_runTD_top_modes or @ref initialize_reduced_model

    Example
    -------
    >>> jam = JAMfit('model.xml', 'thincurr_mesh.h5', 'tokamaker_mesh.h5', B0=2.0, R0=1.7, nthreads=4)
    >>> jam.initialize_reduced_model('reduced_model.h5')
    >>> prep = jam.prepare_lstsq(sigma, n_coils, reg_factor_fil=1e-6, reg_factor_wall=1e-5)
    >>> fil, wall, Ax, diag = jam.solve_lstsq(prep, psi, ip, coil_currs, ip_weight)
    '''

    def __init__(self, xml_file, thincurr_meshfile, tokamaker_meshfile, B0, R0, nthreads=None, oft_env=None):
        '''! Initialize the JAMfit object and its TokaMaker instance.

        TokaMaker is set up here (order 2, `F0 = R0*B0` `[T m]`, the vacuum toroidal field function) and is
        reused by @ref post_process_tidx. The ThinCurr model is only built later, in @ref setup_JAMfit.
        Only one OFT_env may exist per Python process, so pass `oft_env` if one already exists.

        @param xml_file str, path to the ThinCurr XML configuration file (Icoils: non-filament coils first, then filaments)
        @param thincurr_meshfile str, path to the ThinCurr mesh file
        @param tokamaker_meshfile str, path to the TokaMaker mesh file
        @param B0 float, vacuum toroidal field at R0 `[T]`
        @param R0 float, reference major radius `[m]`
        @param nthreads int, number of threads for a new OFT_env (required if oft_env is None) (default: None)
        @param oft_env OFT_env, an existing OFT environment instance to reuse (default: None)
        '''
        if oft_env is not None:
            self.myOFT = oft_env
        else:
            if nthreads is None:
                raise ValueError("Either oft_env or nthreads must be provided")
            self.myOFT = OFT_env(nthreads=nthreads)

        self.xml_file = xml_file
        self.thincurr_meshfile = thincurr_meshfile
        self.tokamaker_meshfile = tokamaker_meshfile
        self.reduced_created_flag = False
        self._lap_cache = None  # last filament Laplacian built, see _get_filament_laplacian

        # --- TokaMaker setup ---
        self.mygs = TokaMaker(self.myOFT)
        mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(self.tokamaker_meshfile)
        self.mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
        self.mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
        self.mygs.setup(order=2, F0=R0*B0)  # F0 = self.F0

    # ==========================================================
    # Model setup and reduced-model creation methods
    # ==========================================================
    def set_xml(self, xml_file):
        '''! Update the XML configuration file path.
        Only the stored path changes, so call this before @ref setup_JAMfit for it to take effect.

        @param xml_file str, path to the new XML configuration file
        '''
        self.xml_file = xml_file

    def set_sensors(self, sensor_array, floops_path='floops.loc'):
        '''! Save a sensor array to file and return the file path.

        @param sensor_array list, ThinCurr sensor objects (e.g. Mirnov, flux loops); real sensors first, then any flux probes
        @param floops_path str, output file path for sensor locations (default: 'floops.loc')
        @result str, path to the saved sensor file, to pass to @ref setup_JAMfit
        '''
        save_sensors(sensor_array, filename=floops_path)
        return floops_path

    def setup_JAMfit(self, floops_path, plot_files=None, use_legacy_io=False, mcoil_cache='HODLR_Mcoil.save', lmat_cache='HODLR_L.save'):
        '''! Set up the ThinCurr model, I/O, and sensor/coil mutual inductance matrices.
        Sets `self.torus` (ThinCurr), `self.Msensor`, `self.Msc` and `self.sensor_obj` (from compute_Msensor),
        `self.Mcoil` `(n_icoils, nelems)` `[Wb/A]`, and computes the wall inductance (HODLR) and resistance
        matrices on `self.torus`.

        @param floops_path str, path to the sensor locations file (e.g. from @ref set_sensors)
        @param plot_files str, base path for ThinCurr output/plot files (default: None -> ThinCurr default)
        @param use_legacy_io bool, if True, uses legacy HDF5 I/O (default: False)
        @param mcoil_cache str, cache file to save/load the coil-to-wall mutual inductance `Mcoil` (default: 'HODLR_Mcoil.save')
        @param lmat_cache str, cache file to save/load the HODLR wall inductance matrix `L` (default: 'HODLR_L.save')
        '''
        self.torus = ThinCurr(self.myOFT)
        self.torus.setup_model(mesh_file=self.thincurr_meshfile, xml_filename=self.xml_file)
        if plot_files is not None:
            self.torus.setup_io(basepath=plot_files, legacy_hdf5=use_legacy_io)
        else:
            self.torus.setup_io()
        self.Msensor, self.Msc, self.sensor_obj = self.torus.compute_Msensor(floops_path)
        self.Mcoil = self.torus.compute_Mcoil(cache_file=mcoil_cache)  # (n_icoils, nelems) `[Wb/A]`, used by compute_freq_response_basis
        self.torus.compute_Lmat(use_hodlr=True, cache_file=lmat_cache)
        self.torus.compute_Rmat(copy_out=True)
        print('JAMfit setup complete.')

    def setup_fil_timeseries(self, time_array, totalip, coil_currs, r_list, z_list, sigma_r, sigma_z, rgrid, zgrid):
        '''! Build a combined coil + plasma current array for a time-dependent run.
        Generates a Gaussian plasma current distribution at each time step and
        appends it to the coil currents for use with run_td.

        @param time_array numpy.ndarray, time points `[s]`, shape (ntimes, 1) (a column, so it can be hstacked)
        @param totalip numpy.ndarray, total plasma current at each time point `[A]`, shape (ntimes,)
        @param coil_currs numpy.ndarray, non-filament coil currents `[A]`, shape (ntimes, num_non_fil_coils)
        @param r_list list, R position of the plasma current centre at each time point `[m]`
        @param z_list list, Z position of the plasma current centre at each time point `[m]`
        @param sigma_r float, Gaussian width in the R direction `[m]`
        @param sigma_z float, Gaussian width in the Z direction `[m]`
        @param rgrid numpy.ndarray, R coordinates of the filament grid `[m]`, shape (n_fil,)
        @param zgrid numpy.ndarray, Z coordinates of the filament grid `[m]`, shape (n_fil,)
        @result numpy.ndarray, shape (ntimes, 1 + num_non_fil_coils + n_fil), columns
            `[time | non-filament coils | filaments]`, the layout run_td expects
        '''
        coil_currs_with_time = numpy.hstack([time_array, coil_currs])                  # [time | coils]
        fil_currs_with_time = setup_synthetic_current(time_array, totalip, sigma_r, sigma_z, r_list, z_list, rgrid, zgrid)
        fil_currs = fil_currs_with_time[:, 1:]                                          # drop the time column
        all_currs = numpy.hstack((coil_currs_with_time, fil_currs))                     # [time | coils | filaments]
        return all_currs

    def gen_synthetic_data(self, coil_currs, dt, nsteps, verbose=False, s_freq=10, p_freq=10):
        '''! Run a synthetic time-dependent simulation with the full ThinCurr model and record sensor signals.
        Uses the L and R matrices from @ref setup_JAMfit. The sensor history is written to
        `<io_basepath>/floops.hist` and plot files are written every `p_freq` steps.

        @param coil_currs numpy.ndarray, current waveforms `[A]` with a leading time column `[s]`, shape
            (ntimes, 1 + num_non_fil_coils + n_fil), e.g. from @ref setup_fil_timeseries
        @param dt float, time step size `[s]`
        @param nsteps int, number of time steps
        @param verbose bool, if True, also builds the XDMF plot files (for viewing in VisIt/ParaView) (default: False)
        @param s_freq int, print run status every s_freq steps (default: 10)
        @param p_freq int, save plot output every p_freq steps (default: 10)
        @result histfile, sensor signal history read from `<io_basepath>/floops.hist`
        '''

        self.torus.run_td(dt, nsteps, coil_currs=coil_currs, sensor_obj=self.sensor_obj, status_freq= s_freq, plot_freq=p_freq)
        self.torus.plot_td(nsteps, plot_freq = p_freq, sensor_obj=self.sensor_obj)
        if verbose:
            self.torus.build_XDMF()
        hist_file = histfile(os.path.join(self.torus._io_basepath, 'floops.hist'))
        return hist_file

    def create_from_runTD_top_modes(self, num_modes, reduced_filename, coil_currs, dt, nsteps, initial_num_eigs=50, plot_modes=False, verbose=False,
                                    s_freq=10, p_freq=10, axisymmetric_modes=None, freq_response=None):
        '''! Build a reduced model using the dominant eigenmodes from a full run.
        Computes eigenvalues, runs a preliminary reduced model to identify the
        most active modes by current amplitude, then builds a final reduced model
        using only those dominant modes.

        The weight of mode i is `sum_t |I_i(t)|`, the summed absolute amplitude of that eigenmode over the
        preliminary run with all `initial_num_eigs` modes. The preliminary model is written to
        `first_reduced_model_temp.h5` in the working directory. The selected modes are kept in ascending
        mode-index order.

        If axisymmetric_modes is provided, only modes whose indices are in
        axisymmetric_modes are eligible to be selected as "top modes". If a
        high-weight mode is not in axisymmetric_modes, it is skipped and the next
        eligible mode is used instead.

        @param num_modes int, number of dominant modes to retain in the reduced model
        @param reduced_filename str, output filename for the reduced model (HDF5)
        @param coil_currs numpy.ndarray, current waveforms `[A]` with a leading time column `[s]`, shape
            (ntimes, 1 + n_icoils), e.g. from @ref setup_fil_timeseries
        @param dt float, time step size `[s]`
        @param nsteps int, number of time steps for the preliminary run
        @param initial_num_eigs int, number of eigenmodes to compute initially (default: 50)
        @param plot_modes bool, if True, plots |amplitude| of every mode over the preliminary run, selected modes in
               colour and the rest in grey, on a new figure left open for the caller to label/show (default: False)
        @param verbose bool, if True, prints the weight of each selected mode and passes verbose on to
               @ref compute_freq_response_basis (default: False)
        @param s_freq int, print run status every s_freq steps (default: 10)
        @param p_freq int, sample the preliminary run's mode currents every p_freq steps; the mode weights are summed
               over these samples (default: 10)
        @param axisymmetric_modes list or array-like, indices of the axisymmetric modes to allow
                                (default: None -> allow all computed modes)
        @param freq_response dict, keyword arguments for @ref compute_freq_response_basis (at least
               `fil_indices` and `num_non_fil_coils`). If given, the frequency-response vectors are
               appended after the selected eigenmodes and always kept. The vectors and their labels
               are stored in `self.freq_vecs` and `self.freq_labels` (default: None -> eigenmodes only)

        @result ThinCurr_reduced, the constructed reduced model object
        @result dict, details of the mode selection, always returned:
            - `'sensors'`: dict, sensor signals from the preliminary run (`'time'` `[s]`, `'sensors'` shape (nsteps+1, n_sensors))
            - `'currents'`: dict, mode currents from the preliminary run (`'time'` `[s]`, `'curr'` shape (ntimes, initial_num_eigs))
            - `'basis'`: numpy.ndarray, basis of the reduced model, shape (n_basis, nelems): the selected eigenmodes,
              followed by any frequency-response vectors
            - `'eig_inds'`: list, indices of the selected eigenmodes, ascending
            - `'weights'`: list, weight `sum_t |I_i(t)|` of each selected eigenmode, in `eig_inds` order

        Example
        -------
        >>> reduced, info = jam.create_from_runTD_top_modes(4, 'reduced.h5', coil_currs, 2.E-4, 150, initial_num_eigs=35, plot_modes=True)
        >>> print(info['eig_inds'])
        '''
        # --- Eigenmodes of the full model ---
        self.eig_vals, self.eig_vecs = self.torus.get_eigs(initial_num_eigs, False)

        # --- Preliminary run with all initial_num_eigs modes ---
        torus_first_reduced = self.torus.build_reduced_model(
            self.eig_vecs,
            filename='first_reduced_model_temp.h5',
            sensor_obj=self.sensor_obj
        )

        sensors_measurement, currents = torus_first_reduced.run_td(
            dt, nsteps, coil_currs, status_freq=s_freq, plot_freq=p_freq
        )

        # --- Rank the modes by their weight sum_t |I_i(t)| ---
        temp_curr = currents['curr']
        temp_curr = temp_curr[:, 0:initial_num_eigs]  # Only consider the modes we computed
        max_weights = [abs(temp_curr[:, i]).sum() for i in range(temp_curr.shape[1])]  # total sum over time

        # axisymmetric filtering
        if axisymmetric_modes is None:
            eligible = list(range(temp_curr.shape[1]))
        else:
            axis_set = set(axisymmetric_modes)
            # Only consider indices that are within the computed range and allowed
            eligible = [i for i in range(temp_curr.shape[1]) if i in axis_set]

        if len(eligible) < num_modes:
            raise ValueError(
                f"Requested num_modes={num_modes}, but only {len(eligible)} eligible modes are available "
                f"within axisymmetric_modes (and < initial_num_eigs={initial_num_eigs})."
            )

        # Pick the top modes *among eligible ones*, then keep them in ascending mode-index order
        top_modenum_indices = sorted(eligible, key=lambda i: max_weights[i], reverse=True)[:num_modes]
        eig_inds = sorted(top_modenum_indices)
        weight_amplitude = [max_weights[i] for i in eig_inds]

        if verbose:
            for i in eig_inds:
                print(f'Saved mode {i} has weight {max_weights[i]}')

        # --- Plot all mode amplitudes, selected modes in colour and the rest in grey ---
        if plot_modes:
            fig, ax = plt.subplots(1, 1)
            for i in range(temp_curr.shape[1]):
                if i in eig_inds:
                    ax.semilogy(currents['time'], abs(currents['curr'][:, i]), label=f'Mode {i}')
                else:
                    ax.semilogy(currents['time'], abs(currents['curr'][:, i]), color='gray', alpha=0.3)

        # --- Build the basis: selected eigenmodes, then any frequency-response vectors (always kept, not ranked) ---
        basis = self.eig_vecs[eig_inds, :]
        self.freq_vecs = None
        self.freq_labels = None
        if freq_response is not None:
            self.freq_vecs, self.freq_labels = self.compute_freq_response_basis(
                verbose=verbose, **freq_response
            )
            basis = numpy.vstack([basis, self.freq_vecs])

        # --- Build the final reduced model ---
        self.torus_reduced = self.create_reduced_model(basis, reduced_filename, compute_B=False)
        self.reduced_created_flag = True

        info = {
            'sensors': sensors_measurement,
            'currents': currents,
            'basis': basis,
            'eig_inds': eig_inds,
            'weights': weight_amplitude,
        }
        return self.torus_reduced, info

    def create_reduced_model(self, eig_vecs, reduced_filename, compute_B=False):
        '''! Build a reduced model from a given set of wall-current basis vectors.
        The basis is usually a subset of the L/R eigenvectors, optionally with frequency-response vectors
        appended (see @ref compute_freq_response_basis). Sets `self.torus_reduced`.

        @param eig_vecs numpy.ndarray, basis vectors for the reduced model, shape (n_basis, nelems)
        @param reduced_filename str, output filename for the reduced model (HDF5)
        @param compute_B bool, if True, computes the B-field reconstruction matrix for the reduced model (default: False)
        @result ThinCurr_reduced, the constructed reduced model object
        '''
        self.torus_reduced = self.torus.build_reduced_model(
            eig_vecs, filename=reduced_filename, compute_B=compute_B, sensor_obj=self.sensor_obj
        )
        print(f"Reduced model created with {eig_vecs.shape[0]} modes")
        self.reduced_created_flag = True
        return self.torus_reduced

    def compute_freq_response_basis(self, fil_indices, num_non_fil_coils, freqs=None, verbose=False):
        '''! Compute wall-current basis vectors from the frequency response to chosen plasma filaments.
        Fast events (e.g. a current quench) drive a flux-conserving wall response that can be spread
        over many short-time-constant eigenmodes, so a truncated eigenbasis may miss it. This method
        computes that response directly and returns it as extra basis vectors for the reduced model.
        The vectors are forced-response patterns, not eigenmodes (they have no L/R time of their own).

        For a filament carrying 1 A at angular frequency w, the wall currents solve
        `(R + i w L) I = -i w M_fil`, i.e. `I = -i w (R + i w L)^-1 M_fil`, where `M_fil` is the
        filament row of `Mcoil` (flux on each wall element per amp, `[Wb/A]`). Limits:
          - freq = inf (inductive limit, `fr_limit=1`): `I = -L^-1 M_fil`, purely real. This is the
            wall current immediately after an instantaneous change of the filament current.
          - freq = 0 (resistive limit, `fr_limit=2`): `I = -i R^-1 M_fil` (per unit w), purely imaginary.
          - finite freq: both the in-phase (real) and quadrature (imaginary) parts are used.

        Each response is kept as computed (it is not orthogonalized against the eigenmodes), only scaled so
        `v^T R v = 1`, the same normalization as ThinCurr eigenvectors. The overlap with the eigenmodes is
        part of the physical response and is kept: removing it would change the current pattern the vector
        carries. The combined basis is therefore not R-orthonormal and the reduced `R` is not the identity;
        `ThinCurr_reduced` uses the full reduced `L` and `R`, so this is fine for time-dependent runs.

        Suggested use: one or two filaments (e.g. plasma centre and the end position of a VDE) at a single
        frequency on the disruption time scale, `f ~ 1/(2 pi tau_CQ)`. The inductive limit alone tends to put too
        much current on the wall surface closest to the plasma. Avoid pairs of drives with nearly identical
        responses (e.g. `numpy.inf` together with a frequency far above the wall rates, where `w tau >> 1`): they
        make the reduced inductance matrix nearly singular.

        Requires @ref setup_JAMfit (needs `self.torus` with L, R and Mcoil).

        @param fil_indices list of int, filaments to drive, counted from the first filament (0 = first filament Icoil)
        @param num_non_fil_coils int, number of non-filament Icoils listed before the filaments in the XML
        @param freqs float or list of float, drive frequencies `[Hz]`; `numpy.inf` = inductive limit,
               `0` = resistive limit (default: None -> `[numpy.inf]`)
        @param verbose bool, if True, prints each kept vector and its R-norm before scaling (default: False)
        @result numpy.ndarray, basis vectors `(nvecs, nelems)`, each scaled to `v^T R v = 1`
        @result list of tuple, `(filament index, freq [Hz], 'real' or 'imag')` source of each vector

        Example
        -------
        >>> # Disruption with the plasma staying near the centre: one filament, inductive limit
        >>> vecs, labels = jam.compute_freq_response_basis([center_fil], num_non_fil_coils=n_coils)
        >>> # VDE then quench with tau_CQ ~ 19 ms: centre and end-position filaments at f_CQ = 1/(2 pi tau_CQ) ~ 8.3 Hz
        >>> vecs, labels = jam.compute_freq_response_basis([center_fil, vde_end_fil], n_coils, freqs=[8.3])
        '''
        if freqs is None:
            freqs = [numpy.inf]
        freqs = numpy.atleast_1d(numpy.asarray(freqs, dtype=numpy.float64))
        if not hasattr(self, 'Mcoil'):
            self.Mcoil = self.torus.compute_Mcoil()
        Rmat = self.torus.Rmat  # scipy sparse (nelems, nelems)

        # --- Frequency response for each filament and frequency ---
        vecs = []
        labels = []
        for fil in fil_indices:
            fdriver = numpy.zeros((2, self.torus.nelems)) #one row for real, one row for img part 
            fdriver[0, :] = self.Mcoil[num_non_fil_coils + fil, :]  # real flux driver for 1 A in the filament, fdriver [0,:] is real, fdriver [1,:] is img part (set to 0)
            for freq in freqs:
                if numpy.isinf(freq):
                    result = self.torus.compute_freq_response(fdriver=fdriver, fr_limit=1)
                    candidates = [(result[0, :], 'real')]
                elif freq == 0.0:
                    result = self.torus.compute_freq_response(fdriver=fdriver, fr_limit=2)
                    candidates = [(result[1, :], 'imag')]
                else:
                    result = self.torus.compute_freq_response(fdriver=fdriver, freq=freq)
                    candidates = [(result[0, :], 'real'), (result[1, :], 'imag')]

                # --- scale each response to v^T R v = 1 (shape unchanged) ---
                # matching magnitudes of the frequency-response vectors to the R-normalized eigenvectors,
                # since we are appending to the reduced-model basis. 
                for v, part in candidates:
                    v_norm = numpy.sqrt(numpy.dot(v, Rmat @ v))
                    if v_norm == 0.0:
                        print(f"Skipping filament {fil} at {freq} Hz ({part}): zero response")
                        continue
                    if verbose:
                        print(f"Keeping filament {fil} at {freq} Hz ({part}): ||v||_R = {v_norm:.3e} before scaling")
                    vecs.append(v/v_norm)
                    labels.append((fil, float(freq), part))

        new_vecs = numpy.array(vecs).reshape((-1, self.torus.nelems))
        print(f"Frequency response: {new_vecs.shape[0]} vectors from {len(fil_indices)} filament(s) x {len(freqs)} frequency(ies)")
        return new_vecs, labels

    # ==========================================================
    # Reconstruction methods
    # ==========================================================

    def initialize_reduced_model(self, reduced_filename):
        '''! Load a previously saved reduced model from file.
        Use this instead of the `create_*` methods when the reduced model already exists (e.g. in a
        reconstruction notebook). Sets `self.torus_reduced`. Does not need @ref setup_JAMfit.

        @param reduced_filename str, path to the HDF5 reduced model file
        '''
        self.torus_reduced = ThinCurr_reduced(reduced_filename)
        self.reduced_created_flag = True
        print('Reduced model initialized from file: {0}'.format(reduced_filename))

    # ----------------------------------------------------------
    # How the three solvers are organized
    # ----------------------------------------------------------
    # Every solver fits the unknowns x = [I_wall | filament unknowns] to the rows
    #   - weighted sensors:   (A_sensor x - B_weighted)       one row per sensor
    #   - regularization:     (A_reg x - 0)                   Tikhonov and/or Laplacian rows
    #   - Ip constraint:      c * (r . x - Ip)                one row, r . x = sum of the filament currents
    # in the least-squares sense. Everything except the Ip row only depends on sigma, the grid and the
    # regularization weights, so prepare_<solver> builds it once: a QR factorization of the fixed rows (no A^T A is
    # formed, which would square the condition number) plus the vector g = (A0^T A0)^-1 r^T that adds the Ip row.
    # solve_<solver> then only builds B_weighted for the slice and adds the Ip row with weight c (which may change
    # every slice: it depends on Ip, and on the signals when lam is automatic) through the exact one-row update in
    # _prepare_fixed_rows / _solve_with_ip_row. run_reconstruction_<solver> is prepare + solve in one call, for
    # sigma or weights that change every slice.

    def _get_filament_laplacian(self, rgrid, zgrid, gaussian=False, sigma_r=1e-1, sigma_z=1e-1, verbose=False):
        '''! Filament Laplacian for the grid, cached on the object.
        Only depends on the grid (and on sigma_r, sigma_z for the Gaussian version), so it is rebuilt only when one
        of those changes. The cache holds the last grid used.

        @param rgrid numpy.ndarray, R coordinates of the filament grid `[m]`, shape (n_fil,)
        @param zgrid numpy.ndarray, Z coordinates of the filament grid `[m]`, shape (n_fil,)
        @param gaussian bool, if True uses @ref get_gaussian_lap_mat, else @ref get_laplace_matrix (default: False)
        @param sigma_r float, R width of the Gaussian weighting `[m]`, only used if gaussian is True (default: 1e-1)
        @param sigma_z float, Z width of the Gaussian weighting `[m]`, only used if gaussian is True (default: 1e-1)
        @param verbose bool, passed to @ref get_laplace_matrix (default: False)
        @result numpy.ndarray, Laplacian matrix, shape (n_fil, n_fil)
        '''
        rgrid = numpy.asarray(rgrid, dtype=float)
        zgrid = numpy.asarray(zgrid, dtype=float)
        widths = (float(sigma_r), float(sigma_z)) if gaussian else None
        cache = getattr(self, '_lap_cache', None)
        if cache is not None:
            same_grid = (cache['gaussian'] == gaussian and cache['widths'] == widths
                         and cache['rgrid'].shape == rgrid.shape
                         and numpy.array_equal(cache['rgrid'], rgrid) and numpy.array_equal(cache['zgrid'], zgrid))
            if same_grid:
                return cache['lap_mat']
        if gaussian:
            lap_mat, _ = get_gaussian_lap_mat(rgrid, zgrid, sigma_r, sigma_z)
        else:
            lap_mat, _ = get_laplace_matrix(rgrid, zgrid, verbose=verbose)
        self._lap_cache = {'rgrid': rgrid.copy(), 'zgrid': zgrid.copy(), 'gaussian': gaussian, 'widths': widths,
                           'lap_mat': lap_mat}
        return lap_mat

    def _weighted_sensor_matrices(self, sigma, num_non_fil_coils, num_sensors):
        '''! Sensor matrices of the reduced model divided by the sensor noise.
        Dividing by sigma puts every sensor in units of its own noise, so no single sensor dominates the fit.

        @param sigma numpy.ndarray, noise on each sensor, same units as the sensor signals, shape (num_sensors,)
        @param num_non_fil_coils int, number of non-filament coils before the filaments in the Icoil list
        @param num_sensors int, number of (real) sensors to use (the first num_sensors sensors)
        @result numpy.ndarray, `Ms / sigma`, shape (n_wall, num_sensors)
        @result numpy.ndarray, `Msc_coils / sigma` (non-filament coils), shape (num_non_fil_coils, num_sensors)
        @result numpy.ndarray, `Msc_fil / sigma` (filaments), shape (n_fil, num_sensors)
        '''
        if not self.reduced_created_flag:
            raise ValueError("Reduced model has not been created yet. Please create or initialize a reduced model before running reconstruction.")
        sigma = numpy.asarray(sigma, dtype=float)
        if sigma.shape != (num_sensors,):
            raise ValueError("sigma has shape {0}, expected ({1},) (one value per sensor)".format(sigma.shape, num_sensors))
        Ms_weighted = self.torus_reduced.Ms[:, :num_sensors] / sigma                      # (n_wall, n_sensors)
        Msc_weighted = self.torus_reduced.Msc[:, :num_sensors] / sigma                    # (n_coils_total, n_sensors)
        return Ms_weighted, Msc_weighted[:num_non_fil_coils, :], Msc_weighted[num_non_fil_coils:, :]

    def prepare_lstsq(self, sigma, num_non_fil_coils, reg_factor_fil, reg_factor_wall, num_sensors=None):
        '''! Build everything for @ref solve_lstsq that does not change between time slices.
        Call once for a fixed `sigma` and fixed regularization; call again if any of them changes (or use
        @ref run_reconstruction_lstsq, which does both steps in one call).

        Unknowns `x = [I_wall | I_fil]` (wall basis amplitudes, then filament currents). Fixed rows:
          - sensors: `(Ms^T I_wall + Msc_fil^T I_fil) / sigma = Psi / sigma - (Msc_coils^T I_coils) / sigma`
          - Tikhonov: `reg_factor_wall * I_wall = 0` and `reg_factor_fil * I_fil = 0`
        The Ip row is added per slice in @ref solve_lstsq.

        @param sigma numpy.ndarray, noise on each sensor, same units as the sensor signals, shape (num_sensors,)
        @param num_non_fil_coils int, number of non-filament coils before the filaments in the Icoil list
        @param reg_factor_fil float, Tikhonov regularization factor for the filament currents
        @param reg_factor_wall float, Tikhonov regularization factor for the wall basis amplitudes
        @param num_sensors int or None, number of (real) sensors to use; None uses all sensors in the reduced model (default: None)
        @result dict, prepared matrices to pass to @ref solve_lstsq

        Example
        -------
        >>> prep = jam.prepare_lstsq(sigma, num_non_fil_coils, reg_factor_fil=4e-4, reg_factor_wall=1e-6, num_sensors=117)
        >>> for t in range(len(times)):
        ...     fil, wall, Ax, diag = jam.solve_lstsq(prep, psi[t], ip[t], coils[t], ip_weight=4e2)
        '''
        if num_sensors is None:
            num_sensors = self.torus_reduced.Ms.shape[1]
        Ms_weighted, Msc_coils_weighted, Msc_fil_weighted = self._weighted_sensor_matrices(sigma, num_non_fil_coils, num_sensors)
        num_Ms = Ms_weighted.shape[0]
        n_fil = Msc_fil_weighted.shape[0]

        # --- Fixed rows: sensors, then Tikhonov on the wall and on the filaments ---
        A_sensor = numpy.hstack([Ms_weighted.T, Msc_fil_weighted.T])                       # (n_sensors, n_wall + n_fil)
        A_reg = numpy.diag(numpy.concatenate([numpy.full(num_Ms, float(reg_factor_wall)),
                                              numpy.full(n_fil, float(reg_factor_fil))]))  # (n_wall + n_fil, n_wall + n_fil)

        # --- Ip row: r . x = sum(I_fil) ---
        ip_row = numpy.concatenate([numpy.zeros(num_Ms), numpy.ones(n_fil)])

        prep = _prepare_fixed_rows(A_sensor, A_reg, ip_row)
        prep.update({
            'method': 'lstsq',
            'sigma': numpy.array(sigma, dtype=float),
            'num_sensors': num_sensors,
            'num_non_fil_coils': num_non_fil_coils,
            'num_Ms': num_Ms,
            'Msc_coils_weighted': Msc_coils_weighted,
            'reg_factor_fil': reg_factor_fil,
            'reg_factor_wall': reg_factor_wall,
        })
        return prep

    def solve_lstsq(self, prep, psi_at_time, ip_at_time, coil_curr_at_time, ip_weight, ip_norm_floor=1.E3):
        '''! Solve one time slice with the lstsq method, using matrices from @ref prepare_lstsq.

        The Ip row is `(ip_weight / ip_norm) * (sum(I_fil) - ip_at_time)`, with `ip_norm = max(|ip_at_time|, ip_norm_floor)`.
        For normal currents this is the relative Ip error; the floor stops a divide-by-zero when Ip passes through
        zero (e.g. at the end of a current quench), where the row becomes the absolute Ip error in units of
        `ip_norm_floor`.

        @param prep dict, from @ref prepare_lstsq
        @param psi_at_time numpy.ndarray, sensor signals at this time (e.g. flux `[Wb]`), shape (num_sensors,)
        @param ip_at_time float, total plasma current measurement at time `[A]`
        @param coil_curr_at_time numpy.ndarray, currents in the non-filament coils at time `[A]`, shape (num_non_fil_coils,)
        @param ip_weight float, weight of the Ip constraint row relative to the (sigma-weighted) sensor rows
        @param ip_norm_floor float, smallest current used to normalize the Ip constraint row and `ip_error_pct` `[A]` (default: 1.E3)
        @result numpy.ndarray, filament currents `[A]`, shape (n_fil,)
        @result numpy.ndarray, wall basis amplitudes, shape (n_wall,)
        @result numpy.ndarray, `A x` for the sensor rows only (fitted weighted signals), shape (num_sensors,)
        @result dict, diagnostics, see @ref _reconstruction_diagnostics; also `ip_weight`, `reg_factor_fil`, `reg_factor_wall`
        '''
        _check_prep(prep, 'lstsq', 'prepare_lstsq')
        num_Ms = prep['num_Ms']
        B_weighted = psi_at_time / prep['sigma'] - coil_curr_at_time @ prep['Msc_coils_weighted']

        ip_norm = max(abs(ip_at_time), ip_norm_floor)
        solution, Ax, total_residual = _solve_with_ip_row(prep, B_weighted, ip_weight / ip_norm, ip_at_time)
        curr_wall = solution[:num_Ms]
        curr_fil = solution[num_Ms:]

        diagnostics = _reconstruction_diagnostics(curr_fil, ip_at_time, ip_norm_floor, Ax, B_weighted, total_residual, solution)
        diagnostics.update({'ip_weight': ip_weight, 'reg_factor_fil': prep['reg_factor_fil'],
                            'reg_factor_wall': prep['reg_factor_wall']})
        return curr_fil, curr_wall, Ax, diagnostics

    def run_reconstruction_lstsq(self, psi_at_time, ip_at_time, num_non_fil_coils, coil_curr_at_time, ip_weight, sigma, reg_factor_fil, reg_factor_wall, num_sensors=None,
                                 ip_norm_floor=1.E3):
        '''! Run the filament current reconstruction with the lstsq method for one time slice.
        Same as @ref prepare_lstsq followed by @ref solve_lstsq. Use this when `sigma` or the regularization change
        every slice; for fixed ones, call @ref prepare_lstsq once and @ref solve_lstsq in the time loop (faster).
        See those two methods for the rows of the fit.

        @param psi_at_time numpy.ndarray, sensor signals at this time (e.g. flux `[Wb]`), shape (num_sensors,)
        @param ip_at_time float, total plasma current measurement at time `[A]`
        @param num_non_fil_coils int, number of non-filament coils before the filaments in the Icoil list
        @param coil_curr_at_time numpy.ndarray, currents in the non-filament coils at time `[A]`, shape (num_non_fil_coils,)
        @param ip_weight float, weight of the Ip constraint row relative to the (sigma-weighted) sensor rows
        @param sigma numpy.ndarray, noise on each sensor, same units as psi_at_time, shape (num_sensors,)
        @param reg_factor_fil float, Tikhonov regularization factor for the filament currents
        @param reg_factor_wall float, Tikhonov regularization factor for the wall basis amplitudes
        @param num_sensors int or None, number of (real) sensors to use; None uses all sensors in the reduced model (default: None)
        @param ip_norm_floor float, smallest current used to normalize the Ip constraint row `[A]` (default: 1.E3)
        @result numpy.ndarray, filament currents `[A]`, shape (n_fil,)
        @result numpy.ndarray, wall basis amplitudes, shape (n_wall,)
        @result numpy.ndarray, `A x` for the sensor rows only (fitted weighted signals), shape (num_sensors,)
        @result dict, diagnostics, see @ref solve_lstsq
        '''
        prep = self.prepare_lstsq(sigma, num_non_fil_coils, reg_factor_fil, reg_factor_wall, num_sensors)
        return self.solve_lstsq(prep, psi_at_time, ip_at_time, coil_curr_at_time, ip_weight, ip_norm_floor)

    def prepare_tsvd_laplace(self, sigma, num_non_fil_coils, num_sensors, rgrid, zgrid, num_modes, lap_lam=1e-8, reg_wall=1e-5, verbose=False):
        '''! Build everything for @ref solve_tsvd_laplace that does not change between time slices.
        Call once for a fixed `sigma` and fixed regularization; call again if any of them changes (or use
        @ref run_reconstruction_tsvd_laplace, which does both steps in one call).

        The SVD is taken of the weighted filament-to-sensor matrix `Msc_fil / sigma` `(n_fil, n_sensors)`. The
        filament currents are then written as `I_fil = U_trun @ w`, with `U_trun = U[:, :num_modes]` and `w` the
        `num_modes` unknown TSVD weights ("TSVD space"). This keeps only the filament patterns the sensors can see.
        Unknowns are `[I_wall | w]`. Fixed rows: weighted sensors, Tikhonov on the wall (`reg_wall`) and the
        filament Laplacian acting on `U_trun @ w` (`lap_lam`). The Ip row is added per slice in @ref solve_tsvd_laplace.

        @param sigma numpy.ndarray, noise on each sensor, same units as the sensor signals, shape (num_sensors,)
        @param num_non_fil_coils int, number of non-filament coils before the filaments in the Icoil list
        @param num_sensors int, number of real sensors in the system
        @param rgrid numpy.ndarray, R coordinates of the filament grid `[m]`, shape (n_fil,)
        @param zgrid numpy.ndarray, Z coordinates of the filament grid `[m]`, shape (n_fil,)
        @param num_modes int, number of SVD modes to truncate to
        @param lap_lam float, Laplacian regularization weight (default: 1e-8)
        @param reg_wall float, Tikhonov regularization weight on the wall amplitudes (default: 1e-5)
        @param verbose bool, if True, plots singular values (default: False)
        @result dict, prepared matrices to pass to @ref solve_tsvd_laplace (includes `U_trun`, shape (n_fil, num_modes))

        Example
        -------
        >>> prep = jam.prepare_tsvd_laplace(sigma, num_non_fil_coils, 117, rgrid, zgrid, num_modes=38, lap_lam=1e-6, reg_wall=1e-3)
        >>> for t in range(len(times)):
        ...     fil, wall, Ax, diag = jam.solve_tsvd_laplace(prep, psi[t], ip[t], coils[t], lam=1e-3)
        '''
        Ms_weighted, Msc_coils_weighted, Msc_fil_weighted = self._weighted_sensor_matrices(sigma, num_non_fil_coils, num_sensors)
        num_Ms = Ms_weighted.shape[0]
        lap_mat = self._get_filament_laplacian(rgrid, zgrid, verbose=verbose)

        # --- Truncated SVD of the filament part ---
        # only the filaments are decomposed (the shaping coil currents are known, not solved for)
        U, S, Vh = numpy.linalg.svd(Msc_fil_weighted, full_matrices=False)
        U_trun = U[:, :num_modes]                                                          # (n_fil, num_modes)

        # --- Fixed rows: sensors [wall | TSVD weights], wall Tikhonov, filament Laplacian in TSVD space ---
        A_sensor = numpy.hstack([Ms_weighted.T, Msc_fil_weighted.T @ U_trun])              # (n_sensors, n_wall + num_modes)
        A_reg = numpy.vstack([
            numpy.hstack([reg_wall * numpy.eye(num_Ms), numpy.zeros((num_Ms, num_modes))]),  # wall Tikhonov
            numpy.hstack([numpy.zeros((lap_mat.shape[0], num_Ms)), lap_lam * (lap_mat @ U_trun)]),  # fil Laplacian
        ])

        # --- Ip row: sum(I_fil) = sum(U_trun @ w) = (column sums of U_trun) . w ---
        ip_row = numpy.concatenate([numpy.zeros(num_Ms), U_trun.sum(axis=0)])

        if verbose:
            plt.figure(figsize=(8, 5))
            plt.semilogy(S/S[0], marker='o')
            plt.title('Normalized Singular Values of plasma filament contribution to sensor signals')
            plt.xlabel('Mode Index')
            plt.ylabel('Singular Value (log scale)')
            plt.grid(True)
            plt.show()

        prep = _prepare_fixed_rows(A_sensor, A_reg, ip_row)
        prep.update({
            'method': 'tsvd_laplace',
            'sigma': numpy.array(sigma, dtype=float),
            'num_sensors': num_sensors,
            'num_non_fil_coils': num_non_fil_coils,
            'num_Ms': num_Ms,
            'Msc_coils_weighted': Msc_coils_weighted,
            'U_trun': U_trun,
            # norm of the Ip row; lam is defined for the unit-norm row, so the row weight is lam / ip_row_norm
            'ip_row_norm': numpy.linalg.norm(ip_row),
            'lap_lam': lap_lam,
            'reg_wall': reg_wall,
        })
        return prep

    def solve_tsvd_laplace(self, prep, psi_at_time, ip_at_time, coil_curr_at_time, lam=None, ip_norm_floor=1.E3):
        '''! Solve one time slice with the svd + laplacian method, using matrices from @ref prepare_tsvd_laplace.
        The filament currents are `U_trun @ w`.

        The Ip row is `lam * (r . x - Ip) / |r|`, with `r . x = sum(I_fil)`, i.e. weight `lam` on the unit-norm row.
        If `lam` is None it is set from the ratio of the typical weighted sensor signal to Ip,
        `lam = 100 * 10**floor(log10(mean(|B|) / ip_norm))` with `ip_norm = max(|Ip|, ip_norm_floor)`, so the
        Ip row is on the same order as the sensor rows (the factor 100 gives the single Ip row more weight
        against the many sensor rows). The floor keeps `lam` finite when Ip passes through zero (e.g. at the
        end of a current quench). `ip_error_pct` is normalized by the same `ip_norm`, so below the floor it is
        the absolute Ip error in percent of `ip_norm_floor`.

        @param prep dict, from @ref prepare_tsvd_laplace
        @param psi_at_time numpy.ndarray, sensor signals at this time (e.g. flux `[Wb]`), shape (num_sensors,)
        @param ip_at_time float, total plasma current measurement at time `[A]`
        @param coil_curr_at_time numpy.ndarray, non-filament coil currents at time `[A]`, shape (num_non_fil_coils,)
        @param lam float, Ip constraint weight; None sets it automatically (see above) (default: None)
        @param ip_norm_floor float, smallest current used for the automatic `lam` and for `ip_error_pct` `[A]` (default: 1.E3)
        @result numpy.ndarray, filament currents `[A]`, shape (n_fil,)
        @result numpy.ndarray, wall basis amplitudes, shape (n_wall,)
        @result numpy.ndarray, `A x` for the sensor rows only (fitted weighted signals), shape (num_sensors,)
        @result dict, diagnostics, see @ref _reconstruction_diagnostics (`solution` is `[I_wall | w]`); also `lam`, `lap_lam`, `reg_wall`
        '''
        _check_prep(prep, 'tsvd_laplace', 'prepare_tsvd_laplace')
        num_Ms = prep['num_Ms']
        B_weighted = psi_at_time / prep['sigma'] - coil_curr_at_time @ prep['Msc_coils_weighted']
        if lam is None:
            lam = _auto_ip_weight(B_weighted, ip_at_time, ip_norm_floor)

        solution, Ax, total_residual = _solve_with_ip_row(prep, B_weighted, lam / prep['ip_row_norm'], ip_at_time)
        curr_wall = solution[:num_Ms]
        curr_fil = prep['U_trun'] @ solution[num_Ms:]   # map the TSVD weights back to filament currents

        diagnostics = _reconstruction_diagnostics(curr_fil, ip_at_time, ip_norm_floor, Ax, B_weighted, total_residual, solution)
        diagnostics.update({'lam': lam, 'lap_lam': prep['lap_lam'], 'reg_wall': prep['reg_wall']})
        return curr_fil, curr_wall, Ax, diagnostics

    def run_reconstruction_tsvd_laplace(self, psi_at_time, ip_at_time, coil_curr_at_time, sigma, num_non_fil_coils, num_sensors, rgrid, zgrid, num_modes,
                                        lam=None, lap_lam=1e-8, reg_wall=1e-5, ip_norm_floor=1.E3):
        '''! Run the filament current reconstruction with the svd + laplacian method for one time slice.
        Same as @ref prepare_tsvd_laplace followed by @ref solve_tsvd_laplace, so the SVD is redone on every call.
        Use this when `sigma` or the regularization change every slice; for fixed ones, call
        @ref prepare_tsvd_laplace once and @ref solve_tsvd_laplace in the time loop (much faster).

        @param psi_at_time numpy.ndarray, sensor signals at this time (e.g. flux `[Wb]`), shape (num_sensors,)
        @param ip_at_time float, total plasma current measurement at time `[A]`
        @param coil_curr_at_time numpy.ndarray, non-filament coil currents at time `[A]`, shape (num_non_fil_coils,)
        @param sigma numpy.ndarray, noise on each sensor, same units as psi_at_time, shape (num_sensors,)
        @param num_non_fil_coils int, number of non-filament coils before the filaments in the Icoil list
        @param num_sensors int, number of real sensors in the system
        @param rgrid numpy.ndarray, R coordinates of the filament grid `[m]`, shape (n_fil,)
        @param zgrid numpy.ndarray, Z coordinates of the filament grid `[m]`, shape (n_fil,)
        @param num_modes int, number of SVD modes to truncate to
        @param lam float, Ip constraint weight; None sets it automatically (see @ref solve_tsvd_laplace) (default: None)
        @param lap_lam float, Laplacian regularization weight (default: 1e-8)
        @param reg_wall float, Tikhonov regularization weight on the wall amplitudes (default: 1e-5)
        @param ip_norm_floor float, smallest current used for the automatic `lam` and for `ip_error_pct` `[A]` (default: 1.E3)
        @result numpy.ndarray, filament currents `[A]`, shape (n_fil,)
        @result numpy.ndarray, wall basis amplitudes, shape (n_wall,)
        @result numpy.ndarray, `A x` for the sensor rows only (fitted weighted signals), shape (num_sensors,)
        @result dict, diagnostics, see @ref solve_tsvd_laplace
        '''
        prep = self.prepare_tsvd_laplace(sigma, num_non_fil_coils, num_sensors, rgrid, zgrid, num_modes, lap_lam, reg_wall)
        return self.solve_tsvd_laplace(prep, psi_at_time, ip_at_time, coil_curr_at_time, lam, ip_norm_floor)

    def prepare_laplace(self, sigma, num_non_fil_coils, num_real_sensors, rgrid, zgrid, lap_lam=1e-6, reg_wall=1e-5,
                        sigma_r=1e-1, sigma_z=1e-1, gaussian=False, verbose=False):
        '''! Build everything for @ref solve_laplace that does not change between time slices.
        Call once for a fixed `sigma` and fixed regularization; call again if any of them changes (or use
        @ref run_reconstruction_laplace, which does both steps in one call). The filament Laplacian itself is
        cached on the object and only rebuilt when the grid changes.

        Solves directly in the physical filament space (no SVD projection). Unknowns are `[I_wall | I_fil]`.
        Fixed rows: weighted sensors, Tikhonov on the wall (`reg_wall`) and the filament Laplacian (`lap_lam`).
        The Ip row is added per slice in @ref solve_laplace.

        @param sigma numpy.ndarray, noise on each sensor, same units as the sensor signals, shape (num_real_sensors,)
        @param num_non_fil_coils int, number of non-filament coils before the filaments in the Icoil list
        @param num_real_sensors int, number of real sensors (the first num_real_sensors sensors are used)
        @param rgrid numpy.ndarray, R coordinates of the filament grid `[m]`, shape (n_fil,)
        @param zgrid numpy.ndarray, Z coordinates of the filament grid `[m]`, shape (n_fil,)
        @param lap_lam float, Laplacian regularization weight (default: 1e-6)
        @param reg_wall float, Tikhonov regularization weight on the wall amplitudes (default: 1e-5)
        @param sigma_r float, R width of the Gaussian weighting `[m]`, only used if gaussian is True (default: 1e-1)
        @param sigma_z float, Z width of the Gaussian weighting `[m]`, only used if gaussian is True (default: 1e-1)
        @param gaussian bool, if True uses the experimental Gaussian-weighted Laplacian (@ref get_gaussian_lap_mat)
            instead of the finite-difference one (@ref get_laplace_matrix) (default: False)
        @param verbose bool, if True prints grid information while building the Laplacian (default: False)
        @result dict, prepared matrices to pass to @ref solve_laplace

        Example
        -------
        >>> prep = jam.prepare_laplace(sigma, num_non_fil_coils, 117, rgrid, zgrid, lap_lam=1e-6, reg_wall=1e-5)
        >>> for t in range(len(times)):
        ...     fil, wall, Ax, diag = jam.solve_laplace(prep, psi[t], ip[t], coils[t], lam=None)
        '''
        Ms_weighted, Msc_coils_weighted, Msc_fil_weighted = self._weighted_sensor_matrices(sigma, num_non_fil_coils, num_real_sensors)
        num_Ms = Ms_weighted.shape[0]
        n_fil = Msc_fil_weighted.shape[0]
        lap_mat = self._get_filament_laplacian(rgrid, zgrid, gaussian, sigma_r, sigma_z, verbose)

        # --- Fixed rows: sensors [wall | filament], wall Tikhonov, filament Laplacian ---
        A_sensor = numpy.hstack([Ms_weighted.T, Msc_fil_weighted.T])                       # (n_sensors, n_wall + n_fil)
        A_reg = numpy.vstack([
            numpy.hstack([reg_wall * numpy.eye(num_Ms), numpy.zeros((num_Ms, n_fil))]),    # wall Tikhonov
            numpy.hstack([numpy.zeros((n_fil, num_Ms)), lap_lam * lap_mat]),               # fil Laplacian
        ])

        # --- Ip row: r . x = sum(I_fil) ---
        ip_row = numpy.concatenate([numpy.zeros(num_Ms), numpy.ones(n_fil)])

        prep = _prepare_fixed_rows(A_sensor, A_reg, ip_row)
        prep.update({
            'method': 'laplace',
            'sigma': numpy.array(sigma, dtype=float),
            'num_sensors': num_real_sensors,
            'num_non_fil_coils': num_non_fil_coils,
            'num_Ms': num_Ms,
            'Msc_coils_weighted': Msc_coils_weighted,
            # norm of the Ip row (sqrt(n_fil)); lam is defined for the unit-norm row, so the row weight is lam / ip_row_norm
            'ip_row_norm': numpy.linalg.norm(ip_row),
            'lap_lam': lap_lam,
            'reg_wall': reg_wall,
        })
        return prep

    def solve_laplace(self, prep, psi_at_time, ip_at_time, coil_curr_at_time, lam=None, verbose=False, ip_norm_floor=1.E3):
        '''! Solve one time slice with direct Laplacian regularization, using matrices from @ref prepare_laplace.
        The Ip row and the automatic `lam` (when `lam` is None) are as in @ref solve_tsvd_laplace, including the
        `ip_norm_floor` used there for the automatic `lam` and for `ip_error_pct`.

        @param prep dict, from @ref prepare_laplace
        @param psi_at_time numpy.ndarray, sensor signals at this time (e.g. flux `[Wb]`), shape (num_real_sensors,)
        @param ip_at_time float, total plasma current measurement at time `[A]`
        @param coil_curr_at_time numpy.ndarray, non-filament coil currents at time `[A]`, shape (num_non_fil_coils,)
        @param lam float, Ip constraint weight; None sets it automatically (default: None)
        @param verbose bool, if True prints the weights, Ip error and residual (default: False)
        @param ip_norm_floor float, smallest current used for the automatic `lam` and for `ip_error_pct` `[A]` (default: 1.E3)
        @result numpy.ndarray, filament currents `[A]`, shape (n_fil,)
        @result numpy.ndarray, wall basis amplitudes, shape (n_wall,)
        @result numpy.ndarray, `A x` for the sensor rows only (fitted weighted signals), shape (num_real_sensors,)
        @result dict, diagnostics, see @ref _reconstruction_diagnostics; also `lam`, `lap_lam`, `reg_wall`
        '''
        _check_prep(prep, 'laplace', 'prepare_laplace')
        num_Ms = prep['num_Ms']
        B_weighted = psi_at_time / prep['sigma'] - coil_curr_at_time @ prep['Msc_coils_weighted']
        if lam is None:
            lam = _auto_ip_weight(B_weighted, ip_at_time, ip_norm_floor)

        solution, Ax, total_residual = _solve_with_ip_row(prep, B_weighted, lam / prep['ip_row_norm'], ip_at_time)
        curr_wall = solution[:num_Ms]
        curr_fil = solution[num_Ms:]

        diagnostics = _reconstruction_diagnostics(curr_fil, ip_at_time, ip_norm_floor, Ax, B_weighted, total_residual, solution)
        diagnostics.update({'lam': lam, 'lap_lam': prep['lap_lam'], 'reg_wall': prep['reg_wall']})

        if verbose:
            print(f"lam = {lam:.3e}, lap_lam = {prep['lap_lam']:.3e}, reg_wall = {prep['reg_wall']:.3e}")
            print(f"Ip reconstructed: {diagnostics['ip_reconstructed']:.3f} A | Ip actual: {ip_at_time:.3f} A | Error: {diagnostics['ip_error_pct']:.2f}%")
            print(f"Fit residual: {diagnostics['fit_residual']:.4e}")

        return curr_fil, curr_wall, Ax, diagnostics

    def run_reconstruction_laplace(self, psi_at_time, ip_at_time, coil_curr_at_time, sigma,
                                    num_non_fil_coils, num_real_sensors, rgrid, zgrid,
                                    lam=None, lap_lam=1e-6, reg_wall=1e-5, sigma_r=1e-1, sigma_z=1e-1, gaussian=False, verbose=False,
                                    ip_norm_floor=1.E3):
        '''! Run filament current reconstruction using direct Laplacian regularization for one time slice.
        Same as @ref prepare_laplace followed by @ref solve_laplace. Use this when `sigma` or the regularization
        change every slice; for fixed ones, call @ref prepare_laplace once and @ref solve_laplace in the time loop
        (faster). The filament Laplacian is cached either way.

        @param psi_at_time numpy.ndarray, sensor signals at this time (e.g. flux `[Wb]`), shape (num_real_sensors,)
        @param ip_at_time float, total plasma current measurement at time `[A]`
        @param coil_curr_at_time numpy.ndarray, non-filament coil currents at time `[A]`, shape (num_non_fil_coils,)
        @param sigma numpy.ndarray, noise on each sensor, same units as psi_at_time, shape (num_real_sensors,)
        @param num_non_fil_coils int, number of non-filament coils before the filaments in the Icoil list
        @param num_real_sensors int, number of real sensors (the first num_real_sensors sensors are used)
        @param rgrid numpy.ndarray, R coordinates of the filament grid `[m]`, shape (n_fil,)
        @param zgrid numpy.ndarray, Z coordinates of the filament grid `[m]`, shape (n_fil,)
        @param lam float, Ip constraint weight; None sets it automatically (default: None)
        @param lap_lam float, Laplacian regularization weight (default: 1e-6)
        @param reg_wall float, Tikhonov regularization weight on the wall amplitudes (default: 1e-5)
        @param sigma_r float, R width of the Gaussian weighting `[m]`, only used if gaussian is True (default: 1e-1)
        @param sigma_z float, Z width of the Gaussian weighting `[m]`, only used if gaussian is True (default: 1e-1)
        @param gaussian bool, if True uses the experimental Gaussian-weighted Laplacian (@ref get_gaussian_lap_mat)
            instead of the finite-difference one (@ref get_laplace_matrix) (default: False)
        @param verbose bool, if True prints the weights, Ip error and residual (default: False)
        @param ip_norm_floor float, smallest current used for the automatic `lam` and for `ip_error_pct` `[A]` (default: 1.E3)
        @result numpy.ndarray, filament currents `[A]`, shape (n_fil,)
        @result numpy.ndarray, wall basis amplitudes, shape (n_wall,)
        @result numpy.ndarray, `A x` for the sensor rows only (fitted weighted signals), shape (num_real_sensors,)
        @result dict, diagnostics, see @ref solve_laplace
        '''
        prep = self.prepare_laplace(sigma, num_non_fil_coils, num_real_sensors, rgrid, zgrid, lap_lam, reg_wall,
                                    sigma_r, sigma_z, gaussian)
        return self.solve_laplace(prep, psi_at_time, ip_at_time, coil_curr_at_time, lam, verbose, ip_norm_floor)

    # ==========================================================
    # Flux probe placement (positions on the TokaMaker mesh of self.mygs)
    # ==========================================================

    def get_bbox_grid_pts(self, n_r=10, n_z=10, verbose=False):
        '''! Uniform rectangular grid over the bounding box of the TokaMaker mesh nodes (`self.mygs`).
        Used to place the flux probes whose wall psi is later interpolated onto every TokaMaker node
        (see @ref get_wall_psi_tidx).

        @param n_r int, number of points along R (default: 10)
        @param n_z int, number of points along Z (default: 10)
        @param verbose bool, if True, plots the grid, the TokaMaker nodes and the limiter contour (default: False)
        @result numpy.ndarray, bounding-box grid points `[m]`, shape (n_r*n_z, 2) with columns [R, Z]
        @result numpy.ndarray, R meshgrid `[m]`, shape (n_z, n_r)
        @result numpy.ndarray, Z meshgrid `[m]`, shape (n_z, n_r)

        Example
        -------
        >>> jam_obj = JAMfit(xml_file, thincurr_meshfile, tokamaker_meshfile, B0, R0, oft_env=myOFT)
        >>> probe_pts, R, Z = jam_obj.get_bbox_grid_pts(n_r=10, n_z=20)
        >>> probes = [circular_flux_loop(r, z, 'psi_grid_{0}'.format(i)) for i, (r, z) in enumerate(probe_pts)]
        '''
        r_pts_grid = self.mygs.r[:, 0]
        z_pts_grid = self.mygs.r[:, 1]

        r_vals = numpy.linspace(r_pts_grid.min(), r_pts_grid.max(), n_r)
        z_vals = numpy.linspace(z_pts_grid.min(), z_pts_grid.max(), n_z)
        R, Z = numpy.meshgrid(r_vals, z_vals, indexing='xy')
        bbox_pts = numpy.column_stack([R.ravel(), Z.ravel()])   # row-major: R varies fastest

        if verbose:
            fig, ax = plt.subplots()
            ax.scatter(r_pts_grid, z_pts_grid, c='lightgray', s=5, alpha=0.5, label='TokaMaker grid')
            ax.scatter(bbox_pts[:, 0], bbox_pts[:, 1], c='tab:blue', s=15, label='Bbox grid pts')
            ax.plot(self.mygs.lim_contour[:, 0], self.mygs.lim_contour[:, 1], 'k-', label='Limiter Contour')
            ax.set_xlabel('R')
            ax.set_ylabel('Z')
            ax.set_title(f'Bounding Box Grid ({n_r} x {n_z} points)')
            ax.legend()
            plt.tight_layout()
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

        return bbox_pts, R, Z

    def get_inside_limiter_pts(self, stride=1, verbose=False):
        '''! TokaMaker mesh nodes (`self.mygs`) inside the limiter contour, to place flux probes for the wall psi.
        Every stride-th inside node gets a probe; see @ref get_wall_psi_tidx.

        @param stride int, subsampling stride for the nodes inside the limiter (default: 1)
        @param verbose bool, if True, plots the nodes with and without a probe and the limiter contour (default: False)
        @result numpy.ndarray, probe points `[m]`, shape (n_probes, 2) with columns [R, Z]
        @result numpy.ndarray, bool mask over all nodes, True where a probe is placed, shape (n_nodes,)
        @result numpy.ndarray, bool mask over all nodes, True inside the limiter, shape (n_nodes,)

        Example
        -------
        >>> probe_pts, probe_mask, inside_mask = jam_obj.get_inside_limiter_pts(stride=4)
        '''
        node_pts = self.mygs.r[:, :2]
        lim_path = Path(self.mygs.lim_contour[:, :2])

        inside_mask = lim_path.contains_points(node_pts)   # full inside-limiter mask
        inside_idx = numpy.where(inside_mask)[0]

        # subsample the indices that actually get a flux loop
        sparse_idx = inside_idx[::stride]
        sparse_mask = numpy.zeros_like(inside_mask)
        sparse_mask[sparse_idx] = True

        inside_lim_pts = node_pts[sparse_idx, :].copy()   # copy: mygs.r is a view on TokaMaker memory

        if verbose:
            fig, ax = plt.subplots()
            ax.scatter(node_pts[:, 0], node_pts[:, 1], c=sparse_mask, cmap='RdYlGn', s=10, alpha=0.7)
            ax.plot(self.mygs.lim_contour[:, 0], self.mygs.lim_contour[:, 1], 'k-', label='Limiter Contour')
            ax.set_xlabel('R')
            ax.set_ylabel('Z')
            ax.set_title(f'Grid Points With Flux Loops (stride={stride})')
            ax.legend()
            plt.colorbar(ax.collections[0], ax=ax, label='Has flux loop')
            plt.tight_layout()
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

        return inside_lim_pts, sparse_mask, inside_mask

    # ==========================================================
    # Post-processing and visualization methods
    # ==========================================================

    def get_wall_psi_tidx(self, num_sensors, solution_wall_at_time):
        '''! Compute the wall contribution to the flux at a given time index.
        Uses the fitted wall basis amplitudes and the wall-to-sensor matrix `Ms` to evaluate the wall flux on
        the extra sensors after the first num_sensors (the flux probes, placed with @ref get_bbox_grid_pts
        or @ref get_inside_limiter_pts). The division by `2 pi` converts flux `[Wb]` to TokaMaker psi `[Wb/rad]`.
        Before @ref post_process_tidx the result must be put on every TokaMaker node (the reconstruct example
        fills the nodes without a probe by interpolation).

        @param num_sensors int, number of real sensors to skip (the probes start at index num_sensors)
        @param solution_wall_at_time numpy.ndarray, wall basis amplitudes at the given time index, shape (n_wall,)
        @result numpy.ndarray, wall psi on the flux probes `[Wb/rad]`, shape (n_probes,)'''
        wall_psi_probes = self.torus_reduced.Ms[:, num_sensors:]/(2*numpy.pi)
        wall_psi_at_time = solution_wall_at_time @ wall_psi_probes
        return wall_psi_at_time

    def compute_psi_basis(self, rgrid, zgrid, disk_radius=0.1):
        '''! Plasma psi per unit current of each filament, for reuse in post_process_tidx (psi_basis).

        Each filament is represented by a TokaMaker init_psi disk of radius disk_radius centred on (r, z), exactly
        as in post_process_tidx. init_psi is linear in Ip (checked to 2e-16 relative), so the plasma psi of any
        filament current vector is psi_basis @ I_fil. This avoids repeating the init_psi loop for every time slice,
        which matters for grids with thousands of filaments.
        
        @param rgrid numpy.ndarray, R coordinates of the filament mesh `[m]`, shape (n_fil,)
        @param zgrid numpy.ndarray, Z coordinates of the filament mesh `[m]`, shape (n_fil,)
        @param disk_radius float, radius of the init_psi disk per filament `[m]` (default: 0.1, see @ref post_process_tidx)
        @result numpy.ndarray, psi per amp on the TokaMaker nodes `[Wb/rad per A]`, shape (n_nodes, n_fil)

        Example
        -------
        >>> psi_basis = jamfit.compute_psi_basis(rgrid, zgrid, disk_radius=0.1)
        >>> out = jamfit.post_process_tidx(fil, wall_psi, coils, rgrid, zgrid, disk_radius=0.1, psi_basis=psi_basis)
        '''
        psi_basis = []
        for r, z in zip(rgrid, zgrid):
            self.mygs.set_coil_currents()
            self.mygs.set_targets(Ip=1.0)
            self.mygs.init_psi(r, z, disk_radius, 1.0, 0.0)
            psi_basis.append(self.mygs.get_psi(False))
        return numpy.array(psi_basis).T

    def post_process_tidx(self, filaments_at_time, wall_psi, coil_curr_dict, rgrid, zgrid, verbose=False, include_negative_encl=True,
                          disk_radius=0.1, psi_basis=None, touch_tol=0.0005):
        '''! Post-process the results at a given time index to compute plasma parameters and visualize.

        The plasma psi is the sum over positive filaments of a TokaMaker init_psi disk of radius disk_radius centred
        on each filament. Negative filaments are left out on purpose: TokaMaker init_psi cannot represent a negative
        plasma current, and negative toroidal current is unphysical here (a well-tuned reconstruction should have
        essentially none). Where the total current is needed (li), it is the current enclosed by the LCFS, see
        include_negative_encl; the current centroid uses all filaments, negatives included.

        Choosing disk_radius: each disk spreads its filament's current over the area around it, so the natural size is
        about one filament grid cell (radius = spacing / sqrt(pi), 0.085 m for the 0.15 m DIII-D example grid; the
        default 0.1 m is close to that). Much larger disks stick out past the limiter; that part is dropped, which moves
        the filament's current inward (about 100 mm at 0.3 m). Much smaller disks behave like point filaments and leave
        saddle points between neighbours that are taken as X-points (seen at 0.05 m). Scale disk_radius with the grid spacing.

        The internal inductance uses the 'std' normalisation of TokaMaker get_stats,
        li = (int B_p^2 dV / V) / (mu0 Ip / L_LCFS)^2, where the B_p^2 integral and V only cover the region inside
        the LCFS. Ip is therefore the current enclosed by the traced LCFS: the sum of the filaments whose centre
        (rgrid, zgrid) lies inside lcfs_points, not the sum of all filaments.

        @param filaments_at_time numpy.ndarray, filament currents at the given time index `[A]`, shape (n_fil,)
        @param wall_psi numpy.ndarray, wall psi on every TokaMaker node `[Wb/rad]`, shape (n_nodes,)
            (see @ref get_wall_psi_tidx)
        @param coil_curr_dict dict, TokaMaker coil name -> coil current `[A]` at the given time index
        @param rgrid numpy.ndarray, R coordinates of the filament mesh `[m]`, shape (n_fil,)
        @param zgrid numpy.ndarray, Z coordinates of the filament mesh `[m]`, shape (n_fil,)
        @param verbose bool, if True, plots the equilibrium and LCFS (default: False)
        @param include_negative_encl bool, if True, negative filaments inside the LCFS count towards the enclosed Ip
            used for li; if False only positive filaments are summed, consistent with the psi used here (negative
            filaments are left out of psi) (default: True)
        @param disk_radius float, radius of the init_psi disk per filament `[m]`; about one filament grid cell, see above (default: 0.1)
        @param psi_basis numpy.ndarray or None, psi per amp of each filament, shape (n_nodes, n_fil), from
            compute_psi_basis with the same disk_radius; if None the init_psi loop is run here (default: None)
        @param touch_tol float, largest LCFS-limiter distance that counts as touching, passed to
            @ref find_limiting_point `[m]` (default: 0.0005)
        @result numpy.ndarray or None, LCFS points `[m]`, shape (n_pts, 2), traced at psi_N = 0.99 if the plasma is diverted and
            at psi_N = 1 if it is limited (None if the trace fails)
        @result tuple or None, `(R, Z)` `[m]` of the point where the LCFS touches the limiter (see
            @ref find_limiting_point); None if the LCFS does not touch the limiter or there is no LCFS
        @result float or None, psi at the traced LCFS `[Wb/rad]`
        @result numpy.ndarray, total psi (filaments + coils + wall) on the TokaMaker nodes `[Wb/rad]`, shape (n_nodes,)
        @result numpy.ndarray or None, safety factor q on TokaMaker's default psi_N sampling
        @result float, q95 (nan if there is no LCFS)
        @result float, internal inductance li, 'std' normalisation (nan if there is no LCFS or Ip_encl <= 0)
        @result tuple, current centroid `(R_c, Z_c)` `[m]` from all filaments (see @ref calc_current_centroid);
            `(nan, nan)` if the net current is zero or the centroid falls outside the limiter
        @result tuple or float, geometric LCFS centroid `(R_c, Z_c)` `[m]` (nan if there is no LCFS)
        @result float, area enclosed by the LCFS `[m^2]` (nan if there is no LCFS)'''
        # reset mygs in case there are left over coil currents from previous runs
        self.mygs.set_coil_currents()
        limiter = self.mygs.lim_contour

        # --- Current centroid (signed weights), NaN if it is not meaningful ---
        # R_c = sum(I R) / sum(I). With mixed-sign currents and a small net current (e.g. at the end of a current
        # quench) the division by sum(I) can put the centroid far outside the vessel, and it is undefined at
        # sum(I) = 0. Both cases are returned as (nan, nan); a centroid outside the limiter cannot be the centre of
        # a plasma current, so no tuned current threshold is needed.
        net_curr = numpy.sum(filaments_at_time)
        current_cent = (numpy.nan, numpy.nan)
        if net_curr != 0.0:
            centroid = calc_current_centroid(rgrid, zgrid, filaments_at_time)
            if Path(limiter[:, :2]).contains_point(centroid):
                current_cent = centroid

        # --- Filament (plasma) psi ---
        if psi_basis is not None:
            pos = filaments_at_time > 0
            psi_total_fil = psi_basis[:, pos] @ filaments_at_time[pos]
        else:
            fil_points = list(zip(rgrid, zgrid))
            psi_fil = []
            for filcount, (r, z) in enumerate(fil_points):
                self.mygs.set_coil_currents()
                if filaments_at_time[filcount] > 0:
                    self.mygs.set_targets(Ip=filaments_at_time[filcount])
                    self.mygs.init_psi(r, z, disk_radius, 1.0, 0.0)
                    psi_fil.append(self.mygs.get_psi(False))
            psi_fil = numpy.array(psi_fil)
            psi_total_fil = numpy.sum(psi_fil, axis=0)
        # --- Coil psi (vacuum solve with the coil currents) ---
        self.mygs.set_coil_currents(coil_curr_dict)
        psi_coils_obj = self.mygs.vac_solve()
        psi_coils = psi_coils_obj.get_psi(normalized=False)

        # --- Total psi and LCFS ---
        total_psi = psi_total_fil + psi_coils + wall_psi
        self.mygs.set_psi(total_psi, update_bounds=True)

        # --- Trace the LCFS ---
        # Diverted: psi_N = 1 is the separatrix, which passes through the X-point, so tracing exactly there succeeds or
        # fails on round-off (and changes the area by a few % when it does succeed). Trace just inside it, at 0.99, as
        # TokaMaker get_stats does by default (lcfs_pad). Limited: psi_N = 1 touches the limiter and traces cleanly.
        psi_n_lcfs = 0.99 if self.mygs.diverted else 1.0
        lcfs_points = self.mygs.trace_surf(psi_n_lcfs)
        psiatlcfs = self.mygs.psinorm_to_absolute(psi_n_lcfs)

        # --- q, li and LCFS geometry (defaults if there is no LCFS) ---
        limiting_pt = None
        q_vals = None
        q95 = numpy.nan
        internal_inductance = numpy.nan
        area_cent = numpy.nan
        area = numpy.nan

        if lcfs_points is not None:
            lim_R, lim_Z, _ = find_limiting_point(lcfs_points, limiter, touch_tol=touch_tol)
            if lim_R is not None:
                limiting_pt = (lim_R, lim_Z)
            _, q_vals, _, _, _, _= self.mygs.get_q()

            # --- Ip enclosed by the traced LCFS, for li ---
            # get_stats integrates B_p^2 only inside the LCFS, so li must be normalised by the current inside it.
            # Late in a VDE much of the filament current sits outside the LCFS, and using the total (ip) makes li
            # too small by (Ip_encl/Ip_total)^2.
            fil_inside = Path(lcfs_points[:, :2]).contains_points(numpy.column_stack((rgrid, zgrid)))
            if include_negative_encl:
                # Negative filaments are unphysical for this plasma and should vanish for a good enough fit, but
                # they are part of the fitted solution, so they are counted for now (they are still left out of psi).
                ip_encl = numpy.sum(filaments_at_time[fil_inside])
            else:
                ip_encl = numpy.sum(filaments_at_time[fil_inside & (filaments_at_time > 0)])

            if ip_encl > 0:
                internal_inductance = self.mygs.get_stats(beta_Ip = ip_encl)['l_i']
            else:
                # get_stats divides by Ip^2 with numpy (inf, no exception), so catch it here
                print(f"Couldn't compute li, Ip enclosed by the LCFS was {ip_encl}, skipping this one")
                internal_inductance = numpy.nan

            _, q95_arr, _, _, _, _ = self.mygs.get_q(psi=0.95)
            q95 = float(q95_arr[0])  # get_q returns q at the sampled psi as an array, here of length 1
            area, area_cent = calc_lcfs_geo(lcfs_points)
        else:
            psiatlcfs = None

        # --- Plot ---
        if verbose:
            fig, ax = plt.subplots(1,1)
            self.mygs.plot_machine(fig,ax, cond_color='blue')
            self.mygs.plot_psi(fig,ax ,plasma_nlevels = 75, normalized=False)
            plt.gca().set_aspect('equal', adjustable='box')
            if lcfs_points is not None:
                ax.plot(lcfs_points[:,0], lcfs_points[:,1], 'r--', label='LCFS')

        return lcfs_points, limiting_pt, psiatlcfs, total_psi, q_vals, q95, internal_inductance, current_cent, area_cent, area

# ==========================================================
# JAMfit Helper Functions
# ==========================================================

def _prepare_fixed_rows(A_sensor, A_reg, ip_row):
    '''! Factor the rows of a reconstruction that do not change between time slices (see the solver notes in JAMfit).
    The fixed system is `A0 = [A_sensor; A_reg]` with right-hand side `[B_weighted; 0]`. A0 is QR factorized
    (`A0 = Q R`), so the least-squares solution without the Ip row is `x0 = R^-1 Q_sensor^T B_weighted`; the
    matrix `P = R^-1 Q_sensor^T` is stored so each slice is one matrix-vector product. No `A0^T A0` is formed,
    which would square the condition number.

    The Ip row `c (r . x - Ip)` is added per slice with the exact one-row least-squares update (Sherman-Morrison):
    `x = x0 + g c^2 (Ip - r . x0) / (1 + c^2 r . g)`, with `g = (A0^T A0)^-1 r^T = R^-1 R^-T r^T` stored here.

    @param A_sensor numpy.ndarray, sigma-weighted sensor rows, shape (n_sensors, n_unknowns)
    @param A_reg numpy.ndarray, regularization rows (right-hand side zero), shape (n_reg, n_unknowns)
    @param ip_row numpy.ndarray, Ip row `r` with `r . x = sum(I_fil)`, shape (n_unknowns,)
    @result dict, `A_sensor`, `A_reg`, `ip_row`, `P` (n_unknowns, n_sensors), `g` (n_unknowns,), `rg` (`r . g`)
    '''
    n_sensors = A_sensor.shape[0]
    Q, R = numpy.linalg.qr(numpy.vstack([A_sensor, A_reg]))                                # reduced QR
    diag_R = numpy.abs(numpy.diag(R))
    if diag_R.min() <= 1.E-14 * diag_R.max():
        raise ValueError("The sensor and regularization rows do not determine every unknown "
                         "(rank-deficient system); increase the regularization weights.")
    P = numpy.linalg.solve(R, Q[:n_sensors, :].T)                                          # (n_unknowns, n_sensors)
    g = numpy.linalg.solve(R, numpy.linalg.solve(R.T, ip_row))                             # (n_unknowns,)
    return {'A_sensor': A_sensor, 'A_reg': A_reg, 'ip_row': ip_row, 'P': P, 'g': g, 'rg': float(ip_row @ g)}


def _solve_with_ip_row(prep, B_weighted, ip_row_weight, ip_at_time):
    '''! Least-squares solution of the prepared fixed rows plus the Ip row `ip_row_weight * (r . x - ip_at_time)`.
    See @ref _prepare_fixed_rows for the method.

    @param prep dict, from @ref _prepare_fixed_rows (or a `prepare_*` method)
    @param B_weighted numpy.ndarray, weighted sensor right-hand side, shape (n_sensors,)
    @param ip_row_weight float, weight `c` of the Ip row for this slice
    @param ip_at_time float, measured plasma current `[A]`
    @result numpy.ndarray, solution `x`, shape (n_unknowns,)
    @result numpy.ndarray, `A_sensor @ x` (fitted weighted signals), shape (n_sensors,)
    @result float, norm of the misfit over all rows (sensors, Ip row and regularization rows)
    '''
    c2 = ip_row_weight**2
    x0 = prep['P'] @ B_weighted                                                            # fit without the Ip row
    solution = x0 + prep['g'] * (c2 * (ip_at_time - prep['ip_row'] @ x0) / (1.0 + c2 * prep['rg']))
    Ax = prep['A_sensor'] @ solution
    total_residual = numpy.sqrt(numpy.sum((Ax - B_weighted)**2)
                                + numpy.sum((prep['A_reg'] @ solution)**2)
                                + c2 * (prep['ip_row'] @ solution - ip_at_time)**2)
    return solution, Ax, total_residual


def _check_prep(prep, method, prepare_name):
    '''! Raise a clear error if `prep` did not come from the matching prepare method.

    @param prep dict, prepared matrices passed to a `solve_*` method
    @param method str, expected value of `prep['method']`
    @param prepare_name str, name of the matching prepare method, for the error message
    '''
    if not isinstance(prep, dict) or prep.get('method') != method:
        raise ValueError("prep must come from {0}()".format(prepare_name))


def _reconstruction_diagnostics(curr_fil, ip_at_time, ip_norm_floor, Ax, B_weighted, total_residual, solution):
    '''! Diagnostics shared by all `solve_*` methods.

    @param curr_fil numpy.ndarray, reconstructed filament currents `[A]`, shape (n_fil,)
    @param ip_at_time float, measured plasma current `[A]`
    @param ip_norm_floor float, smallest current used to normalize `ip_error_pct` `[A]`
    @param Ax numpy.ndarray, fitted weighted sensor signals, shape (n_sensors,)
    @param B_weighted numpy.ndarray, weighted sensor right-hand side, shape (n_sensors,)
    @param total_residual float, misfit norm over all rows (sensors, Ip and regularization)
    @param solution numpy.ndarray, full solution vector of the solver's unknowns
    @result dict, `ip_reconstructed` `[A]`, `ip_actual` `[A]`, `ip_error_pct` (`100 |sum(I_fil) - Ip| / max(|Ip|, ip_norm_floor)`),
        `fit_residual` (norm of the weighted sensor misfit), `fit_residual_nonorm` (the misfit vector `Ax - B`),
        `total_residual`, `solution`
    '''
    ip_reconstructed = curr_fil.sum()
    fit_residual_nonorm = Ax - B_weighted
    return {
        'ip_reconstructed':     ip_reconstructed,
        'ip_actual':            ip_at_time,
        'ip_error_pct':         100 * numpy.abs(ip_reconstructed - ip_at_time) / max(abs(ip_at_time), ip_norm_floor),
        'fit_residual':         numpy.linalg.norm(fit_residual_nonorm),
        'fit_residual_nonorm':  fit_residual_nonorm,
        'total_residual':       total_residual,
        'solution':             solution,
    }


def _auto_ip_weight(B_weighted, ip_at_time, ip_norm_floor=1.E3):
    '''! Automatic weight of the Ip constraint row in the Laplacian reconstructions.
    Sets the weight from the order-of-magnitude difference between the typical weighted sensor signal and
    the total plasma current, so the Ip row is on the same order as the sensor rows:
    `lam = 100 * 10**floor(log10(mean(|B_weighted|) / ip_norm))`, with `ip_norm = max(|Ip|, ip_norm_floor)`.
    The factor 100 gives the single Ip row more weight against the many sensor rows. The floor (as in
    lstsq) keeps `lam` finite when Ip passes through zero; without it, Ip = 0 gave `lam ~ 1e29`, which made
    the normal equations numerically singular.

    @param B_weighted numpy.ndarray, sigma-weighted sensor signals with the coil contribution removed, shape (num_sensors,)
    @param ip_at_time float, total plasma current measurement at time `[A]`
    @param ip_norm_floor float, smallest current used in place of |Ip| `[A]` (default: 1.E3)
    @result float, Ip constraint weight `lam` (100 if the mean signal is zero or not finite)

    Example
    -------
    >>> _auto_ip_weight(numpy.array([2.0, 4.0]), 1.E6)    # mean |B| = 3, ratio 3e-6 -> 100 * 1e-6 = 1e-4
    >>> _auto_ip_weight(numpy.zeros(3), 1.E6)              # no signal -> 100
    '''
    ip_norm = max(abs(ip_at_time), ip_norm_floor)
    mean_signal = numpy.mean(numpy.abs(B_weighted))
    # log10 is undefined if the signals are all zero (e.g. t = 0 of a synthetic run with no plasma or wall
    # current, or every sensor masked); fall back to a magnitude difference of 0, i.e. lam = 100
    if (not numpy.isfinite(mean_signal)) or (mean_signal <= 0.0):
        magnitude_diff = 0
    else:
        magnitude_diff = math.floor(math.log10(mean_signal / ip_norm))
    return 100 * 10**magnitude_diff


def setup_synthetic_current(timepoints, ip_list, sigma_r, sigma_z, r0, z0, rgrid, zgrid):
    '''! Generate synthetic filament currents using a Gaussian plasma distribution.
    At each time step, spreads the total plasma current across the filament grid
    using a 2D Gaussian centered at (r0, z0) with widths (sigma_r, sigma_z).

    The Gaussian is normalised over the grid, so the filament currents sum exactly to `ip_list[i]`:
    `I_k = Ip * g_k / sum(g)`, `g_k = exp(-((R_k - r0)^2 / (2 sigma_r^2) + (Z_k - z0)^2 / (2 sigma_z^2)))`.

    @param timepoints list or numpy.ndarray, time values `[s]`, length ntimes
    @param ip_list list or numpy.ndarray, total plasma current at each time step `[A]`
    @param sigma_r float, Gaussian width in the R direction `[m]`
    @param sigma_z float, Gaussian width in the Z direction `[m]`
    @param r0 list or numpy.ndarray, R position of the Gaussian centre at each time step `[m]`
    @param z0 list or numpy.ndarray, Z position of the Gaussian centre at each time step `[m]`
    @param rgrid numpy.ndarray, R coordinates of the filament mesh `[m]`, shape (n_fil,)
    @param zgrid numpy.ndarray, Z coordinates of the filament mesh `[m]`, shape (n_fil,)
    @result numpy.ndarray, shape (ntimes, 1 + n_fil), time column `[s]` followed by filament currents `[A]`
    '''
    fil_currs = []
    for i in range(len(timepoints)):
        gaussian_raw = numpy.exp(
            -((rgrid - r0[i])**2 / (2 * sigma_r**2) + (zgrid - z0[i])**2 / (2 * sigma_z**2))
        )
        gaussian_values = ip_list[i] * (gaussian_raw / numpy.sum(gaussian_raw))
        fil_currs.append(gaussian_values)

    fil_currs = numpy.array(fil_currs)                       # (ntimes, n_fil)
    time_column = numpy.array(timepoints).reshape(-1, 1)
    fil_currs = numpy.hstack((time_column, fil_currs))       # (ntimes, 1 + n_fil)
    return fil_currs

def interpolate_total_current(coil_currs, nsteps, verbose=False):
    '''! Interpolate total plasma current to a higher-resolution time grid.
    Linearly interpolates each current column onto a finer time grid, then sums the columns
    at each time to get the total current.

    @param coil_currs numpy.ndarray, shape (ntimes, 1 + ncols), first column is time `[s]`, the rest are
        currents `[A]` (e.g. filament currents) to be summed
    @param nsteps int, number of high-resolution steps between first and last time
    @param verbose bool, if True, plots the original and interpolated total current (default: False)
    @result numpy.ndarray, high-resolution time grid `[s]`, shape (nsteps + 1,)
    @result numpy.ndarray, total current on that grid `[A]`, shape (nsteps + 1,)
    '''
    times = coil_currs[:, 0]
    currents = coil_currs[:, 1:]

    high_res_time = numpy.linspace(times[0], times[-1], nsteps + 1)

    interpolated_currents = numpy.array([
        numpy.interp(high_res_time, times, col)
        for col in currents.T
    ]).T  # shape: (nsteps+1, ncols)

    total_current_high_res = numpy.sum(interpolated_currents, axis=1)

    if verbose:
        plt.figure(figsize=(8, 5))
        plt.scatter(times, numpy.sum(currents, axis=1), color='red', label='Original Data')
        plt.plot(high_res_time, total_current_high_res, color='blue', label='Interpolated Total Current')
        plt.xlabel('Time [s]')
        plt.ylabel('Total Current [A]')
        plt.title('Total Current vs Time with Higher Resolution')
        plt.legend()
        plt.grid(True)
        plt.show()

    return high_res_time, total_current_high_res


def get_laplace_matrix(rgrid, zgrid, verbose=False):
    '''! Finite-difference (5-point) Laplacian on a filament grid with holes (e.g. a D-shaped grid).
    The filaments must lie on a uniform rectangular grid (spacing dr, dz), with points outside the limiter
    removed. Coordinates are matched to grid indices after rounding to 8 decimals.

    Row for grid point (i, j): `(I[i+1,j] + I[i-1,j]) / dr^2 + (I[i,j+1] + I[i,j-1]) / dz^2 - 2 (1/dr^2 + 1/dz^2) I[i,j]`, with
    off-grid neighbours left out. The diagonal always has the full `-2 (1/dr^2 + 1/dz^2)`, also at the edge
    of the grid, which acts like a zero-current boundary just outside the grid.

    @param rgrid numpy.ndarray, R coordinates of the filament grid `[m]`, shape (n_fil,)
    @param zgrid numpy.ndarray, Z coordinates of the filament grid `[m]`, shape (n_fil,)
    @param verbose bool, if True, prints debug information (default: False)
    @result numpy.ndarray, Laplacian matrix `[1/m^2]`, shape (n_fil, n_fil)
    @result int, number of filaments n_fil
    '''
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

        lap_mat[k, k] = -2.0*(1.0/(dr**2) + 1.0/(dz**2)) 

    return lap_mat, N


def find_limiting_point(lcfs_points, limiter, touch_tol=0.0005):
    '''! Find where the LCFS touches the limiter (the limiting point), if it does.
    Takes the LCFS point closest to any limiter point; it counts as touching if that distance is at most
    touch_tol.

    @param lcfs_points numpy.ndarray, LCFS points `[m]`, shape (n_pts, 2) with columns [R, Z]
    @param limiter numpy.ndarray, limiter contour points `[m]`, shape (n_lim, 2) with columns [R, Z]
    @param touch_tol float, largest LCFS-limiter distance that counts as touching `[m]` (default: 0.0005)
    @result float or None, R of the limiting point `[m]` (None if the LCFS does not touch)
    @result float or None, Z of the limiting point `[m]` (None if the LCFS does not touch)
    @result float or None, LCFS-limiter distance at that point `[m]` (None if the LCFS does not touch)
    '''

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
    '''! Calculate area and geometric centroid of the LCFS with the shoelace formula.
    With `c_i = R_i Z_i+1 - R_i+1 Z_i` and signed area `A = sum(c_i) / 2`:
    `R_c = sum((R_i + R_i+1) c_i) / (6 A)`, `Z_c = sum((Z_i + Z_i+1) c_i) / (6 A)`. The contour is closed
    automatically if the last point differs from the first; either point ordering works.

    @param lcfs_points numpy.ndarray, LCFS points `[m]`, shape (N, 2) with columns [R, Z]
    @result float, area enclosed by the LCFS `[m^2]` (None if lcfs_points is None)
    @result tuple, geometric centroid `(R_c, Z_c)` `[m]` (None if lcfs_points is None)
    '''
    if lcfs_points is None:
        return None, None

    R = numpy.array(lcfs_points[:, 0])
    Z = numpy.array(lcfs_points[:, 1])

    # Ensure contour is closed
    if not (numpy.isclose(R[0], R[-1]) and numpy.isclose(Z[0], Z[-1])):
        R = numpy.append(R, R[0])
        Z = numpy.append(Z, Z[0])

    cross = R[:-1] * Z[1:] - R[1:] * Z[:-1]    # Calculates area using Shoelace method, Shoelace terms

    # Signed area: positive for counter-clockwise points in (R, Z), negative for clockwise
    signed_area = 0.5 * numpy.sum(cross)
    area = numpy.abs(signed_area)

    # Centroid: R_c = sum((R_i + R_i+1) * cross_i) / (6 A), same for Z, with A the SIGNED area so the point
    # ordering cancels. Do not take abs() of the sums: that flips Z_c > 0 for a plasma below the midplane.
    R_c = numpy.sum((R[:-1] + R[1:]) * cross) / (6 * signed_area)
    Z_c = numpy.sum((Z[:-1] + Z[1:]) * cross) / (6 * signed_area)

    return area, (R_c, Z_c)

def calc_current_centroid(R_fil, Z_fil, I_fil):
    '''! Calculate the current centroid of a filament grid using signed current weights
        This variant uses signed currents as weights, so opposing currents can partially
        cancel out the centroid position. Useful for computing the net moment of the
        current distribution.

        Calculation: R_c = sum(I_i * R_i) / sum(I_i)
                    Z_c = sum(I_i * Z_i) / sum(I_i)

        Raises ValueError if the net current is exactly zero.

        @param R_fil numpy.ndarray, 1D or 2D array of R (radial) coordinates of filament points `[m]`
        @param Z_fil numpy.ndarray, 1D or 2D array of Z (vertical) coordinates of filament points `[m]`
        @param I_fil numpy.ndarray, 1D or 2D array of currents at each filament point `[A]`, must match shape of R and Z
        @result tuple, current centroid `(R_c, Z_c)` `[m]`
    '''
    R = numpy.asarray(R_fil)  # Ensure inputs are numpy arrays
    Z = numpy.asarray(Z_fil)
    I_fil = numpy.asarray(I_fil)

    # Verify shapes match
    if not (R.shape == Z.shape == I_fil.shape):
        raise ValueError(
            f"R, Z, and I must have the same shape. "
            f"Got R: {R.shape}, Z: {Z.shape}, I: {I_fil.shape}"
        )

    R_flat = R.ravel()
    Z_flat = Z.ravel()
    I_flat = I_fil.ravel()

    I_net = numpy.sum(I_flat)
    if I_net == 0:
        raise ValueError("Net current is zero - cannot calculate centroid")

    R_centroid = numpy.sum(I_flat * R_flat) / I_net
    Z_centroid = numpy.sum(I_flat * Z_flat) / I_net

    return (R_centroid, Z_centroid)

def get_gaussian_lap_mat(rgrid, zgrid, sigma_r=8e-2, sigma_z=8e-2):
    '''! Construct a Gaussian-weighted Laplacian regularization matrix for the filament grid (Experimental).
    Each filament point is connected to every other point with a weight that is a Gaussian function of the
    distance: `W_ij = exp(-0.5 ((dR_ij / sigma_r)^2 + (dZ_ij / sigma_z)^2))`, `W_ii = 0`, and the matrix is the
    graph Laplacian `L = D - W` with `D = diag(sum_j W_ij)`. Note the sign convention is opposite to
    @ref get_laplace_matrix (positive diagonal here, negative there); this does not matter for a
    least-squares penalty `|lap_lam L I|^2`.

    @param rgrid numpy.ndarray, R coordinates of the filament grid `[m]`, shape (n_fil,)
    @param zgrid numpy.ndarray, Z coordinates of the filament grid `[m]`, shape (n_fil,)
    @param sigma_r float, Gaussian width in the R direction `[m]` (default: 8e-2)
    @param sigma_z float, Gaussian width in the Z direction `[m]` (default: 8e-2)
    @result numpy.ndarray, Gaussian-weighted Laplacian matrix (dimensionless), shape (n_fil, n_fil)
    @result int, number of filaments n_fil
    '''
    N = len(rgrid)
    dr = rgrid[:, None] - rgrid[None, :]  # (N, N)
    dz = zgrid[:, None] - zgrid[None, :]  # (N, N)

    W = numpy.exp(-0.5 * ((dr / sigma_r)**2 + (dz / sigma_z)**2))
    numpy.fill_diagonal(W, 0)  # no self-coupling

    D = numpy.diag(W.sum(axis=1))
    lap_mat = D - W

    return lap_mat, N
