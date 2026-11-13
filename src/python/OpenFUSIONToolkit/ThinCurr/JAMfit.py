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
    >>> fil, wall, res, Ax, B = jam.run_reconstruction_lstsq(psi, ip, n_coils, coil_currs, ip_weight, sigma, 1e-6, 1e-5)
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

    def create_from_runTD_top_modes(self, num_modes, reduced_filename, coil_currs, dt, nsteps, initial_num_eigs=50, verbose=False, s_freq=10, p_freq=10, axisymmetric_modes=None, freq_response=None):
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
        @param verbose bool, if True, plots mode amplitudes and prints mode info (default: False)
        @param s_freq int, print run status every s_freq steps (default: 10)
        @param p_freq int, save plot output every p_freq steps (default: 10)
        @param axisymmetric_modes list or array-like, indices of the axisymmetric modes to allow
                                (default: None -> allow all computed modes)
        @param freq_response dict, keyword arguments for @ref compute_freq_response_basis (at least
               `fil_indices` and `num_non_fil_coils`). If given, the frequency-response vectors are
               appended after the selected eigenmodes and always kept. The vectors and their labels
               are stored in `self.freq_vecs` and `self.freq_labels` (default: None -> eigenmodes only)

        @result ThinCurr_reduced, the constructed reduced model object
        @result numpy.ndarray, sensor signals from the preliminary run
        @result dict, mode currents from the preliminary run (`'time'` `[s]`, `'curr'` shape (ntimes, initial_num_eigs))
        @result numpy.ndarray, basis of the reduced model `(n_basis, nelems)`: the selected eigenmodes, followed by any frequency-response vectors (if verbose is True)
        @result list, indices of the selected dominant modes (if verbose is True)
        @result list, maximum weight amplitudes of the selected dominant modes (if verbose is True)
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

        # --- Plot all mode amplitudes, selected modes in colour and the rest in grey ---
        if verbose:
            fig, ax = plt.subplots(1, 1)
            self.torus.build_XDMF()
            for i in range(temp_curr.shape[1]):
                if i in eig_inds:
                    ax.semilogy(currents['time'], abs(currents['curr'][:, i]), label=f'Mode {i}')
                    print(f'Saved mode {i} has max weight {max_weights[i]}')
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

        if verbose:
            return self.torus_reduced, sensors_measurement, currents, basis, eig_inds, weight_amplitude
        else:
            return self.torus_reduced, sensors_measurement, currents

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

    def run_reconstruction_lstsq(self, psi_at_time, ip_at_time, num_non_fil_coils, coil_curr_at_time, ip_weight, sigma, reg_factor_fil, reg_factor_wall, num_sensors=None,
                                 ip_norm_floor=1.E3):
        '''! Run the filament current reconstruction with the lstsq method.

        Solves, in the least-squares sense, for `x = [I_wall | I_fil]` (wall basis amplitudes, then filament
        currents) with the rows
          - sensors: `(Ms^T I_wall + Msc_fil^T I_fil) / sigma = Psi / sigma - (Msc_coils^T I_coils) / sigma`
          - Ip:      `(ip_weight / ip_norm) * sum(I_fil) = (ip_weight / ip_norm) * ip_at_time`
          - Tikhonov: `reg_factor_wall * I_wall = 0` and `reg_factor_fil * I_fil = 0`
        through the normal equations `A^T A x = A^T B`.

        The Ip constraint row is `(ip_weight / ip_norm) * (sum(I_fil) - ip_at_time)`, with
        `ip_norm = max(|ip_at_time|, ip_norm_floor)`. For normal currents this is the relative Ip error; the floor
        stops a divide-by-zero when Ip passes through zero (e.g. at the end of a current quench), where the row
        becomes the absolute Ip error in units of `ip_norm_floor`.

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
        @result float, norm of `B - A x` over all rows (sensor, Ip and regularization)
        @result numpy.ndarray, `A x`, all rows
        @result numpy.ndarray, `B`, all rows (weighted sensor signals, Ip row, zeros for the regularization rows)
        '''
        if not self.reduced_created_flag:
            raise ValueError("Reduced model has not been created yet. Please create or initialize a reduced model before running reconstruction.")
        if num_sensors is None:
            num_sensors = self.torus_reduced.Ms.shape[1]

        # --- Weighted sensor matrices ---
        # dividing by sigma normalizes each sensor by its noise, so no single sensor dominates the fit
        Ms_weighted = self.torus_reduced.Ms[:, :num_sensors]/sigma[:]               # (n_wall, n_sensors)
        Msc_weighted = self.torus_reduced.Msc[:, :num_sensors]/sigma[:]             # (n_coils_total, n_sensors)
        Msc_fil_weighted = Msc_weighted[num_non_fil_coils:, :]                      # (n_fil, n_sensors)
        Msc_coils_weighted = Msc_weighted[:num_non_fil_coils, :]                    # (n_shaping_coils, n_sensors)
        num_Ms = Ms_weighted.shape[0]

        # --- Ip constraint row ---
        # scales the total Ip row so it has a comparable influence on the least-squares solution as the
        # magnetic measurements, based on ip_weight and the magnitude of ip_at_time;
        # the floor avoids a divide-by-zero when Ip crosses zero (see docstring)
        ip_row_norm = max(abs(ip_at_time), ip_norm_floor)
        ip_col_wall = numpy.zeros((num_Ms, 1))
        ip_col_fil = numpy.ones((Msc_fil_weighted.shape[0], 1)) * (ip_weight / ip_row_norm)

        # appending the Ip constraint as an extra sensor column, then transposing to [sensors + Ip] x [wall | fil]
        Ms_with_ip = numpy.append(Ms_weighted, ip_col_wall, axis=1)                 # (n_wall, n_sensors + 1)
        Msc_fil_with_ip = numpy.append(Msc_fil_weighted, ip_col_fil, axis=1)        # (n_fil, n_sensors + 1)
        A_meas = numpy.vstack((Ms_with_ip, Msc_fil_with_ip)).T                      # (n_sensors + 1, n_wall + n_fil)

        # --- Right-hand side: sensors minus the non-filament coil contribution, then the Ip row ---
        B_weighted = (psi_at_time)/sigma[:] - coil_curr_at_time @ Msc_coils_weighted
        # for |Ip| >= ip_norm_floor, ip_at_time/ip_row_norm reduces to sign(ip_at_time); the magnitude is in the A-side column
        B_weighted = numpy.append(B_weighted, [ip_at_time * ip_weight/ip_row_norm])

        # --- Tikhonov regularization rows ---
        # one identity row per unknown: the first num_Ms rows are the wall amplitudes, the rest the filament currents
        # matching zero rows are added to B
        n_unknowns = Ms_with_ip.shape[0] + Msc_fil_with_ip.shape[0]
        reg_identity = numpy.eye(n_unknowns)
        A = numpy.vstack([
            A_meas,                                                                 # sensors + Ip
            reg_factor_wall * reg_identity[0:num_Ms, :],                            # wall Tikhonov
            reg_factor_fil * reg_identity[num_Ms:, :],                              # filament Tikhonov
        ])
        B = numpy.concatenate([B_weighted, numpy.zeros(A.shape[1])])

        # --- Solve the normal equations A^T A x = A^T B ---
        AtA = A.T @ A
        AtB = A.T @ B
        solution = numpy.linalg.solve(AtA, AtB)
        Ax = numpy.dot(A, solution)
        residual = numpy.sqrt(numpy.sum((B - Ax)**2))
        solution_wall = solution[:num_Ms]
        solution_fil = solution[num_Ms:]

        return solution_fil, solution_wall, residual, Ax, B

    def prepare_tsvd_laplace(self, sigma, num_non_fil_coils, num_sensors, rgrid, zgrid, num_modes, verbose=False):
        '''! Prepare matrices and projections for the svd + laplacian reconstruction method.
        These only depend on the model and sensors, not on the time slice, so compute them once and pass them
        to @ref run_reconstruction_tsvd_laplace for every time slice.

        The SVD is taken of the weighted filament-to-sensor matrix `Msc_fil / sigma` `(n_fil, n_sensors)`. The
        filament currents are then written as `I_fil = U_trun @ w`, with `U_trun = U[:, :num_modes]` and `w` the
        `num_modes` unknown TSVD weights ("TSVD space"). This keeps only the filament patterns the sensors can see.

        @param sigma numpy.ndarray, noise on each sensor, same units as the sensor signals, shape (num_sensors,)
        @param num_non_fil_coils int, number of non-filament coils before the filaments in the Icoil list
        @param num_sensors int, number of real sensors in the system
        @param rgrid numpy.ndarray, R coordinates of the filament grid `[m]`, shape (n_fil,)
        @param zgrid numpy.ndarray, Z coordinates of the filament grid `[m]`, shape (n_fil,)
        @param num_modes int, number of SVD modes to truncate to
        @param verbose bool, if True, plots singular values (default: False)
        @result numpy.ndarray, `Msc_coils_weighted`, weighted non-filament coil matrix, shape (num_non_fil_coils, num_sensors)
        @result numpy.ndarray, `Ms_weighted`, weighted wall matrix, shape (n_wall, num_sensors)
        @result numpy.ndarray, `U_trun`, truncated left singular vectors, shape (n_fil, num_modes)
        @result numpy.ndarray, `ls_mat`, sensor rows of the system matrix `[wall | TSVD]`, shape (num_sensors, n_wall + num_modes)
        @result numpy.ndarray, `ls_mat_fil`, TSVD part of ls_mat, shape (num_sensors, num_modes)
        @result numpy.ndarray, `lap_proj`, filament Laplacian acting on the TSVD weights, shape (n_fil, num_modes)
        @result numpy.ndarray, `ip_row`, unit-norm Ip constraint row, shape (1, n_wall + num_modes)
        @result int, `N`, number of filaments
        '''

        # --- Filament Laplacian (for smoothing) ---
        lap_mat, N = get_laplace_matrix(rgrid, zgrid, verbose)
        num_Ms = self.torus_reduced.Ms.shape[0]

        # --- Weighted sensor matrices ---
        # dividing by sigma normalizes each sensor by its noise, so no single sensor dominates the fit
        Msc_fil_weighted = self.torus_reduced.Msc[num_non_fil_coils:, :num_sensors]/sigma[:]     # (n_fil, n_sensors)
        Msc_coils_weighted = self.torus_reduced.Msc[:num_non_fil_coils, :num_sensors]/sigma[:]   # (n_shaping_coils, n_sensors)
        Ms_weighted = self.torus_reduced.Ms[:, :num_sensors]/sigma[:]                            # (n_wall, n_sensors)

        # --- Truncated SVD of the filament part ---
        # only the filaments are decomposed (the shaping coil currents are known, not solved for)
        U, S, Vh = numpy.linalg.svd(Msc_fil_weighted, full_matrices=False)
        U_trun = U[:, :num_modes]                                                   # (n_fil, num_modes)

        # --- Sensor rows of the system matrix: [wall | TSVD weights] ---
        ls_mat_fil = Msc_fil_weighted.T @ U_trun                                    # (n_sensors, num_modes)
        ls_mat_wall = Ms_weighted.T                                                 # (n_sensors, n_wall)
        ls_mat = numpy.hstack([ls_mat_wall, ls_mat_fil])

        # --- Laplacian projected onto the TSVD space ---
        # a regularization term in the TSVD weights that smooths in the physical filament space
        lap_proj = lap_mat @ U_trun                                                 # (n_fil, num_modes)

        # --- Ip constraint row in the TSVD space ---
        # sum(I_fil) = sum(U_trun @ w) = (column sums of U_trun) . w ; the row is scaled to unit norm
        ip_row_fil = U_trun.sum(axis=0, keepdims=True)
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
        return Msc_coils_weighted, Ms_weighted, U_trun, ls_mat, ls_mat_fil, lap_proj, ip_row, N

    def run_reconstruction_tsvd_laplace(self, psi_at_time, ip_at_time, coil_curr_at_time, Msc_coils, Ms, U_trun, ls_mat, ls_mat_fil, lap_proj, ip_row, N, sigma, lam=None, lap_lam=1e-8, reg_wall=1e-5,
                                        ip_norm_floor=1.E3):
        '''! Run the filament current reconstruction with the svd + laplacian method.
        Unknowns are `[I_wall | w]`, wall basis amplitudes and TSVD weights (see @ref prepare_tsvd_laplace);
        the filament currents are `U_trun @ w`. Rows: weighted sensors, the Ip constraint (weight `lam`),
        Tikhonov on the wall (`reg_wall`) and the filament Laplacian (`lap_lam`), solved through the normal
        equations.

        If `lam` is None it is set from the ratio of the typical weighted sensor signal to Ip,
        `lam = 100 * 10**floor(log10(mean(|B|) / ip_norm))` with `ip_norm = max(|Ip|, ip_norm_floor)`, so the
        Ip row is on the same order as the sensor rows (the factor 100 gives the single Ip row more weight
        against the many sensor rows). The floor keeps `lam` finite when Ip passes through zero (e.g. at the
        end of a current quench). `ip_error_pct` is normalized by the same `ip_norm`, so below the floor it is
        the absolute Ip error in percent of `ip_norm_floor`.

        @param psi_at_time numpy.ndarray, sensor signals at this time (e.g. flux `[Wb]`), shape (num_sensors,)
        @param ip_at_time float, total plasma current measurement at time `[A]`
        @param coil_curr_at_time numpy.ndarray, non-filament coil currents at time `[A]`, shape (num_non_fil_coils,)
        @param Msc_coils numpy.ndarray, weighted non-filament coil matrix from @ref prepare_tsvd_laplace
        @param Ms numpy.ndarray, weighted wall matrix from @ref prepare_tsvd_laplace, shape (n_wall, num_sensors)
        @param U_trun numpy.ndarray, truncated singular vectors from @ref prepare_tsvd_laplace, shape (n_fil, num_modes)
        @param ls_mat numpy.ndarray, sensor rows of the system matrix from @ref prepare_tsvd_laplace
        @param ls_mat_fil numpy.ndarray, TSVD part of ls_mat from @ref prepare_tsvd_laplace
        @param lap_proj numpy.ndarray, Laplacian in TSVD space from @ref prepare_tsvd_laplace
        @param ip_row numpy.ndarray, unit-norm Ip constraint row from @ref prepare_tsvd_laplace
        @param N int, number of filaments
        @param sigma numpy.ndarray, noise on each sensor, same units as psi_at_time, shape (num_sensors,)
        @param lam float, Ip constraint weight; None sets it automatically (see above) (default: None)
        @param lap_lam float, Laplacian regularization weight (default: 1e-8)
        @param reg_wall float, Tikhonov regularization weight on the wall amplitudes (default: 1e-5)
        @param ip_norm_floor float, smallest current used for the automatic `lam` and for `ip_error_pct` `[A]` (default: 1.E3)
        @result numpy.ndarray, filament currents `[A]`, shape (n_fil,)
        @result numpy.ndarray, wall basis amplitudes, shape (n_wall,)
        @result numpy.ndarray, `A x` for the sensor rows only, shape (num_sensors,)
        @result dict, diagnostics: `lam`, `lap_lam`, `ip_reconstructed` `[A]`, `ip_actual` `[A]`, `ip_error_pct`,
            `fit_residual` (norm of the weighted sensor misfit), `fit_residual_nonorm` (the misfit vector),
            `curr_weights` (full solution `[I_wall | w]`), `ip_row`
        '''
        num_Ms = Ms.shape[0]

        # --- Right-hand side: sensors minus the non-filament coil contribution ---
        B_weighted = (psi_at_time)/sigma[:] - coil_curr_at_time @ Msc_coils

        # --- Ip constraint weight ---
        if lam is None:
            lam = _auto_ip_weight(B_weighted, ip_at_time, ip_norm_floor)

        # ip_row has unit norm, so the B-side Ip value is scaled by the same norm
        ip_scale = numpy.linalg.norm(U_trun.sum(axis=0))

        # --- System matrix: measurements, Ip row, and the regularization rows ---
        reg_mat = numpy.vstack([
            ls_mat,                                                                                    # measurements
            lam * ip_row,                                                                              # Ip constraint
            numpy.hstack([reg_wall * numpy.eye(num_Ms), numpy.zeros((num_Ms, ls_mat_fil.shape[1]))]),  # wall Tikhonov
            numpy.hstack([numpy.zeros((N, num_Ms)), lap_lam * lap_proj])                               # fil Laplacian
        ])

        # --- Right-hand side for all rows: measurements, Ip constraint, zeros for the regularization rows ---
        psi_reg = numpy.concatenate([
            B_weighted,
            numpy.array([lam * ip_at_time / ip_scale]),
            numpy.zeros(num_Ms),
            numpy.zeros(N)
        ])

        # --- Solve the normal equations A^T A x = A^T B ---
        AtA = reg_mat.T @ reg_mat
        AtB = reg_mat.T @ psi_reg
        curr_weights = numpy.linalg.solve(AtA, AtB)

        # --- Map the TSVD weights back to filament currents ---
        curr_expand = U_trun @ curr_weights[num_Ms:]  # only the filament part contributes to the expansion
        wall_expand = curr_weights[:num_Ms]  # wall coefficients

        # --- Diagnostics ---
        Ax = ls_mat @ curr_weights
        ip_reconstructed = curr_expand.sum()
        ip_actual        = ip_at_time
        ip_error_pct     = 100 * numpy.abs(ip_reconstructed - ip_actual) / max(abs(ip_actual), ip_norm_floor)  # floor: see docstring
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

    def run_reconstruction_laplace(self, psi_at_time, ip_at_time, coil_curr_at_time, sigma,
                                    num_non_fil_coils, num_real_sensors, rgrid, zgrid,
                                    lam=None, lap_lam=1e-6, reg_wall=1e-5, sigma_r=1e-1, sigma_z=1e-1, gaussian=False, verbose=False,
                                    ip_norm_floor=1.E3):
        '''! Run filament current reconstruction using direct Laplacian regularization.
        Solves directly in the physical filament space (no SVD projection), applying
        Laplacian smoothing to filaments and Tikhonov regularization to wall currents.
        Unknowns are `[I_wall | I_fil]`; rows: weighted sensors, the Ip constraint (weight `lam`), Tikhonov on
        the wall (`reg_wall`) and the filament Laplacian (`lap_lam`), solved through the normal equations.
        If `lam` is None it is set as in @ref run_reconstruction_tsvd_laplace, including the `ip_norm_floor`
        used there for the automatic `lam` and for `ip_error_pct`.

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
        @result numpy.ndarray, `A x` for the sensor rows only, shape (num_real_sensors,)
        @result dict, diagnostics: `lam`, `lap_lam`, `reg_wall`, `ip_reconstructed` `[A]`, `ip_actual` `[A]`,
            `ip_error_pct`, `fit_residual` (norm of the weighted sensor misfit), `fit_residual_nonorm` (the misfit
            vector), `solution` (full solution `[I_wall | I_fil]`)
        '''

        # --- Weighted sensor matrices ---
        # dividing by sigma normalizes each sensor by its noise, so no single sensor dominates the fit
        Ms_weighted = self.torus_reduced.Ms[:, :num_real_sensors] / sigma[:]        # (n_wall, n_sensors)
        Msc_weighted = self.torus_reduced.Msc[:, :num_real_sensors] / sigma[:]      # (n_coils_total, n_sensors)
        Msc_coils_weighted = Msc_weighted[:num_non_fil_coils, :num_real_sensors]    # (n_shaping_coils, n_sensors)
        Msc_fil_weighted = Msc_weighted[num_non_fil_coils:, :num_real_sensors]      # (n_fil, n_sensors)
        n_fil = Msc_fil_weighted.shape[0]
        num_Ms = Ms_weighted.shape[0]

        # --- Right-hand side: sensors minus the non-filament coil contribution ---
        B_weighted = (psi_at_time)/sigma[:] - coil_curr_at_time @ Msc_coils_weighted

        # --- Ip constraint weight and row ---
        if lam is None:
            lam = _auto_ip_weight(B_weighted, ip_at_time, ip_norm_floor)

        # each filament contributes equally to Ip, so the row is all-ones for filaments, zeros for wall
        ip_row_fil = numpy.ones((1, n_fil))
        ip_row_norm = numpy.linalg.norm(ip_row_fil)
        ip_row = numpy.hstack([numpy.zeros((1, num_Ms)), ip_row_fil / ip_row_norm])  # (1, n_wall + n_fil)

        # --- Filament Laplacian (for smoothing the filament solution in space) ---
        if gaussian:
            lap_mat, N = get_gaussian_lap_mat(rgrid, zgrid, sigma_r, sigma_z)        # (n_fil, n_fil)
        else:
            lap_mat, N = get_laplace_matrix(rgrid, zgrid, verbose=verbose)          # (n_fil, n_fil)

        # --- System matrix: measurements, Ip row, and the regularization rows ---
        # measurements corresponding to sensors: [wall | filament] columns
        meas_block = numpy.hstack([Ms_weighted.T, Msc_fil_weighted.T])              # (n_sensors, n_wall + n_fil)

        # Ip constraint row
        ip_block = lam * ip_row                                                     # (1, n_wall + n_fil)

        # Tikhonov on wall
        wall_reg = reg_wall * numpy.hstack([
            numpy.eye(num_Ms),
            numpy.zeros((num_Ms, n_fil))
        ])                                                                          # (n_wall, n_wall + n_fil)

        # Laplacian on filaments
        fil_lap = lap_lam * numpy.hstack([
            numpy.zeros((N, num_Ms)),
            lap_mat
        ])                                                                          # (n_fil, n_wall + n_fil)
        A = numpy.vstack([meas_block, ip_block, wall_reg, fil_lap])

        # --- Right-hand side for all rows: measurements, Ip constraint, zeros for the regularization rows ---
        psi_reg = numpy.concatenate([
            B_weighted,
            numpy.array([lam * ip_at_time / ip_row_norm]),
            numpy.zeros(num_Ms),
            numpy.zeros(N)
        ])

        # --- Solve the normal equations A^T A x = A^T B ---
        AtA = A.T @ A
        AtB = A.T @ psi_reg
        solution = numpy.linalg.solve(AtA, AtB)
        # extracting wall and filament currents from the solution vector
        curr_wall = solution[:num_Ms]
        curr_fil = solution[num_Ms:]

        # --- Diagnostics ---
        Ax = meas_block @ solution
        ip_reconstructed  = curr_fil.sum()
        ip_error_pct      = 100 * numpy.abs(ip_reconstructed - ip_at_time) / max(abs(ip_at_time), ip_norm_floor)  # floor: see docstring
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

    # ==========================================================
    # Post-processing and visualization methods
    # ==========================================================

    def get_wall_psi_tidx(self, num_sensors, solution_wall_at_time):
        '''! Compute the wall contribution to the flux at a given time index.
        Uses the fitted wall basis amplitudes and the wall-to-sensor matrix `Ms` to evaluate the wall flux on
        the extra sensors after the first num_sensors (the in-vessel flux probes, e.g. placed with
        @ref get_inside_limiter_pts). The division by `2 pi` converts flux `[Wb]` to TokaMaker psi `[Wb/rad]`.
        Before @ref post_process_tidx the result must be put on every TokaMaker node (the reconstruct example
        fills the nodes without a probe by interpolation).

        @param num_sensors int, number of real sensors to skip (the probes start at index num_sensors)
        @param solution_wall_at_time numpy.ndarray, wall basis amplitudes at the given time index, shape (n_wall,)
        @result numpy.ndarray, wall psi on the flux probes `[Wb/rad]`, shape (n_probes,)'''
        wall_psi_probes = self.torus_reduced.Ms[:, num_sensors:]/(2*numpy.pi)
        wall_psi_at_time = solution_wall_at_time @ wall_psi_probes
        return wall_psi_at_time

    def compute_psi_basis(self, rgrid, zgrid, disk_radius=0.3):
        '''! Plasma psi per unit current of each filament, for reuse in post_process_tidx (psi_basis).

        Each filament is represented by a TokaMaker init_psi disk of radius disk_radius centred on (r, z), exactly
        as in post_process_tidx. init_psi is linear in Ip (checked to 2e-16 relative), so the plasma psi of any
        filament current vector is psi_basis @ I_fil. This avoids repeating the init_psi loop for every time slice,
        which matters for grids with thousands of filaments (e.g. the M3D-C1 projection).

        @param rgrid numpy.ndarray, R coordinates of the filament mesh `[m]`, shape (n_fil,)
        @param zgrid numpy.ndarray, Z coordinates of the filament mesh `[m]`, shape (n_fil,)
        @param disk_radius float, radius of the init_psi disk per filament `[m]` (default: 0.3)
        @result numpy.ndarray, psi per amp on the TokaMaker nodes `[Wb/rad per A]`, shape (n_nodes, n_fil)

        Example
        -------
        >>> psi_basis = jamfit.compute_psi_basis(rgrid, zgrid, disk_radius=0.2)
        >>> out = jamfit.post_process_tidx(fil, wall_psi, coils, rgrid, zgrid, disk_radius=0.2, psi_basis=psi_basis)
        '''
        psi_basis = []
        for r, z in zip(rgrid, zgrid):
            self.mygs.set_coil_currents()
            self.mygs.set_targets(Ip=1.0)
            self.mygs.init_psi(r, z, disk_radius, 1.0, 0.0)
            psi_basis.append(self.mygs.get_psi(False))
        return numpy.array(psi_basis).T

    def post_process_tidx(self, filaments_at_time, wall_psi, coil_curr_dict, rgrid, zgrid, verbose=False, include_negative_encl=True,
                          disk_radius=0.3, psi_basis=None, touch_tol=0.0005):
        '''! Post-process the results at a given time index to compute plasma parameters and visualize.

        The plasma psi is the sum over positive filaments of a TokaMaker init_psi disk of radius disk_radius centred
        on each filament. Negative filaments are left out on purpose: TokaMaker init_psi cannot represent a negative
        plasma current, and negative toroidal current is unphysical here (a well-tuned reconstruction should have
        essentially none). Where the total current is needed (li), it is the current enclosed by the LCFS, see
        include_negative_encl; the current centroid uses all filaments, negatives included. Smaller disks are more accurate
        but need a fine enough filament grid for trace_surf to find an LCFS.

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
        @param disk_radius float, radius of the init_psi disk per filament `[m]` (default: 0.3, the old hard-coded value)
        @param psi_basis numpy.ndarray or None, psi per amp of each filament, shape (n_nodes, n_fil), from
            compute_psi_basis with the same disk_radius; if None the init_psi loop is run here (default: None)
        @param touch_tol float, largest LCFS-limiter distance that counts as touching, passed to
            @ref find_limiting_point `[m]` (default: 0.0005)
        @result numpy.ndarray or None, LCFS points `[m]`, shape (n_pts, 2) (traced at psi_N = 1, or 0.99 if that fails; None if both fail)
        @result tuple or None, `(R, Z)` `[m]` of the point where the LCFS touches the limiter (see
            @ref find_limiting_point); None if the LCFS does not touch the limiter or there is no LCFS
        @result float or None, psi at the traced LCFS `[Wb/rad]`
        @result numpy.ndarray, total psi (filaments + coils + wall) on the TokaMaker nodes `[Wb/rad]`, shape (n_nodes,)
        @result numpy.ndarray or None, safety factor q on TokaMaker's default psi_N sampling
        @result float, q95 (nan if there is no LCFS)
        @result float, internal inductance li, 'std' normalisation (nan if there is no LCFS or Ip_encl <= 0)
        @result tuple, current centroid `(R_c, Z_c)` `[m]` from all filaments (see @ref calc_current_centroid)
        @result tuple or float, geometric LCFS centroid `(R_c, Z_c)` `[m]` (nan if there is no LCFS)
        @result float, area enclosed by the LCFS `[m^2]` (nan if there is no LCFS)'''
        # reset mygs in case there are left over coil currents from previous runs
        self.mygs.set_coil_currents()
        limiter = self.mygs.lim_contour
        current_cent = calc_current_centroid(rgrid, zgrid, filaments_at_time)

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
        # Throw-away init_psi at the current centroid. If the last init_psi before set_psi has its disk (mostly)
        # outside the plasma region, get_q returns NaN at every psi while trace_surf and li look fine. set_psi below
        # overwrites psi, so this only resets TokaMaker's internal state.
        self.mygs.set_coil_currents()
        self.mygs.set_targets(Ip=1.0)
        self.mygs.init_psi(current_cent[0], current_cent[1], disk_radius, 1.0, 0.0)

        # --- Coil psi (vacuum solve with the coil currents) ---
        self.mygs.set_coil_currents(coil_curr_dict)
        psi_coils_obj = self.mygs.vac_solve()
        psi_coils = psi_coils_obj.get_psi(normalized=False)

        # --- Total psi and LCFS ---
        total_psi = psi_total_fil + psi_coils + wall_psi
        self.mygs.set_psi(total_psi, update_bounds=True)

        # trace the LCFS at psi_N = 1; if that fails, fall back to psi_N = 0.99
        lcfs_points = self.mygs.trace_surf(1)
        if lcfs_points is not None:
            psiatlcfs = self.mygs.psinorm_to_absolute(1)
        else:
            lcfs_points = self.mygs.trace_surf(0.99)
            psiatlcfs = self.mygs.psinorm_to_absolute(0.99)

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
    @result float, Ip constraint weight `lam`

    Example
    -------
    >>> _auto_ip_weight(numpy.array([2.0, 4.0]), 1.E6)    # mean |B| = 3, ratio 3e-6 -> 100 * 1e-6 = 1e-4
    '''
    ip_norm = max(abs(ip_at_time), ip_norm_floor)
    magnitude_diff = math.floor(math.log10(numpy.mean(numpy.abs(B_weighted)) / ip_norm))
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

def get_bbox_grid_pts(meshfile_tokamaker, myOFT, n_r=10, n_z=10, verbose=False):
    '''! Build a uniform rectangular grid spanning the bounding box of the TokaMaker
    evaluation grid, for later interpolation of measured wall values onto the grid.
    Creates its own TokaMaker instance (order 2, `F0 = 1`) on myOFT just to read the mesh nodes.

    @param meshfile_tokamaker str, path to the TokaMaker mesh file
    @param myOFT OFT_env, the Open FUSION Toolkit environment instance
    @param n_r int, number of points along R (default: 10)
    @param n_z int, number of points along Z (default: 10)
    @param verbose bool, if True, plots the bounding box grid and limiter contour (default: False)
    @result numpy.ndarray, bounding-box grid points `[m]`, shape (n_r*n_z, 2) with columns [R, Z]
    @result numpy.ndarray, R meshgrid `[m]`, shape (n_z, n_r)
    @result numpy.ndarray, Z meshgrid `[m]`, shape (n_z, n_r)
    @result numpy.ndarray, copy of all TokaMaker node coordinates `mygs.r` `[m]`, shape (n_nodes, 3)
    '''
    mygs = TokaMaker(myOFT)
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(meshfile_tokamaker)
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.setup(order=2, F0=1 * 1)

    r_pts_grid = mygs.r[:, 0]
    z_pts_grid = mygs.r[:, 1]

    r_min, r_max = r_pts_grid.min(), r_pts_grid.max()
    z_min, z_max = z_pts_grid.min(), z_pts_grid.max()

    r_vals = numpy.linspace(r_min, r_max, n_r)
    z_vals = numpy.linspace(z_min, z_max, n_z)
    R, Z = numpy.meshgrid(r_vals, z_vals, indexing='xy')

    bbox_pts = numpy.column_stack([R.ravel(), Z.ravel()])

    if verbose:
        fig, ax = plt.subplots()
        ax.scatter(r_pts_grid, z_pts_grid, c='lightgray', s=5, alpha=0.5, label='TokaMaker grid')
        ax.scatter(bbox_pts[:, 0], bbox_pts[:, 1], c='tab:blue', s=15, label='Bbox grid pts')
        r_pts_lim = mygs.lim_contour[:, 0]
        z_pts_lim = mygs.lim_contour[:, 1]
        ax.plot(r_pts_lim, z_pts_lim, 'k-', label='Limiter Contour')
        ax.set_xlabel('R')
        ax.set_ylabel('Z')
        ax.set_title(f'Bounding Box Grid ({n_r} x {n_z} points)')
        ax.legend()
        plt.tight_layout()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    return bbox_pts, R, Z, mygs.r.copy()


def get_inside_limiter_pts(meshfile_tokamaker, myOFT, stride=1, verbose=False):
    '''! Get the TokaMaker mesh nodes inside the limiter contour, to place flux probes for the wall psi.
    Every stride-th inside node gets a probe; see @ref get_wall_psi_tidx. Creates its own TokaMaker instance
    (order 2, `F0 = 1`) on myOFT just to read the mesh nodes and limiter.

    @param meshfile_tokamaker str, path to the TokaMaker mesh file
    @param myOFT OFT_env, the Open FUSION Toolkit environment instance
    @param stride int, subsampling stride for the points inside the limiter (default: 1)
    @param verbose bool, if True, plots the grid points and limiter contour (default: False)
    @result numpy.ndarray, probe points `[m]`, shape (n_probes, 2) with columns [R, Z]
    @result numpy.ndarray, bool mask over all nodes, True where a probe is placed, shape (n_nodes,)
    @result numpy.ndarray, bool mask over all nodes, True inside the limiter, shape (n_nodes,)
    @result numpy.ndarray, copy of all TokaMaker node coordinates `mygs.r` `[m]`, shape (n_nodes, 3)
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

    return inside_lim_pts, sparse_mask, inside_mask, mygs.r.copy()  # copy: mygs.r is a view on TokaMaker memory, freed with this local mygs


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
