#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Solvers and helper functions for TokaMaker bootstrap current functionality

@authors Daniel Burgess
@date January 2026
@ingroup doxy_oft_python
'''
import numpy
from ._interface import *
from OpenFUSIONToolkit.TokaMaker.util import get_jphi_from_GS

def parameterize_edge_jBS(psi, amp, center, width, offset, sk, y_sep=0.0, blend_width=0.03, tail_alpha=1.5):
    r'''! Generates a parameterized bootstrap current profile with a tunable 'concave' fall-off
    @param psi Normalized psi profile
    @param amp Amplitude of edge skewed Gaussian spike
    @param center Center location of edge skewed Gaussian spike
    @param width Width of edge skewed Gaussian spike
    @param offset Height of core flat profile
    @param sk Skew parameter of Gaussian
    @param y_sep Separatrix value of parameterized profile
    @param blend_width Normalized psi profile
    @param tail_alpha Controls the shape of the right-side fall-off.
                    1.0 = Standard Cosine (smooth, rounded).
                    >1.0 = Sharper peak, straighter/concave fall-off (closer to typical bootstrap).
    '''
    from scipy.optimize import root_scalar
    from scipy.stats import skewnorm

    def generate_baseline_prof(psi, amp, center, width, offset, sk, y_sep, blend_width, tail_alpha):
        # --- 1. Find Exact Peak of Underlying Skew Gaussian ---
        def raw_shape(x):
            return skewnorm.pdf(x, sk, loc=center, scale=width)

        # Fine grid scan for precision
        x_grid = numpy.linspace(max(0, center - 3*width), min(1, center + 3*width), 10000)
        y_grid = raw_shape(x_grid)
        idx_peak = numpy.argmax(y_grid)
        val_peak_raw = y_grid[idx_peak]
        x_peak = x_grid[idx_peak]

        # --- 2. Calculate Internal Amplitude & Left Side ---

        # Smoothing factor for the blend
        k_smooth = (amp / width) * (blend_width / 4.0) 
        if k_smooth < 1e-5: k_smooth = 1e-5

        # Case A: Standard Spike (Amp > Offset)
        if amp > offset:
            # Solve for internal_amp to compensate for SoftMax boost
            diff = amp - offset
            if diff > 1e-10:
                argument = 1.0 - numpy.exp((offset - amp) / k_smooth)
                if argument <= 0: argument = 1e-16
                internal_amp = amp + k_smooth * numpy.log(argument)
            else:
                internal_amp = amp 

            spike_profile = raw_shape(psi) / val_peak_raw * internal_amp
            profile_left = k_smooth * numpy.logaddexp(offset / k_smooth, spike_profile / k_smooth)
            stitch_height = amp

        # Case B: Dominant Offset
        else:
            profile_left = numpy.full_like(psi, offset)
            stitch_height = offset

        # --- 3. Construct Right Side (Powered Cosine) ---
        # Model: y = stitch_height * cos(omega * (x - x_peak)) ^ tail_alpha

        u = (psi - x_peak) / (1.0 - x_peak)
        dist_to_edge = 1.0 - x_peak
        if dist_to_edge < 1e-5: dist_to_edge = 1e-5

        # Determine omega to hit y_sep exactly
        # y_sep = H * cos(w * L)^alpha  =>  (y_sep/H)^(1/alpha) = cos(w * L)

        # Safety: Base for power must be positive
        # If y_sep > stitch_height, we use Cosh instead of Cos

        if y_sep < stitch_height:
            # Cosine Decay
            # Calculate target cosine value (inverse power)
            target_cos = (y_sep / stitch_height) ** (1.0 / tail_alpha)
            # Clip to valid domain for arccos [-1, 1]
            target_cos = numpy.clip(target_cos, -1.0, 1.0)

            omega = numpy.arccos(target_cos) / dist_to_edge

            # Calculate profile
            # We must protect the cosine argument so it doesn't go negative before power
            # (though strictly defined, it stays in [0, pi/2] for this range)
            arg = omega * (psi - x_peak)
            base_cos = numpy.cos(arg)
            # Ensure positive base for float power
            base_cos = numpy.maximum(base_cos, 0) 

            profile_right = stitch_height * (base_cos ** tail_alpha)

        else:
            # Cosh Rise (Rare)
            target_cosh = (y_sep / stitch_height) ** (1.0 / tail_alpha)
            omega = numpy.arccosh(target_cosh) / dist_to_edge
            profile_right = stitch_height * (numpy.cosh(omega * (psi - x_peak)) ** tail_alpha)

        # --- 4. Stitching ---
        temp_profile = numpy.where(psi <= x_peak, profile_left, profile_right)

        if y_sep >= 0:
            temp_profile = numpy.maximum(temp_profile, 0)
        return temp_profile
    
    def objective_function(alpha, psi, amp, center, width, offset, sk, y_sep, blend_width, tail_alpha):
        r'''! Compute difference between integrated a*j_tor+j_spike profile and Ip_target

        @param alpha Scaling factor to solve for
        @param jtor_prof Input j_inductive profile
        @param spike_profile Isolated j_bootstrap spike (a Gaussian), 0.0 everywhere else
        @param my_psi_N Local psi_N grid
        @param my_Ip_target Ip target
        '''
        temp_profile = generate_baseline_prof(psi, amp, center, width, alpha*offset, sk, y_sep, blend_width, tail_alpha)
        return temp_profile[0] - offset

    if offset == 0.0:
        final_profile = generate_baseline_prof(psi, amp, center, width, -1e10, sk, y_sep, blend_width, tail_alpha)
    else:
        # Find scalar "alpha" that solves 
        result = root_scalar(objective_function,
                                args=(psi, amp, center, width, offset, sk, y_sep, blend_width, tail_alpha),
                                bracket=[-1e+10, 10*amp],
                                method='brentq',
                                rtol=1e-6)

        # Extract the solution
        a_optimal = result.root

        final_profile = generate_baseline_prof(psi, amp, center, width, a_optimal*offset, sk, y_sep, blend_width, tail_alpha)

    # find correct offset to match user input value
    return final_profile

def analyze_bootstrap_edge_spike(psi_N, j_bootstrap, parameterize_jBS, diagnostic_plots=False):
    r'''! Analyze bootstrap edge spike location, width, and height

    @param psi_N Normalized psi profile
    @param j_bootstrap Bootstrap current profile
    Returns:
    dict: Dictionary containing pedestal properties and spike model
    '''

    from scipy.signal import find_peaks, peak_widths
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt

    # Focus on the edge region (psi_N > 0.7)
    edge_mask = psi_N >= 0.7
    psi_edge = psi_N[edge_mask]
    j_edge = j_bootstrap[edge_mask]

    # Find peak in the edge region
    peaks, properties = find_peaks(j_edge, height=0.)#0.5*numpy.max(j_edge))

    if len(peaks) == 0:
        print("No clear peak found in the edge region")
        return None

    # Choose peak closest to psi_N = 1 if multiple peaks exist
    peak_idx = peaks[numpy.argmax(psi_edge[peaks])]
    peak_psi = psi_edge[peak_idx]
    peak_height = j_edge[peak_idx]

    # Calculate initial FWHM (full width at half maximum)
    widths = peak_widths(j_edge, [peak_idx], rel_height=0.5)
    left_idx, right_idx = int(widths[2][0]), int(widths[3][0])

    # Convert to psi_N coordinates
    fwhm = psi_edge[right_idx] - psi_edge[left_idx]

    # Find minimum of j_BS current in core/mantle region and use as values for flat core profile
    masked_j_BS = j_bootstrap[(psi_N > 0.5) & (psi_N < peak_psi)] # minimum between psi_N = 0.5 and peak of j_BS
    masked_psi_N = psi_N[(psi_N > 0.5) & (psi_N < peak_psi)]
    lmin_j_BS = min(masked_j_BS)
    lmin_arg = numpy.argmin(masked_j_BS)
    fit_mask = (psi_N >= masked_psi_N[lmin_arg])
    
    # Splice j_BS profile edge spike onto flat profile with value set to min(j_BS)
    masked_spike = numpy.ones_like(psi_N)*lmin_j_BS
    masked_spike[fit_mask] = j_bootstrap[fit_mask]
    
    # Define the blend window to smooth out junction of flat core profile and edge j_BS spike
    # Set the total width to be half the distance between splice location and peak
    jBS_min_loc = masked_psi_N[lmin_arg]
    dist_to_peak = peak_psi - jBS_min_loc
    blend_width = min(0.5 * dist_to_peak, 0.2)

    x_start = jBS_min_loc - (blend_width / 2.0)
    x_end   = jBS_min_loc + (blend_width / 2.0)

    # Find indices for the start and end of the window
    idx_start = (numpy.abs(psi_N - x_start)).argmin()
    idx_end   = (numpy.abs(psi_N - x_end)).argmin()

    # Validation to ensure indices are ordered and within bounds
    idx_start = max(0, idx_start)
    idx_end   = min(len(psi_N) - 2, idx_end) # -2 to safe-guard derivative calc later

    # Determine Boundary Conditions
    # Left Boundary (Start): strictly flat, so slope is 0
    y_start  = masked_spike[idx_start]
    dy_start = 0.0  

    # Right Boundary (End): located on the peaked profile
    # Estimate the slope (dy/dx) using a central difference at x_end
    y_end  = masked_spike[idx_end]
    dy_end = (masked_spike[idx_end + 1] - masked_spike[idx_end - 1]) / \
                (psi_N[idx_end + 1] - psi_N[idx_end - 1])

    # Generate Cubic Hermite Spline
    # Extract the x-values in the window to interpolate over
    x_window = psi_N[idx_start : idx_end + 1]

    # Normalized coordinate t (goes from 0 to 1 across the window)
    t = (x_window - x_window[0]) / (x_window[-1] - x_window[0])

    # Hermite Basis Functions
    h00 = 2*t**3 - 3*t**2 + 1           # Weight for Start Value
    h10 = t**3 - 2*t**2 + t             # Weight for Start Slope
    h01 = -2*t**3 + 3*t**2              # Weight for End Value
    h11 = t**3 - t**2                   # Weight for End Slope

    # The physical length of the window (needed to scale the derivative terms)
    dx_window = x_window[-1] - x_window[0]

    # Calculate the smooth patch
    y_patch = (h00 * y_start) + \
                (h10 * dx_window * dy_start) + \
                (h01 * y_end) + \
                (h11 * dx_window * dy_end)

    # Apply the patch to the main masked j_BS profile
    masked_spike[idx_start : idx_end + 1] = y_patch
    
    if parameterize_jBS:
        # Attempt to fit bootstrap parameterization to calculated profile (deprecated)
        # Initial parameter guess
        sigma_init = fwhm/2.355  # Convert FWHM to sigma
        p0 = [peak_height, peak_psi, sigma_init, lmin_j_BS, 1.0, j_bootstrap[-1], 0.05]

        print(f'Fitted peak height = {peak_height:.3f}')
        print(f'Fitted peak center in psi_N = {peak_psi:.3f}')
        print(f'Fitted peak sigma = {sigma_init:.3f}')
        print(f'Measured j_BS local min. = {lmin_j_BS:.3f}')

        lower_bounds = [0.99*peak_height, # fix to measured spike height
                        0.8*peak_psi, # don't fix to measured spike location
                        0.0, # width
                        0.99*lmin_j_BS,
                        -50,
                        0.0,
                        0.001]
        upper_bounds = [1.01*peak_height, # fix to measured spike height
                        1.2*peak_psi, # don't fix to measured spike location
                        0.33, # no spikes wider than 0.33 units of psi_N 
                        1.01*lmin_j_BS,
                        50, # sk
                        2*j_bootstrap[-1],
                        0.2]

        # Perform the fit
        popt, pcov = curve_fit(parameterize_edge_jBS, psi_N[fit_mask], j_bootstrap[fit_mask], p0=p0,
                                bounds=(lower_bounds, upper_bounds),
                                maxfev=10000)

        amp, center, width, offset, sk, y_sep, blend_width = popt

        spike_only = parameterize_edge_jBS(psi_N, amp, center, width, offset, sk, y_sep, blend_width)
    else:
        spike_only = numpy.zeros_like(psi_N)
        width = fwhm/2.355
        offset = (jBS_min_loc,lmin_j_BS)
        popt = None

    if diagnostic_plots:
        plt.figure()
        plt.plot(psi_N[fit_mask],j_bootstrap[fit_mask]/1e6,label='Input j_BS')
        plt.plot(psi_N,spike_only/1e6,label='Fitted parameterization')
        plt.plot(psi_N,masked_spike/1e6,label='Isolated')
        plt.xlabel(r'$\hat \psi$')
        plt.ylabel(r'$j_\phi$ [MA]')
        plt.xlim(0.5,1.02)
        plt.legend(loc='best')
        plt.show()
    
    results = {
        'sigma': width,                 # Gaussian width (sigma)
        'background': offset,           # Tuple of masked j_BS blend location OR parameterized offset level
        'gaussian_params': popt,        # Raw parameters [amp, center, width, offset]
        'parameterized_spike': spike_only,    # Array of spike component values
        'masked_spike': masked_spike
    }

    return results

def solve_jphi(mygs,ffp_prof,pp_prof,Ip_target,pax_target):
    r'''! Take input P' and j_phi profiles and solve equilibrium
    @param ffp_prof j_phi profile, with 'type' set to 'jphi-linterp'
    @param pp_prof P' profile
    @param Ip_target Target Plasma Current [A]
    @param pax_target Target pressure-on-axis [Pa]
    '''

    ffp_prof['type'] = 'jphi-linterp'

    mygs.set_targets(Ip=Ip_target, pax=pax_target)
    mygs.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof)

    # Solve Grad-Shafranov
    mygs.solve()
    
def find_optimal_scale(mygs, psi_N, pressure, ffp_prof, pp_prof, j_inductive,  
                            Ip_target, psi_pad, spike_prof=None, find_j0=True, scale_j0=1.0, 
                            tolerance=0.01, max_iter=5, diagnostic_plots=False):
    '''
    Optimizes scale_j0 OR Ip_target to match either input and output current densities at the axis (index 0)
    or the input and output Ip.
    Uses the Secant Method to minimize equilibrium solves calls.

    @param psi_N Normalized poloidal flux profile
    @param pressure Plasma pressure profile [Pa]
    @param ffp_prof FF' profile
    @param pp_prof P' profile
    @param j_inductive Inductive OR total toroidal current profile to be rescaled [A/m^2]
    @param Ip_target Target Plasma Current [A]
    @param psi_pad Padding for flux surface calculations
    @param spike_prof Optional bootstrap current profile, will NOT be rescaled. Can be set to None [A/m^2]
    @param find_j0 Essential logic, will switch between either rescaling core j_phi profile or Ip_target
    @param scale_j0 Only used if find_j0 = False and Ip_target is being rescaled. Useful if find_j0 = True has 
    already been run to find the correct rescaling factor
    @param tolerance Relative error tolerance of secant method search. Possible defaults: 1% for find_j0 = True,
    0.1% for find_j0 = False
    @param max_iter Maximum iterations of secant method
    @param diagnostic_plots Plot input and output j_phi profiles to check alignment
    '''
    import matplotlib.pyplot as plt

    n_psi = len(psi_N)
    
    if spike_prof is None:
        spike_prof = numpy.zeros_like(j_inductive)
    
    # We want: Input_J0 - Output_J0 ~ 0
    def get_j0_error(scale_val, n):
        print(f"\n--- Checking scale_j0 = {scale_val:.4f} ---")
        
        matched_input_jphi = scale_val*j_inductive + spike_prof
        ffp_prof['type'] = 'jphi-linterp'
        ffp_prof['y'] = matched_input_jphi

        pax_target = pressure[0]

        solve_jphi(mygs,ffp_prof,pp_prof,Ip_target,pax_target)
        
        # Check Convergence
        _, f, fp, _, pp = mygs.get_profiles(npsi=n_psi, psi_pad=psi_pad)
        _, _, ravgs, _, _, _ = mygs.get_q(npsi=n_psi, psi_pad=psi_pad)

        tmp_jphi = get_jphi_from_GS(f*fp, pp, ravgs[0], ravgs[1])

        if diagnostic_plots:
            plt.figure()
            plt.plot(psi_N, matched_input_jphi/1e6, linestyle='--', label=r'Input $j_\phi$')
            plt.plot(psi_N, tmp_jphi/1e6, label=r'Output $j_\phi$')
            plt.title(f'Iteration {n}')#', Ip error: {Ip_err:.3f} %')
            plt.legend()
            plt.xlabel(r'$\hat \psi')
            plt.xlabel(r'$j_\phi$ [MA]')
            plt.grid(ls=':')
            plt.show()
            
        # Input j_0 at this scale
        input_j0 = scale_val * j_inductive[0] + spike_prof[0]
        
        # Output j_0 from solver
        output_j0 = tmp_jphi[0]
        
        # Calculate residual (difference) and relative error
        diff = input_j0 - output_j0
        rel_err = abs(diff) / output_j0
        
        print(f"   Input j_0:  {input_j0:.4e}")
        print(f"   Output j_0: {output_j0:.4e}")
        print(f"   Mismatch:  {rel_err*100:.3f}%")
        
        return diff, rel_err, tmp_jphi
    
    # We want: Input_Ip - Output_Ip ~ 0
    def get_Ip_error(scale_val, scale_j0, n):
        print(f"\n--- Checking scale_Ip = {scale_val:.4f} ---")
        
        ffp_prof['type'] = 'jphi-linterp'
        ffp_prof['y'] = scale_j0*j_inductive + spike_prof

        scaled_Ip_target = Ip_target*scale_val
        pax_target = pressure[0]

        solve_jphi(mygs,ffp_prof,pp_prof,scaled_Ip_target,pax_target)
    
        eq_stats = mygs.get_stats(lcfs_pad=psi_pad)
        output_Ip = eq_stats['Ip']
        
        # Calculate residual (difference) and relative error
        diff = output_Ip - Ip_target
        rel_err = abs(diff) / Ip_target
        
        print(f"   Input Ip target:  {Ip_target/1e6:.4e}")
        print(f"   Trial Ip target:  {Ip_target*scale_val/1e6:.4e}")
        print(f"   Output Ip: {output_Ip/1e6:.4e}")
        print(f"   Mismatch:  {rel_err*100:.4f}%")
        
        return diff, rel_err, None

    # --- Step 1: Initial Guess (1.0) ---
    p0 = 1.0
    if find_j0:
        err0, rel_err0, res_jphi = get_j0_error(p0, 0)
    else:
        err0, rel_err0, res_jphi = get_Ip_error(p0, scale_j0, 0)

    # Check if we got lucky immediately
    if rel_err0 < tolerance:
        return p0, res_jphi

    # --- Step 2: Second Guess (Directional) ---
    # If Input < Output (err0 < 0), we need more current -> Try 1.2
    # If Input > Output (err0 > 0), we have too much -> Try 0.8
    
    if find_j0:
        if err0 < 0:
            p1 = 1.2
        else:
            p1 = 0.8
    else:   
        if err0 < 0:
            p1 = 1.1
        else:
            p1 = 0.9
    
    if find_j0:
        err1, rel_err1, res_jphi = get_j0_error(p1, 1)
    else:
        err1, rel_err1, res_jphi = get_Ip_error(p1, scale_j0, 1)

    if rel_err1 < tolerance:
        return p1, res_jphi

    # --- Step 3: Secant Method Loop ---
    # Iterate to find the root where Input - Output = 0
    for i in range(max_iter):
        print(f"--- Optimization Iteration {i+1} ---")
        
        # Avoid division by zero
        if abs(err1 - err0) < 1e-9:
            print("Error difference too small, stopping.")
            break

        # Secant formula: x_new = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        # This projects the line between (p0, err0) and (p1, err1) to zero.
        p_new = p1 - err1 * (p1 - p0) / (err1 - err0)
        
        # Safety clamp: If the projection is wild (negative or huge), constrain it
        # This acts as a loose bracket to prevent the solver from crashing
        if p_new < 0.1: p_new = 0.1
        if p_new > 5.0: p_new = 5.0
        
        # Execute Solver at new point
        if find_j0:
            err_new, rel_err_new, res_jphi = get_j0_error(p_new, i+2)
        else:
            err_new, rel_err_new, res_jphi = get_Ip_error(p_new, scale_j0, i+2)

        if rel_err_new < tolerance:
            print(f"Converged! Optimal scale factor: {p_new:.4f}")
            return p_new, res_jphi
        
        # Update points for next step (move window forward)
        p0, err0 = p1, err1
        p1, err1 = p_new, err_new

    print("Max iterations reached. Returning best last effort.")
    return p1, res_jphi

def solve_with_bootstrap(mygs,
                            ne,
                            Te,
                            ni,
                            Ti,
                            Zeff,
                            Ip_target,
                            inductive_jphi=None,
                            Zis=None,
                            scale_jBS=1.0,
                            isolate_edge_jBS=False,
                            psi_pad=1e-3,
                            iterations=3,
                            diagnostic_plots=False,
                            parameterize_jBS = False):
    '''
    Self-consistently compute bootstrap contribution from H-mode profiles.
    
    If inductive_jphi is set, solve for FF' using Grad-Shafranov equation and 
    iterate solution until all functions of Psi converge.

    @param ne Electron density profile [m^-3]
    @param Te Electron temperature profile [eV]
    @param ni Ion density profile [m^-3]
    @param Ti Ion temperature profile [eV]
    @param Zeff Effective Z profile
    @param Ip_target Target Plasma Current [A]
    @param inductive_jphi Inductive toroidal current profile
    @param Zis List of impurity atomic numbers (default: [1.0])
    @param scale_jBS Factor by which to scale bootstrap current fraction
    @param isolate_edge_jBS If True, isolates edge spike in bootstrap current
    @param psi_pad Padding for flux surface calculations
    @param iterations Maximum number of solver iterations
    @param diagnostic_plots Plot iteration target and output j_tor profiles
    @param initialize_eq Initialize equilibrium solve with flattened pedestal
    '''
    from scipy.optimize import root_scalar
    import matplotlib.pyplot as plt

    try:
        from omfit_classes.utils_fusion import sauter_bootstrap
    except ImportError:
        raise ImportError('omfit_classes.utils_fusion not installed')

    EC = 1.602176634e-19
    
    # Handle mutable default argument
    if Zis is None:
        Zis = [1.]

    ne = numpy.asarray(ne)
    Te = numpy.asarray(Te)
    ni = numpy.asarray(ni)
    Ti = numpy.asarray(Ti)
    Zeff = numpy.asarray(Zeff)
    
    if inductive_jphi is not None:
        inductive_jphi = numpy.asarray(inductive_jphi).copy()
    
    # Calculate Pressure [Pa]
    # p = n * T * k_B. Since T is in eV, k_B is essentially elementary charge e
    pressure = (EC * ne * Te) + (EC * ni * Ti)
    
    # Reconstruct normalized psi grid based on input pressure length
    # Note: Assumes inputs are evenly sampled in psi_norm 0..1
    n_psi = len(pressure)
    psi_N = numpy.linspace(0., 1., n_psi)

    def current_scaling_objective(alpha, j_inductive, j_spike, psi_N, target_ip):
        '''Objective function to match total Ip.'''
        j_total = (alpha * j_inductive) + j_spike
        ip_computed = mygs.flux_integral(psi_N, j_total)
        return ip_computed - target_ip

    def calculate_profiles_and_bootstrap(psi_N, include_jBS):
        '''
        Main physics calculation:
        1. Gets geometry from current equilibrium (self).
        2. Calculates gradients.
        3. Calls Sauter bootstrap model.
        4. Scales inductive current to match Ip_target.
        '''
        # Get geometry and flux functions
        _, f, _, _, _ = mygs.get_profiles(npsi=n_psi, psi_pad=psi_pad)
        _, fc, r_avgs, _ = mygs.sauter_fc(npsi=n_psi, psi_pad=psi_pad)
        
        # Geometry terms
        ft = 1 - fc 
        eps = r_avgs[2] / r_avgs[0]
        _, qvals, ravgs_q, _, _, _ = mygs.get_q(npsi=n_psi, psi_pad=psi_pad)
        R_avg = ravgs_q[0]
        one_over_R_avg = ravgs_q[1]

        # Gradients (using raw psi for derivatives)
        psi_range = mygs.psi_bounds[1] - mygs.psi_bounds[0]
        d_psi = numpy.gradient(psi_N)
        
        # Avoid division by zero in gradients
        d_psi_eff = d_psi * psi_range
        d_psi_eff[d_psi_eff == 0] = 1e-9

        pprime_local = numpy.gradient(pressure) / d_psi_eff
        
        j_BS_final = numpy.zeros_like(pressure)
        
        if include_jBS:
            dn_e_dpsi = numpy.gradient(ne) / d_psi_eff
            dT_e_dpsi = numpy.gradient(Te) / d_psi_eff
            dn_i_dpsi = numpy.gradient(ni) / d_psi_eff
            dT_i_dpsi = numpy.gradient(Ti) / d_psi_eff

            j_BS_neo = sauter_bootstrap(
                psi_N=psi_N, Te=Te, Ti=Ti, ne=ne, p=pressure,
                nis=[ni], Zis=Zis, Zeff=Zeff, gEQDSKs=[None],
                psiraw=psi_N * psi_range + mygs.psi_bounds[0],
                R=R_avg, eps=eps, q=qvals, fT=ft, I_psi=f,
                nt=1, version='neo_2021', debug_plots=False,
                return_units=True, return_package=False,
                charge_number_to_use_in_ion_collisionality='Koh',
                charge_number_to_use_in_ion_lnLambda='Zavg',
                dT_e_dpsi=dT_e_dpsi, dT_i_dpsi=dT_i_dpsi,
                dn_e_dpsi=dn_e_dpsi, dnis_dpsi=[dn_i_dpsi]
            )[0]
            
            # Convert to A/m^2
            j_BS_final = j_BS_neo * (R_avg / f)
            j_BS_final = numpy.nan_to_num(j_BS_final, nan=0.0)

        # Scale Currents to match Ip
        current_jphi_target = inductive_jphi if inductive_jphi is not None else numpy.zeros_like(pressure)

        spike_prof = numpy.zeros_like(j_BS_final)
        
        if include_jBS:
            if isolate_edge_jBS:
                if parameterize_jBS:
                    res = analyze_bootstrap_edge_spike(psi_N, j_BS_final, parameterize_jBS, diagnostic_plots=True)
                    spike_prof = res['parameterized_spike'] * scale_jBS
                else:
                    res = analyze_bootstrap_edge_spike(psi_N, j_BS_final, parameterize_jBS)
                    spike_prof = res['masked_spike'] * scale_jBS
            else:
                spike_prof = j_BS_final * scale_jBS
        
        # Solve for alpha: integral(alpha * j_ind + j_spike) = Ip_target
        try:
            sol = root_scalar(current_scaling_objective,
                                args=(current_jphi_target, spike_prof, psi_N, Ip_target),
                                bracket=[1e-4 * Ip_target, 10 * Ip_target],
                                method='brentq', rtol=1e-6)
            alpha_opt = sol.root
        except ValueError:
            print("WARNING: Root scalar failed to bracket. Defaulting to alpha=1.0")
            alpha_opt = 1.0

        matched_j_inductive = alpha_opt * current_jphi_target
        matched_jphi = matched_j_inductive + spike_prof
        
        # Package results
        pp_dict = {'type': 'linterp', 'x': psi_N, 'y': pprime_local / pprime_local[0]}
        ffp_dict = {'type': 'jphi-linterp', 'x': psi_N, 'y': numpy.nan_to_num(matched_jphi)}
        
        return pp_dict, ffp_dict, j_BS_final, matched_j_inductive, spike_prof

    # --- Main Execution Flow ---
    mygs.set_targets(Ip=Ip_target, pax=pressure[0])

    if inductive_jphi is not None:
        
        print(f'\n >>> Matching input core j_phi with G-S solution')

        # Calculate new profiles
        pp_prof, ffp_prof, j_bs_curr, matched_j_inductive, spike_prof = calculate_profiles_and_bootstrap(
            psi_N, include_jBS=True
        )

        # Enforce P' edge condition
        pp_prof['y'][-1] = 0.

        ffp_prof['type'] = 'jphi-linterp'
        ffp_prof['y'] = matched_j_inductive + spike_prof

        pax_target = pressure[0]

        # Run through once for better profile convergence
        solve_jphi(mygs,ffp_prof,pp_prof,Ip_target,pax_target)

        # Calculate new profiles
        pp_prof, ffp_prof, j_bs_curr, matched_j_inductive, spike_prof = calculate_profiles_and_bootstrap(
            psi_N, include_jBS=True
        )

        # Enforce P' edge condition
        pp_prof['y'][-1] = 0.
        
        print(f'\n >>> Finding optimal j_phi scale factor')
        # Find optimal jphi scale
        final_scale_j0, final_jphi = find_optimal_scale(mygs,
            psi_N, pressure, ffp_prof, pp_prof, matched_j_inductive, 
            Ip_target, psi_pad, spike_prof=spike_prof, find_j0=True, 
            diagnostic_plots=diagnostic_plots
        )
        #  final_scale_j0 = 1.0
        print(f'\n >>> Finding optimal Ip scale factor')
        # Find optimal Ip_target scale
        final_scale_Ip, _ = find_optimal_scale(mygs,
            psi_N, pressure, ffp_prof, pp_prof, matched_j_inductive, 
            Ip_target, psi_pad, spike_prof=spike_prof, find_j0=False, 
            scale_j0=final_scale_j0, 
            tolerance=0.001, diagnostic_plots=diagnostic_plots
        )
        
        print(f'\n >>> Iterating on H-mode equilibrium solution')

        for n in range(iterations):            
            # Calculate new profiles
            pp_prof, ffp_prof, j_bs_curr, matched_j_inductive, spike_prof = calculate_profiles_and_bootstrap(
                psi_N, include_jBS=True
            )

            # Enforce P' edge condition
            pp_prof['y'][-1] = 0.

            matched_input_jphi = final_scale_j0*matched_j_inductive + spike_prof
            ffp_prof['type'] = 'jphi-linterp'
            ffp_prof['y'] = matched_input_jphi
            
            scaled_Ip_target = Ip_target*final_scale_Ip
            pax_target = pressure[0]

            solve_jphi(mygs,ffp_prof,pp_prof,scaled_Ip_target,pax_target)

            # Check Convergence
            _, f, fp, _, pp = mygs.get_profiles(npsi=n_psi, psi_pad=psi_pad)
            _, _, ravgs, _, _, _ = mygs.get_q(npsi=n_psi, psi_pad=psi_pad)

            tmp_jphi = get_jphi_from_GS(f*fp, pp, ravgs[0], ravgs[1])

            if diagnostic_plots:
                plt.figure()
                plt.plot(psi_N, matched_input_jphi/1e6, linestyle='--', label=r'Input $j_\phi$')
                plt.plot(psi_N, tmp_jphi/1e6, label=r'Output $j_\phi$')
                plt.title(f'Iteration {n}')#', Ip error: {Ip_err:.3f} %')
                plt.legend()
                plt.xlabel(r'$\hat \psi')
                plt.xlabel(r'$j_\phi$ [MA]')
                plt.grid(ls=':')
                plt.show()
    else:
        print("Error: inductive_jphi must be specified.")

    results = {'total_j_phi' : tmp_jphi,
                'j_BS' : j_bs_curr,
                'j_inductive' : matched_j_inductive,
                'isolated_j_BS' : spike_prof}
    
    return results
