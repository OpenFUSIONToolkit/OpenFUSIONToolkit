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
    r'''! Generate parameterized edge bootstrap current profile with tunable concave fall-off

    @param psi Normalized poloidal flux profile \f$\hat{\psi}\f$
    @param amp Amplitude of edge skewed Gaussian spike
    @param center Center location of edge skewed Gaussian spike
    @param width Width of edge skewed Gaussian spike
    @param offset Height of core flat profile
    @param sk Skewness parameter of Gaussian
    @param y_sep Value at separatrix (\f$\hat{\psi}=1\f$) for profile
    @param blend_width Blending width for smooth transition in \f$\hat{\psi}\f$
    @param tail_alpha Controls right-side fall-off shape (1.0 = cosine, >1.0 = concave)
    @result Parameterized bootstrap current profile as array over \f$\hat{\psi}\f$
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
        if k_smooth < 1e-5: 
            k_smooth = 1e-5

        # Case A: Standard Spike (Amp > Offset)
        if amp > offset:
            # Solve for internal_amp to compensate for SoftMax boost
            diff = amp - offset
            if diff > 1e-10:
                argument = 1.0 - numpy.exp((offset - amp) / k_smooth)
                if argument <= 0: 
                    argument = 1e-16
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

        dist_to_edge = 1.0 - x_peak
        if dist_to_edge < 1e-5: 
            dist_to_edge = 1e-5

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
        r'''! Compute difference between integrated profile and \f$I_p^{target}\f$

        @param alpha Scaling factor to solve for
        @param psi Normalized poloidal flux profile \f$\hat{\psi}\f$
        @param amp Amplitude of edge spike
        @param center Center of edge spike
        @param width Width of edge spike
        @param offset Core offset
        @param sk Skewness parameter
        @param y_sep Value at separatrix
        @param blend_width Blending width
        @param tail_alpha Right-side fall-off shape
        @result Difference between computed and target offset
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

def Hmode_profiles(edge=0.08, ped=0.4, core=2.5, rgrid=201, expin=1.5, expout=1.5, widthp=0.04, xphalf=None):
    r'''! This function generates H-mode density and temperature profiles evenly
    spaced in your favorite radial coordinate. Copied from https://omfit.io/_modules/omfit_classes/utils_fusion.html

    @param edge Separatrix height (float)
    @param ped Pedestal height (float)
    @param core On-axis profile height (float)
    @param rgrid Number of radial grid points (int)
    @param expin Inner core exponent for H-mode pedestal profile (float)
    @param expout Outer core exponent for H-mode pedestal profile (float)
    @param widthp Width of pedestal (float)
    @param xphalf Position of tanh (float, optional)
    @result H-mode profile array over radial grid
    '''

    w_E1 = 0.5 * widthp  # width as defined in eped
    if xphalf is None:
        xphalf = 1.0 - w_E1

    xped = xphalf - w_E1

    pconst = 1.0 - numpy.tanh((1.0 - xphalf) / w_E1)
    a_t = 2.0 * (ped - edge) / (1.0 + numpy.tanh(1.0) - pconst)

    coretanh = 0.5 * a_t * (1.0 - numpy.tanh(-xphalf / w_E1) - pconst) + edge

    xpsi = numpy.linspace(0, 1, rgrid)
    ones = numpy.ones(rgrid)

    val = 0.5 * a_t * (1.0 - numpy.tanh((xpsi - xphalf) / w_E1) - pconst) + edge * ones

    xtoped = xpsi / xped
    for i in range(0, rgrid):
        if xtoped[i] ** expin < 1.0:
            val[i] = val[i] + (core - coretanh) * (1.0 - xtoped[i] ** expin) ** expout

    return val

def analyze_bootstrap_edge_spike(psi_N, j_bootstrap, diagnostic_plots=False):
    r'''! Analyze bootstrap edge spike location, width, and height

    @param psi_N Normalized poloidal flux profile \f$\hat{\psi}\f$
    @param j_bootstrap Bootstrap current profile \f$j_{BS}(\hat{\psi})\f$
    @param diagnostic_plots If True, plot diagnostic figures
    @result Dictionary with spike properties and parameterized spike model
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
    
    # Attempt to fit bootstrap parameterization to calculated profile (deprecated)
    # Initial parameter guess
    sigma_init = fwhm/2.355  # Convert FWHM to sigma
    p0 = [peak_height, peak_psi, sigma_init, lmin_j_BS, 1.0, j_bootstrap[-1], 0.05]

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
    
    if diagnostic_plots:
        plt.figure()
        plt.plot(psi_N[fit_mask],j_bootstrap[fit_mask]/1e6,label='Input j_BS')
        plt.plot(psi_N,spike_only/1e6,label='Fitted parameterization')
        plt.plot(psi_N,masked_spike/1e6,label='Isolated')
        plt.xlabel(r'$\hat \psi$')
        plt.ylabel(r'$j_\phi$ [MA/m$^2$]')
        plt.xlim(0.5,1.02)
        plt.legend(loc='best')
        plt.show()
    
    results = {
        'sigma': width,                 # Gaussian width (sigma)
        'background': offset,           # background level
        'gaussian_params': popt,        # Raw parameters [amp, center, width, offset]
        'parameterized_spike': spike_only,    # Array of spike component values
        'masked_spike': masked_spike
    }

    return results

def solve_jphi(mygs,ffp_prof,pp_prof,Ip_target,pax_target):
    r'''! Solve Grad-Shafranov equilibrium for given profiles

    @param mygs Grad-Shafranov solver object
    @param ffp_prof Toroidal current profile (\f$j_\phi(\hat{\psi})\f$), type 'jphi-linterp'
    @param pp_prof Pressure gradient profile (\f{P'}\f$)
    @param Ip_target Target plasma current \f$I_p\f$ [A]
    @param pax_target Target on-axis pressure \f$p_{axis}\f$ [Pa]
    @result None (updates equilibrium in-place)
    '''

    ffp_prof['type'] = 'jphi-linterp'

    mygs.set_targets(Ip=Ip_target, pax=pax_target)
    mygs.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof)

    # Solve Grad-Shafranov
    mygs.solve()
    
def find_optimal_scale(mygs, psi_N, pressure, ffp_prof, pp_prof, j_inductive,  
                            Ip_target, psi_pad, spike_prof=None, find_j0=True, scale_j0=1.0, 
                            tolerance=0.01, max_iter=5, diagnostic_plots=False):
    r'''! Optimize scaling to match input/output \f$j_0\f$ or \f$I_p\f$

    @param mygs Grad-Shafranov solver object
    @param psi_N Normalized poloidal flux profile \f$\hat{\psi}\f$
    @param pressure Plasma pressure profile \f$p(\hat{\psi})\f$ [Pa]
    @param ffp_prof Toroidal current profile (\f$j_\phi(\hat{\psi})\f$)
    @param pp_prof Pressure gradient profile (\f{P'}\f$)
    @param j_inductive Inductive/total toroidal current profile \f$j_{ind}(\hat{\psi})\f$ [A/m^2]
    @param Ip_target Target plasma current \f$I_p\f$ [A]
    @param psi_pad Padding for flux surface calculations
    @param spike_prof Optional bootstrap current profile \f$j_{BS}(\hat{\psi})\f$ [A/m^2]
    @param find_j0 If True, rescale core \f$j_\phi\f$; else, rescale \f$I_p\f$
    @param scale_j0 Used if find_j0=False; rescaling factor for \f$j_0\f$
    @param tolerance Relative error tolerance for secant method
    @param max_iter Maximum number of secant iterations
    @param diagnostic_plots If True, plot input/output \f$j_\phi\f$
    @result Tuple (optimal scale factor, output \f$j_\phi(\hat{\psi})\f$)
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
            plt.xlabel(r'$\hat \psi$')
            plt.ylabel(r'$j_\phi$ [MA/m$^2$]')
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
        if p_new < 0.1: 
            p_new = 0.1
        if p_new > 5.0: 
            p_new = 5.0
        
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

def calculate_ln_lambda(Te, Ti, ne, ni, Zeff=1.0,
                        electron_lnLambda_model='sauter',
                        ion_lnLambda_model='unity'):
    r'''! Calculate Coulomb logarithm for electrons and ions

    @param Te Electron temperature profile \f$T_e\f$ [eV]
    @param Ti Ion temperature profile \f$T_i\f$ [eV]
    @param ne Electron density profile \f$n_e\f$ [m$^{-3}$]
    @param ni Ion density profile \f$n_i\f$ [m$^{-3}$]
    @param Zeff Effective charge \f$Z_{eff}\f$
    @param electron_lnLambda_model Model for electron Coulomb log ('sauter', 'NRL')
    @param ion_lnLambda_model Model for ion Coulomb log ('unity', 'Zavg')
    @result Tuple (\f$\ln\Lambda_e\f$, \f$\ln\Lambda_{ii}\f$)
    '''
    if electron_lnLambda_model == 'NRL':
        ln_lambda_e = (23.5 
                       - numpy.log(numpy.sqrt(ne / 1e6) * Te**(-5.0/4.0))
                       - numpy.sqrt(1e-5 + (numpy.log(Te) - 2)**2 / 16.0))
    else:
        ln_lambda_e = 31.3 - numpy.log(numpy.sqrt(ne) / Te)

    if ion_lnLambda_model == 'Zavg':
        Z_lnLam = numpy.clip(ne / ni, 1.0, None)
    else:
        Z_lnLam = 1.0
    ln_lambda_ii = 30.0 - numpy.log(Z_lnLam**3 * numpy.sqrt(ni) / Ti**1.5)

    ln_lambda_e = numpy.maximum(ln_lambda_e, 10.0)
    ln_lambda_ii = numpy.maximum(ln_lambda_ii, 10.0)
    return ln_lambda_e, ln_lambda_ii


def redl_bootstrap(
    psi_N=None,
    Te=None, Ti=None,
    ne=None, ni=None,
    pe=None, pi=None,
    Zeff=None,
    R=None, q=None, eps=None, fT=None, I_psi=None,
    dT_e_dpsi=None, dT_i_dpsi=None,
    dn_e_dpsi=None, dn_i_dpsi=None,
    dp_dpsi=None,
    ln_lambda_e=17.0, ln_lambda_ii=17.0,
    # --- Legacy-matching toggles ---
    use_legacy_L34=False,
    use_sign_q=False,
    ion_collisionality_model='Zeff',
    formula_form='jB',
    Zeff_override=None,          # ← add this back
    nu_e_star_override=None,      # ← NEW
    nu_i_star_override=None,      # ← NEW
):
    r'''! Bootstrap current via Redl et al., Phys. Plasmas 28, 022502 (2021)

    @param psi_N Normalized poloidal flux profile \f$\hat{\psi}\f$
    @param Te Electron temperature profile \f$T_e\f$ [eV]
    @param Ti Ion temperature profile \f$T_i\f$ [eV]
    @param ne Electron density profile \f$n_e\f$ [m$^{-3}$]
    @param ni Ion density profile \f$n_i\f$ [m$^{-3}$]
    @param pe Electron pressure profile \f$p_e\f$ [Pa]
    @param pi Ion pressure profile \f$p_i\f$ [Pa]
    @param Zeff Effective charge \f$Z_{eff}\f$
    @param R Major radius \f$R\f$ [m]
    @param q Safety factor \f$q\f$
    @param eps Inverse aspect ratio \f$\epsilon\f$
    @param fT Trapped particle fraction \f$f_T\f$
    @param I_psi Toroidal current function \f$I(\psi)\f$
    @param dT_e_dpsi Derivative of \f$T_e\f$ w.r.t. \f$\psi\f$
    @param dT_i_dpsi Derivative of \f$T_i\f$ w.r.t. \f$\psi\f$
    @param dn_e_dpsi Derivative of \f$n_e\f$ w.r.t. \f$\psi\f$
    @param dn_i_dpsi Derivative of \f$n_i\f$ w.r.t. \f$\psi\f$
    @param dp_dpsi Derivative of total pressure w.r.t. \f$\psi\f$
    @param ln_lambda_e Electron Coulomb logarithm \f$\ln\Lambda_e\f$
    @param ln_lambda_ii Ion Coulomb logarithm \f$\ln\Lambda_{ii}\f$
    @param use_legacy_L34 Use legacy L34 formula (see Redl Eq. 19)
    @param use_sign_q Multiply by sign of \f$q\f$
    @param ion_collisionality_model Model for ion collisionality ('Zeff', 'Koh', 'unity')
    @param formula_form Use 'jB' (default) or 'jboot1' equation
    @param Zeff_override Override \f$Z_{eff}\f$ if set
    @param nu_e_star_override Override electron collisionality if set
    @param nu_i_star_override Override ion collisionality if set
    @result Tuple (bootstrap current profile \f$j_{BS}(\hat{\psi})\f$, coefficient dictionary)
    '''
    
    # Apply Zeff override if requested
    if Zeff_override is not None:
        Zeff = (numpy.full_like(Zeff, Zeff_override) 
                if numpy.isscalar(Zeff_override) else Zeff_override)

    inputs = [Te, Ti, ne, ni, pe, pi, Zeff, R, q, eps, fT, I_psi,
              dT_e_dpsi, dT_i_dpsi, dn_e_dpsi]
    if any(x is None for x in inputs):
        raise ValueError("All primary inputs must be provided.")
    if formula_form == 'jboot1' and dn_i_dpsi is None:
        raise ValueError("dn_i_dpsi required for jboot1 form.")

    EC = 1.602176634e-19
    p_total = pe + pi
    R_pe = pe / p_total
        
    # =====================================================================
    # 1. Collisionality
    # =====================================================================

    if nu_e_star_override is not None:
        nu_e_star = nu_e_star_override
    else:
        nu_e_star = (6.921e-18 * q * R * ne * Zeff * ln_lambda_e
                     / (Te**2 * eps**1.5))

    if nu_i_star_override is not None:
        nu_i_star = nu_i_star_override
    else:
        if ion_collisionality_model == 'Koh':
            nu_i_star = (4.90e-18 * q * R * ni * Zeff * ln_lambda_ii
                         / (Ti**2 * eps**1.5))
        elif ion_collisionality_model == 'unity':
            nu_i_star = (4.90e-18 * q * R * ni * 1.0 * ln_lambda_ii
                         / (Ti**2 * eps**1.5))
        else:
            nu_i_star = (4.90e-18 * q * R * ni * (Zeff**4) * ln_lambda_ii
                         / (Ti**2 * eps**1.5))

    # =====================================================================
    # 2. L31 (Redl Eqs. 10-11)
    # =====================================================================

    ft_31_d1 = ((0.67 * (1 - 0.7 * fT) * numpy.sqrt(nu_e_star))
                / (0.56 + 0.44 * Zeff))
    ft_31_d2 = (((0.52 + 0.086 * numpy.sqrt(nu_e_star))
                 * (1 + 0.87 * fT) * nu_e_star)
                / (1 + 1.13 * numpy.sqrt(numpy.maximum(Zeff - 1, 0.0))))
    X31 = fT / (1 + ft_31_d1 + ft_31_d2)

    def F31_poly(X):
        dZ = Zeff**1.2 - 0.71
        return ((1 + 0.15 / dZ) * X
                - (0.22 / dZ) * X**2
                + (0.01 / dZ) * X**3
                + (0.06 / dZ) * X**4)

    L31 = F31_poly(X31)

    # =====================================================================
    # 2b. L34
    # =====================================================================

    if use_legacy_L34:
        f33_d1 = (0.25 * (1 - 0.7 * fT) * numpy.sqrt(nu_e_star)
                  * (1 + 0.45 * numpy.sqrt(numpy.maximum(Zeff - 1, 0.0))))
        f33_d2 = ((0.61 * (1 - 0.41 * fT) * nu_e_star)
                  / numpy.sqrt(Zeff))
        f33teff = fT / (1 + f33_d1 + f33_d2)
        L34 = F31_poly(f33teff)
    else:
        L34 = L31

    # =====================================================================
    # 3. L32 (Redl Eqs. 12-16)
    # =====================================================================

    dee_2 = ((0.23 * (1 - 0.96 * fT) * numpy.sqrt(nu_e_star))
             / numpy.sqrt(Zeff))
    dee_3 = ((0.13 * (1 - 0.38 * fT) * nu_e_star / (Zeff**2))
             * (numpy.sqrt(1 + 2 * numpy.sqrt(numpy.maximum(Zeff - 1, 0.0)))
                + fT**2 * numpy.sqrt(
                    (0.075 + 0.25 * (Zeff - 1)**2) * nu_e_star)))
    X32_ee = fT / (1 + dee_2 + dee_3)
    F32_ee = (
        (0.1 + 0.6 * Zeff)
        / (Zeff * (0.77 + 0.63 * (1 + (Zeff - 1)**1.1)))
        * (X32_ee - X32_ee**4)
        + 0.7 / (1 + 0.2 * Zeff)
        * (X32_ee**2 - X32_ee**4 - 1.2 * (X32_ee**3 - X32_ee**4))
        + 1.3 / (1 + 0.5 * Zeff) * X32_ee**4)

    dei_2 = ((0.87 * (1 + 0.39 * fT) * numpy.sqrt(nu_e_star))
             / (1 + 2.95 * (Zeff - 1)**2))
    dei_3 = (1.53 * (1 - 0.37 * fT) * nu_e_star
             * (2 + 0.375 * (Zeff - 1)))
    X32_ei = fT / (1 + dei_2 + dei_3)
    F32_ei = (
        -(0.4 + 1.93 * Zeff)
        / (Zeff * (0.8 + 0.6 * Zeff))
        * (X32_ei - X32_ei**4)
        + 5.5 / (1.5 + 2 * Zeff)
        * (X32_ei**2 - X32_ei**4 - 0.8 * (X32_ei**3 - X32_ei**4))
        - 1.3 / (1 + 0.5 * Zeff) * X32_ei**4)

    L32 = F32_ee + F32_ei

    # =====================================================================
    # 4. Alpha (Redl Eqs. 20-21)
    # =====================================================================

    alpha0 = (
        -(0.62 + 0.055 * (Zeff - 1))
        / (0.53 + 0.17 * (Zeff - 1))
        * (1 - fT)
        / (1 - (0.31 - 0.065 * (Zeff - 1)) * fT - 0.25 * fT**2))

    alpha = (
        ((alpha0 + 0.7 * Zeff * numpy.sqrt(fT) * numpy.sqrt(nu_i_star))
         / (1 + 0.18 * numpy.sqrt(nu_i_star))
         - 0.002 * nu_i_star**2 * fT**6)
        / (1 + 0.004 * nu_i_star**2 * fT**6))

    # =====================================================================
    # 5. Assemble
    # =====================================================================

    if formula_form == 'jboot1':
        if dp_dpsi is None:
            dp_dpsi = (ne * dT_e_dpsi + Te * dn_e_dpsi
                       + ni * dT_i_dpsi + Ti * dn_i_dpsi) * EC
        bra1 = L31 * dp_dpsi / pe
        bra2 = L32 * dT_e_dpsi / Te
        bra3 = L34 * alpha * (1 - R_pe) / R_pe * dT_i_dpsi / Ti
        j_bootstrap = -I_psi * pe * (bra1 + bra2 + bra3)
        term_decomp = {
            'bra1': -I_psi * pe * bra1,
            'bra2': -I_psi * pe * bra2,
            'bra3': -I_psi * pe * bra3,
        }
    else:
        term_n = p_total * L31 * (1.0 / ne) * dn_e_dpsi
        term_Te = pe * (L31 + L32) * (1.0 / Te) * dT_e_dpsi
        term_Ti = pi * (L31 + L34 * alpha) * (1.0 / Ti) * dT_i_dpsi
        j_bootstrap = -I_psi * (term_n + term_Te + term_Ti)
        term_decomp = {
            'term_n': -I_psi * term_n,
            'term_Te': -I_psi * term_Te,
            'term_Ti': -I_psi * term_Ti,
        }

    if use_sign_q:
        j_bootstrap *= numpy.sign(q)
        for k in term_decomp:
            term_decomp[k] *= numpy.sign(q)
    
    coeffs = {
        'L31': L31, 'L32': L32, 'L34': L34,
        'alpha': alpha, 'alpha0': alpha0,
        'nu_e_star': nu_e_star, 'nu_i_star': nu_i_star,
        'Zeff_used': Zeff,
        'R_pe': R_pe,
        'formula_form': formula_form,
        'F32_ee': F32_ee,
        'F32_ei': F32_ei,
    }
    coeffs.update(term_decomp)

    return j_bootstrap, coeffs

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
                         parameterize_jBS = False,
                         use_OMFIT_sauter = False):
    r'''! Self-consistently compute bootstrap current from H-mode profiles

    @param mygs Grad-Shafranov solver object
    @param ne Electron density profile \f$n_e(\hat{\psi})\f$ [m$^{-3}$]
    @param Te Electron temperature profile \f$T_e(\hat{\psi})\f$ [eV]
    @param ni Ion density profile \f$n_i(\hat{\psi})\f$ [m$^{-3}$]
    @param Ti Ion temperature profile \f$T_i(\hat{\psi})\f$ [eV]
    @param Zeff Effective charge profile \f$Z_{eff}(\hat{\psi})\f$
    @param Ip_target Target plasma current \f$I_p\f$ [A]
    @param inductive_jphi Inductive toroidal current profile \f$j_{ind}(\hat{\psi})\f$
    @param Zis List of impurity atomic numbers (default: [1.0])
    @param scale_jBS Scaling factor for bootstrap current
    @param isolate_edge_jBS If True, isolate edge spike in bootstrap current
    @param psi_pad Padding for flux surface calculations
    @param iterations Number of solver iterations
    @param diagnostic_plots If True, plot diagnostic figures
    @param parameterize_jBS If True, use parameterized edge spike
    @param use_OMFIT_sauter If True, use OMFIT Sauter model
    @result Dictionary with total, bootstrap, inductive, and isolated edge current profiles
    '''
    from scipy.optimize import root_scalar
    import matplotlib.pyplot as plt

    if use_OMFIT_sauter:
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

            if use_OMFIT_sauter:
                j_BS_neo = sauter_bootstrap( # legacy OMFIT implementation
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
            else:
                # Coulomb logarithm: NRL formulary
                ln_le, ln_lii = calculate_ln_lambda(
                    Te, Ti, ne, ni, Zeff,
                    electron_lnLambda_model='NRL',   # more accurate than Sauter
                    ion_lnLambda_model='Zavg',       # consistent multi-species treatment
                )

                # Ion collisionality: Koh multi-species formula
                Zdom = 1.0  # deuterium
                Zavg = ne / ni
                Zion = (Zdom**2 * Zavg * Zeff)**0.25
                nu_i_star = (4.90e-18 * numpy.abs(qvals) * R_avg * ni 
                            * Zion**4 * ln_lii / (Ti**2 * eps**1.5))

                # Electron collisionality
                nu_e_star = (6.921e-18 * numpy.abs(qvals) * R_avg * ne 
                            * Zeff * ln_le / (Te**2 * eps**1.5))
                
                j_BS_neo, _ = redl_bootstrap( # native re-implementation of Redl 2021
                    psi_N=psi_N, Te=Te, Ti=Ti, ne=ne, ni=ni,
                    pe=EC*(ne*Te), pi=EC*(ni*Ti),
                    Zeff=Zeff, R=R_avg, q=qvals, eps=eps, fT=ft, I_psi=f,
                    dT_e_dpsi=dT_e_dpsi, dT_i_dpsi=dT_i_dpsi,
                    dn_e_dpsi=dn_e_dpsi, dn_i_dpsi=dn_i_dpsi,
                    ln_lambda_e=ln_le, ln_lambda_ii=ln_lii,
                    # No Zeff_override — OMFIT preserves input Zeff
                    nu_e_star_override=nu_e_star,
                    nu_i_star_override=nu_i_star,
                    use_legacy_L34=False, # L34 = L31 (correct Redl 2021)
                    use_sign_q=True, # convention-dependent
                    formula_form='jboot1', # no d(ln ne)=d(ln ni) assumption
                )
            
            # Convert to A/m^2
            j_BS_final = j_BS_neo * (R_avg / f)
            j_BS_final = numpy.nan_to_num(j_BS_final, nan=0.0)

            # to-do: project j_BS_parallel to j_phi more accurately?

        # Scale Currents to match Ip
        current_jphi_target = inductive_jphi if inductive_jphi is not None else numpy.zeros_like(pressure)

        spike_prof = numpy.zeros_like(j_BS_final)
        
        if include_jBS:
            if isolate_edge_jBS:
                if parameterize_jBS:
                    res = analyze_bootstrap_edge_spike(psi_N, j_BS_final, diagnostic_plots=diagnostic_plots)
                    spike_prof = res['parameterized_spike'] * scale_jBS
                else:
                    res = analyze_bootstrap_edge_spike(psi_N, j_BS_final)
                    spike_prof = res['masked_spike'] * scale_jBS
            else:
                spike_prof = j_BS_final * scale_jBS
        
        # Solve for alpha: integral(alpha * j_ind + j_spike) = Ip_target
        try:
            sol = root_scalar(current_scaling_objective,
                                args=(current_jphi_target, spike_prof, psi_N, Ip_target),
                                bracket=[1.0, 10. * Ip_target],
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
        
        print('\n >>> Matching input core j_phi with G-S solution')

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
        
        print('\n >>> Finding optimal j_phi scale factor')
        # Find optimal jphi scale
        final_scale_j0, final_jphi = find_optimal_scale(mygs,
            psi_N, pressure, ffp_prof, pp_prof, matched_j_inductive, 
            Ip_target, psi_pad, spike_prof=spike_prof, find_j0=True, 
            diagnostic_plots=diagnostic_plots
        )
        #  final_scale_j0 = 1.0
        print('\n >>> Finding optimal Ip scale factor')
        # Find optimal Ip_target scale
        final_scale_Ip, _ = find_optimal_scale(mygs,
            psi_N, pressure, ffp_prof, pp_prof, matched_j_inductive, 
            Ip_target, psi_pad, spike_prof=spike_prof, find_j0=False, 
            scale_j0=final_scale_j0, 
            tolerance=0.001, diagnostic_plots=diagnostic_plots
        )
        
        print('\n >>> Iterating on H-mode equilibrium solution')

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
                plt.xlabel(r'$\hat \psi$')
                plt.ylabel(r'$j_\phi$ [MA/m$^2$]')
                plt.grid(ls=':')
                plt.show()
    else:
        raise ValueError("inductive_jphi must be specified.")

    results = {'total_j_phi' : tmp_jphi,
                'j_BS' : j_bs_curr,
                'j_inductive' : matched_j_inductive,
                'isolated_j_BS' : spike_prof,
                'scale_j0' : final_scale_j0,
                'scale_Ip' : final_scale_Ip}
    
    return results
