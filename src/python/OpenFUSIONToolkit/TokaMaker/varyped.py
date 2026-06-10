'''! Functions and helpers to run varyped equilibrium scan on TokaMaker instance 

Runs equilibrium scan over range of scaling values using profiles from 
input an g-File and p-File. Pedestal of pressure and edge spike of 
current density profiles are scaled self-consistently, and the results 
are saved to new g-files and p-files generated with TokaMaker.

@authors Kevin Clavijo
@date June 2026
'''

import os
import copy
import numpy as np

from scipy.signal import find_peaks, peak_widths
from scipy.optimize import curve_fit, leastsq
from scipy.stats import skewnorm
from scipy.optimize import root_scalar, least_squares
from scipy.constants import constants

from OpenFUSIONToolkit.TokaMaker import eqdsk

#=============================================================================
#                       PRESSURE SCALING METHODS
#=============================================================================

def locate_pedestal_start(t_object, psi_n, pressure):
    psi_n_6_idx = np.where(np.isclose(psi_n, 0.6, atol=0.01) == True)[0][0]

    pprime = np.gradient(pressure) / (np.gradient(psi_n) * (t_object.psi_bounds[1] - t_object.psi_bounds[0]))
    shifted_pprime_slice = pprime[psi_n_6_idx:] - 0.5 * np.min(pprime[psi_n_6_idx:])
    
    # Find sign changes in shifted_pprime_slice (approximates zero crossings)
    sign_changes = np.where(np.diff(np.sign(shifted_pprime_slice)))[0]
    
    if len(sign_changes) < 2:
        # Fallback: use steepest gradient location
        steep_idx = np.argmin(pprime)
        return max(0.6, psi_n[steep_idx] - 0.08)
    
    # Map sign change indices back to original psi_n array
    roots_idxs = sign_changes + psi_n_6_idx
    
    # Use the last two zero crossings (as in original code)
    pedestal_start = psi_n[roots_idxs[-2]] - 0.5 * (psi_n[roots_idxs[-1]] - psi_n[roots_idxs[-2]])
    
    return pedestal_start

def pressure_scaling_function(psi_n, pedestal_factor, norm_factor, pedestal_start):

    u = np.where( psi_n < pedestal_start, (pedestal_start-psi_n)/pedestal_start, 0.)
    if pedestal_factor <= 1.0:
        psf = pedestal_factor*(np.cosh(norm_factor[0]*u)**0.5) ##used when pedestal is scaled down
    else:
        psf = pedestal_factor/(np.cosh(norm_factor[0]*u)**0.5) ##used when pedestal is scaled up

    return psf

def pressure_objective_function(norm_factor, t_object, psi_n, pressure, pedestal_factor, pedestal_start):
    #calculate initial energy stored (v integral of pressure)
    w_0 = (1.5e3)*t_object.flux_integral(psi_n, pressure)

    #get pressure scaling function
    psf = pressure_scaling_function(psi_n, pedestal_factor, norm_factor, pedestal_start)

    #scale pressure using psf
    scaled_pressure = pressure*psf
    w_1 = (1.5e3)*t_object.flux_integral(psi_n, scaled_pressure)

    #return value to be minimized as numpy array
    return np.array([w_1 - w_0])

def scale_pressure(t_object, psi_n, pressure, pedestal_factor):
    w_0 = (1.5e3)*t_object.flux_integral(psi_n, pressure)
    norm_factor_0 = [0.0]
    pedestal_start = locate_pedestal_start(t_object, psi_n, pressure)

    # Use leastsq with tight tolerances 
    result_lsq = leastsq(pressure_objective_function, 
                         norm_factor_0, 
                         args=(t_object, psi_n, pressure, pedestal_factor, pedestal_start),
                         ftol=1.e-7, xtol=1.e-7, epsfcn=1.e-7, full_output=True)
    
    if result_lsq[4] > 4:  # imerit > 4 indicates poor fit
        # Grid search fallback
        norm_factor_guesses = np.logspace(-2, 1, 100)
        w_differences = []
        for nf in norm_factor_guesses:
            psf = pressure_scaling_function(psi_n, pedestal_factor, [nf], pedestal_start)
            scaled_pressure = pressure*psf
            w_1 = (1.5e3)*t_object.flux_integral(psi_n, scaled_pressure)
            w_differences.append(w_1 - w_0)
        norm_factor = [norm_factor_guesses[np.argmin(np.abs(w_differences))]]
        method = "grid_search"
    else:
        norm_factor = result_lsq[0]
        method = "leastsq"
    
    psf = pressure_scaling_function(psi_n, pedestal_factor, norm_factor, pedestal_start)
    scaled_pressure = pressure*psf
    w_1 = (1.5e3)*t_object.flux_integral(psi_n, scaled_pressure)
    w_diff = w_1 - w_0
    return scaled_pressure, psf, w_diff, w_0, w_1, method


#=============================================================================
#                       CURRENT SCALING METHODS
#=============================================================================

def skew_gaussian(x, amp, center, width, offset, sk):
    skew = skewnorm.pdf(x, sk, center, width)
    skew = skew/max(skew)
    skew *= amp
    skew += offset
    return skew

def analyze_bootstrap_edge_spike(psi_n, j_bootstrap):
    r'''! Analyze bootstrap edge spike location, width, and height

    @param psi_n Normalized psi profile
    @param j_bootstrap Bootstrap current profile
    Returns:
    dict: Dictionary containing pedestal properties and spike model
    '''
    # Focus on the edge region (psi_n > 0.8)
    edge_mask = psi_n >= 0.7
    psi_edge = psi_n[edge_mask]
    j_edge = j_bootstrap[edge_mask]

    # Find peak in the edge region
    peaks, properties = find_peaks(j_edge, height=0.)#0.5*np.max(j_edge))

    if len(peaks) == 0:
        print("No clear peak found in the edge region")
        return None

    # Choose peak closest to psi_n = 1 if multiple peaks exist
    peak_idx = peaks[np.argmax(psi_edge[peaks])]
    peak_psi = psi_edge[peak_idx]
    peak_height = j_edge[peak_idx]

    # Calculate initial FWHM (full width at half maximum)
    widths = peak_widths(j_edge, [peak_idx], rel_height=0.5)
    left_idx, right_idx = int(widths[2][0]), int(widths[3][0])

    # Convert to psi_n coordinates
    fwhm = psi_edge[right_idx] - psi_edge[left_idx]

    # Get good range for fitting - wider than the spike to capture baseline
    fit_range = max(0.15, 3*fwhm)  # At least 0.15 in psi_n or 3x FWHM
    fit_mask = (psi_n >= peak_psi - fit_range) & (psi_n <= min(1.0, peak_psi + fit_range))

    # Initial parameter guess
    # Estimate the background level from points away from the peak
    background_mask = (psi_n >= 0.7) & (psi_n <= 0.75)
    if np.sum(background_mask) > 5:
        background_level = np.median(j_bootstrap[background_mask])
    else:
        background_level = np.min(j_edge)

    sigma_init = fwhm/2.355  # Convert FWHM to sigma
    p0 = [peak_height - background_level, peak_psi, sigma_init, background_level, 1.0]

    # Perform the fit
    popt, pcov = curve_fit(skew_gaussian, psi_n[fit_mask], j_bootstrap[fit_mask], p0=p0, 
                          maxfev=10000)

    amp, center, width, offset, sk = popt
    spike_only = skew_gaussian(psi_n, amp, center, width, offset, sk)

    results = {
        'sigma': width,                 # Gaussian width (sigma)
        'background': offset,           # background level
        'gaussian_params': popt,        # Raw parameters [amp, center, width, offset]
        'spike_profile': spike_only,    # Array of spike component values
    }

    return results

def fit_jphi_profiles(xdata, ped=0.05, core=1.0, widthp=0.02, xphalf=0.98, amp=0.1, center=0.95, width=0.3, sk=0., c=0., alpha=1.5, beta=1.5, gamma=1.5, delta=1.5):
    edge = 0.0
    
    w_E1 = 0.5 * widthp  # width as defined in eped
    if xphalf is None:
        xphalf = 1.0 - w_E1

    xped = xphalf - w_E1

    pconst = 1.0 - np.tanh((1.0 - xphalf) / w_E1)
    a_t = 2.0 * (ped - edge) / (1.0 + np.tanh(1.0) - pconst)

    coretanh = 0.5 * a_t * (1.0 - np.tanh(-xphalf / w_E1) - pconst) + edge
    
    rgrid=len(xdata)
    xpsi = np.linspace(0, 1, rgrid)
    ones = np.ones(rgrid)

    val = 0.5 * a_t * (1.0 - np.tanh((xpsi - xphalf) / w_E1) - pconst) + edge * ones
    
    xtoped = xpsi / xped
    
    for i in range(0, rgrid):
        val[i] = val[i] + (core - coretanh) * (alpha*xpsi[i] + beta*xpsi[i]**2 + gamma*xpsi[i]**3 + delta*xpsi[i]**4) + c

    val += skew_gaussian(xdata, amp, center, width, 0.0, sk)

    return val

def fit_jphi_from_input_jtor(psi, input_jtor, scale_factor=None):
    # get edge spike parameters
    jBS_results = analyze_bootstrap_edge_spike(psi, input_jtor)
    spike_popt = jBS_results['gaussian_params']
    amp, center, width, offset, sk = spike_popt

    if scale_factor is None:
        scale_factor = 1
    
    tmp_jtor = input_jtor
    # set curve_fit initial guesses
    
    p0 = [offset, tmp_jtor[0], 0.02, 0.98, (amp), center, width, sk, 0., -0.1, -0.7, -1.5, 1.5]

    # set reasonable lower and upper bounds on all values
    lower_bounds = [0.8*offset, 0.8*input_jtor[0], 0.0, 0.96, 0.99*(amp), 0.99*center, 0.99*width, sk+(0.01*sk), 0.0, -10., -10., -10., -10.]
    upper_bounds = [1.2*offset, 1.2*input_jtor[0], 0.04, 1.00, 1.01*(amp), 1.01*center, 1.01*width, sk-(0.01*sk), 1.0e+20, 10.,  10.,  10.,  10.]
    bounds = (lower_bounds, upper_bounds)

    popt, pcov = curve_fit(fit_jphi_profiles, psi, input_jtor, p0, bounds=bounds)

    fit_jphi = fit_jphi_profiles(psi, ped=popt[0], core=popt[1], widthp=popt[2], xphalf=popt[3], amp=popt[4]*scale_factor, center=popt[5], width=popt[6], sk=popt[7], c=popt[8], alpha=popt[9], beta=popt[10], gamma=popt[11], delta=popt[12])

    return fit_jphi

##################### #DB
def current_objective_function(alpha, t_object, jtor_prof, psi_n, Ip_target):
    scale_func = alpha*np.tanh(2.0*psi_n)-(alpha-1.0)
    integrand = scale_func*jtor_prof
    Ip_computed = t_object.flux_integral(psi_n, integrand)
    return Ip_computed - Ip_target
#####################


def scale_curr(t_object, psi_n, j_tor, scale_c):

    j_tor_copy = j_tor.copy()
    fit_jphi = fit_jphi_from_input_jtor(psi_n, j_tor_copy, scale_c)
    
    ################ #DB
    Ip_target = t_object.flux_integral(psi_n, j_tor)
    
    # Find scalar "alpha" that solves integral(alpha*inductive_jtor) = Ip_target
    result = root_scalar(current_objective_function,
                         args=(t_object, fit_jphi, psi_n, Ip_target), # Remaining arguments for objective_function()
                         bracket=[-1.0, 1.0],  # Provide a reasonable initial bracket 
                         method='brentq',    # Robust bracketing method
                         rtol=1e-6)          # Relative tolerance

    a_optimal = result.root
    

    scale_func = a_optimal*np.tanh(3.0*psi_n)-(a_optimal-1.0)
    scaled_jtor = scale_func*fit_jphi
    
    Ip_final = t_object.flux_integral(psi_n, scaled_jtor)
    Ip_diff = Ip_final - Ip_target
    #################

    return scaled_jtor, Ip_diff, Ip_final

def scale_current(t_object, psi_n, j_tor, scale_c):
    try:

        if scale_c <= 0:
            raise ValueError(f"scale_c must be positive, got {scale_c}")
        
        if j_tor is None or len(j_tor) == 0:
            raise ValueError("j_tor cannot be None or empty")
        
        if psi_n is None:
            psi_n = np.linspace(0., 1., len(j_tor))
        
        scaled_jtor, Ip_diff, Ip_final = scale_curr(t_object, psi_n, j_tor, scale_c)            
        return scaled_jtor, Ip_diff, Ip_final
    
    except Exception as e:
        print(f"Error during current scaling: {str(e)}")
        raise


#=============================================================================
#                       UPDATE PFILE METHODS
#=============================================================================

#checks if a profile exists within provided pfile
def profile_exists(pfile, key):
    """Check if a profile exists and has valid data.
    
    Special handling for 'N Z A' (ion species metadata block).
    For all other keys, checks for actual profile data.
    """
    if key == "N Z A":
        # Ion species block: check metadata structure
        if "N Z A" not in pfile:
            return False
        nza = pfile["N Z A"]
        if not isinstance(nza, dict):
            return False
        # Require at least 3 species (2 main + 1 impurity minimum for meaningful impurity handling)
        required_keys = ("N", "Z", "A")
        if not all(k in nza for k in required_keys):
            return False
        n_species = len(np.asarray(nza.get("N", []), dtype=float))
        return n_species >= 2  # At least electrons and main ion
    
    # Regular profile: check key and data
    if key not in pfile:
        return False
    profile = pfile[key]
    if not isinstance(profile, dict) or "data" not in profile:
        return False
    data = np.asarray(profile["data"], dtype=float)
    return data.size > 1 and not np.all(data == 0)

#resample the needed gfile profiles onto the normalized psi coordinate being used in equilbrium scan
def resample_gfile_profiles(gfile, psi_n):
    remapped_profiles = {'psi_n': np.asarray(psi_n, dtype=float)}
    psi_n_gfile = np.asarray(gfile.psi_N, dtype=float)
    remapped_profiles['psi_n_gfile'] = psi_n_gfile
    remapped_profiles['Ip'] = float(gfile.Ip)
    
    jtor_tmp = -np.interp(remapped_profiles['psi_n'], psi_n_gfile, gfile.j_tor_averaged_direct) if gfile.j_tor_averaged_direct[0] < 0 else gfile.j_tor_averaged_direct
    remapped_profiles['j_tor_averaged_direct']  = jtor_tmp

    psi_gfile_len = len(gfile.R_grid)
    psi_gfile = np.linspace(gfile.psi_axis, gfile.psi_boundary, psi_gfile_len)

    remapped_profiles['psi'] = np.interp(remapped_profiles['psi_n'], psi_n_gfile, psi_gfile)
    # Resample all g-file profiles used downstream onto unified psi_n.
    for key in ['R', 'Bp', 'Bt', 'Z', 'a', 'kappa', 'delta']:
        if key in gfile.averages:
            profile_data = np.asarray(gfile.averages[key])
        elif key in gfile.geometry:
            profile_data = np.asarray(gfile.geometry[key])
        else:
            raise AttributeError(f"gfile object does not have attribute '{key}'")
        remapped_profiles[key] = np.interp(remapped_profiles['psi_n'], psi_n_gfile, profile_data)

 
            
    remapped_profiles['dpsi'] = np.gradient(remapped_profiles['psi'], remapped_profiles['psi_n'])
    return remapped_profiles

def update_impurity_species(pfile, resampled_gfile_profiles, psi_n):
    """
    Function to update the density, toroidal velocity, and poloidal velocity of a given ion impurity species
    Renormalizes the profiles to maintain quasineutrality and net charge balance
    Profiles updated: [nz{i}, vtor{i}, vpol{i}]
    """

    num_impurity_species = 0
    if profile_exists(pfile, 'N Z A') and len(pfile['N Z A']['N']) > 2:
        num_impurity_species = len(pfile['N Z A']['N']) - 2
    else:
        # No impurities: mark n_imps = 0 and return early
        resampled_gfile_profiles['n_imps'] = 0
        resampled_gfile_profiles['ZIs'] = []
        return

    # Update impurity profiles
    for i in range(num_impurity_species):
        for key in [f'vtor{i+1}', f'vpol{i+1}', f'nz{i+1}']:
            if profile_exists(pfile, key):
                old_psi_n = np.asarray(pfile.psinorm_for(key), dtype=float)
                tmp = np.interp(psi_n, old_psi_n, pfile[key]['data'])
            else:
                tmp = np.zeros_like(psi_n)
            pfile.set_profile(key, psi_n, tmp)
    
    nz_fracs = []
    for i in range(num_impurity_species): 
        nz_frac = pfile[f'nz{i+1}']['data'] / (np.abs(pfile['ni']['data']) + 1e-14)
        nz_fracs.append(nz_frac)
    
    ZIs = []
    for i in range(num_impurity_species): ZIs.append(pfile['N Z A']['Z'][i])
    Zavg = sum(nz_fracs[i]*ZIs[i] for i in range(len(ZIs)))
    
    # Avoid division by zero on profile array
    Zavg = np.where(np.abs(Zavg) < 1e-14, 1.0, Zavg)
    
    nz_tot = (pfile['ne']['data'] - pfile['ni']['data'] - pfile['nb']['data'])/Zavg

    for i in range(num_impurity_species):
        pfile.set_profile(f'nz{i+1}', psi_n, nz_tot * nz_fracs[i])

    resampled_gfile_profiles['n_imps'] = num_impurity_species
    resampled_gfile_profiles['ZIs'] = ZIs

def update_main_profiles(pfile, psf, psi_n):
    """
    Function to update density and temperature profiles
    for ions and electrons
    Profiles updated: [ne,ni,te,ti]
    """
    ne = pfile['ne']['data']*(psf**(2./3.))
    pfile.set_profile('ne',psi_n, ne)

    ni = pfile['ni']['data']*(psf**(2./3.))
    pfile.set_profile('ni', psi_n, ni)

    te = pfile['te']['data']*(psf**(1./3.))
    pfile.set_profile('te', psi_n, te)

    ti = pfile['ti']['data']*(psf**(1./3.))
    pfile.set_profile('ti', psi_n, ti)

def update_fast_ion_profiles(pfile, psf, psi_n):
    """
    Function to update fast ion (beam) density and fast ion (beam) pressure in provided pfile
    Profiles updated: [nb,pb]
    """
    if profile_exists(pfile, 'nb'):
        nb = pfile['nb']['data']*(psf**(2./3.))
    else:
        nb = np.zeros_like(psi_n)
    pfile.set_profile('nb', psi_n, nb)
    
    if profile_exists(pfile, 'pb'):
        pb = pfile['pb']['data']*psf
    else:
        pb = np.zeros_like(psi_n)
    pfile.set_profile('pb', psi_n, pb)

def update_total_pressure(pfile, resampled_gfile_profiles, psi_n):
    """
    Function to update total pressure according to pressure scaling function that maintains
    constant collisionality
    Profiles updated: [ptot] 
    """
    kb = constants.eV
    num_impurity_species = resampled_gfile_profiles['n_imps']
    
    p_electron = pfile['ne']['data'] * pfile['te']['data'] * kb * 1e20
    nzs = sum(pfile[f'nz{i+1}']['data'] for i in range(num_impurity_species))

    p_ions = (pfile['ni']['data'] + nzs) * pfile['ti']['data'] * kb * 1e20
    p_beam = pfile['pb']['data'] 

    pfile.set_profile('ptot', psi_n, (p_electron + p_ions + p_beam))

def update_rotation_profiles(pfile, resampled_gfile_profiles, psi_n):
    """
    Function that updates rotation measurement arrays
    Uses g-file and p-file profiles to calculate kpol, toroidal rotation, and poloidal rotation
    Profiles updated: [kpol, omeg, omegp]
    """
    if profile_exists(pfile, 'vpol1'):
        kpol = pfile['vpol1']['data']/resampled_gfile_profiles['Bp']
    elif profile_exists(pfile, 'kpol'):
        kpol = np.interp(psi_n, np.asarray(pfile.psinorm_for('kpol'), dtype=float), pfile['kpol']['data'])
    else:
        kpol = np.zeros_like(psi_n)
    pfile.set_profile('kpol', psi_n, kpol)

    if profile_exists(pfile, 'vtor1'):
        omeg = pfile['vtor1']['data']/resampled_gfile_profiles['R']
    elif profile_exists(pfile, 'omeg'):
        omeg = np.interp(psi_n, np.asarray(pfile.psinorm_for('omeg'), dtype=float), pfile['omeg']['data'])
    else:
        omeg = np.zeros_like(psi_n)
    pfile.set_profile('omeg', psi_n, omeg)
        
    if profile_exists(pfile, 'vpol1'):
        omegp_tmp = (pfile['vpol1']['data']*resampled_gfile_profiles['Bt']) / (resampled_gfile_profiles['R']*resampled_gfile_profiles['Bp'])
        omegp = omegp_tmp
    elif profile_exists(pfile, 'omegp'):
        omegp = np.interp(psi_n, np.asarray(pfile.psinorm_for('omegp'), dtype=float), pfile['omegp']['data'])
    else:
        omegp = np.zeros_like(psi_n)
    pfile.set_profile('omegp', psi_n, omegp)

def effective_impurity_density(pfile, num_imps, ZIs):
    """Build an effective nz1-equivalent density for pfile compute methods.

    The pfile compute_diamagnetic_rotations method accepts a single impurity
    density nI with charge Z_imp = ZIs[0]. To preserve multi-impurity charge
    content, map total impurity charge density onto that reference charge:

        nI_eff * Z_ref = sum_i(nz_i * Z_i)
    """
    if num_imps <= 0 or len(ZIs) == 0:
        return None

    z_ref = float(ZIs[0])
    if np.isclose(z_ref, 0.0):
        return None

    charge_density = np.zeros_like(pfile['ni']['data'], dtype=float)
    for i in range(num_imps):
        key = f'nz{i+1}'
        if key in pfile:
            charge_density += np.asarray(pfile[key]['data'], dtype=float) * float(ZIs[i])

    return charge_density / z_ref


def update_diamagnetic_and_rotation_decomposition(pfile, resampled_gfile_profiles):
    """Function to update diamagnetic rotation and rotation decomposition terms
    
    Uses the methods built into eqdsk.py to compute the profiles.
    Profiles updated: [omgpp, ommpp, omepp, omgvb, omgeb, ommvb, omevb, er, omghb]
    """
    psi_n = np.asarray(resampled_gfile_profiles['psi_n'], dtype=float)

    # Seed VxB term in the same way as previous varyped implementation.
    if profile_exists(pfile, 'omeg') and profile_exists(pfile, 'omegp'):
        pfile.set_profile('omgvb', psi_n, pfile['omeg']['data'] + pfile['omegp']['data'])
    else:
        pfile.set_profile('omgvb', psi_n, np.zeros_like(psi_n))

    num_imps = int(resampled_gfile_profiles.get('n_imps', 0))
    ZIs = resampled_gfile_profiles.get('ZIs', [])
    nI_eff = effective_impurity_density(pfile, num_imps, ZIs)

    psi = np.asarray(resampled_gfile_profiles['psi'], dtype=float)
    ti = np.asarray(pfile['ti']['data'], dtype=float) if profile_exists(pfile, 'ti') else None

    # Robust diamagnetic terms from pfile.py.
    pfile.compute_diamagnetic_rotations(psi=psi, nI=nI_eff, TI=ti)

    # Robust decomposition + er + omghb from pfile.py.
    pfile.compute_rotation_decomposition(
        R=np.asarray(resampled_gfile_profiles['R'], dtype=float),
        Bp=np.asarray(resampled_gfile_profiles['Bp'], dtype=float),
        Bt=np.asarray(resampled_gfile_profiles['Bt'], dtype=float),
        psi=psi,
    )


def make_updated_pfile(pfile, resampled_gfile_profiles, psf, psi_n):
    '''
    Update pfile_data with profiles from pfile and gfile
    Profile Descriptions w/ Default units:
        ne = 'Electron density' [10^20/m^3] 
        te = 'Electron temperature' [KeV]
        ni = 'Ion density' [10^20/m^3]
        ti = 'Ion temperature' [KeV]
        nb = 'Fast ion density' [10^20/m^3]
        pb = 'Fast ion pressure' [kPa]
        ptot = 'Total pressure'[kPa]
        omeg = 'Toroidal rotation: VTOR/R' [kRad/s]
        omegp = 'Poloidal rotation: Bt * VPOL / (RBp)' [kRad/s]
        omgvb = 'VxB rotation term in the ExB rotation frequency: OMEG + OMEGP' [kRad/s]
        omgpp = 'Diamagnetic term in the ExB rotation frequency: (P_Carbon)/dpsi / (6*n_Carbon)' [kRad/s]
        omgeb = 'ExB rotation frequency: OMGPP + OMGVB = Er/(RBp)' [kRad/s]
        er = 'Radial electric field from force balance: OMGEB * RBp' [kV/m]
        ommvb = 'Main ion VXB term of Er/RBp, considered a flux function' [kRad/s]
        ommpp = 'Main ion pressure term of Er/RBp, considered a flux function' [kRad/s]
        omevb = 'Electron VXB term of Er/RBp, considered a flux function' [kRad/s]
        omepp = 'Electron pressure term of Er/RBp, considered a flux function' [kRad/s]
        kpol = 'KPOL=VPOL/Bp : V_vector = KPOL*B_vector + OMGEB * PHI_Vector' [kRad/s]
        N Z A = 'N Z A of ION SPECIES' [N = Z = Atomic Number, A = Atomic Mass]
        omghb = 'Hahm-Burrell form for the ExB velocity shearing rate: OMGHB = (RBp)**2/Bt * d (Er/RBp)/dpsi' [kRad/s]
        nzi = 'Density of the ith impurity species' [10^20/m^3]
        vtori = 'Toroidal velocity of the ith impurity species' [km/s]
        vpoli = 'Poloidal velocity of the ith impurity species' [km/s]
    '''
    
    pfile_copy = copy.deepcopy(pfile)

    update_main_profiles(pfile_copy, psf, psi_n)
    update_fast_ion_profiles(pfile_copy, psf, psi_n)
    update_impurity_species(pfile_copy, resampled_gfile_profiles, psi_n)
    update_total_pressure(pfile_copy, resampled_gfile_profiles, psi_n)
    update_rotation_profiles(pfile_copy, resampled_gfile_profiles, psi_n)
    update_diamagnetic_and_rotation_decomposition(pfile_copy, resampled_gfile_profiles)

    return pfile_copy


#=============================================================================
#                       RUN VARYPED SCAN METHODS
#=============================================================================

def equilibrium_scan(t_object, pfile, gfile, scaling_values, path_to_output = None, filenames = None):
    r"""! Run equilibrium scan scaling pressure and current density from TokaMaker instance

    @param t_object Grad-Shafranov solver object
    @param pfile Osborne p-File object storing kinetic and rotational profiles
    @param gfile EFIT EQDSK g-File object storing equilibrium solution data
    @param scaling_values array-like object of non-decreasing positive scalar values
    @param path_to_output (str) path to directory where outputs will be created, defaults to cwd
    @param filenames (str) name for output files (typically of the form #####.####, where # are integers), defaults to 10000.0000

    returns:
    List: List object containing output filenames, current targets, work targets, and scaling factors used.
    """
    import pprint

    results = []

    #ensure scaling_values is array-like
    if type(scaling_values) is int or type(scaling_values) is float:
        scaling_values = [scaling_values]

   #Extract psi_n from pfile 
    if not profile_exists(pfile, 'ptot'):
        raise ValueError("pfile must contain 'ptot' (total pressure) profile")
    psi_n_pfile = pfile.psinorm_for('ptot')
    if psi_n_pfile is None:
        raise ValueError("Could not extract psinorm for 'ptot' profile")
    psi_n_pfile = np.asarray(psi_n_pfile, dtype=float)
    
    # Remap all pfile profiles onto ptot's psi_n 
    # ensures all profiles share a common coordinate system
    pfile_copy = pfile.remap(psi_n_pfile)
    
    # Get psi_n from gfile 
    psi_n_gfile = np.asarray(gfile.psi_N, dtype=float)
    
    # Select unified psi_n as one with higher resolution
    # If resolutions are equal, prefer gfile's coordinate system
    if len(psi_n_pfile) > len(psi_n_gfile):
        psi_n = psi_n_pfile
        print(f"[equilibrium_scan] Using pfile psi_n: {len(psi_n)} points (pfile has higher resolution than gfile: {len(psi_n_gfile)})")
    else:
        psi_n = psi_n_gfile
        pfile_copy = pfile_copy.remap(psi_n)
        print(f"[equilibrium_scan] Using gfile psi_n: {len(psi_n)} points (gfile has greater than or equal resolution to pfile: {len(psi_n_pfile)})")
    
    # Resample all needed gfile profiles onto unified psi_n 
    gfile_profiles = resample_gfile_profiles(gfile, psi_n)
    refit_jtor = gfile_profiles['j_tor_averaged_direct']
    # Get Ip target from gfile 
    Ip_target = gfile_profiles['Ip']

    R0 = gfile_profiles['R'][0]
    Z0 = gfile_profiles['Z'][0]
    a = gfile_profiles['a'][-1]
    kappa = gfile_profiles['kappa'][0]
    delta = gfile_profiles['delta'][0]
    
    t_object.init_psi(R0, Z0, a, kappa, delta)

    for scale_p in scaling_values:

        for scale_j in scaling_values:
            
            #get scaled pressure profile, psf and pprime from scaled pressure profile
            ptot_scaled, psf, _, W_target, W_final, _ = scale_pressure(t_object, psi_n, pfile_copy['ptot']['data'], scale_p)
                      
            ptot_scaled_prime = np.gradient(ptot_scaled, psi_n)

            #make tokamaker pprime_profile
            pprime_profile = {'type': 'linterp',
                              'y': ptot_scaled_prime,
                              'x': psi_n
                              }
            
            #get scaled current profile
            jtor_scaled, _, Ip_final = scale_current(t_object, psi_n, refit_jtor, scale_j)
    
            #make tokamaker ffprime_profile 
            ffprime_profile = {'type': 'jphi-linterp',
                               'y': jtor_scaled,
                               'x': psi_n
                               }
            
            #set pax_target (pressure on axis target) to p_scaled[0]
            pax_target = ptot_scaled[0]

            #set t_object targets as pax_target and Ip_target
            t_object.set_targets(Ip=Ip_target, pax=pax_target)

            #set t_object profiles to pprime_profile and ffprime_profile
            t_object.set_profiles(ffp_prof=ffprime_profile, pp_prof=pprime_profile)

            #solve equilibra
            t_object.solve()

            if path_to_output is None:
                path_to_output = os.getcwd()

            if filenames is None:
                filenames = '10000.0000'
            #save output gfile by calling save_eqdsk
            path_to_output_gfile = os.path.join(path_to_output, f'g{filenames}_{scale_p}_{scale_j}' + '.geqdsk')
            t_object.save_eqdsk(path_to_output_gfile, 
                                nr = len(psi_n),
                                nz = len(psi_n),
                                truncate_eq=False,
                                lcfs_pad=1e-2)
            
            #call make_updated_pfile with psf and scaled_ptot
                #make_updated_pfile deepcopies provided pfile and updates every profile within 
                #used newly created gfile to update pfile
            output_gfile = eqdsk.read_geqdsk(path_to_output_gfile)
            output_gfile_profiles = resample_gfile_profiles(output_gfile, psi_n)
            output_pfile = make_updated_pfile(pfile_copy, output_gfile_profiles, psf, psi_n)

            path_to_output_pfile = os.path.join(path_to_output, f'p{filenames}_{scale_p}_{scale_j}')
            output_pfile.save(path_to_output_pfile)
            print(f'Saving pFile: {path_to_output_pfile}' )
            results.append({'scale_p':scale_p, 
                            'scale_j': scale_j,
                            'gfile_name': f'g{filenames}_{scale_p}_{scale_j}' + '.geqdsk',
                            'pfile_name': f'p{filenames}_{scale_p}_{scale_j}',
                            'Ip_target': Ip_target,
                            'Ip_final': Ip_final,
                            'pax_target': f'{float(pax_target)} kPa',
                            'W_target': W_target,
                            'W_final': W_final
                            })

    for i in range(len(results)):
        pprint.pprint(f'Run ({results[i]['scale_p']},{results[i]['scale_j']}) Quick Results: {results[i]}')
    return results


def plot_pressure_results(results, scaling_values, path_to_results=None, pfile_0 = None, gfile_0 = None):
    from matplotlib import pyplot as plt
    from matplotlib.colors import Normalize
    import matplotlib.cm as cm
    scaling_values = np.asarray(scaling_values)
    norm = Normalize(vmin=(scaling_values.min()**2 - scaling_values.max()**2), vmax=(scaling_values.max()))
    cmap = cm.plasma

    if path_to_results is None:
        path_to_results = os.getcwd()    

    for result in results:
        color = cmap(norm(float(result['scale_p'])**2 - float(result['scale_j'])**2))    
        path_to_result_pfile = os.path.join(path_to_results, result['pfile_name'])
        result_pfile = eqdsk.read_pfile(path_to_result_pfile)
        plt.plot(result_pfile.psinorm_for('ptot'), result_pfile.ptot, color = color, label = f'Pres Scaled by: {result['scale_p']}; Curr Scaled by: {result['scale_j']}')
    plt.plot(label = "(Pressure Scale Factor, Current Scale Factor)")
    if pfile_0 is not None:
        plt.plot(pfile_0.psinorm_for('ptot'), pfile_0.ptot, color = 'black', label = 'Input Pressure Profile', linestyle = '--', linewidth = 4)
    else: print("Original pFile Not Provided, Initial Pressure Profile Not Plotted")

    plt.xlabel('Normalized ψ')
    plt.ylabel('Total Pressure [kPa]')
    plt.title('Scaled Pressures')
    plt.legend()
    plt.show()

    for result in results:
        color = cmap(norm(float(result['scale_p'])**2 - float(result['scale_j'])**2))    
        path_to_result_pfile = os.path.join(path_to_results, result['pfile_name'])
        result_pfile = eqdsk.read_pfile(path_to_result_pfile)
        plt.plot(result_pfile.psinorm_for('ptot'), result_pfile.derivative_for('ptot'), color = color, label = f'Pres Scaled by: {result['scale_p']}; Curr Scaled by: {result['scale_j']}')

    if pfile_0 is not None:
        plt.plot(pfile_0.psinorm_for('ptot'), pfile_0.derivative_for('ptot'),  color = 'black', label = 'Input Pressure Derivative Profile', linestyle = '--', linewidth = 4)
    else: print("Original pFile Not Provided, Initial P\' Profile Not Plotted")
    
    plt.xlabel('Normalized ψ')
    plt.ylabel('P\' Profile')
    plt.title('Scaled Pressure Derivatives')
    plt.legend()
    plt.show()

def plot_current_results(results, scaling_values, path_to_results=None, pfile_0 = None, gfile_0 = None):
    
    from matplotlib import pyplot as plt
    from matplotlib.colors import Normalize
    import matplotlib.cm as cm
    scaling_values = np.asarray(scaling_values)
    norm = Normalize(vmin=(scaling_values.min()**2 - scaling_values.max()**2), vmax=(scaling_values.max()))
    cmap = cm.plasma

    if path_to_results is None:
        path_to_results = os.getcwd()    

    for result in results:
        color = cmap(norm(float(result['scale_p'])**2 - float(result['scale_j'])**2))    
        path_to_result_gfile = os.path.join(path_to_results, result['gfile_name'])
        result_gfile = eqdsk.read_geqdsk(path_to_result_gfile)
        plt.plot(result_gfile.psi_N, -result_gfile.j_tor_averaged_direct, color = color, label = f'Pres Scaled by: {result['scale_p']}; Curr Scaled by: {result['scale_j']}')

    if gfile_0 is not None:
        plt.plot(gfile_0.psi_N, -gfile_0.j_tor_averaged_direct, color = 'black', label = 'Input Current Density Profile', linestyle = '--', linewidth = 4)
    else: print("Original gFile Not Provided, Initial Current Density Profile Not Plotted")

    plt.xlabel('Normalized ψ')
    plt.ylabel('Current Density [Ampere/m^2]')
    plt.title('Scaled Current Densities')
    plt.legend()
    plt.show()

def plot_varyped_results(results, scaling_values, path_to_results=None, pfile_0 = None, gfile_0 = None):
    '''! Plots VARYPED results from an equilibrium scan.'''
    
    plot_pressure_results(results, scaling_values, path_to_results, pfile_0, gfile_0)
    plot_current_results(results, scaling_values, path_to_results, pfile_0, gfile_0)