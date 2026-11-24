#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Run the gravitational instability example (MUG_gyro_grav) over a range of beta and compare the measured
# growth rates with the local dispersion relation of Ferraro and Jardin, Phys. Plasmas 13, 092101 (2006),
# Eq. 15, and the result of Roberts and Taylor (Eq. 11).
#
# Usage (from the directory containing MUG_gyro_grav, oft.in, and oft_in.xml):
#   python run_scan.py [--beta 1.E-4 1.E-3 ...] [--no-gyro] [--hall] [--eps EPS] [--plot]
#
# The growth rate is fit over the exponential growth phase of the kinetic energy (see growth_rate). Use
# --exe ./MUG_gyro_grav_lag for the Lagrange variant.
import argparse
import os
import re
import shutil
import subprocess
import numpy as np
from OpenFUSIONToolkit.io import histfile

# Physical constants (as in mhd_utils.F90)
mu0 = 4.e-7*np.pi
proton_mass = 1.672621777E-27
elec_charge = 1.602176565E-19
GAMMA = 5.0/3.0  # Adiabatic index


def read_option(text, name):
    '''Return the value(s) of a namelist entry in an input file as a list of floats'''
    match = re.search(r'^\s*{0}\s*=\s*([^!\n]+)'.format(name), text, re.MULTILINE | re.IGNORECASE)
    if match is None:
        raise ValueError('Option "{0}" not found'.format(name))
    return [float(val.replace('d', 'e').replace('D', 'e')) for val in match.group(1).replace(',', ' ').split()]


def set_option(text, name, value):
    '''Set a namelist entry in an input file, keeping any trailing comment'''
    pattern = re.compile(r'^(\s*{0}\s*=\s*)([^!\n]*?)(\s*(!.*)?)$'.format(name), re.MULTILINE | re.IGNORECASE)
    if pattern.search(text) is None:
        raise ValueError('Option "{0}" not found'.format(name))
    return pattern.sub(lambda m: m.group(1) + value + m.group(3), text, count=1)


class gravity_case:
    '''Normalized parameters of the gravitational instability case defined by an input file

    Quantities are in the paper's units: density n_bar, field B_bar, time 1/Omega_i, length d_i,
    and temperature (pressure over density) B_bar^2/(mu0 n_bar). The midplane density is n0 = 1.
    '''
    def __init__(self, text):
        self.n_bar = read_option(text, 'n_bar')[0]
        self.B_bar = read_option(text, 'B_bar')[0]
        self.g = read_option(text, 'g')[0]
        self.Ln = read_option(text, 'Ln')[0]
        self.I0 = read_option(text, 'I0')[0]
        self.omega_i = elec_charge*self.B_bar/proton_mass
        self.d_i = np.sqrt(proton_mass/(mu0*self.n_bar*elec_charge**2))
        rscale = read_option(text, 'rscale')
        self.kx = 2.0*np.pi*self.d_i/rscale[0]
        self.kz = np.pi*self.d_i/rscale[2]

    def dispersion(self, beta, gyro=True, hall=True):
        '''Growth rate [Omega_i] from the local dispersion relation, Eq. 15 (k rho_i << 1)

        Without gyroviscosity nu = 0; without the Hall term the 1/Omega term is dropped.
        '''
        k = np.sqrt(self.kx**2 + self.kz**2)
        T = beta*self.I0**2/2.0
        gbar = self.g*self.Ln/self.I0**2
        bbar = GAMMA*beta/2.0
        Omega = self.I0/(k*np.sqrt(self.g*self.Ln))
        nu = 0.0
        if gyro:
            nu = (T/self.I0)*k/(2.0*np.sqrt(self.g*self.Ln))
        b = 2.0*nu*(1.0 + gbar + bbar)*(1.0 + bbar) + 2.0*nu*gbar
        if hall:
            b += 1.0/Omega
        roots = np.roots([1.0 + bbar, -b, 1.0 + gbar + bbar])
        return max(0.0, np.max(roots.imag))*self.kx*np.sqrt(self.g/self.Ln)/k

    def roberts_taylor(self, beta, hall=True):
        '''Growth rate [Omega_i] from Roberts and Taylor, Eq. 11 (low beta, evaluated at k = kx)'''
        T = beta*self.I0**2/2.0
        Omega = self.I0/(self.kx*np.sqrt(self.g*self.Ln))
        nu = (T/self.I0)*self.kx/(2.0*np.sqrt(self.g*self.Ln))
        b = 2.0*nu
        if hall:
            b += 1.0/Omega
        roots = np.roots([1.0, b, 1.0])
        return max(0.0, np.max(roots.imag))*np.sqrt(self.g/self.Ln)


def smoothed_slope(t, y, width):
    '''Local slope dy/dt from least-squares line fits over `width` neighboring points'''
    n = len(t)
    half = max(1, width//2)
    slope = np.zeros(n)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        slope[i] = np.polyfit(t[lo:hi], y[lo:hi], 1)[0]
    return slope


def exponential_window(t, y, tol=0.1, min_frac=0.1):
    '''Indices [i, j) of the longest window where the local growth rate d(y)/dt is positive and constant to
    within +/- `tol` (relative), and that spans at least `min_frac` of the run, or None if there is none'''
    n = len(t)
    slope = smoothed_slope(t, y, max(5, n//25))
    min_len = max(5, int(min_frac*n))
    best = None
    for i in range(n):
        smin = slope[i]
        smax = slope[i]
        if smin <= 0.0:
            continue
        j = i + 1
        while j < n:
            smin = min(smin, slope[j])
            smax = max(smax, slope[j])
            if (smin <= 0.0) or (smax - smin > tol*(smax + smin)):
                break
            j += 1
        if (j - i >= min_len) and ((best is None) or (j - i > best[1] - best[0])):
            best = (i, j)
    return best


def growth_rate(hist_file, omega_i, fit_frac=None, tol=0.1):
    '''Growth rate [Omega_i] from the kinetic energy, gamma = (1/2) d ln(E_k)/dt

    By default ln(E_k) is fit over the exponential growth phase, the longest part of the run where its slope
    is constant to within +/- `tol`, which excludes the initial transient and the late saturation or blow-up.
    If `fit_frac` is set, the last `fit_frac` of the run is fit instead (the previous behavior).
    Returns the growth rate [Omega_i] and the start and end of the fit window [1/Omega_i], all NaN if no
    exponential phase is found.
    '''
    hist = histfile(hist_file)
    t = np.asarray(hist['time'], dtype=float)
    ek = np.asarray(hist['ven'], dtype=float)
    mask = ek > 0.0
    t = t[mask]
    y = np.log(ek[mask])
    if fit_frac is None:
        window = exponential_window(t, y, tol=tol)
        if window is None:
            return np.nan, np.nan, np.nan
        i, j = window
    else:
        i = np.searchsorted(t, t[-1] - fit_frac*(t[-1] - t[0]))
        j = len(t)
    slope = np.polyfit(t[i:j], y[i:j], 1)[0]
    return 0.5*slope/omega_i, t[i]*omega_i, t[j-1]*omega_i


def run_case(exe, template, beta, gyro, hall, nprocs, eps=None, run=True):
    '''Run one case in its own directory (unless `run` is False) and return the name of the directory'''
    case_dir = 'beta_{0:.3E}'.format(beta)
    if not gyro:
        case_dir += '_nogyro'
    if hall:
        case_dir += '_hall'
    if not run:
        return case_dir
    os.makedirs(case_dir, exist_ok=True)
    text = set_option(template, 'beta', '{0:.6E}'.format(beta))
    text = set_option(text, 'nu_gyro', '1.' if gyro else '-1.')
    text = set_option(text, 'xmhd_hall', 'T' if hall else 'F')
    if eps is not None:
        text = set_option(text, 'eps', '{0:.6E}'.format(eps))
    with open(os.path.join(case_dir, 'oft.in'), 'w') as fid:
        fid.write(text)
    # Copy the XML file, since long paths are truncated by the executable
    shutil.copy('oft_in.xml', case_dir)
    cmd = [exe, 'oft.in', 'oft_in.xml']
    if nprocs > 1:
        cmd = ['mpirun', '-np', str(nprocs)] + cmd
    with open(os.path.join(case_dir, 'run.log'), 'w') as log:
        result = subprocess.run(cmd, cwd=case_dir, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print('  WARNING: run for beta = {0:.3E} stopped early (see {1})'.format(beta, case_dir))
    return case_dir


def main():
    parser = argparse.ArgumentParser(description='Gravitational instability growth rate vs beta')
    parser.add_argument('--beta', type=float, nargs='+', default=[1.E-4, 5.E-4, 1.E-3, 1.5E-3, 2.E-3, 2.5E-3, 3.E-3])
    parser.add_argument('--no-gyro', action='store_true', help='Disable gyroviscosity (nu_gyro<0)')
    parser.add_argument('--hall', action='store_true',
                        help='Include the Hall term (currently shows a spurious faster instability, see the docs)')
    parser.add_argument('--np', type=int, default=1, help='Number of MPI processes per run')
    parser.add_argument('--exe', default='./MUG_gyro_grav', help='Path to the MUG_gyro_grav executable')
    parser.add_argument('--eps', type=float, default=None,
                        help='Relative amplitude of the density perturbation (default: "eps" in oft.in)')
    parser.add_argument('--fit-frac', type=float, default=None,
                        help='Fit the last FIT_FRAC of each run instead of the exponential growth phase')
    parser.add_argument('--fit-tol', type=float, default=0.1,
                        help='Relative variation of the local growth rate allowed in the exponential phase')
    parser.add_argument('--no-run', action='store_true', help='Only fit the results of a previous scan')
    parser.add_argument('--plot', action='store_true', help='Plot growth rate vs beta (gyro_grav_scan.png)')
    options = parser.parse_args()
    gyro = not options.no_gyro
    hall = options.hall
    with open('oft.in', 'r') as fid:
        template = fid.read()
    case = gravity_case(template)
    exe = os.path.abspath(options.exe)
    print('  beta       gamma (MUG)  gamma (Eq. 15)  gamma (Eq. 11)   fit window   [Omega_i]')
    results = []
    for beta in options.beta:
        case_dir = run_case(exe, template, beta, gyro, hall, options.np, eps=options.eps, run=not options.no_run)
        gamma, t_start, t_end = growth_rate(os.path.join(case_dir, 'xmhd.hist'), case.omega_i,
                                            fit_frac=options.fit_frac, tol=options.fit_tol)
        theory = case.dispersion(beta, gyro=gyro, hall=hall)
        results.append([beta, gamma, theory])
        print('  {0:.3E}  {1:12.4E}  {2:14.4E}  {3:14.4E}   {4:5.0f}-{5:5.0f}'.format(
            beta, gamma, theory, case.roberts_taylor(beta, hall=hall), t_start, t_end))
        if np.isnan(gamma):
            print('  WARNING: no exponential growth phase found for beta = {0:.3E}'.format(beta))
    if options.plot:
        import matplotlib.pyplot as plt
        results = np.asarray(results)
        beta_fine = np.linspace(0.0, 1.2*max(options.beta), 200)
        fig, ax = plt.subplots(1, 1)
        ax.plot(beta_fine, [case.dispersion(b, gyro=gyro, hall=hall) for b in beta_fine], 'k-', label='Eq. 15')
        if gyro:
            ax.plot(beta_fine, [case.roberts_taylor(b, hall=hall) for b in beta_fine], 'k--',
                    label='Roberts and Taylor (Eq. 11)')
        ax.plot(results[:, 0], results[:, 1], 'ro', label='MUG')
        ax.set_xlabel(r'$\beta$')
        ax.set_ylabel(r'$\gamma / \Omega_i$')
        ax.set_ylim(bottom=0.0)
        ax.grid(True)
        ax.legend()
        fig.savefig('gyro_grav_scan.png', dpi=150)


if __name__ == '__main__':
    main()
