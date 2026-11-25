#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Run the linearized 2D Rayleigh-Taylor example (MUG_rt_2d) with the exponential profile over a range of in-plane
# field strengths, B_0/B_c, plus a case with an out-of-plane field, and compare the measured growth rates with the
# incompressible, ideal result of S. Chandrasekhar, Hydrodynamic and Hydromagnetic Stability (1961), sections 92
# and 97. The profile, box (1 m x 1 m), mesh and dissipation are set here; the other options come from oft.in.
#
# Each case starts from the velocity-only m=1 eigenmode, so the velocity grows as cosh(gamma*t). The growth rate
# is obtained by fitting E_k = A*cosh^2(gamma*t) to the whole kinetic energy history, and, as a check, from the
# slope of ln(E_k)/2 over the part of the run with gamma*t > 3 (where tanh(gamma*t) > 0.995).
#
# Usage (from the directory containing MUG_rt_2d, oft.in, oft_in.xml, and oft_in_mpi.xml for --np > 1):
#   python run_scan.py [--b_frac 0 0.25 ...] [--by0 1.E-2] [--np 1] [--plot]
import argparse
import os
import re
import shutil
import subprocess
import numpy as np
from scipy.optimize import curve_fit
from OpenFUSIONToolkit.io import histfile

# Physical constants (as in mhd_utils.F90)
mu0 = 4.e-7*np.pi
proton_mass = 1.672621777E-27
elec_charge = 1.602176565E-19


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


def parallel_settings(text, nprocs):
    '''Adapt an input file for a run on nprocs MPI processes

    Returns the input text and the solver XML file to use. With more than one process the base mesh is distributed
    over the processes (nlevels=2, a transfer level, not a refinement), the overlapping block-Jacobi preconditioner
    in oft_in_mpi.xml is used, and ittarget is raised to at least 200, since that preconditioner needs more linear
    iterations per step (see the example documentation).
    '''
    if nprocs <= 1:
        return text, 'oft_in.xml'
    text = set_option(text, 'nlevels', '2')
    text = set_option(text, 'ittarget', str(max(200, int(read_option(text, 'ittarget')[0]))))
    return text, 'oft_in_mpi.xml'


def run_command(exe, xml_file, nprocs):
    '''Command line for one run of the example'''
    cmd = [exe, 'oft.in', xml_file]
    if nprocs > 1:
        cmd = ['mpirun', '-np', str(nprocs)] + cmd
    return cmd


class rt_case:
    '''Parameters of the Rayleigh-Taylor case defined by an input file (SI units)'''
    def __init__(self, text):
        self.n_0 = read_option(text, 'n_0')[0]
        self.L_n = read_option(text, 'L_n')[0]
        self.g = read_option(text, 'g')[0]
        rscale = read_option(text, 'rscale')
        self.k = 2.0*np.pi/rscale[0]
        self.K2 = self.k**2 + (np.pi/rscale[1])**2 + 1.0/(4.0*self.L_n**2)
        self.B_c = np.sqrt(mu0*proton_mass*self.n_0*abs(self.g/self.L_n)/self.K2)

    def gamma(self, b_frac):
        '''Incompressible, ideal growth rate [1/s] of the m=1 mode for B_0 = b_frac*B_c (0 if stable)'''
        gamma2 = self.k**2*(self.g/self.L_n)/self.K2*(1.0 - b_frac**2)
        return np.sqrt(max(gamma2, 0.0))


def fit_growth(hist_file):
    '''Growth rate [1/s] from the kinetic energy history

    Returns (gamma_cosh, gamma_late): a fit of E_k = A*cosh^2(gamma*t) to the whole history, and the slope of
    ln(E_k)/2 over gamma_cosh*t > 3. Both are None if the kinetic energy does not grow (stable case).
    '''
    hist = histfile(hist_file)
    t = np.asarray(hist['time'], dtype=float)
    ek = np.asarray(hist['ven'], dtype=float)
    if ek[-1] < 10.0*ek[0]:
        return None, None
    log_ek = np.log(ek)
    # Initial guess from the slope over the second half of the run
    gamma0 = 0.5*np.polyfit(t[t.shape[0]//2:], log_ek[t.shape[0]//2:], 1)[0]

    def model(tt, log_a, gamma):
        return log_a + 2.0*np.log(np.cosh(gamma*tt))
    popt, _ = curve_fit(model, t, log_ek, p0=[log_ek[0], gamma0])
    gamma_cosh = popt[1]
    mask = gamma_cosh*t > 3.0
    gamma_late = None
    if np.count_nonzero(mask) > 5:
        gamma_late = 0.5*np.polyfit(t[mask], log_ek[mask], 1)[0]
    return gamma_cosh, gamma_late


def run_case(exe, template, case, b_frac, by0, nprocs):
    '''Run one case in its own directory and return the measured growth rates [1/s]'''
    case_dir = 'b_{0:.3f}'.format(b_frac)
    if by0 != 0.0:
        case_dir += '_by_{0:.2E}'.format(by0)
    os.makedirs(case_dir, exist_ok=True)
    # Run until gamma*t = 6 (theory), or for 0.6 ms if stable, keeping gamma*dt <= 0.025
    dt = read_option(template, 'dt')[0]
    gamma_ref = case.gamma(b_frac)
    t_end = 6.0/gamma_ref if gamma_ref > 0.0 else 6.E-4
    nsteps = int(np.ceil(t_end/dt))
    text = set_option(template, 'linear', 'T')
    text = set_option(text, 'B_0', '{0:.6E}'.format(b_frac*case.B_c))
    text = set_option(text, 'by0', '{0:.6E}'.format(by0))
    text = set_option(text, 'nsteps', str(nsteps))
    text = set_option(text, 'rst_freq', str(nsteps))
    text, xml_file = parallel_settings(text, nprocs)
    with open(os.path.join(case_dir, 'oft.in'), 'w') as fid:
        fid.write(text)
    # Copy the XML file, since long paths are truncated by the executable
    shutil.copy(xml_file, case_dir)
    cmd = run_command(exe, xml_file, nprocs)
    with open(os.path.join(case_dir, 'run.log'), 'w') as log:
        result = subprocess.run(cmd, cwd=case_dir, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print('  WARNING: run for B_0/B_c = {0:.3f} stopped early (see {1})'.format(b_frac, case_dir))
    return fit_growth(os.path.join(case_dir, 'xmhd_2d.hist'))


def main():
    parser = argparse.ArgumentParser(description='Rayleigh-Taylor growth rate vs in-plane field strength')
    parser.add_argument('--b_frac', type=float, nargs='+', default=[0.0, 0.25, 0.5, 0.75, 0.9, 1.1],
                        help='In-plane field strengths, B_0/B_c')
    parser.add_argument('--by0', type=float, default=1.E-2,
                        help='Out-of-plane field [T] for the extra B_0 = 0 case (0 to skip)')
    parser.add_argument('--ni', type=int, default=8, help='Cells in each direction')
    parser.add_argument('--order', type=int, default=3, help='FE order')
    parser.add_argument('--dt', type=float, default=2.E-6, help='Time step [s]')
    parser.add_argument('--V1', type=float, default=1.0, help='Seed amplitude [m/s]')
    parser.add_argument('--nu', type=float, default=1.0, help='Viscosity [m^2/s]')
    parser.add_argument('--eta', type=float, default=0.1,
                        help='Resistivity [m^2/s] (resistive slippage raises gamma near B_c, see the example docs)')
    parser.add_argument('--np', type=int, default=1, help='Number of MPI processes per run')
    parser.add_argument('--exe', default='./MUG_rt_2d', help='Path to the MUG_rt_2d executable')
    parser.add_argument('--plot', action='store_true', help='Plot growth rate vs B_0/B_c (rt_2d_scan.png)')
    options = parser.parse_args()
    with open('oft.in', 'r') as fid:
        template = fid.read()
    # Scan settings: small, well resolved runs with fixed time step and weak dissipation (see the example docs)
    template = set_option(template, 'profile', "'exp'")
    template = set_option(template, 'rscale', '1.,1.,0.')
    template = set_option(template, 'shift', '-0.5,-0.5,0.')
    template = set_option(template, 'ni', '{0},{0},0'.format(options.ni))
    template = set_option(template, 'order', str(options.order))
    template = set_option(template, 'dt', '{0:.6E}'.format(options.dt))
    template = set_option(template, 'V1', '{0:.6E}'.format(options.V1))
    template = set_option(template, 'plot_run', 'F')
    template = set_option(template, 'nu', '{0:.6E}'.format(options.nu))
    template = set_option(template, 'eta', '{0:.6E}'.format(options.eta))
    for name, value in (('chi', '0.'), ('D_diff', '0.'), ('lin_tol', '1.E-10'),
                        ('nl_tol', '1.E-8'), ('use_mfnk', 'F'), ('pre_freq', '1'), ('ittarget', '1000'),
                        ('timestep_cn', 'T')):
        template = set_option(template, name, value)
    case = rt_case(template)
    exe = os.path.abspath(options.exe)
    print('B_c = {0:.4E} T'.format(case.B_c))
    print('  B_0/B_c  by0 [T]    gamma (theory)  gamma (cosh^2 fit)  ratio    gamma (gamma*t>3)  ratio')
    cases = [(b_frac, 0.0) for b_frac in options.b_frac]
    if options.by0 != 0.0:
        cases.append((0.0, options.by0))
    results = []
    for b_frac, by0 in cases:
        gamma_cosh, gamma_late = run_case(exe, template, case, b_frac, by0, options.np)
        theory = case.gamma(b_frac)
        if gamma_cosh is None:
            print('  {0:7.3f}  {1:9.2E}  {2:14.4E}  {3:>18s}'.format(b_frac, by0, theory, 'stable'))
            results.append([b_frac, by0, theory, 0.0])
            continue
        line = '  {0:7.3f}  {1:9.2E}  {2:14.4E}  {3:18.4E}  {4:7.4f}'.format(b_frac, by0, theory, gamma_cosh,
                                                                          gamma_cosh/theory)
        if gamma_late is not None:
            line += '  {0:17.4E}  {1:7.4f}'.format(gamma_late, gamma_late/theory)
        print(line)
        results.append([b_frac, by0, theory, gamma_cosh])
    if options.plot:
        import matplotlib.pyplot as plt
        results = np.asarray(results)
        b_fine = np.linspace(0.0, 1.2*max(1.0, max(options.b_frac)), 200)
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        ax.plot(b_fine, [case.gamma(b) for b in b_fine], 'k-', label='Theory (incompressible, ideal)')
        mask = results[:, 1] == 0.0
        ax.plot(results[mask, 0], results[mask, 3], 'ro', label=r'MUG, $B_0 \hat{x}$')
        if np.any(~mask):
            ax.plot(results[~mask, 0], results[~mask, 3], 'bs', fillstyle='none', markersize=10,
                    label=r'MUG, $b_{y0} \hat{y}$')
        ax.set_xlabel(r'$B_0 / B_c$')
        ax.set_ylabel(r'$\gamma$ [s$^{-1}$]')
        ax.set_ylim(bottom=0.0)
        ax.grid(True)
        ax.legend()
        fig.tight_layout()
        fig.savefig('rt_2d_scan.png', dpi=100)


if __name__ == '__main__':
    main()
