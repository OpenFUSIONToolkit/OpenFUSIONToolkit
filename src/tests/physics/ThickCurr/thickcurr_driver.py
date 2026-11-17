#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Python driver for ThickCurr tests, run as a separate process (optionally using MPI) by the test scripts:
#
#   [mpirun -np N] python thickcurr_driver.py config.json
#
# Results are written to "thickcurr_results.json" by the head process.
import os
import sys
import json
import numpy
import h5py
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
import oft_testing  # noqa: F401 (sets up path to OFT python package)
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.ThickCurr import ThickCurr
from OpenFUSIONToolkit.ThinCurr.sensor import Mirnov, flux_loop
from OpenFUSIONToolkit._interface import oftpy_dump_cov

# Probe locations and directions for the box decay mode (A = sin(pi*x)*sin(pi*y) z_hat)
BOX_MIRNOVS = (
    ((0.3, 0.2, 0.5), (1.0, 0.0, 0.0)),
    ((0.7, 0.35, 0.3), (0.0, 1.0, 0.0)),
    ((0.15, 0.6, 0.8), (1.0, 1.0, 0.0)),
)
# Rectangular flux loop in the x=0.3 plane (normal +x), 0.1 < y < 0.5, 0.2 < z < 0.8
BOX_FLOOP = numpy.array([[0.3, 0.1, 0.2], [0.3, 0.5, 0.2], [0.3, 0.5, 0.8], [0.3, 0.1, 0.8], [0.3, 0.1, 0.2]])


def box_mode(pts):
    '''Vector potential of the slowest diffusion mode of the unit cube, A = sin(pi*x)*sin(pi*y) z_hat'''
    vals = numpy.zeros(pts.shape)
    vals[:, 2] = numpy.sin(numpy.pi*pts[:, 0])*numpy.sin(numpy.pi*pts[:, 1])
    return vals


def main(config):
    myOFT = OFT_env(nthreads=config.get('nthreads', 1), quiet=True)
    tc_model = ThickCurr(myOFT)
    #---Setup model
    setup_args = {
        'order': config.get('order', 2),
        'eta_values': config['eta'],
        'xml_filename': config.get('xml_filename', None),
        'pmg': config.get('pmg', False),
        'dt': config.get('dt', 1.E-3),
        'theta': config.get('theta', 1.0),
    }
    if config.get('in_memory', False):
        with h5py.File(config['mesh_file'], 'r') as h5_file:
            r = numpy.asarray(h5_file['mesh/R'])
            lc = numpy.asarray(h5_file['mesh/LC'])-1
            reg = numpy.asarray(h5_file['mesh/REG'])
        tc_model.setup_model(r=r, lc=lc, reg=reg, **setup_args)
    else:
        tc_model.setup_model(mesh_file=config['mesh_file'], **setup_args)
    results = {'nprocs': myOFT.nprocs}
    #---Eigenvalue solve
    if config['case'] == 'eigs':
        tau = tc_model.get_eigs(config.get('nev', 5), shift=config.get('eig_shift', None), save_plot=True)
        results['tau'] = tau.tolist()
    #---Time advance of the slowest box mode
    elif config['case'] == 'decay':
        sensors = [Mirnov(pt, norm, 'BP{0}'.format(i+1), dx=None) for i, (pt, norm) in enumerate(BOX_MIRNOVS)]
        sensors.append(flux_loop(BOX_FLOOP, 'FLOOP'))
        sensor_obj = tc_model.setup_sensors(sensors)
        e0 = tc_model.set_field(box_mode)
        results['signals0'] = tc_model.eval_sensors(sensor_obj).tolist()
        if config.get('eig_shift', None) is not None:
            # Eigenvalue solve with a different shift (time advance operators must be restored)
            tc_model.get_eigs(1, shift=config['eig_shift'])
        t_final = tc_model.run_td(nsteps=config['nsteps'], plot_freq=config.get('plot_freq', 10), sensor_obj=sensor_obj)
        results['t_final'] = t_final
        results['energy'] = [e0, tc_model.get_energy()]
        results['signals1'] = tc_model.eval_sensors(sensor_obj).tolist()
        tc_model.plot_td(plot_freq=config.get('plot_freq', 10))
        if myOFT.head_proc:
            sensor_hist = tc_model.load_history('floops.hist')
            results['hist_signals1'] = [float(sensor_hist[name][-1]) for name in sensor_obj['names']]
            results['hist_t1'] = float(sensor_hist['time'][-1])
    else:
        raise ValueError('Unknown test case "{0}"'.format(config['case']))
    #---Save results and plot metadata
    if myOFT.head_proc:
        tc_model.build_XDMF()
        with open('thickcurr_results.json', 'w+') as fid:
            json.dump(results, fid)
    oftpy_dump_cov()


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as fid:
        main(json.load(fid))
