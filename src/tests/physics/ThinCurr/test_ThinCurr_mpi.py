#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThinCurr tests: MPI runs (HODLR matrices and solvers distributed by leaf blocks) against serial runs
#
# Each case runs thincurr_mpi_driver.py serially and with mpirun on 2 and 3 processes, and compares HODLR
# products, L/R eigenmodes (parallel ARPACK), frequency response, time-dependent results, and B-field
# reconstruction. HODLR tolerances are tight so results agree closely, even though ACA+ pivots differ between runs.
import os
import sys
import json
import shutil
import pytest
import numpy as np
from thincurr_testing import test_dir
import oft_testing
from oft_testing import run_OFT
from OpenFUSIONToolkit.util import mu0
from OpenFUSIONToolkit.io import write_oft_xml
from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_Icoil, ThinCurr_Vcoil, ThinCurr_XML
from OpenFUSIONToolkit.ThinCurr.sensor import save_sensors, circular_flux_loop, Mirnov

# Tight HODLR tolerances and small leaves (several leaf blocks per process)
hodlr_opts = {'target_size': 200, 'L_svd_tol': 1.E-10, 'B_svd_tol': 1.E-10}


def write_cyl_xml():
    '''Resistivity, coils, and sensors for the high-resolution cylinder

    One circular and one polygon I-coil, one circular V-coil, two flux loops and a point Mirnov probe, covering the
    ring and polygon paths of the coil and sensor coupling matrices.
    '''
    xml = ThinCurr_XML()
    xml.set_eta([10.0*mu0])
    icoil = ThinCurr_Icoil(name="IC")
    icoil.add_circular_subcoil([1.1, 0.25])
    xml.add_Icoil(icoil)
    phi = np.linspace(0.0, 2.0*np.pi, 91)
    icoil = ThinCurr_Icoil(name="IC_POLY")
    icoil.add_subcoil(pts=np.c_[1.2*np.cos(phi), 1.2*np.sin(phi), 0.0*phi+0.4])
    xml.add_Icoil(icoil)
    vcoil = ThinCurr_Vcoil(name="VC", resistivity_per_length=1.256637E-5, radius=1.E-2)
    vcoil.add_circular_subcoil([1.1, -0.25])
    xml.add_Vcoil(vcoil)
    write_oft_xml([xml], "oft_in_mpi.xml")


def write_sensors():
    '''Two flux loops and a point Mirnov probe'''
    sensors = [circular_flux_loop(0.9, 0.5, 'FLOOP_0', npts=180), circular_flux_loop(0.9, 0.0, 'FLOOP_1', npts=180),
               Mirnov([0.9, 0.0, 0.2], [0.0, 0.0, 1.0], 'MIRNOV_0', dx=None)]
    save_sensors(sensors, 'floops_mpi.loc')


def write_junction_mesh():
    '''Junction mesh (two cylinders joined by a flange) with an I-coil and a V-coil'''
    from test_ThinCurr_junctions import ring, write_pair
    r, lc, reg = ring()
    files, _ = write_pair(r, lc, reg, 'mpi_junction', vcoil=True)
    return files[0]


def run_driver(name, nproc, case, cache=None, expect_cache_read=False):
    '''Run the driver for `case` on `nproc` processes and return its results

    With `cache`, HODLR and coil/sensor coupling matrices are cached in files with that prefix.
    '''
    case = dict(case)
    case['out'] = '{0}_np{1}.npz'.format(name, nproc)
    case['basepath'] = '{0}_np{1}/'.format(name, nproc)
    if cache is not None:
        case['cache'] = cache
        case['cm_cache'] = cache
    case_file = '{0}_np{1}.json'.format(name, nproc)
    with open(case_file, 'w+') as fid:
        json.dump(case, fid)
    if os.path.isfile(case['out']):
        os.remove(case['out'])
    result, std_out = run_OFT('{0} thincurr_mpi_driver.py {1}'.format(sys.executable, case_file), nproc, 600,
                              return_stdout=True)
    if os.path.isdir(case['basepath']):
        shutil.rmtree(case['basepath'])
    assert result
    assert std_out.find('DRIVER_COMPLETE') > -1
    if expect_cache_read:
        for msg in ('Reading HODLR matrix from file', 'Reading coil mutual matrices', 'Reading sensor mutual matrices'):
            assert std_out.find(msg) > -1, msg
        assert std_out.find('Ignoring stored matrix') < 0
        assert std_out.find('Error reading matrix from file') < 0
    with np.load(case['out']) as data:
        return {key: data[key] for key in data.files}


def compare(ref, out, rtol=1.E-6):
    '''Compare driver results (relative to the largest value of each quantity)

    Iterative solves (frequency response and time-dependent) are compared to `10*rtol`, as differences in the
    operators at the level of the HODLR tolerance and round-off in the reductions change the solver iterations.
    '''
    assert out['nprocs'] > 1
    for key in ('Msc', 'red_Ms', 'red_Mc', 'red_L', 'red_R', 'Lx', 'eig_vals', 'eig_vec0', 'fr', 'td_U', 'td_sensors', 'B'):
        if key not in ref:
            continue
        scale = np.max(np.abs(ref[key]))
        assert scale > 0.0, key
        err = np.max(np.abs(out[key]-ref[key]))/scale
        print("{0} (nprocs = {1}): relative difference {2:.3E}".format(key, int(out['nprocs']), err))
        tol = rtol
        if key in ('fr', 'td_U', 'td_sensors'):
            tol = 10.0*rtol
        assert err < tol, "{0}: relative error {1:.3E}".format(key, err)


cases = {
    'cyl': {'mesh': 'tw_test-cyl_hr.h5', 'xml': 'oft_in_mpi.xml', 'jumper_start': 2, 'sensors': 'floops_mpi.loc',
            'hodlr': dict(hodlr_opts)},
    'cyl_coarse': {'mesh': 'tw_test-cyl_hr.h5', 'xml': 'oft_in_mpi.xml', 'jumper_start': 2,
                   'hodlr': dict(hodlr_opts, pre_coarse=True, pre_neigs=4), 'compute_B': False},
    'junction': {'mesh': None, 'sensors': 'floops_mpi.loc', 'hodlr': dict(hodlr_opts, target_size=100)},
}


def setup_case(name):
    os.chdir(test_dir)
    case = dict(cases[name])
    if case['mesh'] is None:
        case['mesh'] = write_junction_mesh()
    if 'xml' in case:
        write_cyl_xml()
    write_sensors()
    return case


@pytest.mark.mpi
@pytest.mark.parametrize("name", ('cyl', 'cyl_coarse', 'junction'))
def test_mpi_vs_serial(name):
    '''HODLR products and solves on 2 and 3 processes match a serial run'''
    if not oft_testing.env_flag('OFT_HAVE_MPI'):
        pytest.skip("Not compiled with MPI")
    case = setup_case(name)
    ref = run_driver(name, 1, case)
    for nproc in (2, 3):
        compare(ref, run_driver(name, nproc, case))


@pytest.mark.mpi
def test_mpi_cache():
    '''HODLR and coupling matrix cache files written by any number of processes are read by any other number of
    processes'''
    if not oft_testing.env_flag('OFT_HAVE_MPI'):
        pytest.skip("Not compiled with MPI")
    case = setup_case('cyl')
    case['compute_B'] = False
    cache = 'mpi_hodlr_cache.h5'
    cache_files = (cache, cache + '_coils', cache + '_sensors')
    for file in cache_files:
        if os.path.isfile(file):
            os.remove(file)
    ref = run_driver('cache', 2, case, cache=cache) # Written by 2 processes
    compare(ref, run_driver('cache', 3, case, cache=cache, expect_cache_read=True), rtol=1.E-8)
    serial = run_driver('cache', 1, case, cache=cache, expect_cache_read=True)
    serial['nprocs'] = 2 # Compared as a parallel result
    compare(ref, serial, rtol=1.E-8)
    for file in cache_files:
        os.remove(file)


@pytest.mark.mpi
def test_mpi_dense_error():
    '''Dense inductance matrices are rejected with MPI'''
    if not oft_testing.env_flag('OFT_HAVE_MPI'):
        pytest.skip("Not compiled with MPI")
    os.chdir(test_dir)
    write_cyl_xml()
    case_file = 'dense_error.json'
    with open(case_file, 'w+') as fid:
        json.dump({'mesh': 'tw_test-cyl.h5', 'xml': 'oft_in_mpi.xml', 'jumper_start': 2, 'expect_error': True}, fid)
    result, std_out = run_OFT('{0} thincurr_mpi_driver.py {1}'.format(sys.executable, case_file), 2, 120,
                              return_stdout=True)
    assert result
    assert std_out.find('Expected failure detected') > -1


@pytest.mark.mpi
def test_mpi_too_few_leaves():
    '''A HODLR partition with fewer leaf blocks than processes is rejected with a clear error (not a hang)'''
    if not oft_testing.env_flag('OFT_HAVE_MPI'):
        pytest.skip("Not compiled with MPI")
    os.chdir(test_dir)
    write_cyl_xml()
    case = {'mesh': 'tw_test-cyl.h5', 'xml': 'oft_in_mpi.xml', 'jumper_start': 2, 'basepath': 'few_leaves/',
            'out': 'few_leaves.npz', 'hodlr': dict(hodlr_opts, target_size=5000)}
    with open('few_leaves.json', 'w+') as fid:
        json.dump(case, fid)
    assert run_OFT('{0} thincurr_mpi_driver.py few_leaves.json'.format(sys.executable), 2, 120,
                   expect_abort='fewer leaf blocks')
    if os.path.isdir(case['basepath']):
        shutil.rmtree(case['basepath'])
