#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Driver for ThinCurr MPI tests, run directly (serial) or with mpirun (one copy per rank)
#
# Usage: python thincurr_mpi_driver.py <case.json>
#
# Runs HODLR operations and solves on a ThinCurr model and saves the results (head process only) to the file
# named in the case, for comparison between serial and parallel runs.
import os
import sys
import json
import numpy as np
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
import oft_testing # noqa: F401 (adds the OFT Python package to the path)
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.ThinCurr import ThinCurr


def main(case):
    myOFT = OFT_env(nthreads=1)
    tw = ThinCurr(myOFT)
    #---Expected errors for unsupported operations with MPI
    if case.get('expect_error', False):
        tw.setup_model(mesh_file=case['mesh'], xml_filename=case.get('xml'), jumper_start=case.get('jumper_start', 0))
        tw.compute_Mcoil()
        try:
            tw.compute_Lmat(use_hodlr=False)
        except ValueError as e:
            if myOFT.head_proc:
                print('Expected failure detected: {0}'.format(e))
        return
    tw.setup_model(mesh_file=case['mesh'], xml_filename=case.get('xml'), jumper_start=case.get('jumper_start', 0))
    tw.setup_io(basepath=case['basepath'], legacy_hdf5=True)
    tw.set_hodlr_options(**case['hodlr'])
    results = {'nprocs': myOFT.nprocs}
    #---Coil and sensor coupling matrices (computed in parallel, replicated on all processes)
    Mc = tw.compute_Mcoil()
    if tw.n_icoils > 0:
        results['Mc'] = Mc.copy()
    if case.get('sensors') is not None:
        Ms, Msc, _ = tw.compute_Msensor(case['sensors'])
        results['Ms'] = Ms.copy()
        if tw.n_icoils > 0:
            results['Msc'] = Msc.copy()
    tw.compute_Lmat(use_hodlr=True, cache_file=case.get('cache'))
    tw.compute_Rmat()
    rng = np.random.default_rng(42)
    x = rng.standard_normal(tw.nelems)
    #---Matrix-vector products
    results['Lx'] = tw.apply_Lmat(x)
    #---L/R eigenmodes (parallel ARPACK)
    eig_vals, eig_vecs = tw.get_eigs(4)
    results['eig_vals'] = eig_vals
    results['eig_vec0'] = eig_vecs[0,:]*np.sign(eig_vecs[0,np.argmax(np.abs(eig_vecs[0,:]))])
    #---Frequency response (complex GMRES with block-Jacobi preconditioner)
    vdriver = np.vstack([x, 0.5*x])
    results['fr'] = tw.compute_freq_response(vdriver=vdriver, freq=case.get('freq', 1.E3))
    #---Time-dependent run (CG with block-Jacobi or two-level preconditioner)
    coil_currs = None
    if tw.n_icoils > 0:
        coil_currs = np.array([[0.0, 0.0] + [0.0]*(tw.n_icoils-1), [1.E-4, 1.E3] + [0.0]*(tw.n_icoils-1),
                               [1.0, 1.E3] + [0.0]*(tw.n_icoils-1)])
    full_volts = None
    if coil_currs is None:
        full_volts = np.zeros((2, tw.nelems+1))
        full_volts[:, 0] = [0.0, 1.0]
        full_volts[:, 1:] = x[np.newaxis, :]
    tw.run_td(1.E-4, 10, coil_currs=coil_currs, full_volts=full_volts, lin_tol=1.E-12, lin_rtol=1.E-20,
              status_freq=5, plot_freq=10)
    if myOFT.head_proc:
        import h5py
        with h5py.File(os.path.join(case['basepath'], 'pThinCurr_0010.rst'), 'r') as fid:
            results['td_U'] = np.array(fid['U'])
    #---Magnetic field reconstruction (HODLR B-field operator)
    if case.get('compute_B', True):
        tw.compute_Bmat()
        results['B'] = tw.reconstruct_Bfield(x)
    if myOFT.head_proc:
        np.savez(case['out'], **results)
        print('DRIVER_COMPLETE')


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as fid:
        main(json.load(fid))
