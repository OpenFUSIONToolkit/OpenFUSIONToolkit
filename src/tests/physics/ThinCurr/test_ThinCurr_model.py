#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThinCurr tests: model construction API on the dummy plate model (user point maps, models without Icoils)
import os
import pytest
import numpy as np
from thincurr_testing import mp_run, mu0, write_thincurr_xml, get_oft_env, save_dummy_mesh
from OpenFUSIONToolkit._interface import oftpy_dump_cov
test_dir = os.path.abspath(os.path.dirname(__file__))


def run_bfield_no_icoils(use_hodlr,mp_q):
    '''B-field operators and reconstruction for a model without I-coils, with input checks'''
    result = None
    try:
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        save_dummy_mesh('tw_no_icoils.h5', np.r_[1.E4*mu0])
        tw_model = ThinCurr(get_oft_env())
        tw_model.setup_model(mesh_file='tw_no_icoils.h5')
        tw_model.setup_io(basepath=f'no_icoils_{use_hodlr}/')
        if use_hodlr:
            tw_model.set_hodlr_options(target_size=20)
        tw_model.compute_Lmat(use_hodlr=use_hodlr)
        pot = np.sin(np.arange(tw_model.nelems)*0.37)
        checks = []
        try:
            tw_model.reconstruct_Bfield(pot)  # operators not computed yet
            checks.append(False)
        except ValueError:
            checks.append(True)
        _, Bdr = tw_model.compute_Bmat()
        checks.append(Bdr.shape == (3, 0, tw_model.np))
        try:
            tw_model.reconstruct_Bfield(pot[:-1])  # wrong size
            checks.append(False)
        except ValueError:
            checks.append(True)
        field = np.asarray(tw_model.reconstruct_Bfield(pot))
        checks.append(np.array_equal(field, tw_model.reconstruct_Bfield(pot[:, None])))  # column input
        result = (all(checks), field)
    except BaseException as e:
        print(e)
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_bfield_no_icoils():
    '''Magnetic field operators of a model without I-coils, with dense and compressed inductance matrices: invalid
    input is rejected, column vectors are accepted, and both give the same field'''
    out = {use_hodlr: mp_run(run_bfield_no_icoils, (use_hodlr,)) for use_hodlr in (False, True)}
    for res in out.values():
        assert res  # None (timeout, crash, or exception)
        assert res[0]
        assert np.all(np.isfinite(res[1]))
    assert np.allclose(out[True][1], out[False][1], rtol=0.0, atol=1.E-4*np.abs(out[False][1]).max())


def _default_pmap(lc, npts):
    '''Point map ThinCurr builds by default: interior vertices numbered in order, boundary vertices inactive (0)'''
    edges = np.sort(np.vstack((lc[:,[0,1]], lc[:,[1,2]], lc[:,[2,0]])), axis=1)
    uedges, counts = np.unique(edges, axis=0, return_counts=True)
    active = np.ones((npts,), dtype=bool)
    active[uedges[counts == 1].ravel()] = False
    pmap = np.zeros((npts,), dtype=np.int32)
    pmap[active] = np.arange(1, active.sum()+1)
    return pmap


@pytest.mark.coverage
def test_user_pmap():
    '''A user-supplied point map equal to the default reproduces the default model, a reversed map reverses the
    matrices, and a map of the wrong length is rejected'''
    assert mp_run(run_user_pmap, ())


def run_user_pmap(mp_q):
    '''Matrices for the default point map, the same map passed explicitly, and the map in reverse order'''
    result = True
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        from OpenFUSIONToolkit.ThinCurr.meshing import build_ThinCurr_dummy
        os.chdir(test_dir)
        xml_filename = 'oft_in_pmap.xml'
        write_thincurr_xml(xml_filename, eta_values=np.r_[1.E4*mu0])
        myOFT = OFT_env(nthreads=-1)
        r, lc = build_ThinCurr_dummy([0.0, 0.0, 10.0], size=0.25, nsplit=2)
        pmap = _default_pmap(lc, r.shape[0])
        nactive = pmap.max()
        reverse = pmap.copy()
        reverse[pmap > 0] = nactive+1-pmap[pmap > 0]
        mats = []
        for pmap_in in (None, pmap, reverse):
            model = ThinCurr(myOFT)
            model.setup_model(r=r, lc=lc, xml_filename=xml_filename, pmap=pmap_in)
            model.compute_Lmat()
            model.compute_Rmat()
            mats.append((np.array(model.Lmat), model.Rmat.toarray()))
        if not (np.allclose(mats[1][0], mats[0][0], rtol=1.E-12, atol=0.0)
                and np.allclose(mats[1][1], mats[0][1], rtol=1.E-12, atol=0.0)):
            print("FAILED: default point map does not reproduce the model")
            result = False
        if not (np.allclose(mats[2][0], mats[0][0][::-1,::-1], rtol=0.0, atol=1.E-6*np.abs(mats[0][0]).max())
                and np.allclose(mats[2][1], mats[0][1][::-1,::-1], rtol=0.0, atol=1.E-12*np.abs(mats[0][1]).max())):
            print("FAILED: reordered point map does not reorder the model")
            result = False
        try:
            ThinCurr(myOFT).setup_model(r=r, lc=lc, xml_filename=xml_filename, pmap=pmap[:-1])
            print("FAILED: point map of wrong length accepted")
            result = False
        except ValueError:
            pass
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)
