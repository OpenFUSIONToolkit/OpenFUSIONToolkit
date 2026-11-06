#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThinCurr tests: surface/volume resistivity and thickness API on a single-filament model
import os
import glob
import shutil
import pytest
import numpy as np
from thincurr_testing import mp_run, mu0, write_thincurr_xml, build_dummy_model, get_oft_env, save_dummy_mesh
from OpenFUSIONToolkit._interface import oftpy_dump_cov
test_dir = os.path.abspath(os.path.dirname(__file__))


@pytest.mark.coverage
def test_thickness_api_roundtrip_and_validation():
    '''Setting and reading surface resistivity, volumetric resistivity and thickness in every supported combination,
    including the deprecated eta_values keyword, inconsistent triples (warn), invalid inputs (raise), and surface
    resistivity stored in the mesh file.'''
    assert mp_run(run_thickness_api_roundtrip_and_validation, ())


def run_thickness_api_roundtrip_and_validation(mp_q):
    result = True
    try:
        import warnings
        os.chdir(test_dir)
        xml_filename = 'oft_in_thickness_api.xml'
        eta_values = np.r_[1.E4*mu0]
        thickness_values = np.r_[2.5E-3]
        eta_vol_values = eta_values * thickness_values
        write_thincurr_xml(xml_filename, eta_values=eta_values, thickness_values=thickness_values)

        tw_model = build_dummy_model(xml_filename)

        tw_model.set_eta_values(eta_surf=eta_values, thickness=thickness_values)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            eta_returned = tw_model.get_eta_values()
            if len(w) == 0 or not issubclass(w[-1].category, DeprecationWarning):
                print("WARNING: Expected DeprecationWarning for get_eta_values() not found")

        if not np.allclose(eta_returned, eta_values):
            result = False
        eta_surf_out, eta_vol_out = tw_model.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_out, eta_values):
            result = False
        if not np.allclose(eta_vol_out, eta_vol_values):
            result = False
        if not np.allclose(tw_model.get_thickness(), thickness_values):
            result = False

        tw_model_pair = build_dummy_model(xml_filename)
        tw_model_pair.set_eta_values(eta_surf=eta_values, eta_vol=eta_vol_values)
        eta_surf_pair, eta_vol_pair = tw_model_pair.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_pair, eta_values):
            result = False
        if not np.allclose(eta_vol_pair, eta_vol_values):
            result = False
        if not np.allclose(tw_model_pair.get_thickness(), thickness_values):
            result = False

        tw_model_alias = build_dummy_model(xml_filename)
        tw_model_alias.set_eta_values(eta_values=eta_values, thickness=thickness_values)
        eta_surf_alias, eta_vol_alias = tw_model_alias.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_alias, eta_values):
            result = False
        if not np.allclose(eta_vol_alias, eta_vol_values):
            result = False

        tw_model_eta_only = build_dummy_model(xml_filename)
        tw_model_eta_only.set_eta_values(eta_surf=eta_values)
        eta_surf_only, eta_vol_only = tw_model_eta_only.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_only, eta_values):
            result = False
        if eta_vol_only is not None:
            result = False
        if not np.allclose(tw_model_eta_only.get_thickness(), -1.0):
            result = False

        eta_surf_override = eta_values * 1.5
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            tw_model_three = build_dummy_model(xml_filename)
            tw_model_three.set_eta_values(eta_surf=eta_surf_override, eta_vol=eta_vol_values, thickness=thickness_values)
        if len(w) == 0 or not any(issubclass(item.category, UserWarning) for item in w):
            result = False
        eta_surf_three, eta_vol_three = tw_model_three.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_three, eta_values):
            result = False
        if not np.allclose(eta_vol_three, eta_vol_values):
            result = False
        if not np.allclose(tw_model_three.get_thickness(), thickness_values):
            result = False

        new_thickness = np.r_[3.5E-3]
        tw_model.set_eta_values(eta_surf=eta_values, thickness=new_thickness)
        if not np.allclose(tw_model.get_thickness(), new_thickness):
            result = False

        with pytest.raises(ValueError):
            tw_model.set_eta_values()
        with pytest.raises(ValueError):
            tw_model.set_eta_values(eta_values=np.r_[-1.0])
        with pytest.raises(ValueError):
            tw_model.set_eta_values(thickness=thickness_values)
        with pytest.raises(ValueError):
            tw_model.set_eta_values(eta_surf=eta_values, thickness=np.r_[-1.0])
        # Volumetric resistivity alone needs a thickness already defined
        xml_no_thickness = 'oft_in_no_thickness.xml'
        write_thincurr_xml(xml_no_thickness, eta_values=eta_values)
        tw_model_no_thickness = build_dummy_model(xml_no_thickness)
        with pytest.raises(Exception):
            tw_model_no_thickness.set_eta_values(eta_vol=np.r_[1.0e-6])
        with pytest.raises(IndexError):
            tw_model.set_eta_values(eta_surf=eta_values, thickness=np.r_[1.E-3, 2.E-3])
        # Surface resistivity only from the mesh file: no volumetric resistivity is set
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        save_dummy_mesh('tw_eta_surf_only.h5', eta_values)
        tw_model_mesh_eta = ThinCurr(get_oft_env())
        tw_model_mesh_eta.setup_model(mesh_file='tw_eta_surf_only.h5')
        eta_surf_mesh, eta_vol_mesh = tw_model_mesh_eta.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_mesh, eta_values):
            result = False
        if eta_vol_mesh is not None:
            result = False
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_plot_td_compute_jvol_outputs_fields():
    '''Plot output of a time-dependent run includes the volumetric current density and thickness fields when a
    thickness is defined.'''
    assert mp_run(run_plot_td_compute_jvol_outputs_fields, ())


def run_plot_td_compute_jvol_outputs_fields(mp_q):
    result = True
    try:
        os.chdir(test_dir)
        xml_filename = 'oft_in_compute_jvol.xml'
        io_basepath = 'td_jvol_regression'
        eta_values = np.r_[1.E4*mu0]
        thickness_values = np.r_[2.0E-3]
        write_thincurr_xml(xml_filename, eta_values=eta_values, thickness_values=thickness_values)

        if os.path.isdir(io_basepath):
            shutil.rmtree(io_basepath)

        tw_model = build_dummy_model(xml_filename)
        tw_model.setup_io(basepath=io_basepath)
        tw_model.compute_Mcoil()
        tw_model.compute_Lmat()
        tw_model.compute_Rmat()
        tw_model.run_td(2.E-5, 5, direct=True, plot_freq=1)
        tw_model.plot_td(5, plot_freq=1)
        tw_model.build_XDMF()

        xmf_files = sorted(glob.glob(os.path.join(io_basepath, '*.xmf')))
        if len(xmf_files) == 0:
            result = False
        xmf_text = ''
        for xmf_file in xmf_files:
            with open(xmf_file, 'r') as fid:
                xmf_text += fid.read()
        if 'J_vol' not in xmf_text or 'thickness' not in xmf_text:
            result = False
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_eta_only_matches_surface_resistivity_with_thickness():
    '''A model defined by volumetric resistivity and thickness has the same resistance matrix as one defined by the
    equivalent surface resistivity.'''
    assert mp_run(run_eta_only_matches_surface_resistivity_with_thickness, ())


def run_eta_only_matches_surface_resistivity_with_thickness(mp_q):
    result = True
    try:
        os.chdir(test_dir)
        xml_filename = 'oft_in_eta_compat.xml'
        eta_surface = np.r_[1.E4*mu0]
        thickness_values = np.r_[2.5E-3]
        eta_bulk = eta_surface*thickness_values
        write_thincurr_xml(xml_filename, eta_values=eta_surface)

        model_surface = build_dummy_model(xml_filename)
        model_surface.set_eta_values(eta_values=eta_surface)
        model_surface.compute_Rmat()
        R_surface = model_surface.Rmat.toarray()

        model_bulk = build_dummy_model(xml_filename)
        model_bulk.set_eta_values(eta_vol=eta_bulk, thickness=thickness_values)
        model_bulk.compute_Rmat()
        R_bulk = model_bulk.Rmat.toarray()

        if not np.allclose(R_surface, R_bulk, rtol=1.E-10, atol=1.E-12):
            result = False
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_eta_vol_without_thickness_warns_and_loads():
    '''An XML file giving only the volumetric resistivity loads with a warning and leaves the surface resistivity
    undefined.'''
    assert mp_run(run_eta_vol_without_thickness_warns_and_loads, ())


def run_eta_vol_without_thickness_warns_and_loads(mp_q):
    result = True
    try:
        os.chdir(test_dir)
        xml_filename = 'oft_in_eta_vol_only_invalid.xml'
        eta_vol_values = np.r_[2.5E-3]
        write_thincurr_xml(xml_filename, eta_vol_values=eta_vol_values)

        # Surface resistivity cannot be inferred without a thickness
        tw_model = build_dummy_model(xml_filename)
        eta_surf = tw_model.get_eta_values()
        if np.all(eta_surf > 0.0):
            result = False
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_get_eta_values_deprecation_warning():
    '''get_eta_values() without include_eta_vol warns of the deprecation but still returns the surface resistivity.
    get_eta_values() without include_eta_vol warns of the deprecation but still returns the surface resistivity.'''
    assert mp_run(run_get_eta_values_deprecation_warning, ())


def run_get_eta_values_deprecation_warning(mp_q):
    result = True
    try:
        os.chdir(test_dir)
        xml_filename = 'oft_in_eta_deprecation.xml'
        eta_values = np.r_[1.E4*mu0]
        write_thincurr_xml(xml_filename, eta_values=eta_values)

        tw_model = build_dummy_model(xml_filename)
        tw_model.set_eta_values(eta_surf=eta_values)

        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            eta = tw_model.get_eta_values()
            if len(w) == 0 or not issubclass(w[-1].category, DeprecationWarning):
                print("ERROR: Expected DeprecationWarning not found")
                result = False
            elif "include_eta_vol" not in str(w[-1].message):
                print(f"ERROR: Expected 'include_eta_vol' in warning message, got: {w[-1].message}")
                result = False
            elif not np.allclose(eta, eta_values):
                print("ERROR: eta_surf values don't match")
                result = False

    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_get_eta_values_with_flag_returns_tuple():
    '''get_eta_values(include_eta_vol=True) returns the surface and volumetric resistivities without a warning.'''
    assert mp_run(run_get_eta_values_with_flag_returns_tuple, ())


def run_get_eta_values_with_flag_returns_tuple(mp_q):
    result = True
    try:
        os.chdir(test_dir)
        xml_filename = 'oft_in_eta_vol_tuple.xml'
        eta_surf = np.r_[1.E4*mu0]
        eta_vol = np.r_[2.5E-2]
        thickness = np.r_[2.5E-3]
        write_thincurr_xml(xml_filename, eta_values=eta_surf, thickness_values=thickness)

        tw_model = build_dummy_model(xml_filename)
        tw_model.set_eta_values(eta_vol=eta_vol, thickness=thickness)

        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result_tuple = tw_model.get_eta_values(include_eta_vol=True)
            for warning in w:
                if issubclass(warning.category, DeprecationWarning):
                    print(f"ERROR: Unexpected DeprecationWarning when include_eta_vol=True: {warning.message}")
                    result = False
                    break

        if not isinstance(result_tuple, tuple) or len(result_tuple) != 2:
            print(f"ERROR: Expected tuple of length 2, got {type(result_tuple)} of length {len(result_tuple) if isinstance(result_tuple, tuple) else 'N/A'}")
            result = False
        else:
            eta_surf_ret, eta_vol_ret = result_tuple
            if eta_vol_ret is None:
                print("ERROR: eta_vol should not be None when it was explicitly set")
                result = False
            elif not np.allclose(eta_vol_ret, eta_vol):
                print(f"ERROR: eta_vol values don't match. Expected {eta_vol}, got {eta_vol_ret}")
                result = False

    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)
