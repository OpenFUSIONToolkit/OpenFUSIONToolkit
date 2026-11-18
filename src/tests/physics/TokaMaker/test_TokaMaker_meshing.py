#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: pure-Python mesh domain setup (extents, region definitions, JSON round-trip)
import numpy as np
import pytest
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain


def build_domain(scale, **kwargs):
    '''Plasma rectangle inside a vacuum rectangle, both scaled vertically by `scale`, with an automatic boundary'''
    gs_mesh = gs_Domain(**kwargs)
    gs_mesh.define_region('air',0.1,'boundary')
    gs_mesh.define_region('plasma',0.1,'plasma')
    gs_mesh.define_region('vac',0.1,'vacuum')
    gs_mesh.add_rectangle(1.0,0.0,0.5,0.5*scale,'plasma')
    gs_mesh.add_rectangle(1.0,0.0,0.9,0.9*scale,'vac')
    return gs_mesh


@pytest.mark.coverage
def test_default_zextents_not_shared():
    '''Two domains built in one process with default z-extents must each pad their own geometry.

    Previously the default `zextents` list was shared between instances and filled in by
    `build_mesh()`, so the second domain reused (and could fail against) the first one's extents.
    '''
    for scale in (1.0, 2.0):
        gs_mesh = build_domain(scale)
        gs_mesh.build_mesh(setup_only=True)
        # Default padding is 1.2x the vertical extent of the vacuum region (half-height 0.45*scale)
        assert np.allclose(gs_mesh.zextents, [-1.2*0.45*scale, 1.2*0.45*scale])


@pytest.mark.coverage
def test_user_zextents_validated():
    '''User-specified z-extents that do not enclose all regions must raise, for both the lower and upper bound.'''
    with pytest.raises(ValueError, match=r'zextents\[0\]'):
        build_domain(1.0, zextents=[-0.1, 1.0]).build_mesh(setup_only=True)
    with pytest.raises(ValueError, match=r'zextents\[1\]'):
        build_domain(1.0, zextents=[-1.0, 0.1]).build_mesh(setup_only=True)


@pytest.mark.coverage
def test_mesh_reg_defs_not_shared():
    '''Region definitions must not accumulate across mesh setups, or leak into the domain's own list.'''
    nreg_defs = []
    for _ in range(2):
        gs_mesh = build_domain(1.0)
        gs_mesh.build_mesh(setup_only=True)
        nreg_defs.append(len(gs_mesh.mesh._reg_defs))
        assert len(gs_mesh._extra_reg_defs) == 0
    assert nreg_defs[0] == nreg_defs[1] == 3


@pytest.mark.coverage
def test_json_roundtrip(tmp_path):
    '''A domain saved with `save_json()` (including a region added with `add_enclosed()`) and reloaded
    with `gs_Domain(json_filename=...)` must produce the same extents and region definitions.'''
    json_file = tmp_path / 'domain.json'
    gs_ref = build_domain(1.0)
    gs_ref.define_region('gap',0.1,'vacuum')
    gs_ref.add_enclosed([1.0, 0.4],'gap')
    # Non-default small-feature threshold (default is dx_curve/2) must also survive the round-trip
    gs_ref.define_region('coil',0.1,'coil')
    gs_ref.add_polygon(np.array([[1.6,-0.1],[1.7,-0.1],[1.7,0.1],[1.6,0.1]]),'coil',small_thresh=0.3)
    gs_ref.save_json(str(json_file))
    gs_ref.build_mesh(setup_only=True)

    gs_load = gs_Domain(json_filename=str(json_file))
    gs_load.build_mesh(setup_only=True)
    assert len(gs_load.regions) == len(gs_ref.regions)
    assert np.allclose(gs_load.zextents, gs_ref.zextents)
    assert np.isclose(gs_load.rextent, gs_ref.rextent)
    assert np.allclose(gs_load.mesh._reg_defs, gs_ref.mesh._reg_defs)
    assert [region._small_thresh for region in gs_load.regions] == [region._small_thresh for region in gs_ref.regions]
