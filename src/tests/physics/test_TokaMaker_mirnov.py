#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Regression tests for TokaMaker reconstruction Mirnov normals.'''
import io
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'python')))
from OpenFUSIONToolkit.TokaMaker.reconstruction import Mirnov_con, Mirnov_con_id, reconstruction


@pytest.mark.parametrize('normal', [
    [1.0, 0.0],
    [1.0, 0.0, 0.0, 0.0],
    1.0,
    np.ones((3, 1)),
])
def test_mirnov_rejects_invalid_normal(normal):
    '''! Reject incorrectly shaped normals before adding a constraint.

    @param normal Invalid sensor normal
    '''
    recon = SimpleNamespace(_mirnovs=[])
    with pytest.raises(ValueError, match=r'three components \(R, phi, Z\)'):
        reconstruction.add_Mirnov(recon, [0.3, 0.0], normal, 0.01, 0.001)
    assert recon._mirnovs == []


@pytest.mark.parametrize('normal_type', [list, tuple, np.array])
def test_mirnov_normal_round_trip(normal_type):
    '''! Preserve cylindrical component order through constraint serialization.

    @param normal_type Supported container for the sensor normal
    '''
    normal = normal_type([0.36, 0.48, 0.8])
    recon = SimpleNamespace(_mirnovs=[])
    reconstruction.add_Mirnov(recon, [0.3, 0.1], normal, 0.01, 0.001)
    output = io.StringIO()
    recon._mirnovs[0].write(output)
    output.seek(0)
    assert int(output.readline()) == Mirnov_con_id
    restored = Mirnov_con()
    restored.read(output)
    assert restored.loc == pytest.approx([0.3, 0.1])
    assert restored.norm == pytest.approx([0.36, 0.48, 0.8])
    assert restored.val == pytest.approx(0.01)
    assert restored.err == pytest.approx(0.001)
