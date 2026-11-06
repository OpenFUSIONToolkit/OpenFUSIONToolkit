#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: reconstruction constraint file I/O without a native solver
import os
import sys
from unittest.mock import Mock

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'python')))
from OpenFUSIONToolkit.TokaMaker.reconstruction import Ip_con, Press_con, reconstruction


@pytest.fixture
def recon(tmp_path):
    '''Reconstruction object writing to a temporary directory, with a mocked TokaMaker environment'''
    maker = Mock()
    maker._oft_env.oft_in_groups = {}
    maker.coil_sets = {'PF1': {'id': 0}}
    return reconstruction(maker, in_filename=str(tmp_path / 'fit.in'),
                          out_filename=str(tmp_path / 'fit.out'))


@pytest.mark.parametrize('constraint', ['current', 'pressure'])
def test_reconstruction_positive_constraint_round_trip(recon, constraint):
    '''Plasma current and pressure constraints survive a write/read of the constraint file without tripping the
    positivity check on an uninitialized value.'''
    if constraint == 'current':
        recon.set_Ip(100000.0, 1000.0)
    else:
        recon.add_pressure([0.3, 0.1], 1000.0, 10.0)
    recon.write_fit_in()
    recon.read_fit_in()
    if constraint == 'current':
        assert recon._Ip_con.val == pytest.approx(100000.0)
        assert recon._Ip_con.err == pytest.approx(1000.0)
    else:
        assert len(recon._pressure_cons) == 1
        restored = recon._pressure_cons[0]
        assert restored.loc == pytest.approx([0.3, 0.1])
        assert restored.val == pytest.approx(1000.0)
        assert restored.err == pytest.approx(10.0)


def test_reconstruction_repeated_coil_constraint_read(recon):
    '''Reading the constraint file twice replaces the coil-current constraints instead of duplicating them, and a reset
    writes an empty file.'''
    recon.set_coil_currents({'PF1': 100.0}, {'PF1': 10.0})
    recon.write_fit_in()
    for _ in range(2):
        recon.read_fit_in()
        assert len(recon._coil_current_cons) == 1
        restored = recon._coil_current_cons[0]
        assert restored.ind == 0
        assert restored.val == pytest.approx(100.0)
        assert restored.err == pytest.approx(10.0)
    recon.reset_constraints()
    recon.write_fit_in()
    with open(recon.con_file) as stream:
        assert int(stream.readline()) == 0


@pytest.mark.parametrize('constraint_class', [Ip_con, Press_con])
@pytest.mark.parametrize('value', [0.0, -1.0])
def test_reconstruction_rejects_nonpositive_constraint(constraint_class, value):
    '''Zero or negative current and pressure measurements are rejected when a constraint is created.'''
    with pytest.raises(ValueError, match='must be positive'):
        constraint_class(val=value, err=1.0)


@pytest.mark.parametrize('constraint_type,location', [(2, ''), (9, '0.3 0.1\n')])
@pytest.mark.parametrize('value', [0.0, -1.0])
def test_reconstruction_rejects_nonpositive_file_value(recon, constraint_type, location, value):
    '''Zero or negative current and pressure measurements are rejected when read from the constraint file.'''
    with open(recon.con_file, 'w') as stream:
        stream.write('1\n\n{0}\n{1}{2} 1.0\n\n'.format(constraint_type, location, value))
    with pytest.raises(ValueError, match='must be positive'):
        recon.read_fit_in()
