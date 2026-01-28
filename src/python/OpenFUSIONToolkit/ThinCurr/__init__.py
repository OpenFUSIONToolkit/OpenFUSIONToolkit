#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Python interface for ThinCurr thin-wall E-M functionality

@authors Chris Hansen
@date March 2024
@ingroup doxy_oft_python
'''
from ._core import ThinCurr, ThinCurr_reduced

__all__ = ["ThinCurr", "ThinCurr_reduced"]