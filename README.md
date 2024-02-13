Open Flexible Unstructured Simulation Infrastructure with Open Numerics (FUSION) Toolkit (OFT)
=====================================

![CI status](https://github.com/hansec/OpenFUSIONToolkit/actions/workflows/build_test.yaml/badge.svg?branch=main)
[![DOI](https://zenodo.org/badge/710415041.svg)](https://zenodo.org/doi/10.5281/zenodo.10306801)

The Open FUSION Toolkit (OFT) is a suite of modeling tools, and their underlying finite element
framework, for problems in plasma and fusion research and engineering in arbitrary 2D and 3D geometries.
The underlying framework and its component tools support the use of variable order finite element methods on
unstructured tetrahedral (triangular) or hexahedral (quadralateral) grids.

**NOTE:** This project is current in a beta status as the transition to a fully open source code, including
associated documentation and examples is completed. For more information see [v1.0 release status](https://github.com/hansec/OpenFUSIONToolkit/milestone/1) for more information.

Component tools
------------

The suite currently includes the following tools:

* **TokaMaker:** Axisymmetric statice and time-dependent ideal MHD equilibria

* **ThinCurr:** Inductively-excited currents in the thin-wall limit

* **MUG:** Time-dependent nonlinear (linear) extendend MHD

* **Marklin:** 3D force-free ideal MHD equilibria

Installation
------------

Binaries are available for Linux (x86) and macOS (x86 and arm64) for each [release on GitHub](https://github.com/hansec/OpenFUSIONToolkit/releases). The framework and its components can also be built from source using the instructions provided on [the wiki](https://github.com/hansec/OpenFUSIONToolkit/wiki).

To use the python interfaces add the `python` directory for your installation to your `PYTHONPATH` environment variable.

Documentation
------------

Detailed documentation is included in the binary packages at `/path/to/oft/install/doc/Documentation.html` (note `/path/to/oft/install` is a placeholder for the actual install path). In the future this documentation will also be made available online.

Contributing
-----------

See [CONTRIBUTING.md](https://github.com/hansec/OpenFUSIONToolkit/blob/main/CONTRIBUTING.md) for information on how to contribute.


Copyright
---------

Open FUSION Toolkit code development project, up to version 1.0 Copyright (c) 2023, Open FUSION Toolkit team.

Written by Open FUSION Toolkit team and collaborators with Christopher J. Hansen as principle developer. All rights reserved.
