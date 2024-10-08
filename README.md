Open Flexible Unstructured Simulation Infrastructure with Open Numerics (FUSION) Toolkit (OFT)
=====================================

[![CI status](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/ci_build.yaml/badge.svg?branch=main)](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/ci_build.yaml)
[![CD status](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/cd_nightly.yaml/badge.svg)](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/cd_nightly.yaml)
[![codecov](https://codecov.io/gh/openfusiontoolkit/OpenFUSIONToolkit/graph/badge.svg?token=GG282HKNAO)](https://codecov.io/gh/openfusiontoolkit/OpenFUSIONToolkit)
[![DOI](https://zenodo.org/badge/710415041.svg)](https://zenodo.org/doi/10.5281/zenodo.10306801)

The Open FUSION Toolkit (OFT) is a suite of modeling tools, and their underlying finite element
framework, for problems in plasma and fusion research and engineering in arbitrary 2D and 3D geometries.
The underlying framework and its component tools support the use of variable order finite element methods on
unstructured tetrahedral (triangular) or hexahedral (quadralateral) grids.

**NOTE:** This project is current in a beta status as the transition to a fully open source code, including
associated documentation and examples is completed. For more information see [v1.0 release status](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/milestone/1) for more information.

Component tools
------------

The suite currently includes the following tools:

* **TokaMaker:** Axisymmetric statice and time-dependent ideal MHD equilibria

* **ThinCurr:** Inductively-excited currents in the thin-wall limit

* **MUG:** Time-dependent nonlinear (linear) extendend MHD

* **Marklin:** 3D force-free ideal MHD equilibria

Installation
------------

Binaries are available for Linux (x86) and macOS (x86 and arm64) for each [release on GitHub](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/releases). The framework and its components can also be built from source using the instructions provided on [the wiki](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/wiki).

To use the python interfaces add the `python` directory for your installation to your `PYTHONPATH` environment variable.

Documentation
------------

Detailed documentation is [available online](https://openfusiontoolkit.github.io/OpenFUSIONToolkit/docs/index.html) and is also included in the binary packages at `/path/to/oft/install/doc/Documentation.html` (please note `/path/to/oft/install` is a placeholder for the actual install path).

Contributing
-----------

See [CONTRIBUTING.md](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/blob/main/CONTRIBUTING.md) for information on how to contribute.


Copyright
---------

Open FUSION Toolkit code development project, up to version 1.0 Copyright (c) 2023, Open FUSION Toolkit team.

Written by Open FUSION Toolkit team and collaborators with Christopher J. Hansen as principle developer. All rights reserved.
