Open Flexible Unstructured Simulation Infrastructure with Open Numerics (FUSION) Toolkit (OFT)
=====================================

[![CI status](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/ci_build.yaml/badge.svg?branch=main)](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/ci_build.yaml)
[![CD status](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/cd_combined.yaml/badge.svg)](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/cd_combined.yaml)
[![Container CD status](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/container_cd.yaml/badge.svg)](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/container_cd.yaml)
[![codecov](https://codecov.io/gh/openfusiontoolkit/OpenFUSIONToolkit/graph/badge.svg?token=GG282HKNAO)](https://codecov.io/gh/openfusiontoolkit/OpenFUSIONToolkit)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10306801.svg)](https://doi.org/10.5281/zenodo.10306801)

<img src="https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit/raw/main/logos/oft_logo_bg.png" width="300px" />

The Open FUSION Toolkit (OFT) is a suite of modeling tools, and their underlying finite element
framework, for problems in plasma and fusion research and engineering in arbitrary 2D and 3D geometries.
The underlying framework and its component tools support the use of variable order finite element methods on
unstructured tetrahedral (triangular) or hexahedral (quadralateral) grids.

**NOTE:** This project is under active development, please watch [releases](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/releases) for any new features and breaking changes that may be introduced.

Component tools
------------

The suite currently includes the following tools:

* **TokaMaker:** A time-dependent free-boundary Grad-Shafranov equilibrium code

* **ThinCurr:** A 3D thin-wall electromagnetic simulation package

* **MUG:** A 2D/3D linear/nonlinear extended MHD simulation package

* **Marklin:** A 3D force-free, uniform &lambda; equilibrium solver

Installation
------------

### PyPI packages
Pre-built Python packages are available on [PyPI](https://pypi.org/project/openfusiontoolkit/) and can be installed using `pip install openfusiontoolkit` (using a virtual environment is strongly recommended).

### Binary packages
Binaries are available for Linux (x86 and arm64) and macOS (x86 and arm64) for each [release on GitHub](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/releases) as well as recent commits to `main` via the [CD workflow](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/actions/workflows/cd_combined.yaml) (see artifacts on each run).

### Container images
Pre-built container images capable of running Open FUSION Toolkit via python scripts and/or Jupyter notebooks are available on the [GitHub Package Registry](https://github.com/orgs/OpenFUSIONToolkit/packages) and can be used by [Docker](https://www.docker.com/) or other container runtimes at `ghcr.io/openfusiontoolkit/{base|jupyter|extras}` (see registry for available images/tags).

### Installation from source
The framework and its components can also be built from source using the instructions provided on [the wiki](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/wiki).

**NOTE:** To use the python interfaces for a binary or source install add the `python` directory for your installation to your `PYTHONPATH` environment variable.

Documentation
------------

Detailed documentation is [available online](https://openfusiontoolkit.github.io/OpenFUSIONToolkit/docs/index.html) and is also included in the binary packages at `/path/to/oft/install/doc/Documentation.html` (please note `/path/to/oft/install` is a placeholder for the actual install path).

Contributing
-----------

See [CONTRIBUTING.md](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/blob/main/CONTRIBUTING.md) for information on how to contribute.


Copyright
---------

Open FUSION Toolkit code development project, up to version 26.6 Copyright (c) 2023 - 2026, Open FUSION Toolkit team.

Written by Open FUSION Toolkit team and collaborators with Christopher J. Hansen as principle developer. All rights reserved.
