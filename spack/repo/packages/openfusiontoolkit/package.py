#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *


class Openfusiontoolkit(CMakePackage):
    """The Open FUSION Toolkit (OFT) is a suite of modeling tools, and their underlying finite element
    framework, for problems in plasma and fusion research and engineering in arbitrary 2D and 3D geometries.
    The underlying framework and its component tools support the use of variable order finite element methods on
    unstructured tetrahedral (triangular) or hexahedral (quadralateral) grids.
    
    Component tools:
      - TokaMaker: A time-dependent free-boundary Grad-Shafranov equilibrium code
      - ThinCurr: A 3D thin-wall electromagnetic simulation package
      - MUG: A 3D linear/nonlinear extended MHD simulation package
      - Marklin: A 3D force-free, uniform λ equilibrium solver"""

    # Project and repo information
    homepage = "https://openfusiontoolkit.github.io/OpenFUSIONToolkit/"
    url = "https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit/archive/refs/tags/v1.0.0-beta6.tar.gz"
    git = "https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit.git"
    maintainers("hansec")
    license("LGPL-3.0-only", checked_by="hansec")

    # CMakeLists.txt is in the src/ directory
    root_cmakelists_dir = "src"

    # Version list supported by spack build (needs 1.0.0-beta8 or newer)
    version("main", branch="main")
    # version("1.0.0-beta8", sha256="d687c788f05118e88b3bcb78c14fbf726286b5359b122ecdb9829fc26084f3c4") # FIXME: Uncomment and update hash when available

    # Build variants
    variant("build_type",default="Release",description="The build type to build",values=("Debug", "Release"))
    variant("examples", default=False, description="Whether to build examples")
    variant("tests", default=False, description="Whether to build tests")
    variant("python", default=True, description="Whether to build python interface")
    variant("mpi", default=False, description="Whether to enable MPI")
    variant("openmp", default=True, description="Whether to build with OpenMP support")
    variant("arpack", default=True, description="Whether to build with ARPACK support")
    variant("netcdf", default=False, description="Whether to build with NetCDF support")
    variant("umfpack", default=False, description="Whether to build with UMFPACK support through SuiteSparse")
    variant("mumps", default=False, description="Whether to build with MUMPS support")
    variant("superlu", default=False, description="Whether to build with SuperLU support")
    variant("superlu-dist", default=False, description="Whether to build with SuperLU-DIST support")

    # Language dependencies
    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")
    depends_on("python@3.7:", when="+python")

    # Build system dependencies
    depends_on("cmake@3.27:", type="build")

    # Core dependencies
    depends_on("libxml2")
    depends_on("lapack")
    depends_on("metis")
    depends_on("hdf5+fortran~mpi")

    # MPI support
    depends_on("mpi", when="+mpi")

    # LU solvers
    depends_on("suite-sparse", when="+umfpack")
    depends_on("superlu@7:", when="+superlu")
    depends_on("superlu-dist", when="+superlu-dist")
    conflicts("+superlu-dist", when="~mpi")
    with when("+mpi"):
        depends_on("mumps~complex+mpi", when="+mumps")
    with when("~mpi"):
        depends_on("mumps~complex~mpi", when="+mumps")

    # Other libraries
    depends_on("netcdf-c~mpi", when="+netcdf")
    depends_on("netcdf-fortran", when="+netcdf")
    with when("+mpi"):
        depends_on("arpack-ng+mpi", when="+arpack")
    with when("~mpi"):
        depends_on("arpack-ng~mpi", when="+arpack")

    def cmake_args(self):
        '''Generate a list of CMake arguments at build time based on the variants and dependencies'''
        args = [
            self.define_from_variant("OFT_BUILD_PYTHON", "python"),
            self.define_from_variant("OFT_BUILD_EXAMPLES", "examples"),
            self.define_from_variant("OFT_BUILD_TESTS", "tests"),
            self.define_from_variant("OFT_USE_OpenMP", "openmp"),
            self.define_from_variant("OFT_USE_MPI", "mpi"),
        ]
        if '+hl'  in self.spec["hdf5"]:
            args.append(self.define("OFT_HDF5_HL", True))
        if "+netcdf" in self.spec:
            args.append(self.define("OFT_NETCDF_ROOT", "{0}".format(self.spec["netcdf"].prefix)))
        if "arpack" in self.spec:
            args.append(self.define("OFT_ARPACK_ROOT", "{0}".format(self.spec["arpack-ng"].prefix)))
        if "umfpack" in self.spec:
            args.append(self.define("OFT_UMFPACK_ROOT", "{0}".format(self.spec["suite-sparse"].prefix)))
        if "mumps" in self.spec:
            args.append(self.define("OFT_MUMPS_ROOT", "{0}".format(self.spec["mumps"].prefix)))
        if "superlu" in self.spec:
            args.append(self.define("OFT_SUPERLU_ROOT", "{0}".format(self.spec["superlu"].prefix)))
        if "superlu-dist" in self.spec:
            args.append(self.define("OFT_SUPERLU_DIST_ROOT", "{0}".format(self.spec["superlu-dist"].prefix)))
        return args
    
    def setup_run_environment(self, env: EnvironmentModifications) -> None:
        '''Set environment variables for runtime use of Open FUSION Toolkit in Python'''
        if self.spec.satisfies("+python"):
            env.set("OFT_ROOTPATH", self.prefix)
            env.prepend_path("PYTHONPATH", self.prefix.python)
