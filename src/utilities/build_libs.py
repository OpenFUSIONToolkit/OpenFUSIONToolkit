from __future__ import print_function
import os
import shutil
import sys
import time
import subprocess
import tarfile
import re
import argparse
try:
    import urllib.request
except ImportError:
    import urllib2 as urllib
    URLError = urllib.URLError
    PY3K = False
else:
    from urllib.error import URLError
    PY3K = True


def error_exit(error_str, extra_info=None):
    # Exit build script with error
    print("\n\n============  BUILD FAILED!  ============")
    print("ERROR: {0}".format(error_str))
    if extra_info is not None:
        for info in extra_info:
            print("INFO: {0}".format(info))
    print()
    sys.exit(-1)


def fetch_file(url, file):
    def format_bytes(size):
        if size <= 0:
            return "? MB"
        elif size < 1.E6:
            return "{0:.1F} KB".format(size/1.E3)
        else:
            return "{0:.1F} MB".format(size/1.E6)
    def update_progress(curr_size, total_size, old_len, final_update=False):
        if old_len > 0:
            print(" "*old_len, end="\r")
        if final_update:
            line = "    Total size: {0}".format(format_bytes(total_size))
            print(line)
        else:
            line = "    ({0}/{1})".format(format_bytes(curr_size), format_bytes(total_size))
            print(line, end="\r")
        if not PY3K:
            sys.stdout.flush()
        return len(line)
    # Download file from url
    try:
        if PY3K:
            response = urllib.request.urlopen(url)
            resolved_url = response.geturl()
            if resolved_url != url: # Handle redirects
                response = urllib.request.urlopen(resolved_url)
            try:
                file_size = int(response.headers["content-length"])
            except:
                file_size = -1
        else:
            response = urllib.urlopen(url)
            resolved_url = response.geturl()
            if resolved_url != url: # Handle redirects
                response = urllib.request.urlopen(resolved_url)
            try:
                file_size = int(response.info().getheaders("Content-Length")[0])
            except:
                file_size = -1
    except ValueError:
        error_exit('Invalid download URL: "{0}"'.format(url))
    except:
        error_exit('Download failed for file: "{0}"'.format(url))
    else:
        line = "  Downloading: {0}".format(url)
        print(line)
        if fetch_progress:
            old_len = update_progress(0, file_size, 0)
        else:
            update_progress(file_size, file_size, 0, final_update=True)
        blocksize = max(4096, file_size//100)
        size = 0
        full_buf = b""
        while True:
            buf1 = response.read(blocksize)
            if not buf1:
                break
            full_buf += buf1
            size += len(buf1)
            if fetch_progress:
                old_len = update_progress(size, file_size, old_len)
        with open(file, "wb") as handle:
            handle.write(full_buf)
        if fetch_progress:
            update_progress(file_size, file_size, old_len, final_update=True)


def extract_archive(file):
    # Extract tarball
    try:
        with tarfile.open(file, errorlevel=2) as tar:
            tar.extractall()
    except:
        error_exit('Extraction failed for file: "{0}"'.format(file))


def run_command(command, timeout=10):
    # Run shell command
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # Wait for process to complete or timeout
    try:
        outs, _ = pid.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        pid.kill()
        outs, _ = pid.communicate()
        print("WARNING: Command timeout")
    errcode = pid.poll()
    if PY3K:
        result = outs.decode("utf-8")
    else:
        result = outs
    return result, errcode


def ver_lt(ver_string, ver_test):
    v1 = ver_test.split('.')
    v2 = ver_string.split('.')
    if int(v2[0]) < int(v1[0]):
        return True
    elif int(v2[0]) == int(v1[0]):
        if int(v2[1]) < int(v1[1]):
            return True
    return False


def ver_gt(ver_string, ver_test):
    v1 = ver_test.split('.')
    v2 = ver_string.split('.')
    if int(v2[0]) > int(v1[0]):
        return True
    elif int(v2[0]) == int(v1[0]):
        if int(v2[1]) > int(v1[1]):
            return True
    return False


def ver_range(ver_string, ver_min, ver_max):
    return ver_gt(ver_string, ver_min) and ver_lt(ver_string, ver_max)


def check_c_compiles_and_runs(source, flags, config_dict):
    def cleanup(pwd):
        os.remove('test.c')
        try:
            os.remove('test')
        except OSError:
            pass
        os.chdir(pwd)
    # Write test file
    tmp_dict = config_dict.copy()
    pwd = os.getcwd()
    os.chdir(tmp_dict['base_dir'])
    with open('test.c', "w+") as handle:
        handle.write(source)
    # Compile source
    tmp_dict["TMP_FLAGS"] = flags
    _, errcode = run_command("{CC} {TMP_FLAGS} test.c -o test".format(**tmp_dict))
    if errcode != 0:
        cleanup(pwd)
        return False
    # Run source
    _, errcode = run_command("./test")
    cleanup(pwd)
    if errcode != 0:
        return False
    else:
        return True


def setup_build_env(build_dir="build", build_cmake_ver=None):
    # Setup build environment
    # Set defaults
    config_dict = {"CC": "gcc", "CXX": "g++", "FC": "gfortran", "LD": None,
                   "AR": "ar rv", "RANLIB": "ranlib", "CMAKE": "cmake", "LD_FLAGS": "",
                   "OMP_FLAGS": "", "DEBUG_FLAGS": "", "CHK_FLAGS": "", "OPT_FLAGS": "",
                   "LOPT_FLAGS": "", "BASE_FFLAGS": "", "BASE_CFLAGS": "", "OMP_LIB": "",
                   "MKL_INC_LINES": "", "MAKE_THREADS": 1, "COMP_DEFS": [],
                   "MAKE_INC_CONTENTS": "", "INC_DEFS": [], "CC_VENDOR": "unknown",
                   "CC_VERSION": "unknown", "EXT_LIBS": ["$(OMP_LIB) -lm -lstdc++"]
                  }
    # Check for environment variables
    for key in ('CC', 'CXX', 'FC', 'LD', 'CMAKE'):
        temp = os.environ.get(key)
        if temp is not None:
            config_dict[key] = temp
    if config_dict["LD"] is None:
        config_dict["LD"] = config_dict["FC"]
    # Get build location
    config_dict['base_dir'] = os.getcwd()
    config_dict['build_dir'] = os.path.join(config_dict['base_dir'], build_dir)
    if (not os.path.isdir(config_dict['build_dir'])):
        os.mkdir(config_dict['build_dir'])
    # Check for "patch" command
    _, errcode = run_command("patch --version")
    if errcode != 0:
        error_exit('"patch" not found')
    # Check CMAKE version
    if build_cmake_ver is None:
        result, errcode = run_command("{CMAKE} --version".format(**config_dict))
        cmake_err_string = ''
        if errcode != 0:
            config_dict['CMAKE'] = None
            cmake_err_string = "specified CMAKE does not appear to work"
        else:
            try:
                line = result.split("\n")[0]
                ver_string = line.split("version")[1]
                ver_string = ver_string.split("-")[0]  # Needed if patch release
                if ver_lt(ver_string,"3.12"):
                    config_dict['CMAKE'] = None
                    cmake_err_string = "specified CMAKE version < 3.12"
                else:
                    config_dict['CMAKE_VERSION'] = ver_string
            except:
                config_dict['CMAKE'] = None
                cmake_err_string = "could not determine system CMAKE version"
        if config_dict['CMAKE'] is None:
            error_exit('CMAKE required, but {0}'.format(cmake_err_string), ('Update or retry with "--build_cmake=1" to build a compatible version',))
    else:
        config_dict['CMAKE'] = 'tobuild'
        config_dict['CMAKE_VERSION'] = build_cmake_ver
    # Check FORTRAN compiler
    result, errcode = run_command("{FC} --version".format(**config_dict))
    if errcode != 0:
        error_exit("FORTRAN compiler does not appear to work!")
    line = result.split("\n")[0]
    fc_vendor = 'unknown'
    if line.find('GNU') >= 0:
        fc_vendor = 'gnu'
    elif (line.find('ifort') >= 0) or (line.find('ifx') >= 0):
        fc_vendor = 'intel'
        if line.find('ifx') >= 0:
            if ver_lt(config_dict.get('CMAKE_VERSION','0.0'),"3.20"):
                error_exit('CMAKE >= 3.20 required for Intel "ifx" Fortran compiler', ('Update or retry with "--build_cmake=1" to build a compatible version',))
    # Check C++ compiler
    result, errcode = run_command("{CXX} --version".format(**config_dict))
    if errcode != 0:
        error_exit("C++ compiler does not appear to work!")
    # Check C compiler and determine compiler types
    result, errcode = run_command("{CC} --version".format(**config_dict))
    if errcode != 0:
        error_exit("C compiler does not appear to work!")
    if result.find('LLVM') >= 0:
        error_exit("Apple compiler wrappers detected!", ("Native GCC required on MacOS",))
    elif result.find('oneAPI DPC++/C++') >= 0:
        if ver_lt(config_dict.get('CMAKE_VERSION','0.0'),"3.20"):
            error_exit('CMAKE >= 3.20 required for Intel "icx" C/C++ compiler', ('Update or retry with "--build_cmake=1" to build a compatible version',))
    line = result.split("\n")[0]
    cc_vendor = 'unknown'
    cc_version = 'unknown'
    if line.find('gcc') >= 0:
        cc_vendor = 'gnu'
        cc_version = line.split()[-1]
    elif (line.find('icc') >= 0) or (line.find('oneAPI') >= 0):
        cc_vendor = 'intel'
    # Make sure we are using compaitble C and Fortran compilers
    if cc_vendor != fc_vendor:
        error_exit("C and FORTRAN compilers appear to be from different vendors!",
                   ("FORTRAN vendor = " + fc_vendor, "C vendor = " + cc_vendor))
    # Set default flags based on compiler vendor
    config_dict['CC_VENDOR'] = cc_vendor
    config_dict['CC_VERSION'] = cc_version
    if cc_vendor == 'gnu':
        if int(config_dict['CC_VERSION'].split(".")[0]) > 9:
            config_dict['BASE_FFLAGS'] = "-fallow-argument-mismatch"
        config_dict['OMP_FLAGS'] = "-fopenmp"
        config_dict['DEBUG_FLAGS'] = "-g"
        config_dict['CHK_FLAGS'] = "-O0 -fcheck=all"
        config_dict['OPT_FLAGS'] = "-O2"
        config_dict["LOPT_FLAGS"] = "-O2 -fPIC"
        config_dict['OMP_LIB'] = "-lgomp"
    elif cc_vendor == 'intel':
        config_dict['OMP_FLAGS'] = "-qopenmp"
        config_dict['DEBUG_FLAGS'] = "-g"
        config_dict['CHK_FLAGS'] = "-O0 -check bounds,pointers,shape,uninit"
        config_dict['OPT_FLAGS'] = ""
        config_dict["LOPT_FLAGS"] = "-fPIC"
        config_dict['OMP_LIB'] = "-liomp5"
    # Determine OS type
    import platform
    config_dict['OS_TYPE'] = platform.uname()[0]
    if config_dict['OS_TYPE'] == 'Darwin':
        result, errcode = run_command('sw_vers -productVersion')
        config_dict['OS_VER'] = result
    else:
        config_dict['OS_VER'] = platform.uname()[2]
    # Return dictionary
    return config_dict


def build_make_include(mydict):
    # Create make_libs.mk file containing library and build information
    tmp_dict = mydict.copy()
    tmp_dict['date'] = time.strftime("%c")
    tmp_dict['machine'] = os.uname()[1]
    tmp_dict['script_args'] = ' '.join(sys.argv[1:])
    tmp_dict['CC'] = mydict['HDF5_CC']
    tmp_dict['CXX'] = mydict['CXX']
    tmp_dict['FC'] = mydict['HDF5_FC']
    tmp_dict['LD'] = mydict['HDF5_FC']
    tmp_dict['COMP_DEFS_LINE'] = ' '.join(mydict['COMP_DEFS'])
    tmp_dict['INC_DEFS_LINE'] = ' '.join(reversed(mydict['INC_DEFS']))
    tmp_dict['EXT_LIBS_LINE'] = ' '.join(reversed(mydict['EXT_LIBS']))
    # Template string
    template = """# Auto-Generated on {date}
# using library build at {base_dir}
# on machine {machine}
# settings: {script_args}

# Compilers
CC = {CC}
CXX = {CXX}
FC = {FC}
LD = {LD}
AR = {AR}
RANLIB = {RANLIB}

# Compiler flags
BASE_FFLAGS = {BASE_FFLAGS}
BASE_CFLAGS = {BASE_CFLAGS}
OMP_FLAGS = {OMP_FLAGS}
DEBUG_FLAGS = {DEBUG_FLAGS}
CHK_FLAGS = {CHK_FLAGS}
OPT_FLAGS = {OPT_FLAGS}
LD_FLAGS = {LD_FLAGS}
{MKL_INC_LINES}

# Optimization flag for included 3rd party libraries (LSODE, MINPACK)
LOPT_FLAGS = {LOPT_FLAGS}

# OpenMP library
OMP_LIB = {OMP_LIB}

# Available libraries
COMP_DEFS = {COMP_DEFS_LINE}

# ==== External library specification ==== #
{MAKE_INC_CONTENTS}
# ==== Final external library/include grouping ==== #
INC_DEFS = {INC_DEFS_LINE}
EXT_LIBS = {EXT_LIBS_LINE}
\n"""
    string = template.format(**tmp_dict)
    # Output
    with open('make_libs.mk', 'w+') as fid:
        fid.write(string)


def build_cmake_script(mydict,build_debug=False,use_openmp=False,build_python=False,build_tests=False, 
                       build_examples=False,build_docs=False,build_coverage=False,package_build=False,package_release=False):
    def bool_to_string(val):
        if val:
            return "TRUE"
        else:
            return "FALSE"
    # Create "config_cmake.sh" file containing library and build information
    tmp_dict = mydict.copy()
    tmp_dict['date'] = time.strftime("%c")
    tmp_dict['machine'] = os.uname()[1]
    tmp_dict['script_args'] = ' '.join(sys.argv[1:])
    tmp_dict['LD'] = mydict['FC']
    tmp_dict['cmake_install_dir'] = "install_debug" if build_debug else "install_release"
    tmp_dict['cmake_build_dir'] = "build_debug" if build_debug else "build_release"
    cmake_lines = [
        "{CMAKE}",
        "-DCMAKE_BUILD_TYPE={0}".format("Debug" if build_debug else "Release"),
        "-DCMAKE_INSTALL_PREFIX:PATH={0}".format(os.path.join("$ROOT_PATH",tmp_dict['cmake_install_dir'])),
        "-DOFT_BUILD_TESTS:BOOL={0}".format(bool_to_string(build_tests)),
        "-DOFT_BUILD_EXAMPLES:BOOL={0}".format(bool_to_string(build_examples)),
        "-DOFT_BUILD_PYTHON:BOOL={0}".format(bool_to_string(build_python)),
        "-DOFT_BUILD_DOCS:BOOL={0}".format(bool_to_string(build_docs)),
        "-DOFT_USE_OpenMP:BOOL={0}".format(bool_to_string(use_openmp)),
        "-DOFT_PACKAGE_BUILD:BOOL={0}".format(bool_to_string(package_build)),
        "-DOFT_PACKAGE_NIGHTLY:BOOL={0}".format(bool_to_string(not package_release)),
        "-DOFT_COVERAGE:BOOL={0}".format(bool_to_string(build_coverage)),
        "-DCMAKE_C_COMPILER:FILEPATH={CC}",
        "-DCMAKE_CXX_COMPILER:FILEPATH={CXX}",
        "-DCMAKE_Fortran_COMPILER:FILEPATH={FC}"
    ]
    if mydict['BASE_CFLAGS'] != '':
        cmake_lines.append('-DCMAKE_C_FLAGS:STRING="{BASE_CFLAGS}"')
    if mydict['BASE_FFLAGS'] != '':
        cmake_lines.append('-DCMAKE_Fortran_FLAGS:STRING="{BASE_FFLAGS}"')
    if mydict['LD_FLAGS'] != '':
        cmake_lines.append('-DCMAKE_EXE_LINKER_FLAGS:STRING="{LD_FLAGS}"')
    have_mpi = False
    if "MPI_ROOT" in mydict:
        have_mpi = True
        cmake_lines.append("-DOFT_USE_MPI:BOOL=TRUE")
        cmake_lines.append("-DMPI_HOME:PATH={0}".format(mydict["MPI_ROOT"]))
    else:
        if "MPI_CC" in mydict:
            have_mpi = True
            cmake_lines.append("-DOFT_USE_MPI:BOOL=TRUE")
            cmake_lines.append("-DMPI_C_COMPILER:PATH={0}".format(mydict["MPI_CC"]))
        if "MPI_CXX" in mydict:
            cmake_lines.append("-DMPI_CXX_COMPILER:PATH={0}".format(mydict["MPI_CXX"]))
        if "MPI_FC" in mydict:
            cmake_lines.append("-DMPI_Fortran_COMPILER:PATH={0}".format(mydict["MPI_FC"]))
    if not have_mpi:
        cmake_lines.append("-DOFT_USE_MPI:BOOL=FALSE")
    if "PETSC_ROOT" in mydict:
        cmake_lines.append("-DOFT_PETSc_ROOT:PATH={0}".format(mydict["PETSC_ROOT"]))
    else:
        if "METIS_ROOT" in mydict:
            cmake_lines.append("-DOFT_METIS_ROOT:PATH={0}".format(mydict["METIS_ROOT"]))
    if "HDF5_ROOT" in mydict:
        cmake_lines.append("-DHDF5_ROOT:PATH={0}".format(mydict["HDF5_ROOT"]))
    if "NETCDF_ROOT" in mydict:
        cmake_lines.append("-DOFT_NETCDF_ROOT:PATH={0}".format(mydict["NETCDF_ROOT"]))
    if "BLAS_ROOT" in mydict:
        cmake_lines.append("-DBLAS_ROOT:PATH={0}".format(mydict["BLAS_ROOT"]))
        cmake_lines.append("-DLAPACK_ROOT:PATH={0}".format(mydict["LAPACK_ROOT"]))
        cmake_lines.append("-DBLA_VENDOR:STRING={0}".format(mydict["BLAS_VENDOR"]))
    elif "MKL_FLAGS" in mydict:
        cmake_lines.append("-DBLA_VENDOR:STRING={0}".format(mydict["BLAS_VENDOR"]))
    if "ARPACK_ROOT" in mydict:
        cmake_lines.append("-DOFT_ARPACK_ROOT:PATH={0}".format(mydict["ARPACK_ROOT"]))
    if "FOX_ROOT" in mydict:
        cmake_lines.append("-DOFT_FoX_ROOT:PATH={0}".format(mydict["FOX_ROOT"]))
    if "ONURBS_ROOT" in mydict:
        cmake_lines.append("-DOFT_OpenNURBS_ROOT:PATH={0}".format(mydict["ONURBS_ROOT"]))
    if "PETSC_ROOT" in mydict:
        cmake_lines.append("-DOFT_PETSc_ROOT:PATH={0}".format(mydict["PETSC_ROOT"]))
    else:
        if "SUPERLU_ROOT" in mydict:
            cmake_lines.append("-DOFT_SUPERLU_ROOT:PATH={0}".format(mydict["SUPERLU_ROOT"]))
            if "SUPERLU_VER_MAJOR" in mydict:
                cmake_lines.append('-DOFT_SUPERLU_VER_MAJOR:STRING="{0}"'.format(mydict["SUPERLU_VER_MAJOR"]))
        if "SUPERLU_DIST_ROOT" in mydict:
            cmake_lines.append("-DOFT_SUPERLU_DIST_ROOT:PATH={0}".format(mydict["SUPERLU_DIST_ROOT"]))
        if "UMFPACK_ROOT" in mydict:
            cmake_lines.append("-DOFT_UMFPACK_ROOT:PATH={0}".format(mydict["UMFPACK_ROOT"]))
    if tmp_dict['OS_TYPE'] == 'Darwin':
        try:
            result, errcode = run_command("pkgutil --pkg-info=com.apple.pkg.CLTools_Executables | grep version")
            result = result.split(':')[1]
        except:
            errcode = -1
        if errcode == 0:
            if ver_range(result, '14.2', '15.3'):
                cmake_lines.append('-DCMAKE_EXE_LINKER_FLAGS:STRING="-Wl,-ld_classic"')
                cmake_lines.append('-DCMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,-ld_classic"')
        else:
            print('Error detecting macOS version using "sw_vers --productVersion"')
    cmake_lines += [os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))]
    cmake_lines_str = ' \\\n  '.join(cmake_lines)
    # Template string
    template = """# Auto-Generated on {date}
# using library build at {base_dir}
# on machine {machine}
# settings: {script_args}

# Get current directory
ROOT_PATH=$(pwd)

# Create fresh build directory
rm -rf {cmake_build_dir}
mkdir {cmake_build_dir} && cd {cmake_build_dir}

""" + cmake_lines_str + '\n'
    string = template.format(**tmp_dict)
    # Output
    with open('config_cmake.sh', 'w+') as fid:
        fid.write(string)


class package:
    name = "Unitialized"
    url = None
    version = None
    file = None
    build_dir = None
    install_dir = None
    root_path = None
    root_build_path = None
    skip = False
    build_timeout = 15
    extra_fetch = []
    patch_files = []
    install_chk_files = []
    config_dict = {}
    children = []

    def setup(self, config_dict):
        error_exit("No setup method defined for package!")
        return config_dict

    def post_child_setup(self, config_dict):
        return config_dict

    def build(self):
        error_exit("No build method defined for package!")

    def install(self, config_dict, force=False):
        def handle_children(config_dict):
            for child in self.children:
                config_dict = child.install(config_dict, force)
            return self.post_child_setup(config_dict)
        #
        if self.url is not None:
            # error_exit("No download URL specified for package: " + self.name)
            self.file = self.url.split("/")[-1]
            if self.build_dir is None:
                self.build_dir = os.path.splitext(self.file)[0]
                if os.path.splitext(self.build_dir)[-1] == '.tar':
                    self.build_dir = os.path.splitext(self.build_dir)[0]
            #
            if self.install_dir is None:
                self.install_dir = re.sub(r"[\.]", "_", self.build_dir)
            self.root_path = config_dict['base_dir']
            self.root_build_path = config_dict['build_dir']
        config_dict = self.setup(config_dict)
        if self.install_dir is not None:
            if os.path.isdir(self.install_dir) and (not force):
                config_dict = self.post_install(config_dict)
                if self.check_install(do_abort=False):
                    return handle_children(config_dict)
        if self.skip:
            return handle_children(config_dict)
        print("=========================================================")
        print("Building library: {0}".format(self.name))
        #
        os.chdir(self.root_build_path)
        self.fetch(force)
        if config_dict['DOWN_ONLY']:
            return handle_children(config_dict)
        if (not os.path.isdir(self.build_dir)) or force:
            self.extract()
        else:
            print("  Using existing folder: {0}".format(self.build_dir))
        #
        os.chdir(self.build_dir)
        print("  Excecuting build (this may take a few minutes)")
        self.build()
        os.chdir(self.root_path)
        # Check to make sure Installation succeeded
        config_dict = self.post_install(config_dict)
        self.check_install()
        return handle_children(config_dict)

    def post_install(self, config_dict):
        return config_dict

    def fetch(self, force=False):
        # Download library archive
        if self.skip:
            return
        if (not os.path.isfile(self.file)) or force:
            fetch_file(self.url, self.file)
        else:
            print("  Using existing file: {0}".format(self.file))
        for args in self.extra_fetch + self.patch_files:
            url = args[0]
            if len(args) == 1:
                tmp_file = url.split("/")[-1]
            elif len(args) == 2:
                tmp_file = args[1]
            else:
                error_exit('Invalid "extra_fetch" or "patch_file" object!')
            if (not os.path.isfile(tmp_file)) or force:
                fetch_file(url, tmp_file)
            else:
                print("  Using existing file: {0}".format(tmp_file))

    def extract(self):
        # Extract archive containing library
        if self.skip:
            return
        print("  Extracting file: {0}".format(self.file))
        extract_archive(self.file)
        # Apply patches
        for patch in self.patch_files:
            if len(patch) == 1:
                tmp_file = patch[0].split("/")[-1]
            elif len(patch) == 2:
                tmp_file = patch[1]
            os.chdir(self.build_dir)
            result, errcode = run_command("patch -N -p0 < {0}".format(os.path.join("..", tmp_file)))
            os.chdir("..")
            if errcode == 0:
                print('  Applied patch "{0}"'.format(tmp_file))
            else:
                shutil.rmtree(self.build_dir)
                with open("{0}.log".format(tmp_file), "w+") as fid:
                    fid.write(result)
                error_exit('Failed to apply patch "{0}"'.format(tmp_file),
                           ('See "{0}.log" for more information'.format(os.path.join(self.root_build_path, tmp_file)),))

    def setup_root_struct(self, lib_path="lib", inc_path="include", bin_path="bin"):
        # Setup default directory structure
        self.config_dict[self.name + '_ROOT'] = os.path.join(self.root_path, self.install_dir)
        self.config_dict[self.name + '_LIB'] = os.path.join(self.config_dict[self.name + '_ROOT'], lib_path)
        self.config_dict[self.name + '_INCLUDE'] = os.path.join(self.config_dict[self.name + '_ROOT'], inc_path)
        self.config_dict[self.name + '_BIN'] = os.path.join(self.config_dict[self.name + '_ROOT'], bin_path)

    def run_build(self, build_lines, config_dict):
        # Run build script
        # Write build script to file and add step status checks
        script_lines = ['build_base=$(pwd)', 'rm -f $build_base/build_tmp.stat']
        for line in build_lines:
            line_stripped = line.strip()
            if line_stripped != "":
                if line_stripped.startswith("export ") or line_stripped.endswith("\\"):
                    script_lines += ['{0}'.format(line.format(**config_dict))]
                else:
                    script_lines += ['{0}; echo $? >> $build_base/build_tmp.stat'.format(line.format(**config_dict))]
        with open("build_tmp.sh", "w+") as fid:
            fid.write("\n".join(script_lines))
        if self.config_dict['SETUP_ONLY']:
            return
        result, _ = run_command("bash build_tmp.sh", timeout=self.build_timeout*60)
        with open("build_tmp.log", "w+") as fid:
            fid.write(result)
        # Check for build success
        build_success = True
        with open("build_tmp.stat", "r") as fid:
            for line in fid:
                build_success = build_success and (int(line) == 0)
        if not build_success:
            log_path = os.path.join(os.path.join(self.root_build_path, self.build_dir),'build_tmp.log')
            shutil.copy(log_path,os.path.join(self.root_build_path, 'build_error.log'))
            error_exit("Build failed for package: " + self.name,
                       ("See '{0}' for more information".format(log_path),))

    def check_install(self, do_abort=True):
        # Check to make sure Installation succeeded
        # Look for each expected install files
        if self.config_dict['SETUP_ONLY']:
            return True
        for chk_file in self.install_chk_files:
            if not os.path.isfile(chk_file):
                if do_abort:
                    log_path = os.path.join(os.path.join(self.root_build_path, self.build_dir),'build_tmp.log')
                    shutil.copy(log_path,os.path.join(self.root_build_path, 'build_error.log'))
                    error_exit("Installation failed for package: " + self.name,
                               ("Did not find expected install file at path '{0}'".format(chk_file),
                                "See '{0}' for more information".format(log_path)
                                )
                               )
                else:
                    return False
        return True


class CMAKE(package):
    def __init__(self):
        self.name = "CMAKE"
        self.url = "https://github.com/Kitware/CMake/releases/download/v3.27.9/cmake-3.27.9.tar.gz"
        self.version = '3.27'

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        self.setup_root_struct()
        # Installation check files
        self.install_chk_files = [os.path.join(self.config_dict['CMAKE_BIN'], 'cmake')]
        # Replace CMAKE executable
        self.config_dict['CMAKE'] = os.path.join(self.config_dict['CMAKE_BIN'], 'cmake')
        self.config_dict['CMAKE_VERSION'] = self.version
        return self.config_dict

    def build(self):
        build_lines = [
            "rm -rf build",
            "mkdir build",
            "cd build",
            "unset CC",
            "unset CXX",
            "unset FC",
            "../bootstrap -- -DCMAKE_USE_OPENSSL=OFF -DCMAKE_INSTALL_PREFIX={CMAKE_ROOT}",
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, self.config_dict)


class METIS(package):
    def __init__(self, comp_wrapper=False):
        self.name = "METIS"
        self.url = "http://faculty.washington.edu/hansec/libs/metis-5.1.0-mod.tar.gz"
        self.build_dir = "metis-5.1.0"
        self.comp_wrapper = comp_wrapper

    def detect_sizes(self):
        print("  Testing METIS sizes")
        with open('tmp.cc', 'w+') as fid:
            fid.write('#include "metis.h"\n\nTEST_INDWIDTH IDXTYPEWIDTH\nTEST_REALWIDTH REALTYPEWIDTH')
        result, _ = run_command("{CC} -E tmp.cc | tail -2".format(**config_dict))
        os.remove('tmp.cc')
        values = result.split()
        idx_size = 4
        real_size = 4
        if int(values[1]) == 64:
            idx_size = 8
        if int(values[3]) == 64:
            real_size = 8
        print("    idx_width  = {0}".format(values[1]))
        print("    real_width = {0}".format(values[3]))
        return idx_size, real_size

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        if self.comp_wrapper:
            self.skip = True
            print("Detecting METIS compiler wrapper:")
            try:
                self.detect_sizes()
            except:
                self.skip = False
                print("  Failed, building")
            else:
                print("  METIS found, skiping build")
        if self.skip:
            return self.config_dict
        self.setup_root_struct()
        self.config_dict['METIS_LIBS'] = "-lmetis"
        # Installation check files
        self.install_chk_files = [os.path.join(self.config_dict['METIS_LIB'], 'libmetis.a')]
        #
        self.config_dict['INC_DEFS'].append('$(METIS_INC)')
        self.config_dict['EXT_LIBS'].append('$(METIS_LIB)')
        template = '\n# METIS library\nMETIS_INC = -I{METIS_INCLUDE}\nMETIS_LIB = -L{METIS_LIB} {METIS_LIBS}\n'
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        return self.config_dict

    def build(self):
        build_dir = os.getcwd()
        build_lines = [
            "GKLIB_PATH={0}".format(os.path.join(build_dir,"GKlib")),
            "rm -rf build",
            "mkdir -p build",
            "cd build",
            "{CMAKE} -DCMAKE_INSTALL_PREFIX={METIS_ROOT} -DCMAKE_C_COMPILER={CC} -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON -DGKLIB_PATH=$GKLIB_PATH .. ",
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, self.config_dict)


class MPI(package):
    def __init__(self):
        self.name = "MPI"
        self.url = "http://faculty.washington.edu/hansec/libs/mpich-3.3.2.tar.gz"
        self.build_timeout = 20

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        self.config_dict['COMP_DEFS'].append("-DHAVE_MPI")
        if "MPI_CC" in config_dict:
            self.skip = True
            print("MPI provided by compiler wrappers: Skipping build")
            return self.config_dict
        if "MPI_LIBS" in config_dict:
            template = '\n'.join([
                '',
                '# MPI library',
                'MPI_INC = -I{MPI_INCLUDE}',
                'MPI_LIB = -L{MPI_LIB} {MPI_LIBS}',
                ''
            ])
            self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
            self.config_dict['INC_DEFS'].append('$(MPI_INC)')
            self.config_dict['EXT_LIBS'].append('$(MPI_LIB)')
            self.skip = True
            print("MPI libraries specified: Skipping build")
            return self.config_dict
        self.setup_root_struct()
        bin_dir = self.config_dict['MPI_BIN']
        self.config_dict['MPI_CC'] = os.path.join(bin_dir, "mpicc")
        self.config_dict['MPI_FC'] = os.path.join(bin_dir, "mpif90")
        print('To use MPI please add "{0}" to your path'.format(bin_dir))
        # Installation check files
        self.install_chk_files = [self.config_dict['MPI_CC'], self.config_dict['MPI_FC']]
        return self.config_dict

    def build(self):
        build_lines = [
            "rm -rf build",
            "mkdir build",
            "cd build",
            "export CC={CC}",
            "export FC={FC}"]
        if config_dict['CC_VENDOR'] == 'gnu' and int(config_dict['CC_VERSION'].split(".")[0]) > 9:
            build_lines.append('export FFLAGS=-fallow-argument-mismatch')
        build_lines += [
            "../configure --prefix={MPI_ROOT} --enable-fortran=yes --enable-shared=no --with-pic --disable-opencl --disable-nvml --disable-cuda",
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, self.config_dict)

    def post_install(self, config_dict):
        # Add MPI bin directory to path for following builds
        os.environ["PATH"] = "{0}:".format(self.config_dict['MPI_BIN']) + os.environ.get("PATH", "")
        return config_dict


class HDF5(package):
    def __init__(self, parallel=False):
        self.name = "HDF5"
        # self.url = "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.2/src/hdf5-1.14.2.tar.gz"
        self.url = "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.10/src/hdf5-1.10.10.tar.gz"
        self.parallel = parallel

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        if "HDF5_CC" in config_dict:
            if self.parallel:
                self.config_dict['COMP_DEFS'].append("-DHAVE_PHDF5")
            self.skip = True
            print("HDF5 provided by compiler wrappers: Skipping build")
            return self.config_dict
        self.setup_root_struct()
        bin_dir = self.config_dict['HDF5_BIN']
        if self.parallel:
            self.config_dict['HDF5_CC'] = os.path.join(bin_dir, "h5pcc")
            self.config_dict['HDF5_FC'] = os.path.join(bin_dir, "h5pfc")
            self.config_dict['COMP_DEFS'].append("-DHAVE_PHDF5")
        else:
            self.config_dict['HDF5_CC'] = os.path.join(bin_dir, "h5cc")
            self.config_dict['HDF5_FC'] = os.path.join(bin_dir, "h5fc")
        # Installation check files
        self.install_chk_files = [self.config_dict['HDF5_CC'], self.config_dict['HDF5_FC']]
        return self.config_dict

    def build(self):
        build_lines = [
            "rm -rf build",
            "mkdir build",
            "cd build"
        ]
        # cmake_options = [
        #     '-DCMAKE_INSTALL_PREFIX:PATH="{HDF5_ROOT}"',
        #     '-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON',
        #     '-DHDF5_BUILD_FORTRAN:BOOL=ON',
        #     '-DBUILD_SHARED_LIBS:BOOL=OFF',
        #     '-DHDF5_BUILD_CPP_LIB=OFF',
        #     '-DBUILD_TESTING=OFF',
        #     '-DHDF5_BUILD_EXAMPLES=OFF'
        # ]
        options = [
            "--prefix={HDF5_ROOT}",
            "--enable-fortran",
            "--enable-shared=no",
            "--disable-tests",
            "--disable-examples",
            "--with-pic"
        ]
        if "MPI_CC" in self.config_dict:
            build_lines += [
                "export CC={MPI_CC}",
                "export FC={MPI_FC}"
            ]
        else:
            build_lines += [
                "export CC={CC}",
                "export FC={FC}"
            ]
        if self.parallel:
            # cmake_options += ["-DHDF5_ENABLE_PARALLEL:BOOL=ON"]
            options += ["--enable-parallel"]
            if "MPI_LIBS" in self.config_dict:
                build_lines += [
                    'export CPPFLAGS="-I{MPI_INCLUDE}"',
                    'export FCFLAGS="-I{MPI_INCLUDE}"',
                    'export LDFLAGS="-L{MPI_LIB}"',
                    'export LIBS="{MPI_LIBS}"'
                ]
        #
        build_lines += [
            # "cmake " + " ".join(cmake_options) + " ..",
            "../configure " + " ".join(options),
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, self.config_dict)


class NETCDF(package):
    def __init__(self, comp_wrapper=False):
        self.name = "NETCDF"
        self.url = "http://faculty.washington.edu/hansec/libs/netcdf-c-4.6.2.tar.gz"
        self.install_dir = "netcdf-4_6_2"
        self.comp_wrapper = comp_wrapper
        self.children = [NETCDF_Fortran(comp_wrapper)]

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        if self.comp_wrapper:
            self.config_dict['COMP_DEFS'].append("-DHAVE_NCDF")
            self.skip = True
            print("NETCDF provided by compiler wrappers: Skipping build")
            return self.config_dict
        self.setup_root_struct()
        self.config_dict['NETCDF_C_LIBS'] = "-lnetcdf"
        self.config_dict['COMP_DEFS'].append("-DHAVE_NCDF")
        # Installation check files
        self.install_chk_files = [os.path.join(self.config_dict['NETCDF_LIB'], 'libnetcdf.a')]
        #
        self.config_dict['INC_DEFS'].append('$(NETCDF_INC)')
        self.config_dict['EXT_LIBS'].append('$(NETCDF_LIB)')
        return self.config_dict

    def post_child_setup(self, config_dict):
        if self.comp_wrapper:
            return self.config_dict
        self.config_dict['NETCDF_LIBS'] = ' '.join([
            self.children[0].config_dict['NETCDF_Fortran_LIBS'],
            self.config_dict['NETCDF_C_LIBS']
        ])
        template = '\n'.join([
            '\n# NETCDF library',
            'NETCDF_INC = -I{NETCDF_INCLUDE}',
            'NETCDF_LIB = -L{NETCDF_LIB} {NETCDF_LIBS}\n'
        ])
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        return self.config_dict

    def build(self):
#         build_script = """rm -rf build
# mkdir build && cd build
# export CC={CC}
# export FC={FC}
# {CMAKE} \
#   -DCMAKE_INSTALL_PREFIX:PATH={NETCDF_ROOT} \
#   -DHDF5_ROOT:PATH={HDF5_ROOT} \
#   -DBUILD_SHARED_LIBS:BOOL=FALSE \
#   -DENABLE_DAP:BOOL=FALSE \
#   ..
# make -j{MAKE_THREADS}
# make install""".format(**self.config_dict)
#         self.run_build(build_script)
        build_lines = [
            "rm -rf build",
            "mkdir build",
            "cd build",
            "export CC={HDF5_CC}",
            "export FC={HDF5_FC}",
            "../configure --prefix={NETCDF_ROOT} --enable-netcdf-4 --enable-shared=no --with-pic --disable-dap",
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, self.config_dict)


class NETCDF_Fortran(package):
    def __init__(self, comp_wrapper=False):
        self.name = "NETCDF_Fortran"
        self.url = "http://faculty.washington.edu/hansec/libs/netcdf-fortran-4.4.4.tar.gz"
        self.install_dir = "netcdf-4_6_2"
        self.comp_wrapper = comp_wrapper

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        if self.comp_wrapper:
            self.skip = True
            return self.config_dict
        self.setup_root_struct()
        self.config_dict['NETCDF_Fortran_LIBS'] = "-lnetcdff"
        # Installation check files
        self.install_chk_files = [os.path.join(self.config_dict['NETCDF_Fortran_LIB'], 'libnetcdff.a')]
        return self.config_dict

    def build(self):
#         build_script = """rm -rf build
# mkdir build && cd build
# export CC={CC}
# export FC={FC}
# {CMAKE} \
#   -DCMAKE_INSTALL_PREFIX:PATH={NETCDF_ROOT} \
#   -DBUILD_SHARED_LIBS:BOOL=FALSE \
#   ..
# make -j{MAKE_THREADS}
# make install""".format(**self.config_dict)
#         self.run_build(build_script)
        build_lines = [
            "rm -rf build",
            "mkdir build",
            "cd build",
            "export CC={CC}",
            "export FC={FC}",
            'export CPPFLAGS="-I{NETCDF_INCLUDE}"',
            'export LDFLAGS="-L{NETCDF_LIB}"',
            'export LIBS="{NETCDF_C_LIBS}"']
        if config_dict['CC_VENDOR'] == 'gnu' and int(config_dict['CC_VERSION'].split(".")[0]) > 9:
            build_lines.append('export FCFLAGS="-fallow-argument-mismatch -I{NETCDF_INCLUDE}"')
        else:
            build_lines.append('export FCFLAGS="-I{NETCDF_INCLUDE}"')
        build_lines += [
            "../configure --prefix={NETCDF_Fortran_ROOT} --enable-shared=no --with-pic",
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, self.config_dict)


class OpenBLAS(package):
    def __init__(self, build_threaded=False, dynamic_arch=False):
        self.name = "OpenBLAS"
        self.url = "https://github.com/xianyi/OpenBLAS/archive/refs/tags/v0.3.23.tar.gz"
        self.build_dir = "OpenBLAS-0.3.23"
        self.install_dir = "OpenBLAS-0_3_23"
        self.threaded = build_threaded
        self.dynamic_arch = dynamic_arch

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        self.setup_root_struct()
        self.config_dict['BLAS_ROOT'] = self.config_dict['OpenBLAS_ROOT']
        self.config_dict['LAPACK_ROOT'] = self.config_dict['OpenBLAS_ROOT']
        self.config_dict['BLAS_VENDOR'] = "OpenBLAS"
        self.config_dict['LAPACK_LIB'] = self.config_dict[self.name + '_LIB']
        self.config_dict['LAPACK_LIBS'] = "-lopenblas"
        self.config_dict['BLAS_LIB_PATH'] = os.path.join(self.config_dict['LAPACK_LIB'], 'libopenblas.a')
        self.config_dict['LAPACK_LIB_PATH'] = self.config_dict['BLAS_LIB_PATH']
        if self.threaded:
            self.config_dict['COMP_DEFS'].append("-DBLAS_THREADED")
            self.config_dict['COMP_DEFS'].append("-DLAPACK_THREADED")
            self.config_dict['OpenBLAS_THREADS'] = True
        else:
            self.config_dict['OpenBLAS_THREADS'] = False
        # Installation check files
        self.install_chk_files = [self.config_dict['BLAS_LIB_PATH']]
        #
        self.config_dict['EXT_LIBS'].append('$(LAPACK_LIB)')
        template = '\n# LAPACK library\nLAPACK_LIB = -L{LAPACK_LIB} {LAPACK_LIBS}\n'
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        return self.config_dict

    def build(self):
        avx_test_source = """
#include <immintrin.h>
int main(int argc, char** argv) {
  double a[4] = {1.0, 2.0, 3.0, 4.0};
  __m256d t = _mm256_loadu_pd(a);
  return 0;
}
        """
        avx2_test_source = """
#include <immintrin.h>
int main(int argc, char** argv) {
  double a[4] = {1.0, 2.0, 3.0, 4.0};
  __m128i vindex = _mm_set_epi32(0, 2, 1, 3);
  __m256d t = _mm256_i32gather_pd(a, vindex, 8);
  return 0;
}
        """
        tmp_dict = self.config_dict.copy()
        oblas_options = ['NO_CBLAS=1', 'NO_LAPACKE=1', 'NO_SHARED=1']
        if self.config_dict['MAKE_THREADS'] == 1:
            oblas_options += ['NO_PARALLEL_MAKE=1']
        else:
            oblas_options += ['MAKE_NB_JOBS={MAKE_THREADS}']
        if self.threaded:
            oblas_options += ['USE_THREAD=1', 'USE_OPENMP=1', 'FCOMMON_OPT="-frecursive {OMP_FLAGS} -fPIC"']
        else:
            oblas_options += ['USE_THREAD=0', 'FCOMMON_OPT="-frecursive -fPIC"']
        if self.config_dict['OS_TYPE'] == 'Darwin':
            if not check_c_compiles_and_runs(avx_test_source, "-mavx", tmp_dict):
                oblas_options += ['NO_AVX=1']
            if not check_c_compiles_and_runs(avx2_test_source, "-mavx2", tmp_dict):
                oblas_options += ['NO_AVX2=1']
        if self.dynamic_arch:
            oblas_options += ['DYNAMIC_ARCH=1']
        build_lines = [
            'export CC={CC}',
            'export FC={FC}',
            'make clean',
            'make {0}'.format(' '.join(oblas_options)),
            'make NO_PARALLEL_MAKE=1 NO_SHARED=1 PREFIX={OpenBLAS_ROOT} install'
        ]
        self.run_build(build_lines, tmp_dict)


class MKL(package):
    def __init__(self, mkl_root=None):
        self.name = "MKL"
        self.skip = True
        if mkl_root is None:
            error_exit("MKL root path required when using MKL")
        self.mkl_root = mkl_root

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        print("Using MKL library:")
        self.config_dict['BLAS_VENDOR'] = "Intel10_64lp"
        if self.config_dict['CC_VENDOR'] == 'intel':
            print("  Intel compiler flag = {0}".format("-mkl"))
            self.config_dict['MKL_FLAGS'] = "-mkl"
        else:
            iomp_path = os.path.normpath(os.path.join(self.mkl_root, "../lib"))
            mkl_libs = ("libmkl_intel_lp64.a", "libmkl_intel_thread.a", "libmkl_core.a")
            flag_string = ""
            for mkl_lib in mkl_libs:
                flag_string += "$(MKL_ROOT)/lib/{0} ".format(mkl_lib)
            self.config_dict['MKL_FLAGS'] = flag_string[:-1] #"$(MKL_ROOT)/lib/libmkl_intel_lp64.a $(MKL_ROOT)/lib/libmkl_intel_thread.a $(MKL_ROOT)/lib/libmkl_core.a"
            self.config_dict['OMP_LIB'] = "-L{0} -liomp5 -lpthread -ldl".format(iomp_path)
            self.config_dict['BLAS_ROOT'] = self.mkl_root
            self.config_dict['LAPACK_ROOT'] = self.mkl_root
            print("  MKL libraries    = {0}".format(" ".join(mkl_libs)))
            print("  Intel OpenMP lib = {0}".format(self.config_dict['OMP_LIB']))
        self.config_dict['COMP_DEFS'].append("-DHAVE_MKL")
        self.config_dict['MKL_ROOT'] = self.mkl_root
        self.config_dict['MKL_INC_LINES'] = "MKL_ROOT = {0}\n".format(self.mkl_root)
        self.config_dict['MKL_INC_LINES'] += "MKL_FLAGS = {0}".format(self.config_dict['MKL_FLAGS'])
        self.config_dict['INC_DEFS'].append('-I$(MKL_ROOT)/include')
        return self.config_dict


class BLAS_LAPACK(package):
    def __init__(self, comp_wrapper=False, blas_lib_path=None, lapack_lib_path=None):
        self.name = "BLAS_LAPACK"
        self.url = "http://faculty.washington.edu/hansec/libs/lapack-3.5.0.tgz"
        self.comp_wrapper = comp_wrapper
        self.blas_lib_path = blas_lib_path
        self.lapack_lib_path = lapack_lib_path
        if (self.blas_lib_path is not None) and (self.lapack_lib_path is not None):
            if not os.path.isfile(self.blas_lib_path):
                error_exit('Specified BLAS library "{0}" does not exist'.format(self.blas_lib_path))
            if not os.path.isfile(self.lapack_lib_path):
                error_exit('Specified LAPACK library "{0}" does not exist'.format(self.lapack_lib_path))

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        if self.comp_wrapper:
            self.config_dict['BLAS_LIB_PATH'] = "-lm"
            self.config_dict['LAPACK_LIB_PATH'] = "-lm"
            self.skip = True
            print("BLAS/LAPACK provided by compiler wrappers: Skipping build")
            return self.config_dict
        if (self.blas_lib_path is not None) and (self.lapack_lib_path is not None):
            self.config_dict['BLAS_LIB_PATH'] = self.blas_lib_path
            self.config_dict['LAPACK_LIB_PATH'] = self.lapack_lib_path
            self.config_dict['EXT_LIBS'].append('$(LAPACK_LIB)')
            template = '\n# LAPACK library\nLAPACK_LIB = {LAPACK_LIB_PATH} {BLAS_LIB_PATH}\n'
            self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
            self.skip = True
            print("Pre-built BLAS/LAPACK provided: Skipping build")
            return self.config_dict
        self.setup_root_struct()
        self.config_dict['BLAS_ROOT'] = self.config_dict['BLAS_LAPACK_ROOT']
        self.config_dict['LAPACK_ROOT'] = self.config_dict['BLAS_LAPACK_ROOT']
        self.config_dict['BLAS_VENDOR'] = "Generic"
        self.config_dict['LAPACK_LIB'] = self.config_dict[self.name + '_LIB']
        self.config_dict['LAPACK_LIBS'] = "-llapack -lblas"
        self.config_dict['BLAS_LIB_PATH'] = os.path.join(self.config_dict['LAPACK_LIB'], 'libblas.a')
        self.config_dict['LAPACK_LIB_PATH'] = os.path.join(self.config_dict['LAPACK_LIB'], 'liblapack.a')
        # Installation check files
        self.install_chk_files = [self.config_dict['BLAS_LIB_PATH'], self.config_dict['LAPACK_LIB_PATH']]
        #
        self.config_dict['EXT_LIBS'].append('$(LAPACK_LIB)')
        template = '\n# LAPACK library\nLAPACK_LIB = -L{LAPACK_LIB} {LAPACK_LIBS}\n'
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        return self.config_dict

    def build(self):
        # -DBLAS_LIBRARIES:PATH=
        fflags = []
        if (config_dict['CC_VENDOR'] == 'gnu') and (int(config_dict['CC_VERSION'].split(".")[0]) > 9):
            fflags.append("-fallow-argument-mismatch")
        build_lines = [
            "rm -rf build",
            "mkdir build",
            "cd build",
            "export CC={CC}",
            "export FC={FC}"
        ]
        if len(fflags) > 0:
            build_lines += [
                'export FFLAGS="{0}"'.format(" ".join(fflags)),
                'export FCFLAGS="{0}"'.format(" ".join(fflags))
            ]
        build_lines += [    
            "{CMAKE} -DCMAKE_INSTALL_PREFIX:PATH={BLAS_LAPACK_ROOT} -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON -Wno-dev -DCMAKE_BUILD_TYPE=Release ..",
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, self.config_dict)


class ARPACK(package):
    def __init__(self, parallel=False, link_omp=False):
        self.name = "ARPACK"
        self.url = "http://faculty.washington.edu/hansec/libs/arpack-ng-3.5.0.tar.gz"
        self.parallel = parallel
        self.link_omp = link_omp

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        self.setup_root_struct()
        if self.parallel:
            self.config_dict['ARPACK_LIBS'] = "-lparpack -larpack"
        else:
            self.config_dict['ARPACK_LIBS'] = "-larpack"
        self.config_dict['COMP_DEFS'].append("-DHAVE_ARPACK")
        # Installation check files
        self.install_chk_files = [os.path.join(self.config_dict['ARPACK_LIB'], 'libarpack.a')]
        #
        self.config_dict['EXT_LIBS'].append('$(ARPACK_LIB)')
        template = '\n# ARPACK library\nARPACK_LIB = -L{ARPACK_LIB} {ARPACK_LIBS}\n'
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        return self.config_dict

    def build(self):
        tmp_dict = self.config_dict.copy()
        options = [
            "-DCMAKE_INSTALL_PREFIX:PATH={ARPACK_ROOT}",
            "-DCMAKE_INSTALL_LIBDIR=lib",
            "-DEXAMPLES=OFF",
            "-DBUILD_SHARED_LIBS=OFF",
            "-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
        ]
        if "BLAS_LIB_PATH" in self.config_dict:
            options += [
                '-DBLAS_LIBRARIES:PATH={BLAS_LIB_PATH}',
                '-DLAPACK_LIBRARIES:PATH={LAPACK_LIB_PATH}'
            ]
        build_lines = [
            "rm -rf build",
            "mkdir build",
            "cd build"
        ]
        fflags = []
        if (config_dict['CC_VENDOR'] == 'gnu') and (int(config_dict['CC_VERSION'].split(".")[0]) > 9):
            fflags.append("-fallow-argument-mismatch")
        if self.link_omp:
            fflags.append("{0}".format(tmp_dict['OMP_FLAGS']))
        if self.parallel:
            options.append("-DMPI=ON")
            if "MPI_CC" in tmp_dict:
                build_lines += [
                    "export CC={MPI_CC}",
                    "export F77={MPI_FC}",
                    "export FC={MPI_FC}",
                    "export MPIF77={MPI_FC}"
                ]
            else:
                build_lines += [
                    "export CC={CC}",
                    "export F77={FC}",
                    "export FC={FC}",
                    'export LDFLAGS="-L{MPI_LIB}"',
                    'export LIBS="{MPI_LIBS}"',
                    'export CPPFLAGS="-I{MPI_INCLUDE}"'
                ]
        else:
            options.append("-DMPI=OFF")
            build_lines += [
                "export CC={CC}",
                "export F77={FC}",
                "export FC={FC}"
            ]
        #
        build_lines += [
            'export FFLAGS="{0}"'.format(" ".join(fflags)),
            "{CMAKE} " + " ".join(options) + " ..",
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, tmp_dict)


class SUPERLU(package):
    def __init__(self, comp_wrapper=False):
        self.name = "SUPERLU"
        self.url = "http://faculty.washington.edu/hansec/libs/SuperLU_5.2.0.tar.gz"
        self.libname = '-lsuperlu'
        self.libpath = 'libsuperlu.a'
        self.comp_wrapper = comp_wrapper

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        if self.comp_wrapper:
            self.skip = True
            print("SUPERLU provided by compiler wrappers: Skipping build")
            return self.config_dict
        self.setup_root_struct()
        self.config_dict['SUPERLU_LIB'] = self.config_dict[self.name + '_LIB']
        self.config_dict['SUPERLU_LIBS'] = self.libname
        self.config_dict["SUPERLU_VER_MAJOR"] = 5
        self.config_dict['COMP_DEFS'].append("-DHAVE_SUPERLU")
        self.config_dict['COMP_DEFS'].append("-DSUPERLU_VER_MAJOR=5")
        # Installation check files
        self.install_chk_files = [os.path.join(self.config_dict['SUPERLU_LIB'], self.libpath)]
        #
        self.config_dict['INC_DEFS'].append('$(SUPERLU_INC)')
        self.config_dict['EXT_LIBS'].append('$(SUPERLU_LIB)')
        template = '\n# SuperLU library\nSUPERLU_INC = -I{SUPERLU_INCLUDE}\nSUPERLU_LIB = -L{SUPERLU_LIB} {SUPERLU_LIBS}\n'
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        return self.config_dict

    def build(self):
        build_lines = [
            "rm -rf build",
            "mkdir build",
            "cd build",
            "export CC={CC}",
            "export CXX={CXX}",
            "export FC={FC}",
            "{CMAKE} -DCMAKE_INSTALL_PREFIX:PATH={SUPERLU_ROOT} \\",
            "  -DTPL_BLAS_LIBRARIES:PATH={BLAS_LIB_PATH} \\",
            "  -DCMAKE_POSITION_INDEPENDENT_CODE=TRUE \\",
            "  -Denable_blaslib:BOOL=TRUE ..",
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, self.config_dict)


class SUPERLU_DIST(package):
    def __init__(self, build_openmp, comp_wrapper=False):
        self.name = "SUPERLU_DIST"
        self.url = "https://portal.nersc.gov/project/sparse/superlu/superlu_dist_6.2.0.tar.gz"
        self.build_dir = "SuperLU_DIST_6.2.0"
        self.build_openmp = build_openmp
        self.comp_wrapper = comp_wrapper

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        if self.comp_wrapper:
            self.skip = True
            print("SUPERLU_DIST provided by compiler wrappers: Skipping build")
            return self.config_dict
        self.setup_root_struct()
        self.config_dict['SUPERLU_DIST_LIBS'] = "-lsuperlu_dist"
        self.config_dict['COMP_DEFS'].append("-DHAVE_SUPERLU_DIST")
        # Set makefile fields
        self.config_dict['INC_DEFS'].append('$(SUPERLU_DIST_INC)')
        self.config_dict['EXT_LIBS'].append('$(SUPERLU_DIST_LIB)')
        template = """
# SuperLU-Dist library
SUPERLU_DIST_INC = -I{SUPERLU_DIST_INCLUDE}
SUPERLU_DIST_LIB = -L{SUPERLU_DIST_LIB} {SUPERLU_DIST_LIBS}
"""
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        self.install_chk_files = [os.path.join(self.config_dict['SUPERLU_DIST_LIB'], "libsuperlu_dist.a")]
        return self.config_dict

    def build(self):
        #
        tmp_dict = self.config_dict.copy()
        if self.build_openmp:
            tmp_dict["SUPERLU_DIST_OMP_FLAG"] = "TRUE"
        else:
            tmp_dict["SUPERLU_DIST_OMP_FLAG"] = "FALSE"
        build_lines = [
            "rm -rf build_dir",
            "mkdir build_dir",
            "cd build_dir",
            "export CC={CC}",
            "export CXX={CXX}",
            "export FC={FC}",
            "{CMAKE} -DTPL_ENABLE_PARMETISLIB:BOOL=FALSE \\",
            "  -DCMAKE_INSTALL_PREFIX:PATH={SUPERLU_DIST_ROOT} \\",
            "  -DTPL_BLAS_LIBRARIES:PATH={BLAS_LIB_PATH} \\",
            "  -DTPL_LAPACK_LIBRARIES:PATH={LAPACK_LIB_PATH} \\",
            "  -Denable_openmp:BOOL={SUPERLU_DIST_OMP_FLAG} \\",
            "  -Denable_tests=OFF -Denable_examples=OFF \\",
            "  -DBUILD_SHARED_LIBS:BOOL=FALSE \\",
            "  -DMPI_Fortran_COMPILER={MPI_FC} \\",
            "  -DMPI_C_COMPILER={MPI_CC} \\",
            "  -DCMAKE_POSITION_INDEPENDENT_CODE=TRUE \\",
            "  -DCMAKE_INSTALL_LIBDIR=lib ..",
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, tmp_dict)


class UMFPACK(package):
    def __init__(self, comp_wrapper=False):
        self.name = "UMFPACK"
        self.url = "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.0.1.tar.gz"
        self.build_dir = "SuiteSparse-7.0.1"
        self.install_dir = "UMFPACK-7_0_1"
        self.comp_wrapper = comp_wrapper

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        if self.comp_wrapper:
            self.skip = True
            print("UMFPACK provided by compiler wrappers: Skipping build")
            return self.config_dict
        self.setup_root_struct()
        self.config_dict['UMFPACK_LIBS'] = "-lumfpack -lamd -lsuitesparseconfig"
        self.config_dict['COMP_DEFS'].append("-DHAVE_UMFPACK")
        # Installation check files
        self.install_chk_files = [os.path.join(self.config_dict['UMFPACK_LIB'], 'libumfpack.a')]
        #
        self.config_dict['INC_DEFS'].append('$(UMFPACK_INC)')
        self.config_dict['EXT_LIBS'].append('$(UMFPACK_LIB)')
        template = """
# UMFPACK library
UMFPACK_INC = -I{UMFPACK_INCLUDE}
UMFPACK_LIB = -L{UMFPACK_LIB} {UMFPACK_LIBS}
"""
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        return self.config_dict

    def build(self):
        if ver_lt(self.config_dict.get('CMAKE_VERSION','0.0'),"3.22"):
            error_exit('CMAKE >= 3.22 required for UMFPACK', ('Update or retry with "--build_cmake=1" to build a compatible version',))
        AMD_CMAKE_options = [
            "-DCMAKE_INSTALL_PREFIX={UMFPACK_ROOT}",
            "-DCMAKE_POSITION_INDEPENDENT_CODE=TRUE",
            "-DCMAKE_INSTALL_LIBDIR=lib"
        ]
        if self.config_dict.get('OpenBLAS_THREADS',False):
            AMD_CMAKE_options.append("-DCMAKE_EXE_LINKER_FLAGS={0}".format(self.config_dict['OMP_FLAGS']))
        config_CMAKE_options = AMD_CMAKE_options.copy() + [
            "-DBLAS_ROOT:PATH={BLAS_ROOT}",
            "-DBLA_VENDOR:STRING={BLAS_VENDOR}"
        ]
        UMFPACK_CMAKE_options = config_CMAKE_options.copy()
        if (config_dict['CC_VENDOR'] == 'gnu') and (self.config_dict['OS_TYPE'] == 'Darwin'):
            UMFPACK_CMAKE_options.append("-DCMAKE_SHARED_LINKER_FLAGS:STRING=-lgfortran")
        #
        build_lines = [
            "export CC={CC}",
            "export CX={CXX}",
            "export FC={FC}"
        ]
        if 'CMAKE_BIN' in self.config_dict:
            build_lines += ["export PATH={0}:$PATH".format(self.config_dict['CMAKE_BIN'])]
        build_lines += [
            "export FFLAGS=-fPIC",
            "export JOBS={MAKE_THREADS}",
            "export OPTIMIZATION=-O2",
            "export AUTOCC=no",
            "export UMFPACK_CONFIG=-DNCHOLMOD",
            "make distclean",
            "cd SuiteSparse_config",
            'make library CMAKE_OPTIONS="{0}"'.format(' '.join(config_CMAKE_options)),
            "make install",
            "cd ..",
            "cd AMD",
            'make library CMAKE_OPTIONS="{0}"'.format(' '.join(AMD_CMAKE_options)),
            "make install",
            "cd ..",
            "cd UMFPACK",
            'make library CMAKE_OPTIONS="{0}"'.format(' '.join(UMFPACK_CMAKE_options)),
            "make install"
        ]
        if self.config_dict['OS_TYPE'] == 'Darwin':
            build_lines += ["rm -f {UMFPACK_ROOT}/lib/lib*.dylib*"]
        else:
            build_lines += ["rm -f {UMFPACK_ROOT}/lib/lib*.so*"]
        self.run_build(build_lines, self.config_dict)


class FOX(package):
    def __init__(self):
        self.name = "FOX"
        self.url = "http://faculty.washington.edu/hansec/libs/FoX-4.1.2-dom.tar.gz"
        self.build_dir = "FoX-4.1.2"

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        self.setup_root_struct()
        install_path = os.path.join(self.root_path, self.install_dir)
        self.config_dict["FOX_INCLUDE"] = os.path.join(install_path, "finclude")
        self.config_dict["FOX_LIBS"] = "-lFoX_dom -lFoX_sax -lFoX_fsys -lFoX_utils -lFoX_common"
        self.config_dict['COMP_DEFS'].append("-DHAVE_XML")
        self.config_dict['CROSS_COMPILE_FLAG'] = ""
        if 'CROSS_COMPILE_HOST' in config_dict:
            self.config_dict['CROSS_COMPILE_FLAG'] = '--host="{0}"'.format(config_dict['CROSS_COMPILE_HOST'])
        # Installation check files
        self.install_chk_files = [os.path.join(self.config_dict['FOX_LIB'], 'libFoX_dom.a')]
        #
        self.config_dict['INC_DEFS'].append('$(FOX_INC)')
        self.config_dict['EXT_LIBS'].append('$(FOX_LIB)')
        template = '\n# FOX library\nFOX_INC = -I{FOX_INCLUDE}\nFOX_LIB = -L{FOX_LIB} {FOX_LIBS}\n'
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        return self.config_dict

    def build(self):
        build_lines = [
            "make distclean",
            "export CC={CC}",
            "export FC={FC}",
            "export CFLAGS=-fPIC",
            "export FCFLAGS=-fPIC"]
        if config_dict['OS_TYPE'] == 'Darwin': # Prevent configuration error with GCC
            build_lines.append("export GFORTRAN_UNBUFFERED_ALL=Y")
        build_lines += [
            "./configure --prefix={FOX_ROOT} --enable-dom {CROSS_COMPILE_FLAG}",
            "make -j{MAKE_THREADS}",
            "make install"
        ]
        self.run_build(build_lines, self.config_dict)


class ONURBS(package):
    def __init__(self):
        self.name = "ONURBS"
        self.url = "http://faculty.washington.edu/hansec/libs/opennurbs-5.0.tar.gz"
        self.install_dir = "opennurbs-5_0"

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        install_path = os.path.join(self.root_path, self.install_dir)
        self.config_dict['ONURBS_ROOT'] = install_path
        self.config_dict['ONURBS_LIB'] = os.path.join(install_path, "lib")
        self.config_dict['ONURBS_INCLUDE'] = os.path.join(install_path, "include")
        self.config_dict["ONURBS_LIBS"] = "-lopenNURBS"
        self.config_dict['COMP_DEFS'].append("-DHAVE_ONURBS")
        # Installation check files
        self.install_chk_files = [os.path.join(self.config_dict['ONURBS_LIB'], 'libopenNURBS.a')]
        #
        self.config_dict['INC_DEFS'].append('$(ONURBS_INC)')
        self.config_dict['EXT_LIBS'].append('$(ONURBS_LIB)')
        template = '\n# OpenNURBS library\nONURBS_INC = -I{ONURBS_INCLUDE}\nONURBS_LIB = -L{ONURBS_LIB} {ONURBS_LIBS}\n'
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        return self.config_dict

    def build(self):
        with open("makefile.in", "r") as fid:
            contents = fid.read()
        with open("makefile", "w+") as fid:
            fid.write(contents.format(**self.config_dict))
        build_lines = [
            "make -j{MAKE_THREADS} libopenNURBS.a",
            "mkdir {ONURBS_ROOT}",
            "mkdir {ONURBS_LIB}",
            "mkdir {ONURBS_INCLUDE}",
            "cp libopenNURBS.a {ONURBS_LIB}",
            "cp *.h {ONURBS_INCLUDE}",
            "cp -r zlib {ONURBS_INCLUDE}/"
        ]
        self.run_build(build_lines, self.config_dict)


class PETSC(package):
    def __init__(self, debug=False, with_superlu=True, with_umfpack=True, with_mumps=True, version=3.8,
                 comp_wrapper=False):
        self.name = "PETSC"
        self.version = version
        if self.version == '3.6':
            self.url = "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.6.4.tar.gz"
        elif self.version == '3.7':
            self.url = "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.7.7.tar.gz"
        elif self.version == '3.8':
            self.url = "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.8.4.tar.gz"
        elif self.version == '3.9':
            self.url = "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.9.4.tar.gz"
        elif self.version == '3.10':
            self.url = "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.10.5.tar.gz"
        elif self.version == '3.11':
            self.url = "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.11.4.tar.gz"
        else:
            error_exit('Invalid PETSc version requested')
        self.debug = debug
        self.with_superlu = with_superlu
        self.with_umfpack = with_umfpack
        self.with_mumps = with_mumps
        self.comp_wrapper = comp_wrapper

    def detect_version(self):
        print("  Testing PETSc version")
        with open('tmp.c', 'w+') as fid:
            fid.write('#include "petscversion.h"\n\nPETSC_VERSION_MAJOR PETSC_VERSION_MINOR')
        command = "{CC} -E tmp.c | tail -1".format(**self.config_dict)
        result, errcode = run_command(command)
        os.remove('tmp.c')
        if errcode == 0:
            values = result.split()
            ver_major = int(values[0])
            ver_minor = int(values[1])
            print("    Found PETSc version {0}.{1} ".format(ver_major, ver_minor))
            return ver_major, ver_minor
        else:
            return None, None

    def setup(self, config_dict):
        self.config_dict = config_dict.copy()
        if self.comp_wrapper:
            self.skip = True
            # Detect PETSc version
            ver_major, ver_minor = self.detect_version()
            if ver_major is None:
                print("    PETSc not found!")
                sys.exit(1)
            self.version = '{0}.{1}'.format(ver_major, ver_minor)
            self.config_dict['COMP_DEFS'].append('-DHAVE_PETSC -DPETSC_VERSION_MAJOR={0} -DPETSC_VERSION_MINOR={1}'.format(ver_major, ver_minor))
            try:
                metis_tmp = METIS()
                metis_tmp.detect_sizes()
            except:
                print("    Cannot detect PETSc METIS!")
                sys.exit(1)
            self.config_dict['COMP_DEFS'].append("-DHAVE_MUMPS")
            self.config_dict['COMP_DEFS'].append("-DHAVE_SUPERLU")
            return self.config_dict
        self.setup_root_struct()
        # Setup libraries
        self.config_dict['METIS_ROOT'] = self.config_dict["PETSC_ROOT"]
        self.config_dict['PETSC_LIB'] = self.config_dict[self.name + '_LIB']
        self.config_dict['PETSC_LIBS'] = "-lpetsc"
        if self.with_umfpack:
            self.config_dict['PETSC_LIBS'] += " -lumfpack -lamd"
            self.config_dict['COMP_DEFS'].append("-DHAVE_UMFPACK")
            self.install_chk_files.append(os.path.join(self.config_dict['PETSC_LIB'], 'libumfpack.a'))
        if self.with_superlu:
            if self.version == '3.6':
                self.config_dict['PETSC_LIBS'] += " -lsuperlu_4.3 -lsuperlu_dist_4.1"
                self.install_chk_files.append(os.path.join(self.config_dict['PETSC_LIB'], 'libsuperlu_4.3.a'))
            elif not ver_lt(self.version,'3.7'):
                self.config_dict['PETSC_LIBS'] += " -lsuperlu"
                if ver_lt(self.version,'3.10'):
                    self.config_dict['PETSC_LIBS'] += " -lsuperlu_dist"
                self.config_dict['COMP_DEFS'].append("-DHAVE_SUPERLU")
                self.config_dict['COMP_DEFS'].append("-DSUPERLU_VER_MAJOR=5")
                self.install_chk_files.append(os.path.join(self.config_dict['PETSC_LIB'], 'libsuperlu.a'))
            else:
                self.config_dict['PETSC_LIBS'] += " -lsuperlu_4.3 -lsuperlu_dist_3.3"
                self.install_chk_files.append(os.path.join(self.config_dict['PETSC_LIB'], 'libsuperlu_4.3.a'))
        if self.with_mumps:
            self.config_dict['PETSC_LIBS'] += " -ldmumps -lmumps_common -lpord -lscalapack"
            self.config_dict['COMP_DEFS'].append("-DHAVE_MUMPS")
            self.install_chk_files.append(os.path.join(self.config_dict['PETSC_LIB'], 'libdmumps.a'))
        self.config_dict['PETSC_LIBS'] += " -lparmetis -lmetis"
        # Installation check files
        self.install_chk_files = [os.path.join(self.config_dict['PETSC_LIB'], 'libpetsc.a')]
        #
        self.config_dict['INC_DEFS'].append('$(PETSC_INC)')
        self.config_dict['EXT_LIBS'].append('$(PETSC_LIB)')
        self.config_dict['PETSC_COMP_DEFS'] = ''
        ver_string = self.version.split('.')
        ver_major = ver_string[0]
        ver_minor = ver_string[1]
        self.config_dict['COMP_DEFS'].append('-DHAVE_PETSC -DPETSC_VERSION_MAJOR={0} -DPETSC_VERSION_MINOR={1}'.format(ver_major, ver_minor))
        template = """\n# PETSc library
PETSC_INC = -I{PETSC_INCLUDE}
PETSC_LIB = -L{PETSC_LIB} {PETSC_LIBS}\n"""
        self.config_dict["MAKE_INC_CONTENTS"] += template.format(**self.config_dict)
        return self.config_dict

    def build(self):
        #
        def_lines = []
        options = ['--with-cxx=0']
        if 'MPI_CC' in self.config_dict:
            options += ['--CC={MPI_CC}', '--FC={MPI_FC}']
        else:
            options += ['--with-mpi-dir={MPI_ROOT}']
        if config_dict['CC_VENDOR'] == 'gnu' and int(config_dict['CC_VERSION'].split(".")[0]) > 9:
            options.append('--FFLAGS="-fallow-argument-mismatch"')
        options += ['--download-metis', '--download-parmetis', '--with-x=no', '--with-shared-libraries=0']
        if not ver_lt(self.version,'3.5'):
            options += [' --with-ssl=0']
            if not ver_lt(self.version,'3.8'):
                options += ['--with-cmake-exec={CMAKE}']
            else:
                options += ['--with-cmake={CMAKE}']
        if self.with_superlu:
            # Fix SDK issue on MacOS "Catalina" (10.15)
            if (self.config_dict['OS_TYPE'] == 'Darwin') and (not ver_lt(self.config_dict['OS_VER'], '10.15')):
                result, _ = run_command("{0} -print-sysroot".format(self.config_dict['CC']))
                def_lines.append('export SDKROOT={0}'.format(result.strip()))
            options += ['--download-superlu']
            if ver_lt(self.version,'3.10'):
                options += ['--download-superlu_dist']
        if self.with_umfpack:
            options += ['--download-umfpack']
        if self.with_mumps:
            options += ['--download-mumps', '--download-scalapack', '--download-blacs']
        #
        if 'MKL_ROOT' in self.config_dict:
            options += ['--with-blas-lapack-dir={MKL_ROOT}']
        else:
            options += ['--with-blas-lib={BLAS_LIB_PATH}', '--with-lapack-lib={LAPACK_LIB_PATH}']
        #
        if self.debug:
            options += ['--with-debugging=yes']
        else:
            options += ['--with-debugging=no', '--COPTFLAGS=-O2', '--FOPTFLAGS=-O2']
        build_lines = def_lines + [
            "./configure --prefix={PETSC_ROOT} " + " ".join(options),
            "make MAKE_NP={MAKE_THREADS} all",
            "make install"
        ]
        run_command("make distclean", timeout=30)
        self.run_build(build_lines, self.config_dict)


# Start of main script
parser = argparse.ArgumentParser()
parser.description = "Third-party library build script for the Open FUSION Toolkit"
parser.add_argument("--download_only", action="store_true", default=False, help="Only download packages")
parser.add_argument("--setup_only", action="store_true", default=False, help="Download and setup build, but do not actually build")
parser.add_argument("--nthread", default=1, type=int, help="Number of threads to use for make (default=1)")
parser.add_argument("--opt_flags", default=None, type=str, help="Compiler optimization flags")
parser.add_argument("--ld_flags", default=None, type=str, help="Linker flags")
parser.add_argument("--cross_compile_host", default=None, type=str, help="Host type for cross-compilation")
parser.add_argument("--no_dl_progress", action="store_false", default=True, help="Do not report progress during file download")
#
group = parser.add_argument_group("CMAKE", "CMAKE configure options for the Open FUSION Toolkit")
group.add_argument("--build_cmake", default=0, type=int, choices=(0,1), help="Build CMAKE instead of using system version?")
group.add_argument("--oft_build_debug", default=0, type=int, choices=(0,1), help="Build debug version of OFT?")
group.add_argument("--oft_build_python", default=1, type=int, choices=(0,1), help="Build OFT Python libraries? (default: 1)")
group.add_argument("--oft_use_openmp", default=1, type=int, choices=(0,1), help="Build OFT with OpenMP support? (default)")
group.add_argument("--oft_build_tests", default=0, type=int, choices=(0,1), help="Build OFT tests?")
group.add_argument("--oft_build_examples", default=0, type=int, choices=(0,1), help="Build OFT examples?")
group.add_argument("--oft_build_docs", default=0, type=int, choices=(0,1), help="Build OFT documentation? (requires doxygen)")
group.add_argument("--oft_package", action="store_true", default=False, help="Perform a packaging build of OFT?")
group.add_argument("--oft_package_release", action="store_true", default=False, help="Perform a release package of OFT?")
group.add_argument("--oft_build_coverage", action="store_true", default=False, help="Build OFT with code coverage flags?")
#
group = parser.add_argument_group("MPI", "MPI package options")
group.add_argument("--build_mpi", default=0, type=int, choices=(0,1), help="Build MPI libraries?")
group.add_argument("--mpi_cc", default=None, type=str, help="MPI C compiler wrapper")
group.add_argument("--mpi_fc", default=None, type=str, help="MPI FORTRAN compiler wrapper")
group.add_argument("--mpi_lib_dir", default=None, type=str, help="MPI library directory")
group.add_argument("--mpi_libs", default=None, type=str, help="MPI libraries")
group.add_argument("--mpi_include_dir", default=None, type=str, help="MPI include directory")
#
group = parser.add_argument_group("HDF5", "HDF5 package options")
group.add_argument("--hdf5_cc", default=None, type=str, help="HDF5 C compiler wrapper")
group.add_argument("--hdf5_fc", default=None, type=str, help="HDF5 FORTRAN compiler wrapper")
group.add_argument("--hdf5_parallel", action="store_true", default=False, help="Use parallel HDF5 interface?")
#
group = parser.add_argument_group("BLAS/LAPACK", "BLAS/LAPACK package options")
group.add_argument("--oblas_threads", action="store_true", default=False, help="Build OpenBLAS with thread support (OpenMP)")
group.add_argument("--oblas_dynamic_arch", action="store_true", default=False, help="Build OpenBLAS with multiple architecure support")
group.add_argument("--ref_blas", action="store_true", default=False, help="Use reference BLAS/LAPACK instead of OpenBLAS")
group.add_argument("--blas_lapack_wrapper", action="store_true", default=False, help="BLAS/LAPACK included in compilers")
group.add_argument("--use_mkl", action="store_true", default=False, help="Use MKL BLAS/LAPACK instead of OpenBLAS")
group.add_argument("--mkl_root", default=None, type=str, help="MKL root directory path (required for MKL)")
group.add_argument("--blas_lib_path", default=None, type=str, help="Path to pre-built BLAS library")
group.add_argument("--lapack_lib_path", default=None, type=str, help="Path to pre-built LAPACK library")
#
group = parser.add_argument_group("METIS", "METIS package options")
group.add_argument("--metis_wrapper", action="store_true", default=False, help="METIS included in compilers")
#
group = parser.add_argument_group("FoX XML", "FoX XML package options")
group.add_argument("--build_fox", default=1, type=int, choices=(0,1), help="Build Fox XML library? (default)")
#
group = parser.add_argument_group("OpenNURBS", "OpenNURBS package options")
group.add_argument("--build_onurbs", default=0, type=int, choices=(0,1), help="Build OpenNURBS library?")
#
group = parser.add_argument_group("NETCDF", "NETCDF package options")
group.add_argument("--build_netcdf", default=0, type=int, choices=(0,1), help="Build NETCDF library?")
group.add_argument("--netcdf_wrapper", action="store_true", default=False, help="NETCDF included in compilers")
#
group = parser.add_argument_group("ARPACK", "ARPACK package options")
group.add_argument("--build_arpack", default=0, type=int, choices=(0,1), help="Build ARPACK library?")
#
group = parser.add_argument_group("SuperLU", "SuperLU package options")
group.add_argument("--build_superlu", default=0, type=int, choices=(0,1), help="Build SuperLU library?")
group.add_argument("--superlu_wrapper", action="store_true", default=False, help="SuperLU included in compilers")
#
group = parser.add_argument_group("SuperLU-DIST", "SuperLU-DIST package options")
group.add_argument("--build_superlu_dist", default=0, type=int, choices=(0,1), help="Build SuperLU-DIST library?")
group.add_argument("--superlu_dist_threads", action="store_true", default=False, help="Build SuperLU-DIST with thread support (OpenMP)")
group.add_argument("--superlu_dist_wrapper", action="store_true", default=False, help="SuperLU-DIST included in compilers")
#
group = parser.add_argument_group("UMFPACK", "UMFPACK package options")
group.add_argument("--build_umfpack", default=0, type=int, choices=(0,1), help="Build UMFPACK library?")
group.add_argument("--umfpack_wrapper", action="store_true", default=False, help="UMFPACK included in compilers")
#
group = parser.add_argument_group("PETSc", "PETSc package options")
group.add_argument("--build_petsc", default=0, type=int, choices=(0,1), help="Build PETSc library?")
group.add_argument("--petsc_debug", default=1, type=int, choices=(0,1), help="Build PETSc with debugging information (default)")
group.add_argument("--petsc_superlu", default=1, type=int, choices=(0,1), help="Build PETSc with SuperLU (default)")
group.add_argument("--petsc_mumps", default=1, type=int, choices=(0,1), help="Build PETSc with MUMPS (default)")
group.add_argument("--petsc_umfpack", default=0, type=int, choices=(0,1), help="Build PETSc with UMFPACK")
group.add_argument("--petsc_version", default="3.8", type=str,
    help="Use different version of PETSc [3.6,3.7,3.8,3.9,3.10] (default is 3.8)")
group.add_argument("--petsc_wrapper", action="store_true", default=False, help="PETSc included in compilers")
#
options = parser.parse_args()
fetch_progress = options.no_dl_progress
build_cmake_ver = None
if options.build_cmake == 1:
    build_cmake_ver = CMAKE().version
config_dict = setup_build_env(build_cmake_ver=build_cmake_ver)
config_dict['DOWN_ONLY'] = options.download_only
config_dict['SETUP_ONLY'] = options.setup_only
if options.nthread > 1:
    config_dict['MAKE_THREADS'] = options.nthread
if options.opt_flags is not None:
    config_dict['OPT_FLAGS'] = options.opt_flags
if options.ld_flags is not None:
    config_dict['LD_FLAGS'] = options.ld_flags
if options.cross_compile_host is not None:
    config_dict['CROSS_COMPILE_HOST'] = options.cross_compile_host
# Building with MPI?
use_mpi = False
if (options.mpi_cc is not None) and (options.mpi_fc is not None):
    config_dict['MPI_CC'] = options.mpi_cc
    config_dict['MPI_FC'] = options.mpi_fc
    use_mpi = True
else:
    if options.mpi_lib_dir is not None:
        use_mpi = True
        config_dict['MPI_LIB'] = options.mpi_lib_dir
        if options.mpi_libs is not None:
            config_dict['MPI_LIBS'] = options.mpi_libs
        else:
            parser.exit(-1, '"MPI_LIBS" required when "MPI_LIB" is specified')
        if options.mpi_include_dir is not None:
            config_dict['MPI_INCLUDE'] = options.mpi_include_dir
        else:
            parser.exit(-1, '"MPI_INCLUDE" required when "MPI_LIB" is specified')
    elif options.build_mpi:
        use_mpi = True
# Setup library builds (in order of dependency)
packages = []
# CMAKE
if options.build_cmake == 1:
    packages.append(CMAKE())
# BLAS/LAPACK
if options.use_mkl:
    packages.append(MKL(options.mkl_root))
else:
    if (options.blas_lib_path is not None) or (options.lapack_lib_path is not None):
        if (options.blas_lib_path is None) or (options.lapack_lib_path is None):
            parser.exit(-1, 'Both "blas_lib_path" and "lapack_lib_path" must be specified to use pre-built BLAS/LAPACK')
        packages.append(BLAS_LAPACK(blas_lib_path=options.blas_lib_path, lapack_lib_path=options.lapack_lib_path))
    else:
        if options.ref_blas or options.blas_lapack_wrapper:
            packages.append(BLAS_LAPACK(options.blas_lapack_wrapper))
        else:
            packages.append(OpenBLAS(options.oblas_threads,options.oblas_dynamic_arch))
# MPI
if use_mpi:
    packages.append(MPI())
else:
    if options.hdf5_parallel:
        print('Warning: Reverting to serial HDF5 library without MPI')
    if (options.build_petsc == 1) or options.petsc_wrapper:
        parser.exit(-1, 'PETSc requires MPI')
# HDF5
if (options.hdf5_cc is not None) and (options.hdf5_fc is not None):
    config_dict['HDF5_CC'] = options.hdf5_cc
    config_dict['HDF5_FC'] = options.hdf5_fc
    packages.append(HDF5(parallel=(options.hdf5_parallel and use_mpi)))
else:
    packages.append(HDF5(parallel=(options.hdf5_parallel and use_mpi)))
# Are we building OpenNURBS?
if options.build_onurbs == 1:
    packages.append(ONURBS())
# Are we building FoX?
if options.build_fox == 1:
    packages.append(FOX())
# Are we building ARPACK?
if options.build_arpack == 1:
    packages.append(ARPACK(parallel=use_mpi, link_omp=options.oblas_threads))
# Are we building NETCDF?
if (options.build_netcdf == 1) or options.netcdf_wrapper:
    packages.append(NETCDF(options.netcdf_wrapper))
# Are we building PETSc?
if (options.build_petsc == 1) or options.petsc_wrapper:
    packages.append(PETSC(debug=options.petsc_debug, with_superlu=options.petsc_superlu,
                          with_umfpack=options.petsc_umfpack, with_mumps=options.petsc_mumps,
                          version=options.petsc_version, comp_wrapper=options.petsc_wrapper))
else:
    packages.append(METIS(options.metis_wrapper))
    if (options.build_superlu == 1) or options.superlu_wrapper:
        packages.append(SUPERLU(options.superlu_wrapper))
    if (options.build_superlu_dist == 1) or options.superlu_dist_wrapper:
        if not use_mpi:
            parser.exit(-1, 'SuperLU-DIST requires MPI')
        packages.append(SUPERLU_DIST(options.superlu_dist_threads, options.superlu_dist_wrapper))
    if (options.build_umfpack == 1) or options.umfpack_wrapper:
        packages.append(UMFPACK(options.umfpack_wrapper))
#
for package in packages:
    config_dict = package.install(config_dict)
#
# print(config_dict)
if not (config_dict['DOWN_ONLY'] or config_dict['SETUP_ONLY']):
    build_make_include(config_dict)
    build_cmake_script(config_dict,
        build_debug=(options.oft_build_debug == 1),
        use_openmp=(options.oft_use_openmp == 1),
        build_python=(options.oft_build_python == 1),
        build_tests=(options.oft_build_tests == 1),
        build_examples=(options.oft_build_examples == 1),
        build_docs=(options.oft_build_docs == 1),
        build_coverage=(options.oft_build_coverage == 1),
        package_build=options.oft_package,
        package_release=options.oft_package_release
    )
