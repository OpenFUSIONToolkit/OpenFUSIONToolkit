#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
"""! Build a BLAS/LAPACK shim library for the macOS Accelerate framework

Accelerate exports two interfaces: a legacy one (eg. "dgemm_"), which is LAPACK 3.2.1 and uses
the f2c return convention for single precision and complex functions, and a modern one
(LAPACK 3.11, macOS >= 13.3) with symbols suffixed by "$NEWLAPACK". Fortran codes cannot reference
the modern symbols directly, so this script generates and compiles a small shim library that maps each
standard Fortran symbol to its "$NEWLAPACK" counterpart. Most symbols are pure tail-call trampolines;
complex-valued functions, which return through a hidden argument in Accelerate, are wrapped in C.
Utility routines that Accelerate only provides through the legacy interface ("lsame", "dcabs1") are
implemented directly in the shim.

The symbol list is read from the text-based stub files (.tbd) of the macOS SDK, so the shim matches
the SDK used for the build. The resulting library is linked against the Accelerate framework and is
used in place of a BLAS/LAPACK library (eg. "-L<dir> -loftaccelerate").

Usage: python accelerate_shim.py --output <dir>/liboftaccelerate.dylib [--cc CC] [--arch ARCH] [--sdk_path PATH]

This script is called by build_libs.py ("--use_accelerate"), but can also be run directly.
"""
import argparse
import os
import platform
import re
import shlex
import subprocess
import sys

## Complex valued functions (returned through hidden first argument in Accelerate)
complex_funcs = {
    'cdotc': ('float complex', 'dot'),
    'cdotu': ('float complex', 'dot'),
    'zdotc': ('double complex', 'dot'),
    'zdotu': ('double complex', 'dot'),
    'cladiv': ('float complex', 'div'),
    'zladiv': ('double complex', 'div')
}
## Reference implementations of routines without a "$NEWLAPACK" version
legacy_funcs = {
    'lsame': [
        '#include <ctype.h>',
        '#include <stddef.h>',
        'int lsame_(const char *ca, const char *cb, size_t ca_len, size_t cb_len)',
        '{ return toupper((unsigned char)*ca) == toupper((unsigned char)*cb); }'
    ],
    'dcabs1': [
        '#include <math.h>',
        'double dcabs1_(const double complex *z) { return fabs(creal(*z)) + fabs(cimag(*z)); }'
    ]
}


def shim_error(message):
    '''! Print error message and exit

    @param message Error message
    '''
    print('accelerate_shim: ERROR: {0}'.format(message), file=sys.stderr)
    sys.exit(1)


def run(cmd):
    '''! Run a command, printing it first, and exit on failure

    @param cmd Command as a list of arguments
    '''
    print(' '.join(shlex.quote(arg) for arg in cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        shim_error('Command failed with exit code {0}'.format(result.returncode))


def get_sdk_path(sdk_path=None):
    '''! Get macOS SDK path

    @param sdk_path SDK path (if None, uses "SDKROOT" or "xcrun --show-sdk-path")
    @result Path to macOS SDK
    '''
    if sdk_path is None:
        sdk_path = os.environ.get('SDKROOT', None)
    if sdk_path is None:
        result = subprocess.run(['xcrun', '--show-sdk-path'], stdout=subprocess.PIPE, universal_newlines=True)
        if result.returncode != 0:
            shim_error('Could not determine macOS SDK path (specify with "--sdk_path")')
        sdk_path = result.stdout.strip()
    return sdk_path


def get_symbols(sdk_path):
    '''! Read modern (non-ILP64) Fortran BLAS/LAPACK symbols from Accelerate's stub files

    @param sdk_path Path to macOS SDK
    @result Sorted list of symbol names (without leading/trailing underscores or suffix)
    '''
    veclib_path = os.path.join(sdk_path, 'System', 'Library', 'Frameworks', 'Accelerate.framework',
                               'Versions', 'A', 'Frameworks', 'vecLib.framework', 'Versions', 'A')
    symbols = set()
    for lib_name in ('libBLAS.tbd', 'libLAPACK.tbd'):
        tbd_path = os.path.join(veclib_path, lib_name)
        if not os.path.isfile(tbd_path):
            shim_error('Could not find Accelerate stub file "{0}"'.format(tbd_path))
        with open(tbd_path, 'r') as fid:
            symbols.update(re.findall(r"_([a-z][a-z0-9_]*)\$NEWLAPACK(?!\$)", fid.read()))
    # Remove C-only interfaces
    symbols = sorted([sym for sym in symbols if not sym.startswith(('cblas_', 'catlas_', 'appleblas_'))])
    if ('dgemm' not in symbols) or ('dgetrf' not in symbols):
        shim_error('Accelerate "$NEWLAPACK" interface not found in SDK "{0}" (macOS SDK >= 13.3 required)'.format(sdk_path))
    return symbols


def write_sources(symbols, arch, build_dir):
    '''! Write shim sources (assembly trampolines and C wrappers)

    @param symbols List of symbols to wrap (from @ref get_symbols)
    @param arch Target architecture ("arm64" or "x86_64")
    @param build_dir Directory for generated sources
    @result Paths to assembly and C sources
    '''
    # Trampolines for all symbols except complex-valued functions
    if arch == 'arm64':
        branch_inst = 'b'
    elif arch == 'x86_64':
        branch_inst = 'jmp'
    else:
        shim_error('Unsupported architecture "{0}" (must be "arm64" or "x86_64")'.format(arch))
    asm_lines = ['.text']
    for sym in symbols:
        if sym in complex_funcs:
            continue
        asm_lines += [
            '.globl _{0}_'.format(sym),
            '.p2align 4',
            '_{0}_:'.format(sym),
            '  {0} "_{1}$NEWLAPACK"'.format(branch_inst, sym)
        ]
    asm_file = os.path.join(build_dir, 'oft_accelerate_tramp.S')
    with open(asm_file, 'w+') as fid:
        fid.write('\n'.join(asm_lines) + '\n')
    # C wrappers for complex-valued functions
    c_lines = ['#include <complex.h>']
    for sym, (ctype, func_type) in complex_funcs.items():
        if sym not in symbols:
            continue
        if func_type == 'dot':
            c_lines += [
                'extern void {0}_new({1} *, const int *, const {1} *, const int *, const {1} *, const int *) __asm__("_{0}$NEWLAPACK");'.format(sym, ctype),
                '{1} {0}_(const int *n, const {1} *x, const int *incx, const {1} *y, const int *incy)'.format(sym, ctype),
                '{{ {1} r; {0}_new(&r, n, x, incx, y, incy); return r; }}'.format(sym, ctype)
            ]
        else:
            c_lines += [
                'extern void {0}_new({1} *, const {1} *, const {1} *) __asm__("_{0}$NEWLAPACK");'.format(sym, ctype),
                '{1} {0}_(const {1} *x, const {1} *y)'.format(sym, ctype),
                '{{ {1} r; {0}_new(&r, x, y); return r; }}'.format(sym, ctype)
            ]
    for sym, func_lines in legacy_funcs.items():
        if sym not in symbols:
            c_lines += func_lines
    c_file = os.path.join(build_dir, 'oft_accelerate_cplx.c')
    with open(c_file, 'w+') as fid:
        fid.write('\n'.join(c_lines) + '\n')
    return asm_file, c_file


def build_shim(output, cc='cc', arch=None, sdk_path=None, build_dir=None):
    '''! Generate and compile the shim library

    @param output Path of shim library to create (also used as its install name)
    @param cc C compiler (may include arguments)
    @param arch Target architecture (default: host architecture)
    @param sdk_path Path to macOS SDK (see @ref get_sdk_path)
    @param build_dir Directory for generated sources and objects (default: current directory)
    '''
    if platform.system() != 'Darwin':
        shim_error('Accelerate is only available on macOS')
    if arch is None:
        arch = platform.machine()
    if build_dir is None:
        build_dir = os.getcwd()
    output = os.path.abspath(output)
    symbols = get_symbols(get_sdk_path(sdk_path))
    print('Found {0} Accelerate BLAS/LAPACK symbols'.format(len(symbols)))
    asm_file, c_file = write_sources(symbols, arch, build_dir)
    asm_obj = os.path.splitext(asm_file)[0] + '.o'
    c_obj = os.path.splitext(c_file)[0] + '.o'
    cc_cmd = shlex.split(cc)
    os.makedirs(os.path.dirname(output), exist_ok=True)
    run(cc_cmd + ['-c', asm_file, '-o', asm_obj])
    run(cc_cmd + ['-O2', '-fPIC', '-c', c_file, '-o', c_obj])
    run(cc_cmd + ['-dynamiclib', asm_obj, c_obj, '-o', output, '-install_name', output, '-framework', 'Accelerate'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build a Fortran BLAS/LAPACK shim library for the macOS Accelerate framework')
    parser.add_argument('--output', required=True, type=str, help='Path of shim library to create (eg. lib/liboftaccelerate.dylib)')
    parser.add_argument('--cc', default=os.environ.get('CC', 'cc'), type=str, help='C compiler (default: $CC or "cc")')
    parser.add_argument('--arch', default=None, type=str, help='Target architecture, "arm64" or "x86_64" (default: host)')
    parser.add_argument('--sdk_path', default=None, type=str, help='macOS SDK path (default: $SDKROOT or "xcrun --show-sdk-path")')
    parser.add_argument('--build_dir', default=None, type=str, help='Directory for generated sources and objects (default: current directory)')
    options = parser.parse_args()
    build_shim(options.output, cc=options.cc, arch=options.arch, sdk_path=options.sdk_path, build_dir=options.build_dir)
