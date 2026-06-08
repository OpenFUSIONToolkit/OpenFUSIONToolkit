/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://www.hdfgroup.org/licenses.               *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "H5private.h"

H5_GCC_DIAG_OFF("larger-than=")
H5_CLANG_DIAG_OFF("overlength-strings")

/* clang-format off */
const char H5build_settings[]=
    "        SUMMARY OF THE HDF5 CONFIGURATION\n"
    "        =================================\n"
    "\n"
    "General Information:\n"
    "-------------------\n"
    "                   HDF5 Version: 1.14.6\n"
    "                  Configured on: Fri Feb 13 16:05:50 EST 2026\n"
    "                  Configured by: taliaangles@Talias-MacBook-Air-1549.local\n"
    "                    Host system: aarch64-apple-darwin24.6.0\n"
    "              Uname information: Darwin Talias-MacBook-Air-1549.local 24.6.0 Darwin Kernel Version 24.6.0: Mon Jul 14 11:30:40 PDT 2025; root:xnu-11417.140.69~1/RELEASE_ARM64_T8132 arm64\n"
    "                       Byte sex: little-endian\n"
    "             Installation point: /Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/hdf5-1_14_6\n"
    "\n"
    "Compiling Options:\n"
    "------------------\n"
    "                     Build Mode: production\n"
    "              Debugging Symbols: no\n"
    "                        Asserts: no\n"
    "                      Profiling: no\n"
    "             Optimization Level: high\n"
    "\n"
    "Linking Options:\n"
    "----------------\n"
    "                      Libraries: shared\n"
    "  Statically Linked Executables: \n"
    "                        LDFLAGS: \n"
    "                     H5_LDFLAGS: \n"
    "                     AM_LDFLAGS: \n"
    "                Extra libraries: -lz -ldl -lm \n"
    "                       Archiver: ar\n"
    "                       AR_FLAGS: cr\n"
    "                         Ranlib: ranlib\n"
    "\n"
    "Languages:\n"
    "----------\n"
    "                              C: yes\n"
    "                     C Compiler: /opt/homebrew/bin/gcc-15 ( gcc-15 (Homebrew GCC 15.2.0_1) 15.2.0)\n"
    "                       CPPFLAGS: \n"
    "                    H5_CPPFLAGS:   -DNDEBUG -UH5_DEBUG_API -I/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/build/hdf5-1.14.6/src/H5FDsubfiling\n"
    "                    AM_CPPFLAGS: \n"
    "                        C Flags: \n"
    "                     H5 C Flags:  -std=c99  -Wall -Wcast-qual -Wconversion -Wextra -Wfloat-equal -Wformat=2 -Winit-self -Winvalid-pch -Wmissing-include-dirs -Wshadow -Wundef -Wwrite-strings -pedantic -Wno-c++-compat -Wlarger-than=2560 -Wlogical-op -Wframe-larger-than=16384 -Wpacked-bitfield-compat -Wsync-nand -Wno-unsuffixed-float-constants -Wdouble-promotion -Wtrampolines -Wstack-usage=8192 -Wmaybe-uninitialized -Wdate-time -Warray-bounds=2 -Wc99-c11-compat -Wduplicated-cond -Whsa -Wnormalized -Wnull-dereference -Wunused-const-variable -Walloca -Walloc-zero -Wduplicated-branches -Wformat-overflow=2 -Wformat-truncation=1 -Wattribute-alias -Wshift-overflow=2 -Wattribute-alias=2 -Wmissing-profile -Wc11-c2x-compat -fstdarg-opt -fdiagnostics-urls=never -fno-diagnostics-color -s  -Wbad-function-cast -Wcast-align -Wformat -Wimplicit-function-declaration -Wint-to-pointer-cast -Wmissing-declarations -Wmissing-prototypes -Wnested-externs -Wold-style-definition -Wpacked -Wpointer-sign -Wpointer-to-int-cast -Wredundant-decls -Wstrict-prototypes -Wswitch -Wunused-but-set-variable -Wunused-variable -Wunused-function -Wunused-parameter -Wincompatible-pointer-types -Wint-conversion -Wshadow -Wrestrict -Wcast-function-type -Wmaybe-uninitialized -Wcast-align=strict -Wno-aggregate-return -Wno-inline -Wno-missing-format-attribute -Wno-missing-noreturn -Wno-overlength-strings -Wno-jump-misses-init -Wstrict-overflow=2 -Wno-suggest-attribute=const -Wno-suggest-attribute=noreturn -Wno-suggest-attribute=pure -Wno-suggest-attribute=format -Wno-suggest-attribute=cold -Wno-suggest-attribute=malloc -O3\n"
    "                     AM C Flags: \n"
    "               Shared C Library: yes\n"
    "               Static C Library: no\n"
    "\n"
    "\n"
    "                        Fortran: yes\n"
    "               Fortran Compiler: /opt/homebrew/bin/gfortran-15 ( GNU Fortran (Homebrew GCC 15.2.0_1) 15.2.0)\n"
    "                  Fortran Flags: \n"
    "               H5 Fortran Flags:  -std=f2008  -Waliasing -Wall -Wcharacter-truncation -Wextra -Wimplicit-interface -Wsurprising -Wunderflow -pedantic -Wintrinsics-std -Wimplicit-procedure -Wreal-q-constant -Wfunction-elimination -Wrealloc-lhs -Wrealloc-lhs-all -Wno-c-binding-type -Winteger-division -Wfrontend-loop-interchange  -fdiagnostics-urls=never -fno-diagnostics-color -s  -Wno-unused-dummy-argument -Wno-array-temporaries -O3\n"
    "               AM Fortran Flags: \n"
    "         Shared Fortran Library: yes\n"
    "         Static Fortran Library: no\n"
    "               Module Directory: ${includedir}\n"
    "\n"
    "                            C++: no\n"
    "\n"
    "                           Java: no\n"
    "\n"
    "\n"
    "Features:\n"
    "---------\n"
    "                     Parallel HDF5: no\n"
    "  Parallel Filtered Dataset Writes: no\n"
    "                Large Parallel I/O: no\n"
    "                High-level library: no\n"
    "Dimension scales w/ new references: no\n"
    "                  Build HDF5 Tests: no\n"
    "                  Build HDF5 Tools: yes\n"
    "                   Build GIF Tools: no\n"
    "                      Threadsafety: no\n"
    "               Default API mapping: v114\n"
    "    With deprecated public symbols: yes\n"
    "            I/O filters (external): deflate(zlib)\n"
    "                  _Float16 support: yes\n"
    "                     Map (H5M) API: no\n"
    "                        Direct VFD: no\n"
    "                        Mirror VFD: no\n"
    "                     Subfiling VFD: no\n"
    "                (Read-Only) S3 VFD: no\n"
    "              (Read-Only) HDFS VFD: no\n"
    "    Packages w/ extra debug output: none\n"
    "                       API tracing: no\n"
    "              Using memory checker: no\n"
    "                  Use file locking: best-effort\n"
    "         Strict file format checks: no\n"
    "      Optimization instrumentation: no\n"
;
/* clang-format on */

H5_GCC_DIAG_ON("larger-than=")
H5_CLANG_DIAG_OFF("overlength-strings")
