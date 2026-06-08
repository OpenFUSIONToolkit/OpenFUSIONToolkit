.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90 .m4 .exe .xml

BUILD_TARGETS=    dom_lib

VPATH=
compile_prefix=/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/fox-4.1.2/objs
install_prefix=/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/fox-4_1_2
LIB_DIR=$(compile_prefix)/lib
MOD_DIR=$(compile_prefix)/finclude

FPP=
FC=gfortran-15
RANLIB=ranlib

FFLAGS=-fPIC 
FPPFLAGS= -DFC_HAVE_FLUSH -DFC_HAVE_ABORT -DFC_EOR_LF
LDFLAGS=

FCFLAGS_free_f90=
FPPFLAGS_free_F90=

INC_PREFIX=-I
MOD_PREFIX=-I
LIB_PREFIX=-L
#
MOD_EXT=mod
MKDIR_P=config/install-sh -c -d
INSTALL=/usr/bin/install -c
OBJEXT=o
EXEEXT=
LIBEXT=a
LINK_O_FLAG=-o 

#INCFLAGS must be set by the user makefile

#Dependency rules are created by autoconf according to whether
#discrete preprocessing is necessary or not.
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<

