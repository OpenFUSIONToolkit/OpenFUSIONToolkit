# Open FUSION Toolkit (OFT) libla makefile

include ../make_cmd.mk
include ../make_inc.mk

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .cxx

MOD_DIRS = -I../include -I../base -I../grid

# Specify dependent modules
LA_OBJS = \
	lin_alg_base.o \
	native_la.o \
	petsc_la.o \
	deriv_matrices.o

SOLVER_OBJS = \
	solver_base.o \
	native_solvers.o \
	petsc_solvers.o \
	superlu_bridge.o \
	superlu_dist_bridge.o \
	umfpack_bridge.o \
	mkl_bridge.o \
	oft_arpack.o \
	oft_lu.o

MGMT_OBJS = \
	lin_alg_utils.o \
	solver_utils.o

LIB_OBJS = $(LA_OBJS) $(SOLVER_OBJS) $(MGMT_OBJS)

# Make alias for eclipse
all: liboftla.a

# Remove built elements
clean:
	$(CLEAN_OBJS)
	$(CLEAN_FORT)
	$(CLEAN_LIBS)

#---Generic compile rules
# Generic compile rule for Fortran 90 files
%.o: %.f90
	$(FC) $(FC_FLAGS) $(COMP_DEFS) -c $< $(MOD_DIRS) $(INC_DEFS)
# Generic compile rule for Fortran 90 files (Pre-Processed)
%.o: %.F90
	$(FC) $(FC_FLAGS) $(COMP_DEFS) -c $< $(MOD_DIRS) $(INC_DEFS)
# Generic compile rule for C files
%.o: %.c
	$(CC) $(CC_FLAGS) $(COMP_DEFS) -c $< $(MOD_DIRS) $(INC_DEFS)
# Generic compile rule for C++ files
%.o: %.cxx
	$(CXX) $(CC_FLAGS) $(COMP_DEFS) -c $< $(MOD_DIRS) $(INC_DEFS)

liboftla.a: $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)
	$(RANLIB) $@

# Module Dependencies
native_la.o: lin_alg_base.o
petsc_la.o: lin_alg_base.o
deriv_matrices.o: lin_alg_base.o
# Solver Modules
solver_base.o: lin_alg_base.o
native_solvers.o: solver_base.o native_la.o
petsc_solvers.o: solver_base.o petsc_la.o
oft_arpack.o: solver_base.o
oft_lu.o: solver_base.o native_la.o
# Management Modules
lin_alg_utils.o: native_la.o petsc_la.o
solver_utils.o: native_solvers.o petsc_solvers.o oft_lu.o
