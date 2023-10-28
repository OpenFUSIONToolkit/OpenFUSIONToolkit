# Open FUSION Toolkit (OFT) libgrid makefile

include ../make_cmd.mk
include ../make_inc.mk

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .cxx

MOD_DIRS = -I../include -I../base -I../grid -I../lin_alg

# Specify dependent modules

FEM_OBJS = \
	fem_base.o \
	fem_composite.o \
	fem_utils.o \
	lagrange_basis.o \
	lagrange_fields.o \
	lagrange_operators.o \
	blag_operators.o \
	hcurl_basis.o \
	hcurl_fields.o \
	hcurl_operators.o \
	h0_basis.o \
	h0_fields.o \
	h0_operators.o \
	h1_basis.o \
	h1_fields.o \
	h1_operators.o

LIB_OBJS = $(FEM_OBJS)

# Make alias for eclipse
all: liboftfem.a

#---Library archiver
liboftfem.a: $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)
	$(RANLIB) $@

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

# Module Dependencies
# FEM Modules
fem_composite.o: fem_base.o
fem_utils.o: fem_base.o fem_composite.o
lagrange_basis.o: fem_base.o fem_composite.o
lagrange_fields.o: lagrange_basis.o
lagrange_operators.o: lagrange_fields.o lagrange_basis.o
blag_operators.o: lagrange_operators.o
hcurl_basis.o: fem_base.o fem_composite.o
hcurl_fields.o: hcurl_basis.o
hcurl_operators.o: hcurl_fields.o hcurl_basis.o lagrange_operators.o
h0_basis.o: fem_base.o fem_composite.o
h0_fields.o: h0_basis.o
h0_operators.o: h0_fields.o h0_basis.o
h1_basis.o: h0_basis.o hcurl_basis.o
h1_fields.o: h0_fields.o hcurl_fields.o h1_basis.o
h1_operators.o: h1_fields.o h1_basis.o hcurl_operators.o
