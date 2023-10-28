# Open FUSION Toolkit (OFT) libgrid makefile

include ../make_cmd.mk
include ../make_inc.mk

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .cxx

MOD_DIRS = -I../include -I../base -I../grid -I../lin_alg -I../fem

# Specify dependent modules
OBJS = \
	mhd_utils.o \
	oft_scalar_inits.o \
	oft_vector_inits.o \
	diagnostic.o \
	taylor.o \
	xmhd.o \
	xmhd_lag.o \
	tracing.o \
	tracing_2d.o \
	gs_eq.o \
	axi_green.o \
	grad_shaf.o \
	grad_shaf_profiles.o \
	grad_shaf_mercier.o \
	grad_shaf_util.o \
	grad_shaf_fit.o \
	grad_shaf_td.o \
	thin_wall.o

LIB_OBJS = $(OBJS)

# Make alias for eclipse
all: liboftphysics.a

#---Library archiver
liboftphysics.a: $(LIB_OBJS)
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
oft_scalar_inits.o: mhd_utils.o
oft_vector_inits.o: mhd_utils.o
diagnositic.o: oft_scalar_inits.o oft_vector_inits.o
taylor.o: diagnostic.o
xmhd.o: mhd_utils.o diagnostic.o
xmhd_lag.o: mhd_utils.o diagnostic.o
grad_shaf.o: axi_green.o tracing_2d.o
grad_shaf_profiles.o: grad_shaf.o
grad_shaf_mercier.o: tracing_2d.o grad_shaf_profiles.o
grad_shaf_util.o: grad_shaf.o grad_shaf_profiles.o
grad_shaf_fit.o: grad_shaf_util.o grad_shaf_profiles.o grad_shaf.o
grad_shaf_td.o: grad_shaf.o mhd_utils.o
thin_wall.o: axi_green.o
