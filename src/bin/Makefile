# OpenFUSIONToolkit makefile
# All platform specific customization is done
# in associated make_libs.mk files for each system

include ../make_cmd.mk
include ../make_inc.mk

.IGNORE:

MOD_DIRS = -I../include -I../ext_libs -I../base -I../grid -I../lin_alg -I../fem -I../physics

LIBS = -L../physics -loftphysics -L../fem -loftfem -L../lin_alg -loftla \
	-L../grid -loftgrid -L../base -loftbase -L../ext_libs -loftext $(EXT_LIBS)

OFT_BINS = oft_mesh_check oft_poincare oft_trace

MARKLING_BINS = marklin_eigs

TOKAMAKER_BINS = tokamaker_gs tokamaker_fit tokamaker_wall

THINCURR_BINS = thincurr_eig thincurr_fr thincurr_td thincurr_from_mode

# Make alias for eclipse
all: oft marklin tokamaker thincurr

oft: $(OFT_BINS)

marklin: $(MARKLING_BINS)

tokamaker: $(TOKAMAKER_BINS)

thincurr: $(THINCURR_BINS)

# Remove built elements
clean:
	rm -f $(OFT_BINS) $(MARKLING_BINS) $(TOKAMAKER_BINS) $(THINCURR_BINS)
	$(CLEAN_OBJS)
	$(CLEAN_FORT)
	$(CLEAN_EXE)

# Make executable object file
%.o: %.f90 FORCE
	$(FC) $(COMP_DEFS) $(FC_FLAGS) -c $< $(INC_DEFS) $(MOD_DIRS)
# Make executable object file (preprocessed)
%.o: %.F90 FORCE
	$(FC) $(COMP_DEFS) $(FC_FLAGS) -c $< $(INC_DEFS) $(MOD_DIRS)
# Make executable
%: %.o FORCE
	$(LD) $(LD_FLAGS) $(FC_FLAGS) $< -o $@ $(LIBS) $(MKL_FLAGS)

FORCE:
