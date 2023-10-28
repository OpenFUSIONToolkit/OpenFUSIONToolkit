# Testing makefile

include ../../make_cmd.mk
include ../../make_inc.mk

.IGNORE:

MOD_DIRS = -I$(OFT_BASE_DIR)/ext_libs -I$(OFT_BASE_DIR)/base -I$(OFT_BASE_DIR)/grid -I$(OFT_BASE_DIR)/lin_alg -I$(OFT_BASE_DIR)/fem -I$(OFT_BASE_DIR)/physics

LIBS = -L$(OFT_BASE_DIR)/ext_libs -L$(OFT_BASE_DIR)/base -L$(OFT_BASE_DIR)/grid -L$(OFT_BASE_DIR)/lin_alg -L$(OFT_BASE_DIR)/fem -L$(OFT_BASE_DIR)/physics \
	-loftphysics -loftfem -loftla -loftgrid -loftbase -loftext $(INT_LIBS) $(EXT_LIBS)

TEST_OBJS = test_gmsh test_cubit test_t3d test_mapping_find test_mapping_jac test_quad

all: $(TEST_OBJS)

# Remove built elements
clean:
	rm -f $(TEST_OBJS)
	$(CLEAN_EXE)
	$(CLEAN_OBJS)
	$(CLEAN_FORT)
	$(CLEAN_RESULTS)

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
