# Testing makefile

include ../../make_cmd.mk
include ../../make_inc.mk

.IGNORE:

MOD_DIRS = -I$(OFT_BASE_DIR)/ext_libs -I$(OFT_BASE_DIR)/base -I$(OFT_BASE_DIR)/grid -I$(OFT_BASE_DIR)/lin_alg -I$(OFT_BASE_DIR)/fem -I$(OFT_BASE_DIR)/physics

LIBS = -L$(OFT_BASE_DIR)/ext_libs -L$(OFT_BASE_DIR)/base -L$(OFT_BASE_DIR)/grid -L$(OFT_BASE_DIR)/lin_alg -L$(OFT_BASE_DIR)/fem -L$(OFT_BASE_DIR)/physics \
        -loftphysics -loftfem -loftla -loftgrid -loftbase -loftext $(INT_LIBS) $(EXT_LIBS)

TEST_OBJS = test_taylor test_taylor_inj test_alfven test_sound test_alfven_lag test_sound_lag

HELPER_OBJS = test_phys_helpers.o

all: $(HELPER_OBJS) $(TEST_OBJS)

# Remove built elements
clean:
	rm -f $(TEST_OBJS)
	$(CLEAN_EXE)
	$(CLEAN_OBJS)
	$(CLEAN_FORT)
	$(CLEAN_RESULTS)
	rm -f floops.loc curr.drive volt.drive tw_fr.dat tw_eigs.dat jumpers_orient.dat
	rm -f f_prof.in p_prof.in bc_mat.dat gs.prof fixed_vflux.dat safety_factor.dat

# Make executable object file
%.o: %.f90 FORCE
	$(FC) $(COMP_DEFS) $(FC_FLAGS) -c $< $(INC_DEFS) $(MOD_DIRS)
# Make executable object file (preprocessed)
%.o: %.F90 FORCE
	$(FC) $(COMP_DEFS) $(FC_FLAGS) -c $< $(INC_DEFS) $(MOD_DIRS)
# Make executable
%: %.o FORCE
	$(LD) $(LD_FLAGS) $(FC_FLAGS) $< -o $@ $(LIBS) $(MKL_FLAGS) $(HELPER_OBJS)

# Module Dependencies
test_alfven.o: test_phys_helpers.o
test_sound.o: test_phys_helpers.o
test_alfven_lag.o: test_phys_helpers.o
test_sound_lag.o: test_phys_helpers.o

FORCE:
