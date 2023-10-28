# Testing makefile

include ../../make_cmd.mk
include ../../make_inc.mk

.IGNORE:

MOD_DIRS = -I$(OFT_BASE_DIR)/ext_libs -I$(OFT_BASE_DIR)/base -I$(OFT_BASE_DIR)/grid -I$(OFT_BASE_DIR)/lin_alg -I$(OFT_BASE_DIR)/fem -I$(OFT_BASE_DIR)/physics

LIBS = -L$(OFT_BASE_DIR)/ext_libs -L$(OFT_BASE_DIR)/base -L$(OFT_BASE_DIR)/grid -L$(OFT_BASE_DIR)/lin_alg -L$(OFT_BASE_DIR)/fem -L$(OFT_BASE_DIR)/physics \
	-loftphysics -loftfem -loftla -loftgrid -loftbase -loftext $(INT_LIBS) $(EXT_LIBS)

LIB_OBJS = marklin_f.o tokamaker_f.o oft_base_f.o

all: liboftpy.so

#---Library archiver
liboftpy.so: $(LIB_OBJS)
	$(FC) -shared -o $@ $(LIB_OBJS) $(LIBS)

# Remove built elements
clean:
	rm -f $(LIB_OBJS)
	$(CLEAN_EXE)
	$(CLEAN_OBJS)
	$(CLEAN_FORT)

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
tokamaker_f.o: oft_base_f.o
marklin_f.o: oft_base_f.o

FORCE:
