# Examples makefile

include ../make_cmd.mk
include ../make_inc.mk

EXAMPLES = example1 example2 example3 example4 example5 example6

.IGNORE:

MOD_DIRS = -I../include -I../ext_libs -I../base -I../grid -I../lin_alg -I../fem -I../physics

LIBS = -L../physics -loftphysics -L../fem -loftfem -L../lin_alg -loftla \
	-L../grid -loftgrid -L../base -loftbase -L../ext_libs -loftext $(EXT_LIBS)

all: $(EXAMPLES)

# Remove built elements
clean:
	rm -f $(EXAMPLES)
	$(CLEAN_EXE)
	$(CLEAN_OBJS)
	$(CLEAN_FORT)
	$(CLEAN_RESULTS)

# Make executable object file
%.o: %.f90 FORCE
	$(FC) $(FC_FLAGS) $(COMP_DEFS) -c $< $(INC_DEFS) $(MOD_DIRS)
# Make executable object file (preprocessed)
%.o: %.F90 FORCE
	$(FC) $(FC_FLAGS) $(COMP_DEFS) -c $< $(INC_DEFS) $(MOD_DIRS)
# Make executable
%: %.o FORCE
	$(LD) $(LD_FLAGS) $(FC_FLAGS) $< -o $@ $(LIBS) $(MKL_FLAGS)

FORCE:
