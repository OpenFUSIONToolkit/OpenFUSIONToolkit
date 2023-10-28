# Open FUSION Toolkit (OFT) liboftext makefile

include ../make_cmd.mk
include ../make_inc.mk

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .cxx

QUADPACK_OBJS = \
	dqagse.o \
	dqelg.o \
	dqk21.o \
	dqpsrt.o

LIB_OBJS = dlsode.o lmdif.o bvls.o $(QUADPACK_OBJS)

#---Make alias for eclipse
all: liboftext.a

#---Library archiver
liboftext.a: $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)
	$(RANLIB) $@

#---Remove built elements
clean:
	$(CLEAN_OBJS)
	$(CLEAN_FORT)
	$(CLEAN_LIBS)

# Compile external libraries (LSODE,MINPACK)
dlsode.o: dlsode.f
	$(FC) $(OMP_FLAGS) $(LOPT_FLAGS) -c dlsode.f
lmdif.o: lmdif.f
	$(FC) $(OMP_FLAGS) $(LOPT_FLAGS) -c lmdif.f

#---Generic compile rules
# Generic compile rule for Fortran 90 files
%.o: %.f
	$(FC) $(FC_FLAGS) $(COMP_DEFS) -c $< -I../include $(INC_DEFS)
# Generic compile rule for Fortran 90 files
%.o: %.f90
	$(FC) $(FC_FLAGS) $(COMP_DEFS) -c $< -I../include $(INC_DEFS)
# Generic compile rule for Fortran 90 files (Pre-Processed)
%.o: %.F90
	$(FC) $(FC_FLAGS) $(COMP_DEFS) -c $< -I../include $(INC_DEFS)
# Generic compile rule for C files
%.o: %.c
	$(CC) $(CC_FLAGS) $(COMP_DEFS) -c $< -I../include $(INC_DEFS)
# Generic compile rule for C++ files
%.o: %.cxx
	$(CXX) $(CC_FLAGS) $(COMP_DEFS) -c $< -I../include $(INC_DEFS)
