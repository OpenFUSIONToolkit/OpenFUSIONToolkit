# OFT libbase makefile

include ../make_cmd.mk
include ../make_inc.mk

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .cxx

BASE_OBJS = \
	oft_local_c.o \
	oft_sort_c.o \
	oft_base.o \
	oft_local.o \
	oft_sort.o \
	oft_stitching.o \
	oft_io.o \
	spline.o

LIB_OBJS = $(BASE_OBJS)

#---Make alias for eclipse
all: liboftbase.a

#---Library archiver
liboftbase.a: $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)
	$(RANLIB) $@

#---Remove built elements
clean:
	$(CLEAN_OBJS)
	$(CLEAN_FORT)
	$(CLEAN_LIBS)

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

#---Module Dependencies
oft_base.o: oft_local.o oft_sort.o
oft_sort.o: oft_sort_c.o oft_local.o
oft_stitching.o: oft_base.o oft_sort.o
oft_io.o: oft_base.o
spline.o: oft_base.o