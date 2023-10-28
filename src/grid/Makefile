# Open FUSION Toolkit (OFT) libgrid makefile

include ../make_cmd.mk
include ../make_inc.mk

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .cxx

MOD_DIRS = -I../include -I../base

UTIL_OBJS = \
	bezier_cad.o \
	c_nurbs.o \
	nurbs_cad.o \
	metis_bridge.o \
	quadrature.o \
	lagrange_poly.o

TM_OBJS = \
	mesh_type.o \
	mesh_local.o \
	mesh_local_util.o \
	mesh_global.o \
	mesh_global_util.o \
	trimesh_type.o \
	trimesh_tessellation.o \
	tetmesh_type.o \
	tetmesh_tessellation.o \
	hexmesh_type.o \
	quadmesh_type.o

TM_CAD_OBJS = \
	mesh_native.o \
	mesh_cube.o \
	mesh_sphere.o \
	mesh_t3d.o \
	mesh_cubit.o \
	mesh_gmsh.o

MG_OBJS = \
	multigrid.o \
	multigrid_build.o

LIB_OBJS = $(TM_OBJS) $(UTIL_OBJS) $(TM_CAD_OBJS) $(MG_OBJS)

#---Make alias for eclipse
all: liboftgrid.a

#---Library archiver
liboftgrid.a: $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)
	$(RANLIB) $@

#---Remove built elements
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

#---Module Dependencies
# CAD modules
nurbs_cad.o: c_nurbs.o

mesh_type.o: quadrature.o
# Mesh modules
mesh_local_util.o: mesh_type.o
mesh_local.o: mesh_type.o mesh_local_util.o
mesh_global.o: mesh_local.o mesh_local_util.o
mesh_global_util.o: mesh_global.o

# Tetmesh modules
trimesh_tessellation.o: mesh_type.o
tetmesh_tessellation.o: mesh_type.o lagrange_poly.o
trimesh_type.o: trimesh_tessellation.o quadrature.o lagrange_poly.o
tetmesh_type.o: tetmesh_tessellation.o quadrature.o lagrange_poly.o
# Hex modules
hexmesh_type.o: mesh_type.o quadrature.o lagrange_poly.o
quadmesh_type.o: mesh_type.o quadrature.o lagrange_poly.o

# Mesh type modules
mesh_native.o: multigrid.o
mesh_sphere.o: multigrid.o
mesh_cube.o: multigrid.o
mesh_gmsh.o: multigrid.o
mesh_t3d.o: multigrid.o bezier_cad.o
mesh_cubit.o: multigrid.o nurbs_cad.o

# MG modules
multigrid.o: mesh_global.o mesh_global_util.o tetmesh_type.o trimesh_type.o hexmesh_type.o quadmesh_type.o
multigrid_build.o: $(TM_CAD_OBJS) multigrid.o
