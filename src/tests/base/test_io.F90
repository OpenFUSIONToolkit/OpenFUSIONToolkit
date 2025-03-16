!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file test_io.F90
!
!> Regression tests for some HDF5 I/O
!!
!! @authors Chris Hansen
!! @date March 2025
!! @ingroup testing
!-----------------------------------------------------------------------------
PROGRAM test_io
USE oft_base
USE oft_io, ONLY: hdf5_write, hdf5_read, hdf5_create_file
!--Grid
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector
!---Lagrange FE space
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup
IMPLICIT NONE
LOGICAL :: success
INTEGER(i4) :: iounit,ierr
INTEGER(i4), PARAMETER :: ival=10,iarray(3)=[100,101,102],imat(2,2)=RESHAPE([100,101,102,103],[2,2])
INTEGER(i4) :: ival_chk,iarray_chk(3),imat_chk(2,2)
REAL(r8), PARAMETER :: rval=1.d-1,rarray(3)=[1.0d-2,1.1d-2,1.2d-2],rmat(2,2)=RESHAPE([1.0d-2,1.1d-2,1.2d-2,1.3d-2],[2,2])
REAL(r8) :: rval_chk,rarray_chk(3),rmat_chk(2,2),norm,norm_chk
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
CLASS(oft_vector), POINTER :: vec_tmp,vec_chk
INTEGER(i4) :: test_id=1
NAMELIST/test_io_options/test_id
!---------------------------------------------------------------------------
! Initialize enviroment
!---------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=iounit,FILE=oft_env%ifile)
READ(iounit,test_io_options,IOSTAT=ierr)
CLOSE(iounit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---------------------------------------------------------------------------
! Run tests
!---------------------------------------------------------------------------
oft_env%pm=.FALSE.
SELECT CASE(test_id)
  CASE(1) !< Test saving/reading 1D
    IF(oft_env%head_proc)THEN
      CALL hdf5_create_file("oft_test.rst")
      CALL hdf5_write(ival,"oft_test.rst","ival")
      CALL hdf5_write(iarray,"oft_test.rst","iarray")
      CALL hdf5_write(rval,"oft_test.rst","rval")
      CALL hdf5_write(rarray,"oft_test.rst","rarray")
      !
      CALL hdf5_read(ival_chk,"oft_test.rst","ival")
      IF(ival/=ival_chk)CALL oft_abort("Integer values do not match","test_io",__FILE__)
      CALL hdf5_read(iarray_chk,"oft_test.rst","iarray")
      IF(ANY(iarray/=iarray_chk))CALL oft_abort("Integer arrays do not match","test_io",__FILE__)
      CALL hdf5_read(rval_chk,"oft_test.rst","rval")
      IF(rval/=rval_chk)CALL oft_abort("Real values do not match","test_io",__FILE__)
      CALL hdf5_read(rarray_chk,"oft_test.rst","rarray")
      IF(ANY(rarray/=rarray_chk))CALL oft_abort("Real arrays do not match","test_io",__FILE__)
    END IF
  CASE(2) !< Test saving/reading 2D
    IF(oft_env%head_proc)THEN
      CALL hdf5_create_file("oft_test.rst")
      CALL hdf5_write(imat,"oft_test.rst","iarray")
      CALL hdf5_write(rmat,"oft_test.rst","rarray")
      !
      CALL hdf5_read(imat_chk,"oft_test.rst","iarray")
      IF(ANY(imat(:,1)/=imat_chk(:,1)))CALL oft_abort("Integer matrices do not match","test_io",__FILE__)
      IF(ANY(imat(:,2)/=imat_chk(:,2)))CALL oft_abort("Integer matrices do not match","test_io",__FILE__)
      CALL hdf5_read(rmat_chk,"oft_test.rst","rarray")
      IF(ANY(rmat(:,1)/=rmat_chk(:,1)))CALL oft_abort("Real matrices do not match","test_io",__FILE__)
      IF(ANY(rmat(:,2)/=rmat_chk(:,2)))CALL oft_abort("Real matrices do not match","test_io",__FILE__)
    END IF
  CASE(3) !< Test saving/reading vectors
    CALL oft_lag_setup(mg_mesh,2,ML_oft_lagrange,minlev=-1)
    CALL ML_oft_lagrange%vec_create(vec_tmp)
    CALL vec_tmp%set(0.d0,random=.TRUE.)
    norm = vec_tmp%dot(vec_tmp)
    CALL ML_oft_lagrange%vec_create(vec_chk)
    !
    IF(oft_env%head_proc)CALL hdf5_create_file("oft_test.rst")
    CALL oft_mpi_barrier(ierr)
    CALL ML_oft_lagrange%current_level%vec_save(vec_tmp,"oft_test.rst","vec",append=.TRUE.)
    !
    CALL ML_oft_lagrange%current_level%vec_load(vec_chk,"oft_test.rst","vec")
    CALL vec_chk%add(1.d0,-1.d0,vec_tmp)
    norm_chk = vec_chk%dot(vec_chk)
    IF(ABS(norm_chk/norm)>1.d-14)CALL oft_abort("Integer matrices do not match","test_io",__FILE__)
    !
    CALL vec_tmp%delete()
    CALL vec_chk%delete()
    DEALLOCATE(vec_tmp,vec_chk)
  CASE(4) !< Test read missing fields
    CALL oft_lag_setup(mg_mesh,2,ML_oft_lagrange,minlev=-1)
    CALL ML_oft_lagrange%vec_create(vec_tmp)
    !
    IF(oft_env%head_proc)THEN
      CALL hdf5_create_file("oft_test.rst")
      CALL hdf5_write(ival,"oft_test.rst","field")
      CALL hdf5_read(ival_chk,"oft_test.rst","empty",success=success)
      IF(success)CALL oft_abort("Success reported for integer read of missing field","test_io",__FILE__)
      CALL hdf5_read(iarray_chk,"oft_test.rst","empty",success=success)
      IF(success)CALL oft_abort("Success reported for integer array read of missing field","test_io",__FILE__)
      CALL hdf5_read(imat_chk,"oft_test.rst","empty",success=success)
      IF(success)CALL oft_abort("Success reported for integer matrix read of missing field","test_io",__FILE__)
      CALL hdf5_read(rval_chk,"oft_test.rst","empty",success=success)
      IF(success)CALL oft_abort("Success reported for real read of missing field","test_io",__FILE__)
      CALL hdf5_read(rarray_chk,"oft_test.rst","empty",success=success)
      IF(success)CALL oft_abort("Success reported for real array read of missing field","test_io",__FILE__)
      CALL hdf5_read(rmat_chk,"oft_test.rst","empty",success=success)
      IF(success)CALL oft_abort("Success reported for real matrix read of missing field","test_io",__FILE__)
    END IF
    CALL oft_mpi_barrier(ierr)
    CALL ML_oft_lagrange%current_level%vec_load(vec_tmp,"oft_test.rst","empty",err_flag=ierr)
    IF(ierr==0)CALL oft_abort("Success reported for vector read of missing field","test_io",__FILE__)
    !
    CALL vec_tmp%delete()
    DEALLOCATE(vec_tmp)
  CASE(5) !< Test overwritting fields
    CALL oft_lag_setup(mg_mesh,2,ML_oft_lagrange,minlev=-1)
    CALL ML_oft_lagrange%vec_create(vec_tmp)
    !
    IF(oft_env%head_proc)THEN
      CALL hdf5_create_file("oft_test.rst")
      CALL hdf5_write(ival,"oft_test.rst","field")
      CALL hdf5_write(ival,"oft_test.rst","field")
      CALL hdf5_write(iarray,"oft_test.rst","field")
      CALL hdf5_write(imat,"oft_test.rst","iarray")
      CALL hdf5_write(rval,"oft_test.rst","field")
      CALL hdf5_write(rarray,"oft_test.rst","field")
      CALL hdf5_write(rmat,"oft_test.rst","rarray")
    END IF
    CALL oft_mpi_barrier(ierr)
    CALL ML_oft_lagrange%current_level%vec_save(vec_tmp,"oft_test.rst","field",append=.TRUE.)
    !
    CALL vec_tmp%delete()
    DEALLOCATE(vec_tmp)
END SELECT
!---Finalize enviroment
CALL oft_finalize
END PROGRAM test_io