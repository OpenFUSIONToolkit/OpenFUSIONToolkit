!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file mugtok_f.F90
!
!> Fortran part of Python wrapper for coupled MUG/TokaMaker (mugtok_td) functionality
!!
!! @ingroup doxy_oft_python
!---------------------------------------------------------------------------------
MODULE mugtok_f
USE iso_c_binding, ONLY: c_int, c_double, c_char, c_loc, c_ptr, &
  c_f_pointer, c_bool, c_associated
USE oft_base
USE oft_mesh_type, ONLY: oft_bmesh
USE mhd_utils, ONLY: mu0
USE oft_gs, ONLY: gs_equil
USE mugtok_td, ONLY: oft_mugtok_td, compute_gradf
USE oft_base_f, ONLY: copy_string
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------------
!> MUG/TokaMaker wrapper object for the Python API
!---------------------------------------------------------------------------------
TYPE :: mugtok_instance
  TYPE(oft_mugtok_td), POINTER :: sim => NULL() !< Coupled MUG/TokaMaker simulation object
  TYPE(gs_equil), POINTER :: equil => NULL() !< Active TokaMaker equilibrium (borrowed, not owned)
  REAL(r8), POINTER, DIMENSION(:,:) :: r_pmesh => NULL() !< Pressure (order-1) node coordinates
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lc_pmesh => NULL() !< Pressure (order-1) cell list
  INTEGER(i4), POINTER, DIMENSION(:) :: reg_pmesh => NULL() !< Pressure (order-1) per-cell region IDs
END TYPE mugtok_instance
CONTAINS
!---------------------------------------------------------------------------------
!> Cast a C pointer back to the Fortran MUG/TokaMaker wrapper object
!---------------------------------------------------------------------------------
FUNCTION mugtok_ccast(mugtok_cptr,sim_obj,error_str) RESULT(success)
TYPE(c_ptr), INTENT(in) :: mugtok_cptr !< C pointer to wrapper object
TYPE(mugtok_instance), POINTER, INTENT(out) :: sim_obj !< Fortran wrapper object
CHARACTER(KIND=c_char), OPTIONAL, INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
LOGICAL :: success
IF(PRESENT(error_str))CALL copy_string('',error_str)
IF(.NOT.c_associated(mugtok_cptr))THEN
  IF(PRESENT(error_str))CALL copy_string('MUG/TokaMaker object not associated',error_str)
  success=.FALSE.
  RETURN
END IF
CALL c_f_pointer(mugtok_cptr,sim_obj)
success=.TRUE.
END FUNCTION mugtok_ccast
!---------------------------------------------------------------------------------
!> Allocate the combined MUG/TokaMaker simulation object from an existing equilibrium
!---------------------------------------------------------------------------------
SUBROUTINE mugtok_alloc(mugtok_ptr,equil_ptr,error_str) BIND(C,NAME="mugtok_alloc")
TYPE(c_ptr), INTENT(out) :: mugtok_ptr !< Pointer to new MUG/TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: equil_ptr !< Pointer to active TokaMaker equilibrium (gs_equil)
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(mugtok_instance), POINTER :: sim_obj
TYPE(gs_equil), POINTER :: equil_obj
CALL copy_string('',error_str)
IF(.NOT.c_associated(equil_ptr))THEN
  CALL copy_string('Equilibrium object not associated',error_str)
  RETURN
END IF
CALL c_f_pointer(equil_ptr,equil_obj)
ALLOCATE(sim_obj)
ALLOCATE(sim_obj%sim)
sim_obj%equil=>equil_obj
mugtok_ptr=C_LOC(sim_obj)
END SUBROUTINE mugtok_alloc
!---------------------------------------------------------------------------------
!> Set up the combined simulation (mesh/FE are reused from the equilibrium's device)
!---------------------------------------------------------------------------------
SUBROUTINE mugtok_setup(mugtok_ptr,dt,lin_tol,nl_tol,mhd_flag_ptr,dens_ptr,visc_ptr, &
                        nreg,incomp,toroidal_flow,error_str) BIND(C,NAME="mugtok_setup")
TYPE(c_ptr), VALUE, INTENT(in) :: mugtok_ptr !< Pointer to MUG/TokaMaker object
REAL(c_double), VALUE, INTENT(in) :: dt !< Timestep [s]
REAL(c_double), VALUE, INTENT(in) :: lin_tol !< Linear solver tolerance
REAL(c_double), VALUE, INTENT(in) :: nl_tol !< Non-linear solver tolerance
TYPE(c_ptr), VALUE, INTENT(in) :: mhd_flag_ptr !< Per-region MHD flag (int32, 1=MHD region)
TYPE(c_ptr), VALUE, INTENT(in) :: dens_ptr !< Per-region density (<0 where unused)
TYPE(c_ptr), VALUE, INTENT(in) :: visc_ptr !< Per-region kinematic viscosity (<0 where unused)
INTEGER(c_int), VALUE, INTENT(in) :: nreg !< Number of mesh regions
LOGICAL(c_bool), VALUE, INTENT(in) :: incomp !< Incompressible flow?
LOGICAL(c_bool), VALUE, INTENT(in) :: toroidal_flow !< Allow toroidal (phi) flow in the MHD regions?
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(mugtok_instance), POINTER :: sim_obj
INTEGER(i4), POINTER, DIMENSION(:) :: mhd_int
REAL(r8), POINTER, DIMENSION(:) :: dens_tmp,visc_tmp
LOGICAL, ALLOCATABLE, DIMENSION(:) :: mhd_flag
IF(.NOT.mugtok_ccast(mugtok_ptr,sim_obj,error_str))RETURN
CALL c_f_pointer(mhd_flag_ptr,mhd_int,[nreg])
CALL c_f_pointer(dens_ptr,dens_tmp,[nreg])
CALL c_f_pointer(visc_ptr,visc_tmp,[nreg])
ALLOCATE(mhd_flag(nreg))
mhd_flag=(mhd_int/=0)
CALL sim_obj%sim%setup(sim_obj%equil,dt,lin_tol,nl_tol,mhd_flag,dens_tmp,visc_tmp, &
                       incomp=LOGICAL(incomp),toroidal_flow=LOGICAL(toroidal_flow))
DEALLOCATE(mhd_flag)
END SUBROUTINE mugtok_setup
!---------------------------------------------------------------------------------
!> Advance the coupled solution by one timestep
!---------------------------------------------------------------------------------
SUBROUTINE mugtok_step(mugtok_ptr,curr_ptr,ncoils,time,dt,nl_its,lin_its,nretry,error_str) &
    BIND(C,NAME="mugtok_step")
TYPE(c_ptr), VALUE, INTENT(in) :: mugtok_ptr !< Pointer to MUG/TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: curr_ptr !< Coil currents [A] at end of step
INTEGER(c_int), VALUE, INTENT(in) :: ncoils !< Number of coils
REAL(c_double), INTENT(inout) :: time !< Current simulation time [s]
REAL(c_double), INTENT(inout) :: dt !< Timestep size [s]
INTEGER(c_int), INTENT(out) :: nl_its !< Number of nonlinear iterations
INTEGER(c_int), INTENT(out) :: lin_its !< Number of linear iterations
INTEGER(c_int), INTENT(out) :: nretry !< Number of retries
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(mugtok_instance), POINTER :: sim_obj
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
REAL(8), ALLOCATABLE, DIMENSION(:) :: coil_currents
IF(.NOT.mugtok_ccast(mugtok_ptr,sim_obj,error_str))RETURN
ALLOCATE(coil_currents(ncoils))
IF(ncoils>0)THEN
  CALL c_f_pointer(curr_ptr,vals_tmp,[ncoils])
  coil_currents=vals_tmp*mu0 ! Python passes [A]; Fortran uses [A]*mu0
END IF
CALL sim_obj%sim%step(coil_currents,time,dt,nl_its,lin_its,nretry)
DEALLOCATE(coil_currents)
END SUBROUTINE mugtok_step
!---------------------------------------------------------------------------------
!> Get the current values of a MUG-owned solution field at node points
!---------------------------------------------------------------------------------
SUBROUTINE mugtok_get_field(mugtok_ptr,field_id,vals,error_str) BIND(C,NAME="mugtok_get_field")
TYPE(c_ptr), VALUE, INTENT(in) :: mugtok_ptr !< Pointer to MUG/TokaMaker object
INTEGER(c_int), VALUE, INTENT(in) :: field_id !< Field index (1=n,2/3/4=vel,5=P,7=F)
TYPE(c_ptr), VALUE, INTENT(in) :: vals !< Output array (size = field's FE block DOFs)
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(mugtok_instance), POINTER :: sim_obj
REAL(r8), POINTER, DIMENSION(:) :: vals_tmp
IF(.NOT.mugtok_ccast(mugtok_ptr,sim_obj,error_str))RETURN
CALL c_f_pointer(vals,vals_tmp,[sim_obj%sim%fe_rep%fields(field_id)%fe%ne])
CALL sim_obj%sim%u%get_local(vals_tmp,field_id) ! fill the caller's buffer directly (no temp/copy)
END SUBROUTINE mugtok_get_field
!---------------------------------------------------------------------------------
!> Get the projected poloidal gradient of F at node points (dF/dR, dF/dZ)
!---------------------------------------------------------------------------------
SUBROUTINE mugtok_get_gradf(mugtok_ptr,gr,gz,error_str) BIND(C,NAME="mugtok_get_gradf")
TYPE(c_ptr), VALUE, INTENT(in) :: mugtok_ptr !< Pointer to MUG/TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: gr !< Output dF/dR array 
TYPE(c_ptr), VALUE, INTENT(in) :: gz !< Output dF/dZ array
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(mugtok_instance), POINTER :: sim_obj
REAL(r8), POINTER, DIMENSION(:) :: gr_tmp,gz_tmp,buf
INTEGER(i4) :: n
IF(.NOT.mugtok_ccast(mugtok_ptr,sim_obj,error_str))RETURN
NULLIFY(buf)
CALL sim_obj%sim%u%get_local(buf,7) ! size the output to field 7 (F)
n=SIZE(buf)
DEALLOCATE(buf)
CALL c_f_pointer(gr,gr_tmp,[n])
CALL c_f_pointer(gz,gz_tmp,[n])
CALL compute_gradf(sim_obj%sim,gr_tmp,gz_tmp)
END SUBROUTINE mugtok_get_gradf
!---------------------------------------------------------------------------------
!> Get the order-1 (pressure) plotting mesh: node coordinates and triangle connectivity
!---------------------------------------------------------------------------------
SUBROUTINE mugtok_get_pmesh(mugtok_ptr,np,r_loc,nc,lc_loc,reg_loc,error_str) BIND(C,NAME="mugtok_get_pmesh")
TYPE(c_ptr), VALUE, INTENT(in) :: mugtok_ptr !< Pointer to MUG/TokaMaker object
INTEGER(c_int), INTENT(out) :: np !< Number of pressure node points
TYPE(c_ptr), INTENT(out) :: r_loc !< Node coordinate list pointer
INTEGER(c_int), INTENT(out) :: nc !< Number of cells
TYPE(c_ptr), INTENT(out) :: lc_loc !< Cell (triangle) list pointer
TYPE(c_ptr), INTENT(out) :: reg_loc !< Per-cell region ID list pointer
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(mugtok_instance), POINTER :: sim_obj
CLASS(oft_bmesh), POINTER :: mesh
INTEGER(i4) :: i,j,k,id
IF(.NOT.mugtok_ccast(mugtok_ptr,sim_obj,error_str))RETURN
mesh=>sim_obj%equil%device%mesh
!---Tessellate at the pressure FE order (main FE order - 1) so the plot nodes align with the
!   order-(N-1) pressure DOFs
CALL mesh%tessellate(sim_obj%r_pmesh,sim_obj%lc_pmesh,sim_obj%equil%device%fe_rep%order-1)
np=SIZE(sim_obj%r_pmesh,DIM=2,KIND=c_int)
nc=SIZE(sim_obj%lc_pmesh,DIM=2,KIND=c_int)
r_loc=c_loc(sim_obj%r_pmesh)
lc_loc=c_loc(sim_obj%lc_pmesh)
!---Per-cell region IDs
ALLOCATE(sim_obj%reg_pmesh(nc))
k=nc/mesh%nc
IF(ASSOCIATED(mesh%reg))THEN
  DO i=1,mesh%nc
    id=mesh%reg(i)
    DO j=1,k
      sim_obj%reg_pmesh((i-1)*k+j)=id
    END DO
  END DO
ELSE
  sim_obj%reg_pmesh=0
END IF
reg_loc=c_loc(sim_obj%reg_pmesh)
END SUBROUTINE mugtok_get_pmesh
!---------------------------------------------------------------------------------
!> Destroy the combined simulation object (the borrowed equilibrium/mesh are left intact)
!---------------------------------------------------------------------------------
SUBROUTINE mugtok_destroy(mugtok_ptr,error_str) BIND(C,NAME="mugtok_destroy")
TYPE(c_ptr), VALUE, INTENT(in) :: mugtok_ptr !< Pointer to MUG/TokaMaker object
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(mugtok_instance), POINTER :: sim_obj
IF(.NOT.mugtok_ccast(mugtok_ptr,sim_obj,error_str))RETURN
CALL sim_obj%sim%delete()
DEALLOCATE(sim_obj%sim)
IF(ASSOCIATED(sim_obj%r_pmesh))DEALLOCATE(sim_obj%r_pmesh)
IF(ASSOCIATED(sim_obj%lc_pmesh))DEALLOCATE(sim_obj%lc_pmesh)
IF(ASSOCIATED(sim_obj%reg_pmesh))DEALLOCATE(sim_obj%reg_pmesh)
DEALLOCATE(sim_obj)
END SUBROUTINE mugtok_destroy
END MODULE mugtok_f
