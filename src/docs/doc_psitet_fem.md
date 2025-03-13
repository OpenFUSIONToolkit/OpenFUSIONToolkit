Finite Element Representations      {#doc_oft_fem}
==============================

[TOC]

This page addresses general considerations with the implementation of finite element methods in the Open FUSION Toolkit (OFT).
In particular basis set definition, field creation, and operator construction will be addressed. OFT
currently support C0 Galerkin finite element methods with volume based interaction only, where the discretized
system matrix has the form

\f[ \int_{\Omega} g( u^T ) f( u ) dV + \int_{\partial \Omega} u^T h( u ) dS \f]

where \f$ u \f$ and \f$ u^T \f$ are the basis and test functions respectively, which for Galerkin weighting are
the same set of functions. Matrix elements exist for this system when the row and column correspond elements that share
a common cell. This limitation is due to the linkage and seam construction methods, which currently do no support
connectivity using internal faces.

This document explains the basic methods to define and implement a finite element representation using a
scalar Lagrange representation as an example. For information about a specific finite element representation
please view its documentation. These documents may be accessed through the \ref fem "finite element" group.

\section doc_oft_fem_def Basis Definition

A finite element representation is defined by the basis set on the unit tetrahedron. Each basis function is
tied to a geometric entity (vertex, edge, etc.) according the logical coordinates it depends on. In the case
of a quadratic Lagrange representation scalars are represented by a single basis function associated with each
vertex and edge. These functions are

\f$ f^v_{i} (u_i) = 2 u_i ( u_i - 1/2 ) \f$

\f$ f^e_{i} (u^1_i,u^2_i) = 4 u^1_i u^2_i \f$

where \f$ u_i \f$ is the logical coordinate associated with the i-th vertex and \f$ u^1_i, u^2_i \f$ are the
coordinates associated with the end points of the i-th edge.

In OFT a finite element representation is defined by the \ref fem_base::oft_fem_type "oft_fem_type"
class, which contains the structures required to define geometric connectivity, vector definition, and
the underlying graph for discretized operators. Most of the setup and construction of this structure is
automated and relies only on the mesh and association of the elements to geometric primitives.

Shown below is an example of the finite element setup procedure for the quadratic Lagrange set presented above.
First the mesh is linked to the structure and the order of the representation is set, quadratic in this case. The
number of basis functions per geometric primitive is then set in the `gstruct` field. Where `gstruct` is a 4 value
integer array containing the number of basis functions associated with each vertex, edge, face and cell respectively.
Finally, the remaining setup can be completed using the \ref fem_base::oft_fem_type::setup "setup" method, where the
argument is the desired quadrature order.

~~~~~~~~~{.F90}
USE tetmesh_local, ONLY: mesh
USE fem_base, ONLY: oft_fem_type
USE oft_lag_basis, ONLY: oft_lagrange_id
TYPE(oft_fem_type), POINTER :: quadratic_lagrange
quadratic_lagrange%mesh=>mesh
quadratic_lagrange%order=2
quadratic_lagrange%type=oft_lagrange_id
quadratic_lagrange%gstruct=(/1,1,0,0/)
CALL quadratic_lagrange%setup(5)
~~~~~~~~~

\subsection doc_oft_fem_def_eval Basis Evaluation

Although the basic structure is setup automatically subroutines must be provided to evaluate the basis functions and
their derivatives if desired. For Lagrange elements the basis functions do not require the Jacobian matrix for the spatial
mapping, so instead the gradient evaluation routine will be presented to illustrate all features. The ordinary basis functions
can be evaluated using the same structure, without the Jacobian terms.

Below is annotated source code for the \ref oft_lag_basis::oft_lag_geval "oft_lag_geval" subroutine, which evaluates the gradient
in physical coordinates of a given basis function. The gradient is dependent both on the derivatives of the basis function in
logical coordinates as well as the Jacobian of the logical to physical mapping for the current cell. Where physical derivatives
are evaluated using the chain rule as

\f$ \frac{\partial f}{\partial x_i} = \frac{\partial f}{\partial u_j} \frac{\partial u_j}{\partial x_i} \f$

Spatial derivatives of the logical coordinates are computed using the subroutine \ref tetmesh_mapping::tetmesh_get_jacobian
"tetmesh_get_jacobian" and passed to the routines through the `gop` parameter. In order to evaluate the function the position in
the mesh must be specified, by the cell index `cell` and logical coordinates within the cell `f`, along with the desired basis
function within the cell to evaluate `dof`. Basis functions are numbered by the \ref fem_base::fem_setup "fem_setup" routine
according to their geometric linkage. The resulting gradient is then returned using the `val` parameter. For the existing finite
element representations element specific evaluation routines are called from this driver routine depending on the element linkage,
specified in the \ref fem_base::oft_fem_type::cmap "cmap" structure.

~~~~~~~~~{.F90}
USE oft_base
USE tetmesh_type, ONLY: tet_ed, tet_fc
USE oft_lag_basis, ONLY: oft_lag_gevalp, oft_lag_gevale, oft_lag_gevalf, oft_lag_gevalc
SUBROUTINE oft_lag_geval(self,cell,dof,f,val,gop)
CLASS(oft_fem_type), INTENT(in) :: self
INTEGER(i4), INTENT(in) :: cell,dof
REAL(r8), INTENT(in) :: f(:)
REAL(r8), OPTIONAL, INTENT(in) :: gop(3,4)
REAL(r8), INTENT(out) :: val(3)
REAL(r8) :: grads(3,4),cofs(4)
INTEGER(i4) :: ed,etmp(2),fc,ftmp(3),i
val=0.d0
IF(PRESENT(gop))THEN
  grads=gop
ELSE
  CALL oft_abort('Logical gradients not specified.','oft_lag_geval',__FILE__)
END IF
SELECT CASE(self%cmap(dof)%type)
~~~~~~~~~

If the requested element is a vertex element no additional handling is needed and the evaluation
routine may be called.

~~~~~~~~~{.F90}
  CASE(1)
    CALL oft_lag_gevalp(self%order,self%cmap(dof)%el,f,cofs)
~~~~~~~~~

If the requested element is a edge element the edge must first be orient to be globally consistent.
This is performed by first retrieving the orientation flag for the local edge \ref tetmesh_type::tet_mesh::lce
"lce" and the end point vertices in the local orientation \ref tetmesh_type::tet_ed "tet_ed". The end
points can then be oriented to the global element orientation using \ref local::orient_list2 "orient_list2" and
the evaluation routine called.

~~~~~~~~~{.F90}
  CASE(2)
    ed=self%mesh%lce(self%cmap(dof)%el,cell)
    etmp=(/tet_ed(1,self%cmap(dof)%el),tet_ed(2,self%cmap(dof)%el)/)
    CALL orient_list2(ed,etmp)
    CALL oft_lag_gevale(self%order,etmp,self%cmap(dof)%ind,f,cofs)
~~~~~~~~~

If the requested element is a face element the face must also be orient to be globally consistent.
This is performed by first retrieving the orientation flag for the local face \ref tetmesh_type::tet_mesh::lcfo
"lcfo" and the corner vertices in the local orientation \ref tetmesh_type::tet_fc "tet_fc". The corner points
can then be oriented to the global element orientation using \ref oft_local::orient_listn "orient_listn" and the evaluation routine called.

~~~~~~~~~{.F90}
  CASE(3)
    fc=self%mesh%lcfo(self%cmap(dof)%el,cell)
    ftmp=(/tet_fc(1,self%cmap(dof)%el),tet_fc(2,self%cmap(dof)%el),tet_fc(3,self%cmap(dof)%el)/)
    CALL orient_listn(fc,ftmp,3_i4)
    CALL oft_lag_gevalf(self%order,ftmp,self%cmap(dof)%ind,f,cofs)
~~~~~~~~~

If the requested element is a cell element no additional handling is needed and the evaluation
routine may be called.

~~~~~~~~~{.F90}
  CASE(4)
    CALL oft_lag_gevalc(self%order,self%cmap(dof)%ind,f,cofs)
END SELECT
~~~~~~~~~

Each of the evaluation routines returns and array of values corresponding to \f$ \frac{\partial f}{\partial u_j} \f$.
In order to reconstruct the gradient in physical coordinates these coefficients must be combined with the jacobian terms.

~~~~~~~~~{.F90}
!---Sum contributions
DO i=1,4
  val=val+grads(:,i)*cofs(i)
END DO
END SUBROUTINE oft_lag_geval
~~~~~~~~~

\section doc_oft_fem_fields Field Creation

Once the finite element structure has been constructed it can be used to create a vector representation of the basis
function weights. These vectors can then be used in the solution of linear systems formed using this finite element
space. Vectors can be created from the finite element structure by using the \ref fem_base::oft_fem_type::vec_create
"vec_create" method.

~~~~~~~~~{.F90}
USE oft_vectors, ONLY: oft_field
CLASS(oft_field), POINTER :: lagrange_vec
CALL quadratic_lagrange%vec_create(lagrange_vec)
~~~~~~~~~

\subsection doc_oft_fem_fields_rst Binary Restart I/O

The \ref fem_base::oft_fem_type "oft_fem_type" class also provides functionality for saving and reading FE weight
vectors in binary restart files, using the \ref fem_base::oft_fem_type::vec_save "vec_save" and
\ref fem_base::oft_fem_type::vec_load "vec_load" methods. Restart files utilize the HDF5 library and perform I/O in
parallel using a common file for all processes.

~~~~~~~~~{.F90}
CALL quadratic_lagrange%vec_save(lagrange_vec,'output.rst','V')
~~~~~~~~~

\section doc_oft_fem_ops Operator Construction

The supporting structure presented so far is primarily used to support the construction of discretized operators
defined by a finite element projection of the weak form of the desired PDE. The example below outlines the basics of this
operation for the Laplacian operator, whose weak form is

\f[ \nabla^2 u = - \int_{\Omega} \nabla \phi \cdot \nabla u dV + \int_{\Gamma} \phi \nabla u \cdot dS \f]

For this example we will show the operator construction of the volume term using a Galerkin method where both the solution
and test functions are expanded using a Lagrange basis. A Dirichlet boundary condition for all boundaries will also be imposed
for this example, by eliminating the boundary rows and replacing the corresponding diagonal entry with 1.

\f[ \phi = \sum \phi_i f_i \f]
\f[ u = \sum u_j f_j \f]
\f[ f \in F_{lagrange} \f]

\note This example performs the integration loop in parallel on each domain using OpenMP. Parallelization is done across mesh cells with
local private variables and thread gaurds used to prevent data races. This method has been seen to scale well dispite the blocking
presented by the thread gaurds as the majority of the loop is spent computing local entries allowing efficient staggering of the serialized
memory accesses.

\subsection doc_oft_fem_ops_alloc Allocation and Setup

Before we can compute the entries of the matrix corresponding to this operator the matrix representation must be set up and space
allocated for the all entries. This process is performed automatically for a standard FE representation by the
\ref fem_base::oft_fem_type::mat_create "mat_create" method.

~~~~~~~~~{.F90}
USE tetmesh_global_util, ONLY: tetmesh_global_curved
USE tetmesh_mapping, ONLY: tetmesh_get_jacobian
!---
USE oft_matrices, ONLY: oft_matrix, oft_graph_ptr
USE oft_la_utils, ONLY: create_matrix
!---
USE oft_lag_basis, ONLY: oft_lag_geval
CLASS(oft_matrix), POINTER :: mat
INTEGER(i4) :: i,k,l,jr,jc,jp,jn,m,node
INTEGER(i4), ALLOCATABLE :: j(:),jk(:,:)
REAL(r8) :: goptmp(3,4)
REAL(r8) :: v,f(4),det,vol
REAL(r8), ALLOCATABLE :: gop(:,:,:),dets(:),lop(:,:)
LOGICAL :: curved
!---------------------------------------------------------------------------
! Allocate Operator
!---------------------------------------------------------------------------
CALL quadratic_lagrange%mat_create(mat)
~~~~~~~~~

The \ref fem_base::oft_fem_type::mat_create "mat_create" method returns a bare matrix representation for an operator with full interaction
between all FE weights that share a common cell. It contains all of the required mapping and structural information to support computing
matrix entries and performing matrix-vector operations.

\note At this point the matix is not usable for linear algebra operations, all entries must be set and `assemble` called at least
once before matrix vector products may be computed.

\subsection doc_oft_fem_ops_int Computing Entries

Matrix entries for the desired operator correpsond to the result of the above volume integral for each basis/test function pair

\f[ A_{ij} = \int \nabla f_i \cdot \nabla f_j dV \f]

where \f$ f_i \f$ and \f$ f_j \f$ are independent Lagrange basis functions. As each function only has a non-zero value in cells
which contain that particular element this integral simplfies to sum of contributions from each individual cell. Integration with in each
cell is then carried out using a numerical integration with a given quadrature stencil. The desired order of the quadrature scheme for this
case was set during finite element setup, so the appropriate quadrature object is already bound to the `quadratic_lagrange` structure.

The first section of the operator integration loop sets up local variables which are unique in each cell. These variables are declared private
in the OpenMP clause and must be allocated after parallel execution has begun.

~~~~~~~~~{.F90}
!---------------------------------------------------------------------------
! Setup parallel integration variables
!---------------------------------------------------------------------------
!$omp parallel private(j,gop,dets,lop,curved,goptmp,m,v,jc,jr,f)
ALLOCATE(j(quadratic_lagrange%nce)) ! Local DOF and matrix indices
ALLOCATE(gop(3,quadratic_lagrange%nce,quadratic_lagrange%quad%np)) ! Reconstructed gradient operator
ALLOCATE(dets(quadratic_lagrange%quad%np)) ! Quadrature determinant values
ALLOCATE(lop(quadratic_lagrange%nce,quadratic_lagrange%nce)) ! Local laplacian matrix
~~~~~~~~~

With the local environment setup the main integration begins, looping over all cells in the parent mesh. For each cell in the mesh a test is
performed to see if any curved geometry exists for that cell. If there are no curved entities then the logical to physical Jacobian is constant
over the entire cell and only needs to be computed once, otherwise it must be computed at each integration point. The list of elements in the
current cell is the retrieved to allowing mapping of the local matrix contributions into the full matrix.

~~~~~~~~~{.F90}
!$omp do
DO i=1,quadratic_lagrange%mesh%nc
  !---Straight cell test
  curved=tetmesh_global_curved(quadratic_lagrange%mesh,i)
  !---Get list of element in the current cell
  CALL quadratic_lagrange%ncdofs(i,j)
~~~~~~~~~

The integrand components are then computed at each quadrature point. The physical Jacobian is first computed if necessary at the current quadrature
point. The weight value is then set by computing the product of the quadrature weight and the volume element Jacobian. Finally, the basis function
gradients \f$ \left( \nabla f \right) \f$ are computed for each of the elements in the current cell.

~~~~~~~~~{.F90}
!---------------------------------------------------------------------------
! Get local reconstructed operators
!---------------------------------------------------------------------------
  DO m=1,quadratic_lagrange%quad%np ! Loop over quadrature points
    IF(curved.OR.(m==1))CALL tetmesh_get_jacobian(quadratic_lagrange%mesh,i,quadratic_lagrange%quad%pts(:,m),goptmp,v)
    dets(m)=v*quadratic_lagrange%quad%wts(m)
    DO jc=1,quadratic_lagrange%nce ! Loop over degrees of freedom
      CALL oft_lag_geval(oft_lagrange,i,jc,quadratic_lagrange%quad%pts(:,m),gop(:,jc,m),goptmp)
    END DO
  END DO
~~~~~~~~~

The local matrix entries are now given as

\f[ A_{ij} = \sum_k \nabla f_i (x_k) \cdot \nabla f_j (x_k) \omega_k \left| \mathbf{J}(x_k) \right| \f]

where \f$ x_k \f$ and \f$ \omega_k \f$ are the k-th quadrature point and weight respectively and \f$ \mathbf{J} \f$ is the Jacobian matrix for the
logical to physical mapping. These entries can be readily computed with the existing information by looping over the rows, columns, and quadrature
points. Rows corresponding to boundary elements are skipped consistent with the desired boundary condition.

Finally, the local matrix contributions are added to the full matrix using the \ref oft_matrices::oft_matrix::add_values "add_values" method. This
method takes the local matrix and list of local variables and uses them to add the local contributions into the full matrix. Depending on the linear
algebra backed end being used the operations performed by this method are different but the calling structure is the same. This call must be
surrounded by an <tt>!$omp critical</tt> clause to prevent data races when added contributions into the main matrix. This clause restricts execution
of the given region to a single thread at a time.

~~~~~~~~~{.F90}
!---------------------------------------------------------------------------
! Compute local matrix contributions
!---------------------------------------------------------------------------
  DO jr=1,quadratic_lagrange%nce
    lop(jr,:)=0.d0
    IF(quadratic_lagrange%global%gbe(j(jr)))CYCLE
    DO jc=1,quadratic_lagrange%nce ! Compute quadrature for M and L
      DO m=1,quadratic_lagrange%quad%np
        lop(jr,jc) = lop(jr,jc) + dot_product(gop(:,jr,m),gop(:,jc,m))*dets(m)
      END DO
    END DO
  END DO
!---------------------------------------------------------------------------
! Add local values to global matrix
!---------------------------------------------------------------------------
  !$omp critical
  CALL mat%add_values(j,j,lop,quadratic_lagrange%nce,quadratic_lagrange%nce)
  !$omp end critical
END DO
~~~~~~~~~

Once the loop has completed local variables are deallocated and the matrix is assembled. This assembly call is needed for linear algebra backends such
as PETSc. The diagonal entries corresponding to removed boundary rows are then set to 1 in order to produce a non-singular system and allowing imposing
the boundary value by modifying the RHS vector. Diagonal entries are only set for rows that correspond to boundary elements and are only set on the
domain with ownership of that element, as indicated by \ref tetmesh_stitching::tetmesh_stitch_type::leo "leo". \ref oft_matrices::oft_matrix::assemble
"Assemble" is then called again to finalize the matrix and make it available for use in linear algebra operations.

~~~~~~~~~{.F90}
DEALLOCATE(j)
DEALLOCATE(gop,dets,lop)
!$omp end parallel
CALL mat%assemble
!---------------------------------------------------------------------------
! Set diagonal entries for dirichlet rows
!---------------------------------------------------------------------------
ALLOCATE(lop(1,1),j(1))
lop(1,1)=1.d0
DO i=1,quadratic_lagrange%nbe
  jr=quadratic_lagrange%lbe(i)
  IF(.NOT.quadratic_lagrange%global%gbe(jr))CYCLE
  IF(.NOT.quadratic_lagrange%linkage%leo(i))CYCLE
  j=jr
  CALL mat%add_values(j,j,lop,1,1)
END DO
DEALLOCATE(j,lop)
CALL mat%assemble(oft_lag_vec)
~~~~~~~~~
