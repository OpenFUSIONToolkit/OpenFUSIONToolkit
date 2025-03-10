Composite Finite Element Representations      {#doc_comp_fem}
========================================

[TOC]

The Open FUSION Toolkit supports the construction of composite vectors and matrices from sets of finite element
fields. This is primarily used to support Multi-Physics solution where multiple fields are advanced
simultaneously. It can also be used to represent vector valued fields with scalar finite elements
for each component. A composite FE representation is defined in terms of the sub-field representations
through the \ref fem_composite::oft_fem_comp_type "oft_fem_comp_type" class.

\section doc_comp_fem_const Defining a composite representation

For most use cases composite representations can be assmbled by simply linking to the sub-field representations
as shown in the example below. The \ref fem_composite::oft_fem_comp_type "oft_fem_comp_type" class requires the
number of sub-fields along with a reference to the FE representation and a character tag for each sub-field. The
character tags are used when for restart files I/O using the \ref fem_composite::oft_fem_comp_type::vec_save "vec_save"
and \ref fem_composite::oft_fem_comp_type::vec_load "vec_load" methods.

~~~~~~~~~{.F90}
USE fem_composite, ONLY: oft_fem_comp_type
USE lagrange_basis, ONLY: oft_lagrange
TYPE(oft_fem_comp_type) :: fem_3vec
!---Specify field count
fem_3vec%nfields = 3
ALLOCATE(fem_3vec%fields(fem_3vec%nfields))
ALLOCATE(fem_3vec%field_tags(fem_3vec%nfields))
!---Define field 1 (X-direction)
fem_3vec%fields(1)%fe => oft_lagrange
fem_3vec%field_tags(1) = "x"
!---Define field 2 (Y-direction)
fem_3vec%fields(2)%fe => oft_lagrange
fem_3vec%field_tags(2) = "y"
!---Define field 3 (Z-direction)
fem_3vec%fields(3)%fe => oft_lagrange
fem_3vec%field_tags(3) = "z"
~~~~~~~~~

\subsection doc_comp_fem_const_ml Construction from existing ML representations

If you intend to define and use a multi-level hierarchy for multi-grid this can be done automatically from existing
\ref fem_base::oft_ml_fem_type "ML FE" definitions. The \ref fem_composite::oft_ml_fem_comp_type "oft_ml_fem_comp_type"
class can be constructed in a similar way to the \ref fem_composite::oft_fem_comp_type "oft_fem_comp_type" class, but
using references to \ref fem_base::oft_ml_fem_type "ML FE" objects instead of single level definitions. Once defined
the \ref fem_composite::oft_fem_comp_type::setup "setup" method is called and the sub-levels are created automatically.
An example of this for the 3-field Lagrange vector is shown below.

~~~~~~~~~{.F90}
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE lagrange_basis, ONLY: ml_oft_lagrange
TYPE(oft_ml_fem_comp_type) :: ML_fem_3vec
!---Specify field count
ML_fem_3vec%nlevels = ml_oft_lagrange%nlevels
ML_fem_3vec%nfields = 3
ALLOCATE(ML_fem_3vec%fields(ML_fem_3vec%nfields))
ALLOCATE(ML_fem_3vec%field_tags(ML_fem_3vec%nfields))
!---Define field 1 (X-direction)
ML_fem_3vec%fields(1)%fe => ml_oft_lagrange
ML_fem_3vec%field_tags(1) = "x"
!---Define field 2 (Y-direction)
ML_fem_3vec%fields(2)%fe => ml_oft_lagrange
ML_fem_3vec%field_tags(2) = "y"
!---Define field 3 (Z-direction)
ML_fem_3vec%fields(3)%fe => ml_oft_lagrange
ML_fem_3vec%field_tags(3) = "z"
!---Setup sub-levels
CALL ML_fem_3vec%setup
~~~~~~~~~

\section doc_comp_fem_mconst Matrix construction

Most matrices for composite representations can be defined using the \ref fem_composite::oft_fem_comp_type::mat_create
"mat_create" method as with single field FE representations. For composite representations the
\ref fem_composite::oft_fem_comp_type::mat_create "mat_create" method supports simple masking of sub-matrices to eliminated
unnecessary storage when only certain sub-fields have non-zero interaction. The example below illustrates creating a matrix
for the 3-field Lagrange vector representation shown above, where sub-fields "y" and "z" do not interact with each other.

~~~~~~~~~{.F90}
USE oft_matrices, ONLY: oft_matrix
INTEGER(i4) :: mask(3,3)
CLASS(oft_matrix), POINTER :: mat
!---Setup graph mask
mask=1      ! Include all blocks by default
mask(2,3)=0 ! Skip block (2,3)
mask(3,2)=0 ! Skip block (3,2)
!---Create block matrix
CALL fem_3vec%mat_create(mat,mask)
~~~~~~~~~

\section doc_comp_fem_layout Storage Layout

Currently the code creates composite fields by stacking the sub-fields in memory. Although there are
instances where an interleaved layout would be desired this method is not amenable to cases where the
finite element representations are not homogeneous, ie when Lagrange and Conforming/Nedelec representations are
used. This is due to the fact that the elements do share a common node ordering so that there is no
way to consistently interleave elements.
