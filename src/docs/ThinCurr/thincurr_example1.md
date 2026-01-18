ThinCurr Example: Eigenstates of a square plate {#doc_thincurr_ex1}
==============

[TOC]

# Running an eigenvalue calculation

To run an eigenstate simulation (compute the characteristic L/R-decaying modes) you can issue the following command from within one of the run directories. Note, make sure `plot_run=F` in the `thincurr_eig_options` group.

    /path/to/oft/bin/thincurr_eig oft.in oft_in.xml

Contents of the input files are provided below and in the examples directory under `ThinCurr/plate`, along with the necessary mesh file `thincurr_ex-plate.h5`.

During the run you will see it report a bunch of general information and then it will report the first few eigenvalues from the calculation. These values should all be purely real (no imaginary component) and should be decreasing in value from first to last. For this case the largest values should be around 9.7 ms.

# Post processing

Once complete you can now generate XDMF files suitable for visualization of results using the [VisIt](https://visit-dav.github.io/visit-website/index.html) code. This is a two step process. First, rerun the code as above but with `plot_run=T` in the `thincurr_eig_options` group. Once complete, you need to run the `build_xdmf.py` script, which generates XDMF metadata files that tells VisIt how to read the data. This can be done using the following command

    python /path/to/oft/bin/build_xdmf.py

Next use VisIt to open the `surf_static.xmf` file, which will contain a series of vector fields named as `J_XX` that correspond to the current distributions of the various eigenstates. If you are running this example remotely and using VisIt locally you will need to copy the `mesh.*.h5`, `scalar_dump.*.h5`, `vector_dump.*.h5`, and `*.xmf` files to your local computer for visualization. The first eigenmode `J_01` should look like the figure below.

\image html thincurr_ex1-result.png "Resulting current distribution for the first eigenmode"

# Supporting information

## Input files

**General global settings** (`oft.in`)
\verbatim
&runtime_options
 debug=0
/

&mesh_options
 cad_type=0
/

&native_mesh_options
 filename="thincurr_ex-plate.h5"
/

\verbatim
&thincurr_eig_options
 direct=T
 plot_run=F
 neigs=4
/
\endverbatim

**XML global settings** (`oft_in.xml`)
```xml
<oft>
  <thincurr>
    <eta>1.257E-5</eta>
    <coil_set>
      <coil>0.5, 0.1</coil>
    </coil_set>
  </thincurr>
</oft>
```

## Mesh definition

**Cubit mesh script** (`thincurr_ex-plate.jou`)
```
reset
undo off

#{mesh_size=0.05}

# Create simple 1x1 plate
create surface rectangle width 1 height 1 zplane

# Generate mesh
set trimesher geometry sizing off
surface all scheme trimesh
surface all size {mesh_size}
mesh surface all

# Create mesh blocks
set duplicate block elements off
block 1 add surface 1

# Export grid
set large exodus file on
export Genesis  "thincurr_ex-plate.g" overwrite block 1
```

The file can then be converted to OFT's native mesh format using the `convert_cubit.py` script as

    python /path/to/OFT/bin/convert_cubit.py --in_file=thincurr_ex-plate.g

which will yield the converted file `thincurr_ex-plate.h5`.