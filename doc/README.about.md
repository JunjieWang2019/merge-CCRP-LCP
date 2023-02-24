General Information
===================

Reference software is being made available to provide a reference
implementation of a profile for geometry-based coding of solid dynamic
content for the G-PCC standard being developed by MPEG-3DGC (ISO/IEC SC29 WG7).

One of the main goals of the reference software is to provide a
basis upon which to conduct experiments in order to determine which coding
tools provide desired coding performance. It is not meant to be a
particularly efficient implementation of anything, and one may notice its
apparent unsuitability for a particular use. It should not be construed to
be a reflection of how complex a production-quality implementation of a
profile for geometry-based solid and dynamic content coding with a future
G-PCC standard would be.

This document aims to provide guidance on the usage of the reference
software. It is widely suspected to be incomplete and suggestions for
improvements are welcome. Such suggestions and general inquiries may be
sent to the general MPEG 3DGC email reflector at
<mpeg-3dgc@gti.ssr.upm.es> (registration required).

Releases notes
--------------

The GeS-TM (for **Ge**ometry-based **S**olid content **T**est **M**odel),
is being based on a subset of the G-PCC adoptions for coding dynamic solid
contents. The original implementation started from the experimental model
being studied in EE13.60, and is based on the `release-v20.0-rc1` of G-PCC
test model TMC13.

### `ges-tm-v1.0`

Tag `ges-tm-v1.0` will be released after additional cleanups of the source
code.

### `ges-tm-v1.0-rc1`

Tag `ges-tm-v1.0-rc1` adds all the latest adoptions integrated in the
release 21.0 of TMC13 and related/applicable to dense geometry content,
RAHT attributes and entropy coding.

- `trisoup/m61577: non-cubic nodes`:  
  Add possibility to use non cubic nodes on boundaries of the slices in
  order to avoid reconstruction artefacts and discontinuities between
  neighbouring slices.
  Node cubic nodes are determined based on a bounding box provided for the
  trisoup volume within the slice.

  It can be enabled separately for lowest coordinates in the slice and
  highest ones.

  Encoder decision is also added to determine the bounding box parameters.

- `trisoup/m61561: provide max num reconstructed points to encoder`:  
  Add possibility to provide separately to the encoder a maximum number of
  points reconstructed using trisoup before having to downsample the
  reconstruction.

  This allows to better control the number of input points per slice such
  that better coherence is obtained in the overall reconstructed content.

- `geom/m61583: reduce dynamic OBUF memory footprint`:  
  Reduce the amount of memory being used by dynamic OBUF.

- `raht/m61151: increase buffer precision to remove rounding ops`:  
  The precision is increased in the memory buffer used in RAHT to reconstruct
  attributes at each layer. Now, values are directly stored with fixed
  point precision, avoiding rounding operations during the transform.

- `entropy/m62220: simpler bypass bit with arithmetic coder`:  
  Simplify the operations for arithmetic coding in case of bypassed bits.

- `trisoup/m61982: rasterization with unitary sampling`:  
  Optimization to replace ray tracing by rasterization when no subsampling
  has to be used (trisoup_sampling_value_minus1 is equal to zero), and when
  integer precision rendering shall be used (trisoup_fine_ray_tracing_flag
  is equal to zero).

- `trisoup/m61982: voxelization of vertices on edges`:  
  Replace previous voxelization of trisoup vertices, which was not aligned
  on the edges, by a voxelization aligned on the edges.

  This voxelisation is not needed anymore when no sub-sampling is used:
  the vertices' points are already included in the rendered triangles.

- `trisoup/m61982: use thickness and adaptive halo with rasterization`:  
  Improve trisoup rendering by using thickness and adaptive halo with
  rasterization.


This release includes some extra optimizations.

- `ctc: consider level limit of 5 million points`:  
  For better efficiency, the GeS-TM currently considers hypothetical profile
  level limit of 5 million points per slice so that every frame is coded in
  an single slice.

- `trisoup: optimization of ray tracing`:  
  The ray tracing is performed along a single direction, and reconstructed
  points are restricted to fit inside the node (this could occur because
  of a lack of precision) by clipping the coordinates.

- `trisoup: unique points are selected at block level`:  
  Since points are restricted to be reconstructed locally to a node,
  the selection of the unique points can be (and is) performed at a block
  level rather than at a global point cloud level.

  Sort has ~N*log(N) complexity. Working locally reduces complexity and
  is better for memory access.

Some additionnal code simplifications and optimizations have been made,
as well as some changes in the motion search presets in order to improve
the execution speed:
- `trisoup: optimize determineTriSoupVertices function`,
- `trisoup: code reorganization for faster execution`,
- `trisoup/ctc: change motion search preset`.

Finally `ges-tm-v1.0-rc1` also includes additional code cleanups, minor bug
fixes and modifications to support the common test conditions which have been
agreed during the meeting.

### `ges-tm-v0.1`

Tag `ges-tm-v0.1` corresponds to the software being output in branch
`mpeg140/ee13.60/m61562_ETM_inter_dense` for the 140th MPEG meeting as
experimental model for dynamic dense content coding, including the latest
tools and improvement in `release-v20.0-rc1` adopted in G-PCC for dense
content, and part of earlier exploratory tools studied in EE13.2 for
inter-frame coding in `mpeg-pcc-em13`.


Bug reporting
-------------
Bugs should be reported on the issue tracker set up at
<http://mpegx.int-evry.fr/software/MPEG/PCC/TM/mpeg-pcc-ges-tm/issues>.

