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

### `ges-tm-v6.0`

Tag `ges-tm-v6.0` will be released after reported issues are resolved,
if any.

### `ges-tm-v6.0-rc2`

Tag `ges-tm-v6.0-rc2` includes a few bug fixes and a refactoring of the RAHT.

- `raht: fix issue #2 - wrong condition used to access a table`:  
    Encoder was crashing when using inter frame coding with encoder parameter
    --rahtLayerInterModeSelectionEnabled=0, because it was accessing a non
    allocated table element. The table is not supposed to be used.

    A condition enabling the use of the table was based on the wrong variable
    name. This commit fixes the condition.

- `raht: fix issue #4 - RDO could select intra layer mode on root node`:  
    Encoder could select intra layer mode on root node, leading to mismatch
    between encoder and decoder output in some situations.
    This commit fixes this issue by enforcing infinite RD cost for intra
    layer mode on first layer.

- `raht: major refactoring`:  
    This commit provides a major refactoring of the RAHT code:

    - The code is split between the encoder and decoder into two
      separate files:
      - move RAHT.cpp -> RAHTEncoder.cpp
      - add RAHTDecoder.cpp
    - Many optimization are introduced within the decoder
    - Some are reported into the encoder, but the work on the encoder will be
      continued after next meeting cycle

- `enc: fix issue #5 - do not use std::log2 for RDO`:  
    This commit removes dependencies on std::log2() for the encoder when
    performing Rate Distortion Optimizations (RDO). Instead, a fixed point
    approximation of log2 pcc::fpLog2 is introduced so that the results will
    not depend on the platform and its implementation dependent math library.

    This commit thus changes computations for RAHT RDO, as well as for motion
    search, and should solve issue #5.

- `misc: fix issue #3`:  
    Theres is a bug in gcc-8.1/8.2 (fixed in 8.3). This commit avoids the
    compiler complaining for no reason.

### `ges-tm-v6.0-rc1`

Tag `ges-tm-v6.0-rc1` includes adoptions made during the meeting #15 of the
WG7 held in Rennes (April 2024), and some additional fixes and cleanups.

- `raht: tidy - remove unnecessary operations`:  
    Performs some cleanups of the code for RAHT.

- `octree: refactor - restore a single function for enc/dec`:  
    Refactors octree's encoding and decoding functions.
    A single function is restored to support both octree and trisoup.

- `geom/attr: refactor - payload handling structure for single pass decoding`:  
    Refactors the decoding process to allow expected feature of
    being able to perform geometry and attributes decoding together in a
    single pass (with more local processing of the data).

- `geom/attr: tidy - dual motion parameters derivation and HLS`:  
    Provides general cleanups and fixes for the dual motion's
    parameters and HLS.

- `geom/attr: refactor - simplify parameters passing to enc/dec functions`:  
    Simplifies parameters passing to several Encoding/Decoding
    functions.

- `geom/attr: tidy - minor cleanups`:  
    Provides minor code cleanups.

- `geom/attr: refactor - template octree enc/dec with/without trisoup`:  
    Uses templatized versions for the octree encoding/decoding
    functions, to obtain better execution performances.

- `trisoup: m67002 - versatile quantization`:  
    Introduces versatile geometry quantization of TriSoup to provide
    better rate control with quantization parameter.

- `trisoup: misc - optimize OBUF memory use`:  
    Optimizes OBUF memory usage and management to improve
    execution speed.

- `enc/trisoup: m67002 - update motion presets`:  
    Improves presets for inter-frame motion search and
    prediction with trisoup.

- `trisoup: m67015 - support for non-powers of two node sizes`:  
    Provides finer granularity for quality control of the
    TriSoup model by introducing the support for TriSoup node sizes which
    are not powers of two.

- `trisoup: m67539 - improve trisoup surface reconstruction`:  
    Improves trisoup surface reconstruction around vertices.
    It reduces visual artefacts when the surface is passing through the
    corner using one vertex instead of multiple vertices.

- `trisoup: m68006 - simplify face vertex`:  
    Provides several cleanups and simplifications of the
    face vertex for TriSoup.

- `trisoup: m67544 - reduce OBUF memory usage`:  
    Refines the OBUF scheme of Trisoup vertex presence flag.

- `octree: m67540 - improve inter states contextualization`:  
    EE13.60-Test2a which improved inter-frame
    contextual modeling for occupancy coding.

- `octree: m67004 - improve inter prediction`:  
    Combination 1 of m67004 with m67540.

- `enc/raht: m67627 - RDOQ improvement`:  
    The reconstruction distortion is also taken into account
    within RDOQ for RAHT.

- `enc/raht: m67627 - fix RDO with lossless`:  
    Fixes the RDO issue when determining prediction mode under CW
    condition.

- `raht: m67627 - fix prediction mode signaling`:  
    Modifies the software to match with the spec when signalling
    the prediction mode.

- `raht: m67628 - remove redundant signalling and RDO calculations`:  
    Removes unnecessary signaling for the prediction mode of
    RAHT, and avoids corresponding RDO calculations.

- `raht: m67541 - cross chroma component prediction`:  
    Introduces cross chroma component prediction, for intra
    frame coding. The first chroma component is used to predict the second
    chroma component.

- `raht: m67542 - remove unused contexts`:  
    Removes unused contexts for prediction Mode coding.

- `raht: m67406 - refactor div/mod 3 RAHT level/layer`:  
    Provides an example of implementation to avoid using
    divisions by 3 and modulo by 3 when computing levels and layers in RAHT
    transform.

- `raht: m67406 - use fixed point for attributes`:  
    Refactors RAHT to use fixed point representation for the
    intermediate values of the attributes.

- `misc: m67406 - approximate division for FixedPoint representation`:  
    Introduces a new division approximation function to be used
    for divisions removal and uniformization with other existing division
    approximations.

- `raht: m67406 - share a single approximate division`:  
    Removes two LUTs and increases precision when computing
    weighted intra/inter prediction.

- `raht: m67406 - rounding inter integer haar predictor`:  
    Slightly improves the prediction by rounding the integer
    haar predictor.

- `tuc/raht: m67093 - RAHT per block`:  
    Adding possibility of performing RAHT transform by blocks as a technology
    under consideration.

- `tuc/geom/attr: m67808 - localized attributes`:  
    Add localized attributes processing as a technology under consideration.
    Localized attributes improve the design of the codec to let it be more
    friendly to being embedded in a device.

- `enc/trisoup: m67838 - fix slice/tile alignment`:  
    From discussions about m67838, try fixing slice/tile alignment in case
    multiple slices/tiles are used, to remove artefacts on slice boundaries.

- `enc/trisoup: m67601 - improve vertex determination`:  
    Introduces an encoding parameter to control the threshold
    for vertex determination. This threshold is put to 1 by default
    (instead of 0).

- `enc/trisoup: m67017 - new CTCs with 6 rate points`:  
    New Common Test Conditions for TriSoup.

- `enc/trisoup/raht: m67018 - new CTCs for colors with trisoup`:  
    New Common Test Conditions for colors with RAHT coded on top of
    a TriSoup geometry.

### `ges-tm-v5.0`

Tag `ges-tm-v5.0` has been released after some cleanups and improvements and
a few minor bug fixes have been integrated. The `README.md` file has also been
updated to advise building for Release, for runtime comparisons.

### `ges-tm-v5.0-rc1`

Tag `ges-tm-v5.0-rc1` includes adoptions made during the meeting #14 of the
WG7 held online (January 2024), and some additional fixes and cleanups.

- `tidy: motion - remove unused AttributeInterPredParams::motionVector`:  
    Removes `AttributeInterPredParams::motionVector` and
    `tmc3/attr_tools.h:MotionVector` which are not used.

- `geom/tidy: motion - remove unused variables and parameters`:  
    Removes `LPUnumInAxis` variables and `numLPUPerLine` parameter which are not
    used.

- `geom/tidy: octree - remove unused PCCOctree3Node::siblingOccupancy`:  
    Removes `PCCOctree3Node::siblingOccupancy` which is not used.

- `geom/tidy: octree - disable remaining IDCM code`:  
    - Comments out some IDCM related coding/decoding functions
    - Comments out `PCCOctree3Node::numSiblingsPlus1` and
      `PCCOctree3Node::numSiblingsMispredicted` which are only used by IDCM
    - Comments out IDCM coding contexts
    - Removes remaining angular mode related coding contexts

- `geom/trisoup: m65829 - Single pass geometry coding`:  
    - two octrees
    - one pass geo for octree-trisoup
    - context class modified to handle interleaved octree and TriSoup
      geometry coding

- `geom/trisoup: m65831 - centroid residual quantization`:  
    Improves centroid residual quantization.

- `attr/fix: m66895 - point cloud should not be offset`:  
    Point cloud and compensated point cloud should not be offset
    before and after attributes coding with current design of GeS-TM.

- `opti: m66238 - preliminary optimizations of NN search`:  
    This commit optimizes the nearest neighbour search used for motion
    search. It has been used as preliminary improvements in m66238 for
    EE13.60 Test 1.

- `attr: m66238 - dual motion field for color`:  
    - add dual motion field for color
    - add encoder options and hls flags for dual motion
    - update cfg files for CTCs (active with TriSoup)

- `attr: m66238 - motion compensated attributes projection`:  
    - add motion compensated attributes projection from reference frame to
      reconstructed geometry, to replace inter-frame attributes predictor,
    - use depth first approximate NN for motion compensated attributes
      projection,
    - add encoder options and hls flags for motion compensated attributes
      projection and enable it by default when dual motion is active (CTCs),
    - use depth first constructed MSOctree.

- `geom: m66238 - accelerated motion search`:  
    - add depth first NN search for motion search,
    - add depth first approximate NN search to further accelerate encoding,
      at the cost of a decrease in geometry coding performances,
    - add encoder option to enable approximate NN search for motion search.

- `raht: m66286 - RDO for inter mode selection`:  
    This is adopted EE13.60 - Test6, providing rate distorsion optimization
    of inter mode selection for RAHT (also adds signaling of the selected
    mode).

- `raht: fix - use safer approximation of tree depth`:  
    Method used for estimating raht tree depth in m66238 does not
    look safe: the distance between the first an last point in morton
    order is not representative of the actual boundaries of the point
    cloud. Also the function roundLog2 does not look correct.

    This is fixed by only considering the last point in morton order.

- `raht: fix - intra behavior shall not be changed`:  
    Complexity reduction introduced by prediction_threshold1 should be kept
    as it was not discussed to change the behavior of intra frame coding with
    m66286.

- `raht: m66274 - weighted average prediction`:  
    This is adopted EE13.60 - Test 5, after being merged with Test 6. It
    provides a prediction of raht coefficient made by a weighted average
    between inter and intra predictions.

- `raht: fix - out of bound votesum in divisionPredictionVotesLUT`:  
    `votesum` could be 12, when `getNeighborsMode()` is called and the parent
    node has no neighbors.

    This is solved by also considering the mode of the parent node itself
    when computing the values for vote(Inter/Intra[Layer])Weight.

- `misc: add helper functions to FixedPoint`:  
    - Add `+`, `-`, `*` and `/` operators to `FixedPoint`.
    - Add `static FixedPoint::fromVal(int64_t)` to build a `FixedPoint` value
      from its internal `int64_t` representation.

- `raht: m66152 - sample domain prediction`:  
    This is adopted EE13.60 - Test 2, after being merged with Test 5 and
    Test 6. It uses sample domain prediction within a raht layer to avoid
    computing a forward transform in decoding process.

- `raht/enc: m66308 - transform domain distorsion estimation`:  
    Compute distorsion estimation in transform domain to avoid an extra
    inverse transform.

- `raht: m66373 - fix prediction mode determination`:  
    This commit aligns the software with the text description of m62583
    which introduced inter RAHT in GeS-TM. It changes the prediction mode
    determination for RAHT.

### `ges-tm-v4.0`

Tag `ges-tm-v4.0` has been released on top of `ges-tm-v4.0-rc1` since no other
issue had been reported. Gitlab references have been update to reflect the new
address of mpeg gitlab server and the new organization of the projects.

### `ges-tm-v4.0-rc1`

Tag `ges-tm-v4.0-rc1` includes adoptions made during the meeting in
Hannover (October 2023), and some additional cleanups.

- `attr/tidy: m58315 - remove unused quantization variable`:  
    This is same as m58315 adoption, but quantization variables was
    not used at all in GeS-TM since predlift is not part of GeS-TM.

- `attr/raht: m65386 - fix overflow`:  
    Use of int64_t instead of int to avoid potential overflow with high
    bit depth attributes.

- `attr/raht: m65386 - fix attribute QP parameter restriction`:  
    Increases attributes QP range to match G-PCC specification text.

- `attr/raht: m65387 - update encoding QP table`:  
    Increases bitdepth resolution for encoding QP. It has been shown
    on TMC13 that it is far better for high bitdepth attributes, and
    slightly better for 8 bits attributes.

- `attr/raht: m65387 - extend number of coding contexts`:  
    Increase the number of contexts for coding attributes. It has
    been shown on TMC13 that it is much better for high bitdepth
    attributes, and slightly better for 8 bits attributes.

- `attr/raht: m65089 - Test 4.3, intra_mode_level=4`:  
    In RAHT, every frame is divided into binary layers for all-intra
    configuration.

- `attr/tidy: remove unused referencePointCloud`:  
    AttributeInterPredParams::referencePointCloud is not used anymore.
    Removes it.

- `geom/tidy: use const references instead of local copies`:  
    Replaces some local instances/copies by const references in Vec3 and
    PointType member functions.

- `geom/tidy: motion - remove unneeded copies of LPU window`:  
    Removes unnecessary local copies of LPU windows.

- `geom/tidy: motion - optimize memory usage`:  
    Builds indices array of points, before extracting points in well
    dimentionned arrays, to get more optimal use of the memory.

- `geom/tidy: motion - put non compensated nodes in compensatedPointCloud`:  
    Parts of refPointCloud that were not compensated are now copied
    into compensatedPointCloud.
  
    refPointCloud can be removed from trisoup (not necessary anymore).

- `geom/refactor: motion - use PCCPointSet3 instead of vectors`:  
    Uses PCCPointSet3 instead of vectors of coordinates in motion related
    algorithm (search and compensation) to be able to more easily compensate
    with attributes.

- `geom/motion: m64918 - octree based motion search`:  
    Uses an octree based motion search and compensation to remove LPU windows
    and there constrains.
    It allows to use wider motion search and compensation without increasing
    complexity.
  
    The motion window size is also removed from the gps as is no more needed
    for motion compensation in decoder side.
  
    Parameters for octree are modified to keep similar complexity

- `geom/trisoup: m65235 - FaceVertex for Trisoup`:  
    Adds possibility of using vertices on faces in trisoup coding.  
    Also updates CTCs.

- `geom/tidy: m64911 - remove HLS, dead code and sparse algo`:  
    Cleaning GeS-TM, as proposed in m64911.

- `geom/trisoup: m64912 - CTC changes and conditions on centroid activation`:  
    Centroid can now be used with at least 3 points, and it is activated for
    all rate points. CTCs' thikness parameters are also updated.

- `tidy: removal of unused variables and code`:  
    Provides attitional code cleanups:  
    - Removal of unused variables,
    - Removal of code in trisoup:  
      The code had been duplicated in a previous loop iterating on the
      nodes of a slice, for centroid determination with face vertices,
      and the original loopbecame unnecessary.

### `ges-tm-v3.0`

Tag `ges-tm-v3.0` has been released on top of `ges-tm-v3.0-rc1` since no other
issue had been reported.

### `ges-tm-v3.0-rc1`

Tag `ges-tm-v3.0-rc1` includes adoptions made during the meeting in
Geneva (July 2023), a few new features, cleanups and some refactoring.

- `misc: add parenthesis delimiter for arguments parsing`:  
    Arguments/arguments-list can now be delimited by parenthesis to allow
    nested sequencial types.

- `trisoup: m63661 - inter skip mode`:  
    In non-moving parts of a dynamic point cloud, it has been observed that
    the geometry inter prediction based on (zero) motion compensation may
    not be optimal due to lack of invariance of the compensation process.
  
    A new inter coding mode based on colocated vertices and nodes has been
    introduced to mimic a kind of skip mode for point cloud coding.

- `trisoup/ctc: m63660 - align slice BB to zero`:  
    Aligning slice origin to (0, 0, 0) coordinates by default.

- `trisoup: m63660 - hls for skip mode and grid alignment`:  
    Adds high level syntax to align trisoup slices according to trisoup node
    sizes. Also adds restrictions to enforce this alignment with skip mode.

- `refactor: m64005 - refactoring Attribute{Enc,Dec}oder`:  
    This is extracted from monolitic commit made for m64005  in EE repository.
    It refactors Attributes code to get a single templated implementation for
    both 1 dimention attribute luminance, and 3 dimention attributes color.

- `refactor: m64005 - refactoring of RAHT intra`:  
    This is extracted from monolitic commit made for m64005  in EE repository.
    It provides some refactoring to RAHT.

- `restore: divisionless RAHT prediction`:  
    Reverts normative change made during the refactoring of RAHT. Integer
    division was (re-)introduced for inter layer prediction, but had not been
    presented nor discussed with WG7 during the meeting.  
    The original behavior (with tabulated division approximation) has been
    restored.

- `restore: fixed point log2 estimation for RAHT RDOQ`:  
    Reverts non-normative change made during the refactoring of RAHT.
    Fixed point arithmetic for RAHT RDOQ was replaced by floating points
    operations, but had not been presented nor discussed with WG7 during the
    meeting.  
    The original behavior (with fixed point approximation) has been
    restored.

- `restore: original neighbours retrieval for intra RAHT`:  
    Reverts normative change made during the refactoring of RAHT.
    The refactoring of the neighbours retrieval for RAHT was modifying the
    behaviour of the inter layer prediction, but had not been presented nor
    discussed with WG7 during the meeting.
    The original code and behavior has been restored.

- `raht/inter: m64005 - inter RAHT`:  
    This is extracted from monolitic commit made for m64005  in EE repository.
    In provides the inter RAHT discussed with WG7 and adopted.  
  
    Integration notes:  
    - Some modifications affecting all-intra behaviour have been removed,
      as they were not discussed in WG7. It includes the coding of a mode
      for choosing the enabling or not of the inter layer prediction.
    - The syntax element `mode_level` has been moved as an inter-prediction
      dependent syntax element, since it should not affect all-intra.

- `attr/raht: m64218 - inference of prediction modes`:  
    Provides the inference of the prediction modes at certains levels of the
    RAHT decomposition.  
  
    Integration notes:  
    - The syntax element `upper_mode_level` has been moved and made
      dependent on RAHT `enable_inter_prediction`, since the mode
      should be infered only when inter prediction is enabled.

- `attr/raht: m64112 - fix rootLevel value`:  
    This is a fix on the derivation of the rootLevel value within RAHT.
    The rootLevel should be obtained by rounding to the upper integer value
    the division by 3 of the number of RAHT decomposition layers, to get the
    number of dyadic decomposition levels. Because latest one might be
    incomplete.

- `attr/raht: m64118 - integer Haar` 
    Introduces lossless extension to RAHT.  
  
    Integration notes:  
    - Added support for inter-RAHT with m64005:  
      - A simple rounding of the motion compensated predictor has
        been added,
      - The refactoring from m64005 has been modified to re-introduce
        templated use of RAHT/Haar kernels,
      - Kernel support for already normalized weights has been added.
    - Added cfg file for octree-raht-inter with lossless attributes

### `ges-tm-v2.0`

Tag `ges-tm-v2.0` was released on top of `ges-tm-v2.0-rc3` since no other
issue had been reported.

### `ges-tm-v2.0-rc3`

Tag `ges-tm-v2.0-rc3`,
- fixes an out of bound memory access, causing an assertion to occure in debug
  mode. This was unsucessfuly fixed in `ges-tm-v2.0-rc2`.

### `ges-tm-v2.0-rc2`

Tag `ges-tm-v2.0-rc2`,
- fixes uninitialized memory access to an encoding parameters;
- restores the output points being reordered in decoding order within octree
  encoder, which had been accidentaly removed during refactoring.

### `ges-tm-v2.0-rc1`

Tag `ges-tm-v2.0-rc1` includes a few adoptions made during the meeting in
Antalya (April 2023).

- `octree: m62531 - nodes processing in raster scan order`:  
    Use lexicographic order (a.k.a. raster scan order) in octree geometry
    nodes processing to get same node ordering as used in trisoup.

- `entropy: m62547 - probability bounds for dynamic OBUF coders`

- `raht/enc: m63002 - replace log2() with LUT in RDOQ`:  
    Current RDOQ code for RAHT uses the log2() function to estimate the
    rate.

    The proposal simplifies this by using a LUT with a maximum of 16
    entries.

- `trisoup/enc: m62527 - align slices to trisoup node size grid`

- `trisoup: m62526 - Optimal weighted centroid`

- `trisoup/enc: m62981 - centroid quantization offset`

- `trisoup/enc: m62982 - determine more precise centroid`


Adoption of raster-scan ordered nodes processing for octree provides an
unified processing order between octree and trisoup elements.
Tag `ges-tm-v2.0-rc1` includes important refactoring changes on trisoup and
optimizations that were made possible by this unified raster scan ordering,
as well as some code cleanups and simplifications.

The input contribution **m63611: \[GPCC\] Report on new architecture for GeS TM 
v2.0-rc1 and related integrations**
is intended to provide more details on theses modifications.

- `trisoup: m62531 - simplify with raster scan ordered nodes`:  
    simplify trisoup processing by reusing raster scan ordered nodes to
    build edges and vertex informations for trisoup.

- `trisoup: optimize findDominantAxis`

- `trisoup: unit sampling only, simplify accordingly`:  
    GeS-TM is intended to adress solid content.
    Unit sampling shall always be used for triangle projections.

    Encoding is simplified and code is optimized accordingly.

- `trisoup: refactor centroid operations to functions`:  

    Centroid operations are moved to dedicated functions:
    1. determination of centroid and dominant axis
    2. determination of normal vector and bounds
    3. determination of residual/drift
    4. determination of inter predictor for residual/drift

- `trisoup/tidy: cleaning in loop on triangles`

- `trisoup/tidy: better memory management for recPointCloud`:  
    Using preallocation for recPointCloud.

- `trisoup/tidy: better memory management for reconstructed voxels`:  
    Reserve memory for reconstructed voxels of current TriSoup node
    to improve memory management.

- `trisoup/tidy: clean and comment RasterScanTrisoupEdges`

- `trisoup: refactor vertex determination`:  
    Vertex determination is performed in same loop as neighbours and edges
    determination thanks to raster-scan order.

- `trisoup: optimize vertex encoding/decoding`

- `trisoup: refactor vertex inter prediction determination`:  
    Vertex determination for inter prediction is performed in same loop
    as neighbours, edges and intra vertex determination thanks to raster-scan
    order.

    Obsolete determineTrisoupVertices() function is removed.

- `trisoup: refactor vertices encoding/decoding`:  
    Vertex encoding/decoding is performed in same loop as neighbours, edges
    intra and inter vertex determination thanks to raster-scan order.

    Obsolete functions encodeTrisoupVertices() and decodeTrisoupVertices()
    are removed.

- `trisoup/tidy: localize some buffers`:  
    Some buffer are localized inside of function rather than being passed
    as parameters, since they are not used outside anymore.

- `trisoup: refactor triangles rendering`:  
    Put rendering of triangles of a node in specific function.
    Apply it in the same loop as the other trisoup operations,
    to get a single raster scan ordered nodes traversal.

    Remove obsolete function decodeTrisoupCommon().

- `trisoup/tidy: localize some arrays`:  
    Some array are now only used locally for raster scan traversal.

- `trisoup: optimize memory usage`:  
    Use queue instead of vector to get more localized memory, with
    edgePattern, xForedgeOfVertex, and TriSoupVerticesPred.

- `trisoup: optimize rendering/ray tracing`:  
  - use int64 values to represent in 1D 3 dimensional coordinates to
    accelerate duplicate points removal.
  - use one optimized function for each ray tracing directions.
  - optimize/simplify ray tracing.

- `trisoup: optimize duplicates removal and buffer allocation`

- `trisoup/tidy: various cleanups and optimizations`


On trisoup, all the divisions at decoder have been removed and replaced
by fixed point inverse multiplications (normative).

- `trisoup: remove divisions`:  
    remove all divisions at decoder.


Additional cleanups, simplifications and code optimization have been made,
mainly on inter-frame motion search and on octree coding.

- `tidy: minor cleanups`

- `geometry: optimize LPU inter search window generation`:  
  - avoid using unnecessary unordered_map -> decoder twice faster.
  - use BB for points inside the loop.
  - offset the points to slice origin during the loop (avoid multiple passes).

- `geometry: reduce size of LPU inter search window`

- `octree/tidy: remove IDCM variables in node + moved planar container for QTBT up`:  
  !! this may impact the code in case there are duplicated point removing turned on.

- `octree/tidy: optimizations and cleanups`

- `trisoup/tidy: more reasonable fifo size with trisoup`


Other changes includes the disabling of non-cubic nodes that is currently not
supported properly, and a bug fix for and issue that could happen outside of
CTCs.

- `ctc: disable non-cubic nodes`

- `trisoup: fix - determination of centroid predictor`:  
    When inter prediction is not used for a node, the centroid predictor
    could be wrongly estimated.


Documentation files for `ges-tm-v2.0-rc1` are also updated.

- `doc: update documentation files`


### `ges-tm-v1.0`

Tag `ges-tm-v1.0` adds a few fixes.

- `hls/entropy: add backward compatibility flag for m62220`:  
  Applies a fix for backward compatibility issue made in TMC13 for
  entropy coding related contribution m62220.

- `fix: read access outside of an array`:  
  The encoder could read outside of an array.
  Even if the value was not used, it could cause issues when debugging.
  This is fixed by accessing the array only when necessary, and so,
  when the index is inside of the array.

- `raht: fix RDOQ overflows`:  
  In rare cases overflows could occur during computation of Rate/Distortion,
  if distortion was really huge (ex: DC coefficient).
  This is fixed by limiting the computation to small coefficients only.


Additional cleanups were also made, mainly removing unused tools previously
commented out.

- `QTBT is not yet supported, avoid using it`:  
  QTBT is not yet supported by inter frame coding. When inter frame is
  used the tool cannot be activated.
- `removing global motion`
- `removing predictive geometry`
- `removing angular`
- `removing bytewise coder`
- `removing planar`  
  Note: some planar things are kept for QTBT handling.
- `removing predlift`


And the documentation files for `ges-tm-v1.0-rc1` have been added.

- `doc: update documentation files`


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
<https://git.mpeg.expert/MPEG/3dgh/g-pcc/software/tm/mpeg-pcc-ges-tm/-/issues>.

