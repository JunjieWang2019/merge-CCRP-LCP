/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2018, ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * Neither the name of the ISO/IEC nor the names of its contributors
 *   may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "geometry.h"

#include "OctreeNeighMap.h"
#include "geometry_octree.h"
#include "geometry_trisoup.h"
//#include "geometry_intra_pred.h"
#include "io_hls.h"
#include "tables.h"
#include "quantization.h"
#include "PCCTMC3Decoder.h"
#include "motionWip.h"
#include <unordered_map>

namespace pcc {

//============================================================================

class GeometryOctreeDecoder {
public:
  GeometryOctreeContexts& ctx;

  GeometryOctreeDecoder(
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    GeometryOctreeContexts& ctxMem,
    EntropyDecoder* arithmeticDecoder);

  GeometryOctreeDecoder(const GeometryOctreeDecoder&) = default;
  GeometryOctreeDecoder(GeometryOctreeDecoder&&) = default;
  //GeometryOctreeDecoder& operator=(const GeometryOctreeDecoder&) = default;
  //GeometryOctreeDecoder& operator=(GeometryOctreeDecoder&&) = default;

  // dynamic OBUF
  void clearMap() { ctx.clearMap(); };
  void resetMap(bool forTrisoup) { ctx.resetMap(forTrisoup); }

  int decodePositionLeafNumPoints();

  uint32_t decodeOccupancyFullNeihbourgs(
    const RasterScanContext::occupancy& occ,
    int planarMaskX,
    int planarMaskY,
    int planarMaskZ,
    int predOcc[8],
    int predOUnComp[8],
    bool isInter
  );

  /*Vec3<int32_t> decodePointPosition(
    const Vec3<int>& nodeSizeLog2, Vec3<int32_t>& deltaPlanar);*/

  /*void decodeOrdered2ptPrefix(
    Vec3<bool> directIdcm,
    Vec3<int>& nodeSizeLog2AfterUnordered,
    Vec3<int32_t> deltaUnordered[2]);*/

  bool decodeNodeQpOffsetsPresent();
  int decodeQpOffset();

  // local QU
  void decodeQU(localQU& qu, int baseQP);

  //bool decodeIsIdcm();

  /*template<class OutputIt>
  int decodeDirectPosition(
    bool geom_unique_points_flag,
    bool joint_2pt_idcm_enabled_flag,
    const Vec3<int>& nodeSizeLog2,
    const Vec3<int>& posQuantBitMasks,
    const PCCOctree3Node& node,
    const OctreeNodePlanar& planar,
    OutputIt outputPoints);*/

  const GeometryOctreeContexts& getCtx() const { return ctx; }

public:
  EntropyDecoder* _arithmeticDecoder;
};

//============================================================================

GeometryOctreeDecoder::GeometryOctreeDecoder(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder* arithmeticDecoder)
  : ctx(ctxtMem)
  , _arithmeticDecoder(arithmeticDecoder)
{
}

//============================================================================
// Decode the number of points in a leaf node of the octree.
int
GeometryOctreeDecoder::decodePositionLeafNumPoints()
{
  int val = _arithmeticDecoder->decode(ctx._ctxDupPointCntGt0);
  if (val)
    val += _arithmeticDecoder->decodeExpGolomb(0, ctx._ctxDupPointCntEgl);

  return val + 1;
}

//-------------------------------------------------------------------------
// decode node occupancy bits
//

static void (*pointer2FunctionContext[8])(
  OctreeNeighours&, int, int&, int&, bool&) = {
  makeGeometryAdvancedNeighPattern0,
  makeGeometryAdvancedNeighPattern1,
  makeGeometryAdvancedNeighPattern2,
  makeGeometryAdvancedNeighPattern3,
  makeGeometryAdvancedNeighPattern4,
  makeGeometryAdvancedNeighPattern5,
  makeGeometryAdvancedNeighPattern6,
  makeGeometryAdvancedNeighPattern7};

static const int LUTinitCoded0[27][6] = {
  {0, 0, 0, 0, 0, 0}, {4, 0, 2, 2, 2, 2}, {0, 4, 2, 2, 2, 2},
  {2, 2, 4, 0, 2, 2}, {4, 2, 4, 2, 3, 3}, {2, 4, 4, 2, 3, 3},
  {2, 2, 0, 4, 2, 2}, {4, 2, 2, 4, 3, 3}, {2, 4, 2, 4, 3, 3},
  {2, 2, 2, 2, 4, 0}, {4, 2, 3, 3, 4, 2}, {2, 4, 3, 3, 4, 2},
  {3, 3, 4, 2, 4, 2}, {4, 3, 4, 3, 4, 3}, {3, 4, 4, 3, 4, 3},
  {3, 3, 2, 4, 4, 2}, {4, 3, 3, 4, 4, 3}, {3, 4, 3, 4, 4, 3},
  {2, 2, 2, 2, 0, 4}, {4, 2, 3, 3, 2, 4}, {2, 4, 3, 3, 2, 4},
  {3, 3, 4, 2, 2, 4}, {4, 3, 4, 3, 3, 4}, {3, 4, 4, 3, 3, 4},
  {3, 3, 2, 4, 2, 4}, {4, 3, 3, 4, 3, 4}, {3, 4, 3, 4, 3, 4}};

//-------------------------------------------------------------------------
// decode node occupancy bits
//
uint32_t
GeometryOctreeDecoder::decodeOccupancyFullNeihbourgs(
  const RasterScanContext::occupancy& occ,
  int planarMaskX,
  int planarMaskY,
  int planarMaskZ,
  int predOcc[8],
  int predOccUnComp[8],
  bool isInter)
{
  // decode occupancy pattern
  uint32_t occupancy = 0;

  // 3 planars => single child and we know its position
  if (planarMaskX && planarMaskY && planarMaskZ) {
    uint32_t cnt = (planarMaskZ & 1);
    cnt |= (planarMaskY & 1) << 1;
    cnt |= (planarMaskX & 1) << 2;
    occupancy = 1 << cnt;
    return occupancy;
  }

  // neighbour empty and only one point => decode index, not pattern
  //------ Z occupancy decoding from here ----------------
  bool sure_planarityX = planarMaskX;
  bool sure_planarityY = planarMaskY;
  bool sure_planarityZ = planarMaskZ;

  bool isPred = predOcc[0] || predOcc[1] || predOcc[2] || predOcc[3]
    || predOcc[4] || predOcc[5] || predOcc[6] || predOcc[7]
    || predOccUnComp[0] || predOccUnComp[1] || predOccUnComp[2]
    || predOccUnComp[3] || predOccUnComp[4] || predOccUnComp[5]
    || predOccUnComp[6] || predOccUnComp[7];

  //int  MaskConfig = !planarMaskX ? 0 : planarMaskX == 15 ? 1 : 2;
  //MaskConfig += !planarMaskY ? 0 : planarMaskY == 51 ? 3 : 6;
  //MaskConfig += !planarMaskZ ? 0 : planarMaskZ == 85 ? 9 : 18;
  int MaskConfig = sure_planarityX * (1 + (planarMaskX != 0x0F));
  MaskConfig += sure_planarityY * 3 * (1 + (planarMaskY != 0x33));
  MaskConfig += sure_planarityZ * 9 * (1 + (planarMaskZ != 0x55));

  int coded0[6] = { 0, 0, 0, 0, 0, 0 };  // mask x0 x1 y0 y1 z0 z1
  if (MaskConfig) {
    memcpy(coded0, LUTinitCoded0[MaskConfig], 6 * sizeof(int));
  }

  OctreeNeighours octreeNeighours;
  prepareGeometryAdvancedNeighPattern(occ, octreeNeighours);

  // loop on occupancy bits from occupancy map
  uint32_t partialOccupancy = 0;
  int maskedOccupancy = planarMaskX | planarMaskY | planarMaskZ;
  for (int i = 0; i < 8; i++) {
    if ((maskedOccupancy >> i) & 1) {
      //  bit is 0 because masked by QTBT or planar
      partialOccupancy <<= 1;
      continue;
    }

    int mask0X = (0xf0 >> i) & 1;
    int mask0Y = 2 + ((0xcc >> i) & 1);
    int mask0Z = 4 + ((0xaa >> i) & 1);

    bool bitIsOne =
      (sure_planarityX && coded0[mask0X] >= 3)
      || (coded0[0] + coded0[1] >= 7)
      || (sure_planarityY && coded0[mask0Y] >= 3)
      || (coded0[2] + coded0[3] >= 7)
      || (sure_planarityZ && coded0[mask0Z] >= 3)
      || (coded0[4] + coded0[5] >= 7);

    if (bitIsOne) {
      occupancy += 1 << i;
      partialOccupancy <<= 1;
      partialOccupancy |= 1;
      continue;
    }

    // OBUF contexts
    int ctx1, ctx2;
    bool Sparse;
    (*pointer2FunctionContext[i])(
      octreeNeighours, occupancy, ctx1, ctx2, Sparse);

    bool isInter2 = isInter && isPred;

    int ctxTable = 0;
    if (!isInter2) // INTRA
      ctxTable = Sparse;
    else  //  INTER
      ctxTable = 2 + predOcc[i];

    if (isInter2)
      ctx1 = (ctx1 << 1) | predOccUnComp[i]>0;

    // decode
    int bit;
    if (Sparse) {
      bit = ctx._MapOccupancySparse[isInter2][i].decodeEvolve(
        _arithmeticDecoder, ctx._CtxMapDynamicOBUF[ctxTable], ctx2, ctx1,
        &ctx._OBUFleafNumber, ctx._BufferOBUFleaves);
    }
    else {
      bit = ctx._MapOccupancy[isInter2][i].decodeEvolve(
        _arithmeticDecoder, ctx._CtxMapDynamicOBUF[ctxTable], ctx2, ctx1,
        &ctx._OBUFleafNumber, ctx._BufferOBUFleaves);
    }

    // update partial occupancy of current node
    occupancy += bit << i;
    coded0[mask0X] += !bit;
    coded0[mask0Y] += !bit;
    coded0[mask0Z] += !bit;
    partialOccupancy <<= 1;
    partialOccupancy |= bit;
  }

  return occupancy;
}
//-------------------------------------------------------------------------

bool
GeometryOctreeDecoder::decodeNodeQpOffsetsPresent()
{
  return _arithmeticDecoder->decode();
}

//-------------------------------------------------------------------------

int
GeometryOctreeDecoder::decodeQpOffset()
{
  if (!_arithmeticDecoder->decode(ctx._ctxQpOffsetAbsGt0))
    return 0;

  int dqp = _arithmeticDecoder->decodeExpGolomb(0, ctx._ctxQpOffsetAbsEgl) + 1;
  int dqp_sign = _arithmeticDecoder->decode(ctx._ctxQpOffsetSign);
  return dqp_sign ? -dqp : dqp;
}

//-------------------------------------------------------------------------
void
GeometryOctreeDecoder::decodeQU(localQU& qu, int baseQP)
{
  qu.isBaseParameters = _arithmeticDecoder->decode(ctx._ctxQUflag);
  if (qu.isBaseParameters) {
    qu.localQP = baseQP;
    return;
  }

  // sign
  bool sign = _arithmeticDecoder->decode(ctx._ctxQUSign);

  // magnitude
  int mag = 1
    + _arithmeticDecoder->decodeExpGolomb(1, ctx._ctxQUQPpref, ctx._ctxQUQPsuf);
  qu.localQP = sign ? mag : -mag;
  qu.localQP += baseQP;
}

/*
//-------------------------------------------------------------------------
// Decode a position of a point in a given volume.
Vec3<int32_t>
GeometryOctreeDecoder::decodePointPosition(
  const Vec3<int>& nodeSizeLog2, Vec3<int32_t>& deltaPlanar)
{
  Vec3<int32_t> delta = deltaPlanar;
  for (int k = 0; k < 3; k++) {
    if (nodeSizeLog2[k] <= 0)
      continue;

    for (int i = nodeSizeLog2[k]; i > 0; i--) {
      delta[k] <<= 1;
      delta[k] |= _arithmeticDecoder->decode();
    }
  }

  return delta;
}
*/
/*
//-------------------------------------------------------------------------
// Decode part of the position of two unordred points  point in a given volume.
void
GeometryOctreeDecoder::decodeOrdered2ptPrefix(
  Vec3<bool> directIdcm, Vec3<int>& nodeSizeLog2, Vec3<int32_t> pointPrefix[2])
{
  if (nodeSizeLog2[0] >= 1 && directIdcm[0]) {
    int ctxIdx = 0;
    bool sameBit = true;
    while (nodeSizeLog2[0] && sameBit) {
      pointPrefix[0][0] <<= 1;
      pointPrefix[1][0] <<= 1;
      nodeSizeLog2[0]--;

      sameBit = _arithmeticDecoder->decode(ctx._ctxSameBitHighx[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      if (sameBit) {
        int bit = _arithmeticDecoder->decode();
        pointPrefix[0][0] |= bit;
        pointPrefix[1][0] |= bit;
      } else {
        pointPrefix[1][0] |= 1;
      }
    }
  }

  if (nodeSizeLog2[1] >= 1 && directIdcm[1]) {
    int ctxIdx = 0;
    bool sameBit = true;
    bool sameX = !directIdcm[0] || pointPrefix[0][0] == pointPrefix[1][0];

    while (nodeSizeLog2[1] && sameBit) {
      pointPrefix[0][1] <<= 1;
      pointPrefix[1][1] <<= 1;
      nodeSizeLog2[1]--;

      sameBit = _arithmeticDecoder->decode(ctx._ctxSameBitHighy[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      int bit = 0;
      if (!(sameX && !sameBit))
        bit = _arithmeticDecoder->decode();
      pointPrefix[0][1] |= bit;
      pointPrefix[1][1] |= sameBit ? bit : !bit;
    }
  }

  if (nodeSizeLog2[2] >= 1 && directIdcm[2]) {
    int ctxIdx = 0;
    bool sameBit = true;
    bool sameXy = (!directIdcm[0] || pointPrefix[0][0] == pointPrefix[1][0])
      && (!directIdcm[1] || pointPrefix[0][1] == pointPrefix[1][1]);

    while (nodeSizeLog2[2] && sameBit) {
      pointPrefix[0][2] <<= 1;
      pointPrefix[1][2] <<= 1;
      nodeSizeLog2[2]--;

      sameBit = _arithmeticDecoder->decode(ctx._ctxSameBitHighz[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      int bit = 0;
      if (!(sameXy && !sameBit))
        bit = _arithmeticDecoder->decode();
      pointPrefix[0][2] |= bit;
      pointPrefix[1][2] |= sameBit ? bit : !bit;
    }
  }
}
*/
/*
//-------------------------------------------------------------------------

bool
GeometryOctreeDecoder::decodeIsIdcm()
{
  return _arithmeticDecoder->decode(ctx._ctxBlockSkipTh);
}
*/
/*
//-------------------------------------------------------------------------
// Direct coding of position of points in node (early tree termination).
// Decoded points are written to @outputPoints
// Returns the number of points emitted.

template<class OutputIt>
int
GeometryOctreeDecoder::decodeDirectPosition(
  bool geom_unique_points_flag,
  bool joint_2pt_idcm_enabled_flag,
  const Vec3<int>& nodeSizeLog2,
  const Vec3<int>& posQuantBitMask,
  const PCCOctree3Node& node,
  const OctreeNodePlanar& planar,
  OutputIt outputPoints)
{
  int numPoints = 1;
  bool numPointsGt1 = _arithmeticDecoder->decode(ctx._ctxNumIdcmPointsGt1);
  numPoints += numPointsGt1;

  int numDuplicatePoints = 0;
  if (!geom_unique_points_flag && !numPointsGt1) {
    numDuplicatePoints = _arithmeticDecoder->decode(ctx._ctxDupPointCntGt0);
    if (numDuplicatePoints) {
      numDuplicatePoints += _arithmeticDecoder->decode(ctx._ctxDupPointCntGt1);
      if (numDuplicatePoints == 2)
        numDuplicatePoints +=
          _arithmeticDecoder->decodeExpGolomb(0, ctx._ctxDupPointCntEgl);
    }
  }

  // nodeSizeLog2Rem indicates the number of bits left to decode
  // the first bit may be inferred from the planar information
  Vec3<int32_t> deltaPlanar{0, 0, 0};
  Vec3<int> nodeSizeLog2Rem = nodeSizeLog2;
  for (int k = 0; k < 3; k++)
    if (nodeSizeLog2Rem[k] > 0 && (planar.planarMode & (1 << k))) {
      deltaPlanar[k] |= (planar.planePosBits & (1 << k) ? 1 : 0);
      nodeSizeLog2Rem[k]--;
    }

  // Indicates which components are directly coded
  Vec3<bool> directIdcm = true;

  // decode (ordred) two points
  Vec3<int32_t> deltaPos[2] = {deltaPlanar, deltaPlanar};
  if (numPoints == 2 && joint_2pt_idcm_enabled_flag)
    decodeOrdered2ptPrefix(directIdcm, nodeSizeLog2Rem, deltaPos);

  Vec3<int32_t> pos;
  for (int i = 0; i < numPoints; i++) {
    *(outputPoints++) = pos =
      decodePointPosition(nodeSizeLog2Rem, deltaPos[i]);
  }

  for (int i = 0; i < numDuplicatePoints; i++)
    *(outputPoints++) = pos;

  return numPoints + numDuplicatePoints;
}
*/
//-------------------------------------------------------------------------
// Helper to inverse quantise positions

Vec3<int32_t>
invQuantPosition(int qp, Vec3<uint32_t> quantMasks, const Vec3<int32_t>& pos)
{
  // pos represents the position within the coded tree as follows:
  //     |pppppqqqqqq|00
  //  - p = unquantised bit
  //  - q = quantised bit
  //  - 0 = bits that were not coded (MSBs of q)
  // The reconstruction is:
  //   |ppppp00qqqqqq| <- just prior to scaling
  //   |pppppssssssss| <  after scaling (s = scale(q))

  QuantizerGeom quantizer(qp);
  int shiftBits = QuantizerGeom::qpShift(qp);
  Vec3<int32_t> recon;
  for (int k = 0; k < 3; k++) {
    int lowPart = pos[k] & (quantMasks[k] >> shiftBits);
    int highPart = pos[k] ^ lowPart;
    int lowPartScaled = PCCClip(quantizer.scale(lowPart), 0, quantMasks[k]);
    recon[k] = (highPart << shiftBits) | lowPartScaled;
  }

  return recon;
}

//-------------------------------------------------------------------------

template <bool forTrisoup>
void
decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int skipLastLayers,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  const Vec3<int> minimum_position,
  InterPredParams& interPredParams,
  PCCTMC3Decoder3& rootDecoder,
  RasterScanTrisoupEdgesDecoder* rste)
{
  if (forTrisoup) {
    assert(rste);
  }

  PCCPointSet3& compensatedPointCloud = interPredParams.compensatedPointCloud;

  const bool isInter = gbh.interPredictionEnabledFlag;

  PCCPointSet3& predPointCloud = interPredParams.referencePointCloud;
  MSOctree& mSOctree = interPredParams.mSOctreeRef;
  auto& mvField = interPredParams.mvField;
  if (isInter) {
    int log2MinPUSize = ilog2(uint32_t(gps.motion.motion_min_pu_size - 1)) + 1;
    predPointCloud = refFrame->cloud;
    // for recoloring, need same depth as in encoder
    mSOctree
      = MSOctree(&predPointCloud, -gbh.geomBoxOrigin, std::min(2,log2MinPUSize));
  }


  // local QU
  if (forTrisoup) {
    if (gbh.qu_size_log2 > 0) {
      int quSize = (1 << gbh.qu_size_log2) * gbh.trisoupNodeSize(gps);
      std::cout << "QU size = " << quSize << "\n";
    }
    ctxtMem.quLastIndex = 0;
    ctxtMem.listOfQUs.clear();
  }

  // init main fifo
  //  -- worst case size is the last level containing every input poit
  //     and each point being isolated in the previous level.
  // NB: some trisoup configurations can generate fewer points than
  //     octree nodes.  Blindly trusting the number of points to guide
  //     the ringbuffer size is problematic.
  // todo(df): derive buffer size from level limit
  int ringBufferSize = gbh.footer.geom_num_points_minus1 + 1;
  if (gps.trisoup_enabled_flag && gbh.trisoupNodeSize(gps) > 1)
     ringBufferSize = std::max(1000,ringBufferSize  >> 2 * gbh.trisoupNodeSizeLog2(gps) - 2);

  std::vector<PCCOctree3Node> fifo;
  std::vector<PCCOctree3Node> fifoNext;
  fifo.reserve(ringBufferSize);
  fifoNext.reserve(ringBufferSize);

  RasterScanContext rsc(fifo);

  size_t processedPointCount = 0;
  std::vector<uint32_t> values;

  //// rotating mask used to enable idcm
  //uint32_t idcmEnableMaskInit = /*mkIdcmEnableMask(gps)*/ 0; //NOTE[FT] : set to 0 by construction

  Vec3<uint32_t> posQuantBitMasks = 0xffffffff;
  int idcmQp = 0;
  int sliceQp = gbh.sliceQp(gps);
  int nodeQpOffsetsSignalled = !gps.geom_scaling_enabled_flag;

  // generate the list of the node size for each level in the tree
  //  - starts with the smallest node and works up
  std::vector<Vec3<int>> lvlNodeSizeLog2{gbh.trisoupNodeSizeLog2(gps)};
  for (auto split : inReverse(gbh.tree_lvl_coded_axis_list)) {
    Vec3<int> splitStv = {!!(split & 4), !!(split & 2), !!(split & 1)};
    lvlNodeSizeLog2.push_back(lvlNodeSizeLog2.back() + splitStv);
  }
  std::reverse(lvlNodeSizeLog2.begin(), lvlNodeSizeLog2.end());

  // Derived parameter used by trisoup.
  gbh.maxRootNodeDimLog2 = lvlNodeSizeLog2[0].max();

  // the termination depth of the octree phase
  // NB: minNodeSizeLog2 is only non-zero for partial decoding (not trisoup)
  int maxDepth = lvlNodeSizeLog2.size() - skipLastLayers - 1;

  // append a dummy entry to the list so that depth+2 access is always valid
  lvlNodeSizeLog2.emplace_back(lvlNodeSizeLog2.back());

  // NB: this needs to be after the root node size is determined to
  //     allocate the planar buffer
  GeometryOctreeDecoder decoder(gps, gbh, ctxtMem, &arithmeticDecoder);

  // saved state for use with parallel bistream coding.
  // the saved state is restored at the start of each parallel octree level
  std::unique_ptr<GeometryOctreeContexts> savedState;

  // variables for local motion
  int log2MotionBlockSize = 0;

  // local motion prediction structure -> LPUs from predPointCloud
  if (isInter) {
    log2MotionBlockSize =
      ilog2(uint32_t(gps.motion.motion_block_size - 1)) + 1;
    if (gbh.maxRootNodeDimLog2 < log2MotionBlockSize) { // LPU is bigger than root note, must adjust if possible
      int log2MotionBlockSizeMin =
        ilog2(uint32_t(gps.motion.motion_min_pu_size - 1)) + 1;
      if (log2MotionBlockSizeMin <= gbh.maxRootNodeDimLog2)
        log2MotionBlockSize = gbh.maxRootNodeDimLog2;
    }
  }

  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(0);
  node00.pos = int32_t(0);
  node00.predStart = uint32_t(0);
  node00.predEnd = isInter ? predPointCloud.getPointCount() : uint32_t(0);
  node00.mSONodeIdx = isInter ? 0 : -1;
  //node00.numSiblingsMispredicted = 0;
  //node00.numSiblingsPlus1 = 8;
  node00.qp = 0;
  //node00.idcmEligible = 0; NOTE[FT]: idcmEligible is already set to false at construction
  //node00.isDirectMode = false;

  // local motion
  node00.hasMotion = 0;
  node00.isCompensated = 0;

  // local QU
  node00.quIndex = -1;

  if (!(isInter && gps.gof_geom_entropy_continuation_enabled_flag) && !gbh.entropy_continuation_flag) {
    decoder.clearMap();
    decoder.resetMap(forTrisoup);
  }

  const int S = gbh.trisoupNodeSize(gps);
  const int S2 = gbh.trisoupNodeSizeLog2(gps);
  const int factorS = (1 << S2) - S;
  const bool flagNonPow2 = factorS != 0;
  // localized attributes point indexes
  std::vector<std::vector<point_t> > laPoints;
  int localSlabIdx = 0;
  int xStartLocalSlab = 0;
  int slabThickness;
  bool useLocalAttr = sps.localized_attributes_enabled_flag && !forTrisoup;
  PCCPointSet3 localPointCloud;
  if (useLocalAttr) {
    slabThickness = sps.localized_attributes_slab_thickness_minus1 + 1;
    int numSlabs =
      ((1 << gbh.maxRootNodeDimLog2) + sps.localized_attributes_slab_thickness_minus1)
        / (sps.localized_attributes_slab_thickness_minus1 + 1);
    laPoints.resize(numSlabs);
    localPointCloud.addRemoveAttributes(pointCloud);
    localPointCloud.reserve(ringBufferSize);
  }

  int lastPos0 = INT_MIN; // used to detect slice change and call for TriSoup

  for (int depth = 0; depth < maxDepth; depth++) {
    // setup at the start of each level
    auto fifoCurrLvlEnd = fifo.end();
    int numNodesNextLvl = 0;

    // derive per-level node size related parameters
    auto nodeSizeLog2 = lvlNodeSizeLog2[depth];
    auto childSizeLog2 = lvlNodeSizeLog2[depth + 1];
    //// represents the largest dimension of the current node
    //int nodeMaxDimLog2 = nodeSizeLog2.max();

    const auto nodeSize = S * (1 << nodeSizeLog2[0] - S2);

    const auto pointSortMask = qtBtChildSize(nodeSizeLog2, childSizeLog2);

    // if one dimension is not split, atlasShift[k] = 0
    int codedAxesPrevLvl = depth ? gbh.tree_lvl_coded_axis_list[depth - 1] : 7;
    int codedAxesCurLvl = gbh.tree_lvl_coded_axis_list[depth];

    // Determine if this is the level where node QPs are sent
    bool nodeQpOffsetsPresent =
      !nodeQpOffsetsSignalled && decoder.decodeNodeQpOffsetsPresent();

    // record the node size when quantisation is signalled -- all subsequnt
    // coded occupancy bits are quantised
    // after the qp offset, idcm nodes do not receive special treatment
    if (nodeQpOffsetsPresent) {
      nodeQpOffsetsSignalled = true;
      idcmQp = 0;
      posQuantBitMasks = Vec3<uint32_t>((1 << nodeSizeLog2) - 1);
    }

    // Idcm quantisation applies to child nodes before per node qps
    if (!nodeQpOffsetsSignalled) {
      auto quantNodeSizeLog2 = nodeSizeLog2;

      for (int k = 0; k < 3; k++)
        quantNodeSizeLog2[k] = std::max(0, quantNodeSizeLog2[k]);

      // limit the idcmQp such that it cannot overquantise the node
      auto minNs = quantNodeSizeLog2.min();
      idcmQp = gps.geom_base_qp + gps.geom_idcm_qp_offset;
      idcmQp <<= gps.geom_qp_multiplier_log2;
      idcmQp = std::min(idcmQp, minNs * 8);

      posQuantBitMasks = Vec3<uint32_t>((1 << quantNodeSizeLog2) - 1);
    }

    // save context state for parallel coding
    if (depth == maxDepth - 1 - gbh.geom_stream_cnt_minus1)
      if (gbh.geom_stream_cnt_minus1)
        savedState.reset(new GeometryOctreeContexts(decoder.ctx));

    // a new entropy stream starts one level after the context state is saved.
    // restore the saved state and flush the arithmetic decoder
    if (depth > maxDepth - 1 - gbh.geom_stream_cnt_minus1) {
      decoder.ctx = *savedState;
      arithmeticDecoder.flushAndRestart();
    }

    int currPUIdx = 0;
    if (isInter) {
      if (nodeSizeLog2[0] == log2MotionBlockSize) {
        // TODO: allocate motion field according to number of nodes
        uint32_t numNodes = fifo.end() - fifo.begin();
        mvField.puNodes.reserve(numNodes * 8); // arbitrary value.
        mvField.mvPool.reserve(numNodes * 8); // arbitrary value.
        // allocate all the PU roots
        mvField.numRoots = numNodes;
        mvField.puNodes.resize(numNodes);
      }
    }

    //// reset the idcm eligibility mask at the start of each level to
    //// support multiple streams
    //auto idcmEnableMask = /*rotateRight(idcmEnableMaskInit, depth)*/ 0; //NOTE[FT]: still 0

    rsc.initializeNextDepth();


    // process all nodes within a single level
    bool isLastDepth = depth == maxDepth - 1;

    if (forTrisoup && isLastDepth) {
      rste->leaves.reserve(8 * fifo.size());
    }

    // planar mode as a container for QTBT at depth level
    OctreeNodePlanar planar;
    int codedAxesCurNode = codedAxesCurLvl;
    int planarMask[3] = { 0, 0, 0 };
    maskPlanar(planar, planarMask, codedAxesCurNode);

    IterOneLevelSubnodesRSO(
      fifo.begin(),
      fifoCurrLvlEnd,
      [&](const decltype(fifo.begin())& itNode) -> const point_t& {
        return itNode->pos;
      },
      [&](const decltype(fifo.begin())& fifoCurrNode, int childIdx) {

      if (childIdx & 1)
        return;

      const int tubeIndex = (childIdx >> 1) & 1;
      const int nodeSliceIndex = (childIdx >> 2) & 1;

      PCCOctree3Node& node0 = *fifoCurrNode;

      if (nodeQpOffsetsPresent && !tubeIndex && !nodeSliceIndex) {
        node0.qp = sliceQp;
        node0.qp += decoder.decodeQpOffset() << gps.geom_qp_multiplier_log2;
      }

      int shiftBits = QuantizerGeom::qpShift(node0.qp);
      auto effectiveNodeSizeLog2 = nodeSizeLog2 - shiftBits;
      auto effectiveChildSizeLog2 = childSizeLog2 - shiftBits;

      if (isLeafNode(effectiveNodeSizeLog2))
        return;

      std::array<int32_t, 8> predUnCompCounts = {};
      if(!tubeIndex && !nodeSliceIndex) {
        // decode local motion PU tree
        if (isInter) {

          if (nodeSizeLog2[0] == log2MotionBlockSize) {
            node0.mvFieldNodeIdx = currPUIdx++;
            node0.hasMotion = true;
          }

          // decode LPU/PU/MV
          if (node0.hasMotion && !node0.isCompensated) {
            decode_splitPU_MV_MC<false>(mSOctree,
              &node0, mvField, node0.mvFieldNodeIdx,
              gps.motion, nodeSize,
              &arithmeticDecoder, &compensatedPointCloud,
              flagNonPow2, S, S2);
          }
        }

        // ...for local motion
        if (isInter) {
          if (node0.isCompensated) {
            countingSort(
              PCCPointSet3::iterator(&compensatedPointCloud, node0.predStart),  // Need to update the predStar
              PCCPointSet3::iterator(&compensatedPointCloud, node0.predEnd),
              node0.predCounts, [=](const PCCPointSet3::Proxy& proxy) {
                const auto& point = *proxy;
                return !!(int(point[2]) & pointSortMask[2])
                  | (!!(int(point[1]) & pointSortMask[1]) << 1)
                  | (!!(int(point[0]) & pointSortMask[0]) << 2);
            });
          }
          else {
            if (depth < mSOctree.depth && node0.mSONodeIdx >= 0) {
              const auto& msoNode = mSOctree.nodes[node0.mSONodeIdx];
              for (int i = 0; i < 8; ++i) {
                uint32_t msoChildIdx = msoNode.child[i];
                if (msoChildIdx) {
                  const auto& msoChild = mSOctree.nodes[msoChildIdx];
                  node0.predCounts[i] = msoChild.end - msoChild.start;
                }
              }
            }
          }

          if (depth < mSOctree.depth && node0.mSONodeIdx >= 0) {
            const auto & msoNode = mSOctree.nodes[node0.mSONodeIdx];
            for (int i = 0; i < 8; ++i) {
              uint32_t msoChildIdx = msoNode.child[i];
              if (msoChildIdx) {
                const auto & msoChild = mSOctree.nodes[msoChildIdx];
                predUnCompCounts[i] = msoChild.end - msoChild.start;
              }
            }
          }
        }
        node0.predPointsStartIdx = node0.predStart;
      }

      // inter information
      // generate the bitmap of child occupancy and count
      // the number of occupied children in node0.
      //int predOccupancy = 0;

      // occupancy inter predictor
      int predOccupancy[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      int predOccupancyUnComp[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      const int L = 1 << nodeSizeLog2[0];

      // TODO avoid computing it at each pass?
      for (int i = 0; i < 8; i++) {
        //predOccupancy |= (node0.predCounts[i] > 0) << i;
        predOccupancy[i] += node0.predCounts[i] > std::max(0, L >> 8);
        predOccupancy[i] += node0.predCounts[i] > std::max(2, L >> 4);
        predOccupancy[i] += node0.predCounts[i] >
          (forTrisoup ? std::max(16, L) : std::max(8, L >> 2));
      }

      for (int i = 0; i < 8; i++) {
        predOccupancyUnComp[i] += predUnCompCounts[i] > std::max(0, L >> 8);
        predOccupancyUnComp[i] += predUnCompCounts[i] > std::max(2, L >> 4);
        predOccupancyUnComp[i] += predUnCompCounts[i] > std::max(8, L >> 2);
      }

      //bool occupancyIsPredictable =
      //  predOccupancy && node0.numSiblingsMispredicted <= 5;
      // The predictor may be cleared for the purpose of context
      // selection if the prediction is unlikely to be good.
      // NB: any other tests should use the original prediction.
      //int predOccupancyReal = predOccupancy;

      if (nodeQpOffsetsPresent) {
        node0.qp = sliceQp;
        node0.qp += decoder.decodeQpOffset() << gps.geom_qp_multiplier_log2;
      }

      // local QU
      if (forTrisoup && !tubeIndex && !nodeSliceIndex) {
        if (gbh.qu_size_log2 && nodeSizeLog2[0] == S2 + gbh.qu_size_log2) {

          ctxtMem.listOfQUs.emplace_back();
          localQU& qu = ctxtMem.listOfQUs.back();
          decoder.decodeQU(qu, gbh.trisoup_QP);
          node0.quIndex = ctxtMem.quLastIndex++;
        }
      }

      // make quantisation work with qtbt and planar.
      /*auto codedAxesCurNode = codedAxesCurLvl;
      if (shiftBits != 0) {
        for (int k = 0; k < 3; k++) {
          if (effectiveChildSizeLog2[k] < 0)
            codedAxesCurNode &= ~(4 >> k);
        }
      }*/

      GeometryNeighPattern gnp{};
      // The position of the node in the parent's occupancy map
      int posInParent = 0;
      posInParent |= (node0.pos[0] & 1) << 2;
      posInParent |= (node0.pos[1] & 1) << 1;
      posInParent |= (node0.pos[2] & 1) << 0;
      posInParent &= codedAxesPrevLvl;

      // At the scaling depth, it is possible for a node that has previously
      // been marked as being eligible for idcm to be fully quantised due
      // to the choice of QP.  There is therefore nothing to code with idcm.
      /*if (isLeafNode(effectiveNodeSizeLog2))
        node0.idcmEligible = false;*/ //NOTE[FT]: always false

      // N.B: contextualOccupancy is only valid during first pass on the node
      RasterScanContext::occupancy contextualOccupancy;
      //OctreeNodePlanar planar;
      uint8_t occupancy = 1;
      if (!isLeafNode(effectiveNodeSizeLog2) && !tubeIndex && !nodeSliceIndex) {
        // update contexts
        rsc.nextNode(&*fifoCurrNode, contextualOccupancy);
        node0.neighPattern = gnp.neighPattern = contextualOccupancy.neighPattern;
        gnp.adjNeighOcc[0] = contextualOccupancy.childOccupancyContext[4];
        gnp.adjNeighOcc[1] = contextualOccupancy.childOccupancyContext[10];
        gnp.adjNeighOcc[2] = contextualOccupancy.childOccupancyContext[12];
      }

      //if (node0.idcmEligible && !tubeIndex && !nodeSliceIndex) {
      //  if (/*!gps.geom_planar_disabled_idcm_angular_flag*/ true) //NOTE[FT]: FORCING geom_planar_disabled_idcm_angular_flag to false
      //    node0.isDirectMode = decoder.decodeIsIdcm();
      //  if (node0.isDirectMode) {
      //    auto idcmSize = effectiveNodeSizeLog2;
      //    if (idcmQp) {
      //      node0.qp = idcmQp;
      //      idcmSize = nodeSizeLog2 - QuantizerGeom::qpShift(idcmQp);
      //    }

      //    int numPoints = decoder.decodeDirectPosition(
      //      gps.geom_unique_points_flag, gps.joint_2pt_idcm_enabled_flag,
      //      idcmSize, posQuantBitMasks,
      //      node0, planar,
      //      &pointCloud[processedPointCount]);

      //    for (int j = 0; j < numPoints; j++) {
      //      auto& point = pointCloud[processedPointCount++];
      //      int childIndex = 0;
      //      for (int k = 0; k < 3; k++)
      //        if (idcmSize[k] > 0 && effectiveNodeSizeLog2[k] != effectiveChildSizeLog2[k])
      //          childIndex += (point[k] >> idcmSize[k] - 1) << (2 - k);
      //      node0.childOccupancy |= 1 << childIndex;
      //      for (int k = 0; k < 3; k++)
      //        point[k] += rotateLeft(node0.pos[k], idcmSize[k]);

      //      point = invQuantPosition(node0.qp, posQuantBitMasks, point);
      //    }

      //    // NB: no further siblings to decode by definition of IDCM
      //    if (gps.inferred_direct_coding_mode <= 1)
      //      assert(node0.numSiblingsPlus1 == 1);

      //    continue;
      //  }
      //}
      //if (node0.isDirectMode)
      //  continue;

      if (!isLeafNode(effectiveNodeSizeLog2) && !tubeIndex && !nodeSliceIndex) {
        // planar mode for current node
        // mask to be used for the occupancy coding
        // (bit =1 => occupancy bit not coded due to not belonging to the plane)
        //int planarMask[3] = {0, 0, 0};
        //maskPlanar(planar, planarMask, codedAxesCurNode);

        // decode child occupancy map
        node0.childOccupancy = decoder.decodeOccupancyFullNeihbourgs(
          contextualOccupancy, planarMask[0], planarMask[1], planarMask[2],
          predOccupancy, predOccupancyUnComp, isInter);
      }
      occupancy = node0.childOccupancy;
      assert(occupancy > 0);

      // population count of occupancy for IDCM
      int numOccupied = popcnt(occupancy);

      if(!tubeIndex && !nodeSliceIndex && isInter) {
        // create motion field child nodes if needed
        if (node0.hasMotion && !node0.isCompensated) {
          auto& puNode = mvField.puNodes[node0.mvFieldNodeIdx];
          assert(occupancy);
          puNode._childsMask = occupancy;
          puNode._firstChildIdx = mvField.puNodes.size();
          for (int puChildIdx = 0; puChildIdx < numOccupied; ++puChildIdx) {
            mvField.puNodes.emplace_back();
          }
        }
      }

      //int predFailureCount = popcnt(uint8_t(occupancy ^ predOccupancyReal));
      if (tubeIndex && nodeSliceIndex && isLeafNode(effectiveChildSizeLog2)) {
        int slabIdx;
        if (useLocalAttr) {
          int nodeposX = node0.pos[0] << !!(codedAxesCurLvl & 4) + effectiveChildSizeLog2[0];
          if (isLastDepth) {
            // rendering attributes of a finished Slab
            while (nodeposX >= xStartLocalSlab + slabThickness) {
              int numPoints = laPoints[localSlabIdx].size();
              localPointCloud.resize(numPoints);
              for (int i = 0; i < numPoints; ++i)
                localPointCloud[i] = laPoints[localSlabIdx][i];
              laPoints[localSlabIdx] = std::vector<point_t>(); // just release memory
              if (numPoints) {
                rootDecoder.processNextSlabAttributes(localPointCloud, xStartLocalSlab, false);
                pointCloud.setFromPartition(localPointCloud, 0, numPoints, processedPointCount);
                localPointCloud.clear();
                processedPointCount += numPoints;
              }
              // keep slabs aligned on a regular grid
              xStartLocalSlab += slabThickness;
              ++localSlabIdx;
            }
            slabIdx = localSlabIdx;
          } else {
            // TODO: better
            slabIdx = nodeposX / slabThickness;
            // TODO: shall we and how to handle case where quantization is bigger than thickness
          }
        }
        // nodeSizeLog2 > 1: for each child:
        //  - determine elegibility for IDCM
        //  - directly decode point positions if IDCM allowed and selected
        //  - otherwise, insert split children into fifo while updating neighbour state
        for (int i = 0; i < 8; i++) {
          // TODO: Answer: do we want to also have points added
          // (mostly if IDCM is used) in raster scan order ?
          uint32_t mask = 1 << i;
          if (!(occupancy & mask)) {
            // child is empty: skip
            continue;
          }

          int x = !!(i & 4);
          int y = !!(i & 2);
          int z = !!(i & 1);

          // point counts for leaf nodes are coded immediately upon
          // encountering the leaf node.
          int numPoints = 1;

          if (!gps.geom_unique_points_flag) {
            numPoints = decoder.decodePositionLeafNumPoints();
          }

          // the final bits from the leaf:
          Vec3<int32_t> point{
            (node0.pos[0] << !!(codedAxesCurLvl & 4)) + x,
            (node0.pos[1] << !!(codedAxesCurLvl & 2)) + y,
            (node0.pos[2] << !!(codedAxesCurLvl & 1)) + z};

          // remove any padding bits that were not coded
          for (int k = 0; k < 3; k++)
            point[k] = rotateLeft(point[k], effectiveChildSizeLog2[k]);

          point = invQuantPosition(node0.qp, posQuantBitMasks, point);

          if (useLocalAttr) {
            int idxLaPoint = laPoints[slabIdx].size();
            laPoints[slabIdx].resize(idxLaPoint + numPoints);
            for (int i = 0; i < numPoints; ++i)
              laPoints[slabIdx][idxLaPoint++] = point;
          } else {
            for (int i = 0; i < numPoints; ++i)
              pointCloud[processedPointCount++] = point;
          }
        }

        // do not recurse into leaf nodes
        return;
      }



      if (!isLeafNode(effectiveChildSizeLog2)) {
        // push child nodes to fifo
        for (int i = 0; i < 2; ++i) {
          int childIndex = (nodeSliceIndex << 2) + (tubeIndex << 1) + i;
          uint32_t mask = 1 << childIndex;
          bool occupiedChild = occupancy & mask;
          if (!occupiedChild) {
            // child is empty: skip
            node0.predPointsStartIdx += node0.predCounts[childIndex];
          }
          else {
            // create & enqueue new child.
            fifoNext.emplace_back();
            auto& child = fifoNext.back();

            child.qp = node0.qp;

            int x = nodeSliceIndex;
            int y = tubeIndex;
            int z = i;

            // only shift position if an occupancy bit was coded for the axis
            child.pos[0] = (node0.pos[0] << !!(codedAxesCurLvl & 4)) + x;
            child.pos[1] = (node0.pos[1] << !!(codedAxesCurLvl & 2)) + y;
            child.pos[2] = (node0.pos[2] << !!(codedAxesCurLvl & 1)) + z;

            //child.numSiblingsPlus1 = numOccupied;
            //child.isDirectMode = false;

            child.predStart = node0.predPointsStartIdx;
            node0.predPointsStartIdx += node0.predCounts[childIndex];
            child.predEnd = node0.predPointsStartIdx;
            //child.numSiblingsMispredicted = predFailureCount;
            if (node0.mSONodeIdx >= 0) {
              child.mSONodeIdx = mSOctree.nodes[node0.mSONodeIdx].child[childIndex];
            }

            //local motion PU inheritance
            child.hasMotion = node0.hasMotion;
            child.isCompensated = node0.isCompensated;

            child.quIndex = node0.quIndex;

            if (node0.hasMotion && !node0.isCompensated) {
              auto& puNode = mvField.puNodes[node0.mvFieldNodeIdx];
              assert(puNode._childsMask);
              assert(puNode._firstChildIdx != MVField::kNotSetMVIdx);
              child.mvFieldNodeIdx = puNode._firstChildIdx;
              for (int j = 0; j < childIndex; ++j)
                if (puNode._childsMask & (1 << j))
                  ++child.mvFieldNodeIdx;
            }


            //if (isInter && /*!gps.geom_angular_mode_enabled_flag*/ true)
            //  child.idcmEligible = isDirectModeEligible_Inter(
            //    /*gps.inferred_direct_coding_mode*/ 0, nodeMaxDimLog2, node0.neighPattern,
            //    node0, child, occupancyIsPredictable);
            //else
            //  child.idcmEligible = isDirectModeEligible(
            //    /*gps.inferred_direct_coding_mode*/ 0, nodeMaxDimLog2, node0.neighPattern,
            //    node0, child, occupancyIsPredictable,
            //    gps.geom_angular_mode_enabled_flag);

#ifdef DONE
            //if (child.idcmEligible) {
            //  child.idcmEligible &= idcmEnableMask & 1; //NOTE[FT] : stays at 0, whatever the childidcmEligible status
            //  child.idcmEligible &= idcmEnableMask & 1;
            //  idcmEnableMask = rotateRight(idcmEnableMask, 1);
            //}
#endif

            if (forTrisoup) {
              if (isLastDepth) {
                rste->leaves.emplace_back(child);
                rste->leaves.back().pos *= S;

                if (flagNonPow2) {
                  int maskS = (1 << S2) - 1;
                  for (int np = child.predStart; np < child.predEnd; np++) {
                    compensatedPointCloud[np][0] = ((compensatedPointCloud[np][0] >> S2) * S) + (compensatedPointCloud[np][0] & maskS);
                    compensatedPointCloud[np][1] = ((compensatedPointCloud[np][1] >> S2) * S) + (compensatedPointCloud[np][1] & maskS);
                    compensatedPointCloud[np][2] = ((compensatedPointCloud[np][2] >> S2) * S) + (compensatedPointCloud[np][2] & maskS);
                  }
                }

                if (lastPos0 != INT_MIN && child.pos[0] != lastPos0)
                  rste->callTriSoupSlice(false); // TriSoup unpile slices (not final = false)
                lastPos0 = child.pos[0];
              } else if (0) {
                /* callback local attributes here */
                // This is not OK
                // we will have to handle local QP...
                // some work to be done here...
              }
            }

            numNodesNextLvl++;
          }
        }
      }
    });

    fifo.resize(0);
    fifo.swap(fifoNext);

    // Check that one level hasn't produced too many nodes
    // todo(df): this check is too weak to spot overflowing the fifo
    assert(numNodesNextLvl <= ringBufferSize);
  }
  if (!(gps.interPredictionEnabledFlag
        && gps.gof_geom_entropy_continuation_enabled_flag)
      && !(gps.trisoup_enabled_flag || gbh.entropy_continuation_flag))
    decoder.clearMap();

  // save the context state for re-use by a future slice if required
  //ctxtMem = decoder.getCtx(); // ctxtMem is now directly used

  ////
  // The following is to render the last slab
  if (useLocalAttr) {
    int numPoints = laPoints[localSlabIdx].size();
    localPointCloud.resize(numPoints);
    for (int i = 0; i < numPoints; ++i)
      localPointCloud[i] = laPoints[localSlabIdx][i];
    if (numPoints) {
      rootDecoder.processNextSlabAttributes(localPointCloud, xStartLocalSlab, true);
      pointCloud.setFromPartition(localPointCloud, 0, numPoints, processedPointCount);
      processedPointCount += numPoints;
    }
  }

  // NB: the point cloud needs to be resized if partially decoded
  // OR: if geometry quantisation has changed the number of points
  // BUT: not for trisoup for which default attributes values are to be kept
  //   for not yet decoded points
  if (!forTrisoup)
    pointCloud.resize(processedPointCount);

  if (forTrisoup) {
    // TriSoup final pass (true)
    rste->callTriSoupSlice(true);
    rste->finishSlice();
  } else if (nodesRemaining) {
    // return partial coding result
    //  - add missing levels to node positions and inverse quantise
    auto nodeSizeLog2 = lvlNodeSizeLog2[maxDepth];
    for (auto& node : fifo) {
      node.pos <<= nodeSizeLog2 - QuantizerGeom::qpShift(node.qp);
      node.pos = invQuantPosition(node.qp, posQuantBitMasks, node.pos);
    }
    *nodesRemaining = std::move(fifo);
  }
}

// instanciate for Trisoup
template void
decodeGeometryOctree<true>(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int skipLastLayers,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  const Vec3<int> minimum_position,
  InterPredParams& interPredParams,
  PCCTMC3Decoder3& rootDecoder,
  RasterScanTrisoupEdgesDecoder* rste);

// instanciate for Octree
template void
decodeGeometryOctree<false>(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int skipLastLayers,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  const Vec3<int> minimum_position,
  InterPredParams& interPredParams,
  PCCTMC3Decoder3& rootDecoder,
  RasterScanTrisoupEdgesDecoder* rste);

//-------------------------------------------------------------------------

void
decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  const Vec3<int> minimum_position,
  InterPredParams& interPredParams,
  PCCTMC3Decoder3& decoder
)
{
  decodeGeometryOctree<false>(
    gps, gbh, 0, pointCloud, ctxtMem, arithmeticDecoder, nullptr,
    refFrame, sps, minimum_position, interPredParams, decoder, nullptr);
}

//-------------------------------------------------------------------------

void
decodeGeometryOctreeScalable(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int minGeomNodeSizeLog2,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  PCCTMC3Decoder3& decoder
)
{
  std::vector<PCCOctree3Node> nodes;
  InterPredParams interPredParams;
  decodeGeometryOctree<false>(
    gps, gbh, minGeomNodeSizeLog2, pointCloud, ctxtMem, arithmeticDecoder,
    &nodes, refFrame, sps, { 0, 0, 0 }, interPredParams, decoder);

  if (minGeomNodeSizeLog2 > 0) {
    size_t size =
      pointCloud.removeDuplicatePointInQuantizedPoint(minGeomNodeSizeLog2);

    pointCloud.resize(size + nodes.size());
    size_t processedPointCount = size;

    if (minGeomNodeSizeLog2 > 1) {
      uint32_t mask = uint32_t(-1) << minGeomNodeSizeLog2;
      for (auto node0 : nodes) {
        for (int k = 0; k < 3; k++)
          node0.pos[k] &= mask;
        node0.pos += 1 << (minGeomNodeSizeLog2 - 1);
        pointCloud[processedPointCount++] = node0.pos;
      }
    } else {
      for (const auto& node0 : nodes)
        pointCloud[processedPointCount++] = node0.pos;
    }
  }
}

//============================================================================

}  // namespace pcc
