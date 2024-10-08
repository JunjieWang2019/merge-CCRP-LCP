﻿/* The copyright in this software is being made available under the BSD
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
#include "TMC3.h"
#include "PCCTMC3Encoder.h"
#include "motionWip.h"
#include <unordered_map>
#include <set>
#include <random>

namespace pcc {

//============================================================================

enum class DirectMode
{
  kUnavailable,
  kAllPointSame,
  kTwoPoints
};

//============================================================================

class GeometryOctreeEncoder {
public:
  GeometryOctreeContexts& ctx;

  GeometryOctreeEncoder(
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    GeometryOctreeContexts& ctxtMem,
    EntropyEncoder* arithmeticEncoder);

  GeometryOctreeEncoder(const GeometryOctreeEncoder&) = default;
  GeometryOctreeEncoder(GeometryOctreeEncoder&&) = default;
  //GeometryOctreeEncoder& operator=(const GeometryOctreeEncoder&) = default;
  //GeometryOctreeEncoder& operator=(GeometryOctreeEncoder&&) = default;

  // dynamic OBUF
  void clearMap() { ctx.clearMap(); };
  void resetMap(bool forTrisoup) { ctx.resetMap(forTrisoup); }

  void encodePositionLeafNumPoints(int count);

  void encodeOccupancyFullNeihbourgs(
    const RasterScanContext::occupancy& occ,
    int occupancy,
    int planarMaskX,
    int planarMaskY,
    int planarMaskZ,
    int predOccupancy[8],
    int predOUnComp[8],
    bool isInter);

  void encodeRasterScanOccupancy(
    const RasterScanContext::occupancy& contextualOccupancy,
    uint8_t childOccupancy);

  /*void encodeOrdered2ptPrefix(
    const point_t points[2], Vec3<bool> directIdcm, Vec3<int>& nodeSizeLog2);*/

  /*void encodePointPosition(
    const Vec3<int>& nodeSizeLog2AfterPlanar, const Vec3<int32_t>& pos);*/

  void encodeNodeQpOffetsPresent(bool);
  void encodeQpOffset(int dqp);

  // local QU
  void createQU(localQU& qu, Vec3<int32_t> pos, int size, int baseQP);
  void encodeQU(localQU& qu, int baseQP);

  //void encodeIsIdcm(DirectMode mode);

 /*void encodeDirectPosition(
    bool geom_unique_points_flag,
    bool joint_2pt_idcm_enabled_flag,
    DirectMode mode,
    const Vec3<uint32_t>& quantMasks,
    const Vec3<int>& nodeSizeLog2,
    int shiftBits,
    const PCCOctree3Node& node,
    const OctreeNodePlanar& planar,
    PCCPointSet3& pointCloud);*/

  const GeometryOctreeContexts& getCtx() const { return ctx; }

public:
  EntropyEncoder* _arithmeticEncoder;
};

//============================================================================

GeometryOctreeEncoder::GeometryOctreeEncoder(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  GeometryOctreeContexts& ctxtMem,
  EntropyEncoder* arithmeticEncoder)
  : ctx(ctxtMem)
  , _arithmeticEncoder(arithmeticEncoder)
{
}


//============================================================================
// Encode the number of points in a leaf node of the octree.

void
GeometryOctreeEncoder::encodePositionLeafNumPoints(int count)
{
  int dupPointCnt = count - 1;
  _arithmeticEncoder->encode(dupPointCnt > 0, ctx._ctxDupPointCntGt0);
  if (dupPointCnt <= 0)
    return;

  _arithmeticEncoder->encodeExpGolomb(dupPointCnt - 1, 0, ctx._ctxDupPointCntEgl);
  return;
}

//-------------------------------------------------------------------------
// encode node occupancy bits
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
// encode node occupancy bits
//
void
GeometryOctreeEncoder::encodeOccupancyFullNeihbourgs(
  const RasterScanContext::occupancy& occ,
  int occupancy,
  int planarMaskX,
  int planarMaskY,
  int planarMaskZ,
  int predOcc[8],
  int predOccUnComp[8],
  bool isInter)
{
  // 3 planars => single child and we know its position
  if (planarMaskX && planarMaskY && planarMaskZ)
    return;

  // encode occupancy bits
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
      // bit is 0 because masked by QTBT or planar
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

    if (bitIsOne) {  // bit is 1
      partialOccupancy <<= 1;
      partialOccupancy |= 1;
      continue;
    }

    // OBUF contexts
    int ctx1, ctx2;
    bool Sparse;
    (*pointer2FunctionContext[i])(octreeNeighours, occupancy, ctx1, ctx2, Sparse);

    bool isInter2 = isInter && isPred;

    int ctxTable = 0;
    if (!isInter2) // INTRA
      ctxTable = Sparse;
    else  //  INTER
      ctxTable = 2 + predOcc[i];

    if (isInter2)
      ctx1 = (ctx1 << 1) | predOccUnComp[i]>0;

    // encode
    int bit = (occupancy >> i) & 1;
    if (Sparse) {
      auto obufIdx = ctx._MapOccupancySparse[isInter2][i].getEvolve(
        bit, ctx2, ctx1, &ctx._OBUFleafNumber, ctx._BufferOBUFleaves);
      _arithmeticEncoder->encode(
        bit, obufIdx >> 3, ctx._CtxMapDynamicOBUF[ctxTable][obufIdx],
        ctx._CtxMapDynamicOBUF[ctxTable].obufSingleBound);
    }
    else {
      auto obufIdx = ctx._MapOccupancy[isInter2][i].getEvolve(
        bit, ctx2, ctx1, &ctx._OBUFleafNumber, ctx._BufferOBUFleaves);
      _arithmeticEncoder->encode(
        bit, obufIdx >> 3, ctx._CtxMapDynamicOBUF[ctxTable][obufIdx],
        ctx._CtxMapDynamicOBUF[ctxTable].obufSingleBound);
    }

    // update partial occupancy of current node
    coded0[mask0X] += !bit;
    coded0[mask0Y] += !bit;
    coded0[mask0Z] += !bit;
    partialOccupancy <<= 1;
    partialOccupancy |= bit;
  }
}

/*
//-------------------------------------------------------------------------
// Encode part of the position of two unordred points  point in a given volume.
void
GeometryOctreeEncoder::encodeOrdered2ptPrefix(
  const point_t points[2], Vec3<bool> directIdcm, Vec3<int>& nodeSizeLog2)
{
  if (nodeSizeLog2[0] >= 1 && directIdcm[0]) {
    bool sameBit = true;
    int ctxIdx = 0;
    while (nodeSizeLog2[0] && sameBit) {
      nodeSizeLog2[0]--;
      int mask = 1 << nodeSizeLog2[0];
      auto bit0 = !!(points[0][0] & mask);
      auto bit1 = !!(points[1][0] & mask);
      sameBit = bit0 == bit1;

      _arithmeticEncoder->encode(sameBit, ctx._ctxSameBitHighx[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      if (sameBit)
        _arithmeticEncoder->encode(bit0);
    }
  }

  if (nodeSizeLog2[1] >= 1 && directIdcm[1]) {
    bool sameX = !directIdcm[0] || points[0][0] == points[1][0];
    bool sameBit = true;
    int ctxIdx = 0;
    while (nodeSizeLog2[1] && sameBit) {
      nodeSizeLog2[1]--;
      int mask = 1 << nodeSizeLog2[1];
      auto bit0 = !!(points[0][1] & mask);
      auto bit1 = !!(points[1][1] & mask);
      sameBit = bit0 == bit1;

      _arithmeticEncoder->encode(sameBit, ctx._ctxSameBitHighy[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      if (!(sameX && !sameBit))
        _arithmeticEncoder->encode(bit0);
    }
  }

  if (nodeSizeLog2[2] >= 1 && directIdcm[2]) {
    bool sameBit = true;
    bool sameXy = (!directIdcm[0] || points[0][0] == points[1][0])
      && (!directIdcm[1] || points[0][1] == points[1][1]);
    int ctxIdx = 0;
    while (nodeSizeLog2[2] && sameBit) {
      nodeSizeLog2[2]--;
      int mask = 1 << nodeSizeLog2[2];
      auto bit0 = !!(points[0][2] & mask);
      auto bit1 = !!(points[1][2] & mask);
      sameBit = bit0 == bit1;

      _arithmeticEncoder->encode(sameBit, ctx._ctxSameBitHighz[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      if (!(sameXy && !sameBit))
        _arithmeticEncoder->encode(bit0);
    }
  }
}
*/
/*
//-------------------------------------------------------------------------
// Encode a position of a point in a given volume.
void
GeometryOctreeEncoder::encodePointPosition(
  const Vec3<int>& nodeSizeLog2AfterPlanar, const Vec3<int32_t>& pos)
{
  for (int k = 0; k < 3; k++) {
    if (nodeSizeLog2AfterPlanar[k] <= 0)
      continue;

    for (int mask = 1 << (nodeSizeLog2AfterPlanar[k] - 1); mask; mask >>= 1) {
      _arithmeticEncoder->encode(!!(pos[k] & mask));
    }
  }
}
*/
//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeNodeQpOffetsPresent(bool flag)
{
  _arithmeticEncoder->encode(flag);
}

//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeQpOffset(int dqp)
{
  _arithmeticEncoder->encode(dqp != 0, ctx._ctxQpOffsetAbsGt0);
  if (dqp == 0)
    return;

  _arithmeticEncoder->encodeExpGolomb(abs(dqp) - 1, 0, ctx._ctxQpOffsetAbsEgl);
  _arithmeticEncoder->encode(dqp < 0, ctx._ctxQpOffsetSign);
}

//-------------------------------------------------------------------------
void
GeometryOctreeEncoder::createQU(
  localQU& qu, Vec3<int32_t> pos, int size, int baseQP)
{
  // simple encoder decision (for testing purpose):
  // any node with starting position having a 'y' lower than 800
  // is encoded with a lower quality
  // TODO: add encoding parameters for more versatile quality definition
  if (pos[1]  < 800) {   // [1] is vertical pos
    qu.isBaseParameters = false;
    qu.localQP = baseQP + 6;
  }
  else {
    qu.isBaseParameters = true;
    qu.localQP = baseQP;
  }
}

//-------------------------------------------------------------------------
void
GeometryOctreeEncoder::encodeQU(localQU& qu, int baseQP)
{
  _arithmeticEncoder->encode(qu.isBaseParameters,  ctx._ctxQUflag);
  if (qu.isBaseParameters)
    return;

  int diffQP = qu.localQP - baseQP;
  // sign
  _arithmeticEncoder->encode(diffQP > 0, ctx._ctxQUSign);

  // magnitude
  _arithmeticEncoder->encodeExpGolomb(
    std::abs(diffQP) - 1, 1, ctx._ctxQUQPpref, ctx._ctxQUQPsuf);
}

//-------------------------------------------------------------------------

template<typename It>
void
setNodeQpsUniform(
  Vec3<int> nodeSizeLog2,
  int qp,
  int geom_qp_multiplier_log2,
  It nodesBegin,
  It nodesEnd)
{
  // Conformance: limit the qp such that it cannot overquantize the node
  qp = std::min(qp, nodeSizeLog2.min() * 8);
  assert(qp % (1 << geom_qp_multiplier_log2) == 0);

  for (auto it = nodesBegin; it != nodesEnd; ++it)
    it->qp = qp;
}

//-------------------------------------------------------------------------
// Sets QP randomly

template<typename It>
void
setNodeQpsRandom(
  Vec3<int> nodeSizeLog2,
  int /* qp */,
  int geom_qp_multiplier_log2,
  It nodesBegin,
  It nodesEnd)
{
  // Conformance: limit the qp such that it cannot overquantize the node
  int maxQp = nodeSizeLog2.min() * 8;

  int seed = getenv("SEED") ? atoi(getenv("SEED")) : 0;
  static std::minstd_rand gen(seed);
  std::uniform_int_distribution<> uniform(0, maxQp);

  // pick a random qp, avoiding unrepresentable values
  for (auto it = nodesBegin; it != nodesEnd; ++it)
    it->qp = uniform(gen) & (~0 << geom_qp_multiplier_log2);
}

//-------------------------------------------------------------------------
// determine delta qp for each node based on the point density

template<typename It>
void
setNodeQpsByDensity(
  Vec3<int> nodeSizeLog2,
  int baseQp,
  int geom_qp_multiplier_log2,
  It nodesBegin,
  It nodesEnd)
{
  // Conformance: limit the qp such that it cannot overquantize the node
  int maxQp = nodeSizeLog2.min() * 8;
  int lowQp = PCCClip(baseQp - 8, 0, maxQp);
  int mediumQp = std::min(baseQp, maxQp);
  int highQp = std::min(baseQp + 8, maxQp);

  // NB: node.qp always uses a step size doubling interval of 8 QPs.
  //     the chosen QPs (and conformance limit) must respect the qp multiplier
  assert(lowQp % (1 << geom_qp_multiplier_log2) == 0);
  assert(mediumQp % (1 << geom_qp_multiplier_log2) == 0);
  assert(highQp % (1 << geom_qp_multiplier_log2) == 0);

  std::vector<int> numPointsInNode;
  std::vector<double> cum_prob;
  int32_t numPointsInLvl = 0;
  for (auto it = nodesBegin; it != nodesEnd; ++it) {
    numPointsInNode.push_back(it->end - it->start);
    numPointsInLvl += it->end - it->start;
  }
  std::sort(numPointsInNode.begin(), numPointsInNode.end());
  double cc = 0;
  for (auto num : numPointsInNode) {
    cc += num;
    cum_prob.push_back(cc / numPointsInLvl);
  }
  int th1 = -1, th2 = -1;
  for (int i = 0; i < cum_prob.size(); i++) {
    if (th1 == -1 && cum_prob[i] > 0.05) {
      th1 = numPointsInNode[i];
    } else if (th2 == -1 && cum_prob[i] > 0.6)
      th2 = numPointsInNode[i];
  }
  for (auto it = nodesBegin; it != nodesEnd; ++it) {
    if (it->end - it->start < th1) {
      it->qp = highQp;
    } else if (it->end - it->start < th2)
      it->qp = mediumQp;
    else
      it->qp = lowQp;
  }
}

//-------------------------------------------------------------------------

template<typename It>
void
calculateNodeQps(
  OctreeEncOpts::QpMethod method,
  Vec3<int> nodeSizeLog2,
  int baseQp,
  int geom_qp_multiplier_log2,
  It nodesBegin,
  It nodesEnd)
{
  auto fn = &setNodeQpsUniform<It>;

  switch (method) {
    using Method = OctreeEncOpts::QpMethod;
  default:
  case Method::kUniform: fn = &setNodeQpsUniform<It>; break;
  case Method::kRandom: fn = &setNodeQpsRandom<It>; break;
  case Method::kByDensity: fn = &setNodeQpsByDensity<It>; break;
  }

  fn(nodeSizeLog2, baseQp, geom_qp_multiplier_log2, nodesBegin, nodesEnd);
}

//-------------------------------------------------------------------------

void
geometryQuantization(
  PCCPointSet3& pointCloud, PCCOctree3Node& node, Vec3<int> nodeSizeLog2)
{
  QuantizerGeom quantizer = QuantizerGeom(node.qp);
  int qpShift = QuantizerGeom::qpShift(node.qp);

  for (int k = 0; k < 3; k++) {
    int quantBitsMask = (1 << nodeSizeLog2[k]) - 1;
    int32_t clipMax = quantBitsMask >> qpShift;

    for (int i = node.start; i < node.end; i++) {
      int32_t pos = int32_t(pointCloud[i][k]);
      int32_t quantPos = quantizer.quantize(pos & quantBitsMask);
      quantPos = PCCClip(quantPos, 0, clipMax);

      // NB: this representation is: |ppppppqqq00|, which, except for
      // the zero padding, is the same as the decoder.
      pointCloud[i][k] = (pos & ~quantBitsMask) | (quantPos << qpShift);
    }
  }
}

//-------------------------------------------------------------------------

void
geometryScale(
  PCCPointSet3& pointCloud, PCCOctree3Node& node, Vec3<int> quantNodeSizeLog2)
{
  QuantizerGeom quantizer = QuantizerGeom(node.qp);
  int qpShift = QuantizerGeom::qpShift(node.qp);

  for (int k = 0; k < 3; k++) {
    int quantBitsMask = (1 << quantNodeSizeLog2[k]) - 1;
    for (int i = node.start; i < node.end; i++) {
      int pos = pointCloud[i][k];
      int lowPart = (pos & quantBitsMask) >> qpShift;
      int lowPartScaled = PCCClip(quantizer.scale(lowPart), 0, quantBitsMask);
      int highPartScaled = pos & ~quantBitsMask;
      pointCloud[i][k] = highPartScaled | lowPartScaled;
    }
  }
}

//-------------------------------------------------------------------------

void
checkDuplicatePoints(
  PCCPointSet3& pointCloud,
  PCCOctree3Node& node,
  std::vector<int>& pointIdxToDmIdx)
{
  auto first = PCCPointSet3::iterator(&pointCloud, node.start);
  auto last = PCCPointSet3::iterator(&pointCloud, node.end);

  std::set<Vec3<int32_t>> uniquePointsSet;
  for (auto i = first; i != last;) {
    if (uniquePointsSet.find(**i) == uniquePointsSet.end()) {
      uniquePointsSet.insert(**i);
      i++;
    } else {
      std::iter_swap(i, last - 1);
      last--;
      pointIdxToDmIdx[--node.end] = -2;  // mark as duplicate
    }
  }
}

/*
//-------------------------------------------------------------------------
// Direct coding of position of points in node (early tree termination).

DirectMode
canEncodeDirectPosition(
  bool geom_unique_points_flag,
  const PCCOctree3Node& node,
  const PCCPointSet3& pointCloud)
{
  int numPoints = node.end - node.start;
  // Check for duplicated points only if there are less than 10.
  // NB: this limit is rather arbitrary
  if (numPoints > 10)
    return DirectMode::kUnavailable;

  bool allPointsAreEqual = numPoints > 1 && !geom_unique_points_flag;
  for (auto idx = node.start + 1; allPointsAreEqual && idx < node.end; idx++)
    allPointsAreEqual &= pointCloud[node.start] == pointCloud[idx];

  if (allPointsAreEqual)
    return DirectMode::kAllPointSame;

  if (numPoints > MAX_NUM_DM_LEAF_POINTS)
    return DirectMode::kUnavailable;

  return DirectMode::kTwoPoints;
}
*/
/*
//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeIsIdcm(DirectMode mode)
{
  bool isIdcm = mode != DirectMode::kUnavailable;
  _arithmeticEncoder->encode(isIdcm, ctx._ctxBlockSkipTh);
}
*/
/*
//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeDirectPosition(
  bool geom_unique_points_flag,
  bool joint_2pt_idcm_enabled_flag,
  DirectMode mode,
  const Vec3<uint32_t>& quantMasks,
  const Vec3<int>& effectiveNodeSizeLog2,
  int shiftBits,
  const PCCOctree3Node& node,
  const OctreeNodePlanar& planar,
  PCCPointSet3& pointCloud)
{
  int numPoints = node.end - node.start;

  switch (mode) {
  case DirectMode::kUnavailable: return;

  case DirectMode::kTwoPoints:
    _arithmeticEncoder->encode(numPoints > 1, ctx._ctxNumIdcmPointsGt1);
    if (!geom_unique_points_flag && numPoints == 1)
      _arithmeticEncoder->encode(0, ctx._ctxDupPointCntGt0);
    break;

  case DirectMode::kAllPointSame:
    _arithmeticEncoder->encode(0, ctx._ctxNumIdcmPointsGt1);
    _arithmeticEncoder->encode(1, ctx._ctxDupPointCntGt0);
    _arithmeticEncoder->encode(numPoints - 1 > 1, ctx._ctxDupPointCntGt1);
    if (numPoints - 1 > 1)
      _arithmeticEncoder->encodeExpGolomb(
        numPoints - 3, 0, ctx._ctxDupPointCntEgl);

    // only one actual psoition to code
    numPoints = 1;
  }

  // if the points have been quantised, the following representation is used
  // for point cloud positions:
  //          |---| = nodeSizeLog2 (example)
  //   ppppppqqqq00 = cloud[ptidx]
  //          |-|   = effectiveNodeSizeLog2 (example)
  // where p are unquantised bits, qqq are quantised bits, and 0 are zero bits.
  // nodeSizeLog2 is the size of the current node prior to quantisation.
  // effectiveNodeSizeLog2 is the size of the node after quantisation.
  //
  // NB: while nodeSizeLog2 may be used to access the current position bit
  //     in both quantised and unquantised forms, effectiveNodeSizeLog2 cannot
  //     without taking into account the padding.
  //
  // NB: this contrasts with node.pos, which contains the previously coded
  //     position bits ("ppppppq" in the above example) without any padding.
  //
  // When coding the direct mode, the zero padding is removed to permit
  // indexing by the effective node size instead.
  Vec3<int> points[2];
  for (int i = 0; i < numPoints; i++)
    points[i] = pointCloud[node.start + i] >> shiftBits;

  Vec3<int> nodeSizeLog2Rem = effectiveNodeSizeLog2;
  for (int k = 0; k < 3; k++)
    if (nodeSizeLog2Rem[k] > 0 && (planar.planarMode & (1 << k)))
      nodeSizeLog2Rem[k]--;

  // Indicates which components are directly coded, or coded using angular
  // contextualisation.
  Vec3<bool> directIdcm = true;

  // Jointly code two points
  if (numPoints == 2 && joint_2pt_idcm_enabled_flag) {
    // Apply an implicit ordering to the two points, considering only the
    // directly coded axes
    if (times(points[1], directIdcm) < times(points[0], directIdcm)) {
      std::swap(points[0], points[1]);
      pointCloud.swapPoints(node.start, node.start + 1);
    }

    encodeOrdered2ptPrefix(points, directIdcm, nodeSizeLog2Rem);
  }

  // code points after planar
  for (auto idx = 0; idx < numPoints; idx++) {
    encodePointPosition(nodeSizeLog2Rem, points[idx]);
  }
}
*/
//-------------------------------------------------------------------------

template <bool forTrisoup>
void
encodeGeometryOctree(
  const EncoderParams& encParams,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Encoder3& rootEncoder,
  RasterScanTrisoupEdgesEncoder* rste
  )
{
  const OctreeEncOpts& params = encParams.geom;

  if (forTrisoup) {
    assert(rste);
  }

  PCCPointSet3& compensatedPointCloud = interPredParams.compensatedPointCloud;

  const bool isInter = gbh.interPredictionEnabledFlag;

  PCCPointSet3& refPointCloud = interPredParams.referencePointCloud;
  refPointCloud = refFrame.cloud;
  PCCPointSet3& predPointCloud = refPointCloud;
  auto& mvField = interPredParams.mvField;

  int log2MinPUSize = ilog2(uint32_t(gps.motion.motion_min_pu_size - 1)) + 1;
  MSOctree& mSOctree = interPredParams.mSOctreeRef;
  if (isInter) {
    mSOctree = MSOctree(&predPointCloud, -gbh.geomBoxOrigin, std::min(2,log2MinPUSize));
  }
  MSOctree mSOctreeCurr;

  auto arithmeticEncoderIt = arithmeticEncoders.begin();
  GeometryOctreeEncoder encoder(gps, gbh, ctxtMem, arithmeticEncoderIt->get());

  // local QU
  if (forTrisoup) {
    if (gbh.qu_size_log2 > 0) {
      int quSize = (1 << gbh.qu_size_log2) * gbh.trisoupNodeSize(gps);
      std::cout << "QU size = " << quSize << "\n";
    }
    ctxtMem.quLastIndex = 0;
    ctxtMem.listOfQUs.clear();
  }

  // saved state for use with parallel bistream coding.
  // the saved state is restored at the start of each parallel octree level
  std::unique_ptr<GeometryOctreeContexts> savedState;

  const int S2 = gbh.trisoupNodeSizeLog2(gps);
  const int S = gbh.trisoupNodeSize(gps);
  const int factorS = (1 << S2) - S;
  const bool flagNonPow2 = factorS != 0;
  if (flagNonPow2)
    gbh.rootNodeSizeLog2 += 1;

  // generate the list of the node size for each level in the tree
  auto lvlNodeSizeLog2 = mkQtBtNodeSizeList(gps, params.qtbt, gbh);

  // variables for local motion

  int log2MotionBlockSize = 0;
  int log2MotionBlockSizeMin = 0;

  // local motion prediction structure -> LPUs from predPointCloud
  if (isInter) {
    log2MotionBlockSize =
      ilog2(uint32_t(gps.motion.motion_block_size - 1)) + 1;
    if (gbh.maxRootNodeDimLog2 < log2MotionBlockSize) { // LPU is bigger than root note, must adjust if possible
      log2MotionBlockSizeMin =
        ilog2(uint32_t(gps.motion.motion_min_pu_size - 1)) + 1;
      if (log2MotionBlockSizeMin <= gbh.maxRootNodeDimLog2)
        log2MotionBlockSize = gbh.maxRootNodeDimLog2;
    }
  }

  // init main fifo
  //  -- worst case size is the last level containing every input poit
  //     and each point being isolated in the previous level.
  int reservedBufferSize = pointCloud.getPointCount();
  if (gps.trisoup_enabled_flag && gbh.trisoupNodeSize(gps) > 1)
    reservedBufferSize = std::max(1000, reservedBufferSize >> 2 * gbh.trisoupNodeSizeLog2(gps) - 1);

  std::vector<PCCOctree3Node> fifo;
  std::vector<PCCOctree3Node> fifoNext;
  fifo.reserve(reservedBufferSize);
  fifoNext.reserve(reservedBufferSize);

  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(pointCloud.getPointCount());
  node00.pos = int32_t(0);
  node00.mSOctreeNodeIdx = uint32_t(0);

  //node00.numSiblingsMispredicted = 0;
  node00.predEnd = isInter ? predPointCloud.getPointCount() : uint32_t(0);
  node00.mSONodeIdx = isInter ? 0 : -1;
  node00.predStart = uint32_t(0);
  //node00.idcmEligible = false;
  //node00.isDirectMode = false;
  //node00.numSiblingsPlus1 = 8;
  node00.qp = 0;
  //node00.idcmEligible = 0;NOTE[FT]: idcmEligible is set to false at construction

  // local motion
  node00.hasMotion = 0;
  node00.isCompensated = 0;

  // local QU
  node00.quIndex = -1;

  RasterScanContext rsc(fifo);

  // map of pointCloud idx to DM idx, used to reorder the points
  // after coding.
  std::vector<int> pointIdxToDmIdx(int(pointCloud.getPointCount()), -1);
  int nextDmIdx = 0;

  //// rotating mask used to enable idcm
  //uint32_t idcmEnableMaskInit = /*mkIdcmEnableMask(gps)*/ 0; //NOTE[FT] : set to 0 by construction

  // the minimum node size is ordinarily 2**0, but may be larger due to
  // early termination for trisoup.
  int minNodeSizeLog2 = gbh.trisoupNodeSizeLog2(gps);

  // prune anything smaller than the minimum node size (these won't be coded)
  // NB: this must result in a cubic node at the end of the list
  // NB: precondition: root node size >= minNodeSizeLog2.
  lvlNodeSizeLog2.erase(
    std::remove_if(
      lvlNodeSizeLog2.begin(), lvlNodeSizeLog2.end(),
      [&](Vec3<int>& size) { return size < minNodeSizeLog2; }),
    lvlNodeSizeLog2.end());
  assert(lvlNodeSizeLog2.back() == minNodeSizeLog2);

  // append a dummy entry to the list so that depth+2 access is always valid
  lvlNodeSizeLog2.emplace_back(lvlNodeSizeLog2.back());

  // the termination depth of the octree phase
  // NB: the tree depth may be greater than the maxNodeSizeLog2 due to
  //     perverse qtbt splitting.
  // NB: by definition, the last two elements are minNodeSizeLog2
  int maxDepth = lvlNodeSizeLog2.size() - 2;

  // generate the qtbt splitting list
  //  - start at the leaf, and work up
  std::vector<int8_t> tree_lvl_partition_list;
  for (int lvl = 0; lvl < maxDepth; lvl++) {
    gbh.tree_lvl_coded_axis_list.push_back(
      ~nonSplitQtBtAxes(lvlNodeSizeLog2[lvl], lvlNodeSizeLog2[lvl + 1]));

    // Conformance: at least one axis must attempt to be coded at each level
    assert(gbh.tree_lvl_coded_axis_list.back() != 0);
  }

  // the node size where quantisation is performed
  Vec3<int> quantNodeSizeLog2 = 0;
  Vec3<uint32_t> posQuantBitMasks = 0xffffffff;
  int idcmQp = 0;
  int sliceQp = gbh.sliceQp(gps);
  int numLvlsUntilQuantization = 0;
  if (gps.geom_scaling_enabled_flag) {
    numLvlsUntilQuantization = params.qpOffsetDepth;

    // Determine the desired quantisation depth after qtbt is determined
    if (params.qpOffsetNodeSizeLog2 > 0) {
      // find the first level that matches the scaling node size
      for (int lvl = 0; lvl < maxDepth; lvl++) {
        if (lvlNodeSizeLog2[lvl].min() > params.qpOffsetNodeSizeLog2)
          continue;
        numLvlsUntilQuantization = lvl;
        break;
      }
    }

    // if an invalid depth is set, use tree height instead
    if (numLvlsUntilQuantization < 0)
      numLvlsUntilQuantization = maxDepth;
    numLvlsUntilQuantization++;
  }

  if (gps.octree_point_count_list_present_flag)
    gbh.footer.octree_lvl_num_points_minus1.reserve(maxDepth);

  if (!(isInter && gps.gof_geom_entropy_continuation_enabled_flag) && !gbh.entropy_continuation_flag) {
    encoder.clearMap();
    encoder.resetMap(forTrisoup);
  }

  // localized attributes point indexes
  std::vector<std::vector<int> > laPointIdx;
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
    laPointIdx.resize(numSlabs);
    localPointCloud.reserve(pointCloud.getPointCount());
  }

  PCCPointSet3 recPointCloud;
  recPointCloud.addRemoveAttributes(
    pointCloud.hasColors(), pointCloud.hasReflectances());
  recPointCloud.resize(pointCloud.getPointCount());
  int nRecPoints = 0;

  int lastPos0 = INT_MIN; // used to detect slice change and call for TriSoup

  // process point cloud positions depending on TriSoup node size
  if (flagNonPow2) {
    for (int np = 0; np < pointCloud.getPointCount(); np++) {
      const int temp0 = pointCloud[np][0] / S;
      const int temp1 = pointCloud[np][1] / S;
      const int temp2 = pointCloud[np][2] / S;
      pointCloud[np][0] += temp0 * factorS;
      pointCloud[np][1] += temp1 * factorS;
      pointCloud[np][2] += temp2 * factorS;
    }
  }

  if (isInter)
    mSOctreeCurr = MSOctree(&pointCloud, {}, gbh.trisoupNodeSizeLog2(gps));

  for (int depth = 0; depth < maxDepth; depth++) {
    // The tree terminated early (eg, due to IDCM or quantisation)
    // Delete any unused arithmetic coders
    if (fifo.empty()) {
      ++arithmeticEncoderIt;
      arithmeticEncoders.erase(arithmeticEncoderIt, arithmeticEncoders.end());
      break;
    }

    // setyo at the start of each level
    auto fifoCurrLvlEnd = fifo.end();
    int numNodesNextLvl = 0;

    // derive per-level node size related parameters
    auto nodeSizeLog2 = lvlNodeSizeLog2[depth];
    auto childSizeLog2 = lvlNodeSizeLog2[depth + 1];
    //// represents the largest dimension of the current node
    //int nodeMaxDimLog2 = nodeSizeLog2.max();

    const auto nodeSize = S * (1 << nodeSizeLog2[0] - S2);

    // if one dimension is not split, atlasShift[k] = 0
    int codedAxesPrevLvl = depth ? gbh.tree_lvl_coded_axis_list[depth - 1] : 7;
    int codedAxesCurLvl = gbh.tree_lvl_coded_axis_list[depth];

    const auto pointSortMask = qtBtChildSize(nodeSizeLog2, childSizeLog2);

    // Idcm quantisation applies to child nodes before per node qps
    if (--numLvlsUntilQuantization > 0) {
      // Indicate that the quantisation level has not been reached
      encoder.encodeNodeQpOffetsPresent(false);

      quantNodeSizeLog2 = nodeSizeLog2;

      for (int k = 0; k < 3; k++)
        quantNodeSizeLog2[k] = std::max(0, quantNodeSizeLog2[k]);

      // limit the idcmQp such that it cannot overquantise the node
      auto minNs = quantNodeSizeLog2.min();
      idcmQp = gps.geom_base_qp + gps.geom_idcm_qp_offset;
      idcmQp <<= gps.geom_qp_multiplier_log2;
      idcmQp = std::min(idcmQp, minNs * 8);
      for (int k = 0; k < 3; k++)
        posQuantBitMasks[k] = (1 << quantNodeSizeLog2[k]) - 1;
    }

    // determing a per node QP at the appropriate level
    if (!numLvlsUntilQuantization) {
      // Indicate that this is the level where per-node QPs are signalled.
      encoder.encodeNodeQpOffetsPresent(true);

      // idcm qps are no longer independent
      idcmQp = 0;
      quantNodeSizeLog2 = nodeSizeLog2;
      for (int k = 0; k < 3; k++)
        posQuantBitMasks[k] = (1 << quantNodeSizeLog2[k]) - 1;
      calculateNodeQps(
        params.qpMethod, nodeSizeLog2, sliceQp, gps.geom_qp_multiplier_log2,
        fifo.begin(), fifoCurrLvlEnd);
    }

    // save context state for parallel coding
    if (depth == maxDepth - 1 - gbh.geom_stream_cnt_minus1)
      if (gbh.geom_stream_cnt_minus1)
        savedState.reset(new GeometryOctreeContexts(encoder.ctx));

    // load context state for parallel coding starting one level later
    if (depth > maxDepth - 1 - gbh.geom_stream_cnt_minus1) {
      encoder.ctx = *savedState;
      encoder._arithmeticEncoder = (++arithmeticEncoderIt)->get();
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

      // encode delta qp for each octree block
      if (numLvlsUntilQuantization == 0 && !tubeIndex && !nodeSliceIndex) {
        int qpOffset = (node0.qp - sliceQp) >> gps.geom_qp_multiplier_log2;
        encoder.encodeQpOffset(qpOffset);
      }

      int shiftBits = QuantizerGeom::qpShift(node0.qp);
      auto effectiveNodeSizeLog2 = nodeSizeLog2 - shiftBits;
      auto effectiveChildSizeLog2 = childSizeLog2 - shiftBits;

      if (isLeafNode(effectiveNodeSizeLog2))
        return;

      // make quantisation work with qtbt and planar.
      /*int codedAxesCurNode = codedAxesCurLvl;
      if (shiftBits != 0) {
        for (int k = 0; k < 3; k++) {
          if (effectiveChildSizeLog2[k] < 0)
            codedAxesCurNode &= ~(4 >> k);
        }
      }*/

      if (numLvlsUntilQuantization == 0 && !tubeIndex && !nodeSliceIndex) {
        geometryQuantization(pointCloud, node0, quantNodeSizeLog2);
        if (gps.geom_unique_points_flag)
          checkDuplicatePoints(pointCloud, node0, pointIdxToDmIdx);
      }

      GeometryNeighPattern gnp{};
      // The position of the node in the parent's occupancy map
      int posInParent = 0;
      posInParent |= (node0.pos[0] & 1) << 2;
      posInParent |= (node0.pos[1] & 1) << 1;
      posInParent |= (node0.pos[2] & 1) << 0;
      posInParent &= codedAxesPrevLvl;

      std::array<int32_t, 8> predUnCompCounts = {};

      if (!tubeIndex && !nodeSliceIndex) {
        //local motion : determine PU tree by motion search and RDO
        if (isInter && nodeSizeLog2[0] == log2MotionBlockSize) {
          node0.hasMotion = motionSearchForNode(
            mSOctreeCurr, mSOctree, &node0, encParams.motion,
            encParams.gps.motion, nodeSize, encoder._arithmeticEncoder,
            mvField, currPUIdx, flagNonPow2, S, S2);
          node0.mvFieldNodeIdx = currPUIdx++;
        }

        // code split PU flag. If not split, code  MV and apply MC
        // results of MC are stacked in compensatedPointCloud that starts empty
        if (node0.mvFieldNodeIdx != -1) {
          encode_splitPU_MV_MC<false>(mSOctree,
            &node0, mvField, node0.mvFieldNodeIdx, gps.motion, nodeSize,
            encoder._arithmeticEncoder, &compensatedPointCloud,
            flagNonPow2, S, S2);
        }

        // split the current node into 8 children
        //  - perform an 8-way counting sort of the current node's points
        //  - (later) map to child nodes
        countingSort(
          PCCPointSet3::iterator(&pointCloud, node0.start),
          PCCPointSet3::iterator(&pointCloud, node0.end), node0.childCounts,
          [=](const PCCPointSet3::Proxy& proxy) {
            const auto& point = *proxy;
            return !!(int(point[2]) & pointSortMask[2])
              | (!!(int(point[1]) & pointSortMask[1]) << 1)
              | (!!(int(point[0]) & pointSortMask[0]) << 2);
          });

        /// sort and partition the predictor for local motion
        // TODO: check we can remove
        //node0.predCounts = {};
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

        // generate the bitmap of child occupancy and count
        // the number of occupied children in node0.
        int occupancy = 0;
        for (int i = 0; i < 8; i++) {
          occupancy |= (node0.childCounts[i] > 0) << i;
        }
        node0.childOccupancy = occupancy;
        node0.predPointsStartIdx = node0.predStart;
      }

      // inter information
      //int predOccupancy = 0;
      //int predFailureCount = 0;

      // occupancy inter predictor
      int predOccupancy[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      int predOccupancyUnComp[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      const int L = 1 << nodeSizeLog2[0];

      // TODO avoid computing it at each pass?
      for (int i = 0; i < 8; i++) {
        //bool childOccupiedTmp = !!node0.childCounts[i];
        //predOccupancy |= (node0.predCounts[i] > 0) << i;
        //predFailureCount += childOccupiedTmp != childPredicted;
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

      // local QU
      if (forTrisoup && !tubeIndex && !nodeSliceIndex) {
        if (gbh.qu_size_log2 && nodeSizeLog2[0] == S2 + gbh.qu_size_log2) {

          ctxtMem.listOfQUs.emplace_back();
          localQU& qu = ctxtMem.listOfQUs.back();
          encoder.createQU(
            qu, (node0.pos << nodeSizeLog2[0] - S2) * S,
            (1 << nodeSizeLog2[0] - S2) * S, gbh.trisoup_QP);
          encoder.encodeQU(qu, gbh.trisoup_QP);
          node0.quIndex = ctxtMem.quLastIndex++;
        }
      }

      //DirectMode mode = DirectMode::kUnavailable;
      // At the scaling depth, it is possible for a node that has previously
      // been marked as being eligible for idcm to be fully quantised due
      // to the choice of QP.  There is therefore nothing to code with idcm.
      //if (isLeafNode(effectiveNodeSizeLog2))
      //  node0.idcmEligible = false; //NOTE[FT]: already false
      //if (node0.idcmEligible && !tubeIndex && !nodeSliceIndex) {
      //  // todo(df): this is pessimistic in the presence of idcm quantisation,
      //  // since that is eligible may only meet the point count constraint
      //  // after quantisation, which is performed after the decision is taken.
      //  mode = canEncodeDirectPosition(
      //    gps.geom_unique_points_flag, node0, pointCloud);
      //}

      // N.B: contextualOccupancy is only valid during first pass on the node
      RasterScanContext::occupancy contextualOccupancy;
      //OctreeNodePlanar planar;
      if (!isLeafNode(effectiveNodeSizeLog2) && !tubeIndex && !nodeSliceIndex) {
        // update contexts
        rsc.nextNode(&*fifoCurrNode, contextualOccupancy);
        node0.neighPattern = gnp.neighPattern = contextualOccupancy.neighPattern;
        gnp.adjNeighOcc[0] = contextualOccupancy.childOccupancyContext[4];
        gnp.adjNeighOcc[1] = contextualOccupancy.childOccupancyContext[10];
        gnp.adjNeighOcc[2] = contextualOccupancy.childOccupancyContext[12];
      }

      //if (node0.idcmEligible && !tubeIndex && !nodeSliceIndex
      //    && !gps.geom_planar_disabled_idcm_angular_flag)
      //  encoder.encodeIsIdcm(mode);

      //if (node0.isDirectMode && !tubeIndex && !nodeSliceIndex) {
      //  int idcmShiftBits = shiftBits;
      //  auto idcmSize = effectiveNodeSizeLog2;

      //  if (idcmQp) {
      //    node0.qp = idcmQp;
      //    idcmShiftBits = QuantizerGeom::qpShift(idcmQp);
      //    idcmSize = nodeSizeLog2 - idcmShiftBits;
      //    geometryQuantization(pointCloud, node0, quantNodeSizeLog2);

      //    if (gps.geom_unique_points_flag)
      //      checkDuplicatePoints(pointCloud, node0, pointIdxToDmIdx);
      //  }

      //  encoder.encodeDirectPosition(
      //    gps.geom_unique_points_flag, gps.joint_2pt_idcm_enabled_flag,
      //    mode, posQuantBitMasks, idcmSize,
      //    idcmShiftBits, node0, planar, pointCloud);

      //  // calculating the number of points coded by IDCM for determining
      //  // eligibility of planar mode
      //  if (checkPlanarEligibilityBasedOnOctreeDepth)
      //    numPointsCodedByIdcm += node0.end - node0.start;

      //  // inverse quantise any quantised positions
      //  geometryScale(pointCloud, node0, quantNodeSizeLog2);

      //  // point reordering to match decoder's order
      //  for (auto idx = node0.start; idx < node0.end; idx++)
      //    pointIdxToDmIdx[idx] = nextDmIdx++;

      //  // NB: by definition, this is the only child node present
      //  if (/*gps.inferred_direct_coding_mode*/ 0 <= 1)
      //    assert(node0.numSiblingsPlus1 == 1);

      //  continue;
      //}
      //if (node0.isDirectMode)
      //  continue;

      // when all points are quantized to a single point
      if (!isLeafNode(effectiveNodeSizeLog2) && !tubeIndex && !nodeSliceIndex) {
        // encode child occupancy map
        assert(node0.childOccupancy > 0);

        // planar mode for current node
        // mask to be used for the occupancy coding
        // (bit =1 => occupancy bit not coded due to not belonging to the plane)
        // int planarMask[3] = {0, 0, 0};
        // maskPlanar(planar, planarMask, codedAxesCurNode);

        encoder.encodeOccupancyFullNeihbourgs(
          contextualOccupancy, node0.childOccupancy, planarMask[0], planarMask[1],
          planarMask[2], predOccupancy, predOccupancyUnComp,isInter);
      }

      // calculating the number of subnodes for determining eligibility
      // population count of occupancy for IDCM
      //int numOccupied = popcnt(node0.childOccupancy);

      // Leaf nodes are immediately coded.  No further splitting occurs.
      if (tubeIndex && nodeSliceIndex && isLeafNode(effectiveChildSizeLog2)) {
        // inverse quantise any quantised positions
        geometryScale(pointCloud, node0, quantNodeSizeLog2);

        if (!useLocalAttr) {
          for (auto idx = node0.start; idx < node0.end; idx++)
            pointIdxToDmIdx[idx] = nextDmIdx++;
        } else {
          int slabIdx;
          if (isLastDepth) {
            int nodeposX = node0.pos[0] << !!(codedAxesCurLvl & 4) + effectiveChildSizeLog2[0];

            // rendering attributes of a finished Slab
            while (nodeposX >= xStartLocalSlab + slabThickness) {
              localPointCloud.appendPartition(pointCloud,laPointIdx[localSlabIdx]);
              laPointIdx[localSlabIdx] = std::vector<int>(); // just release memory
              auto nRecPointsLocal = localPointCloud.getPointCount();
              if (nRecPointsLocal) {
                rootEncoder.processNextSlabAttributes(localPointCloud, xStartLocalSlab, false);
                recPointCloud.setFromPartition(localPointCloud, 0, nRecPointsLocal, nRecPoints);
                localPointCloud.clear();
                nRecPoints += nRecPointsLocal;
              }
              // keep slabs aligned on a regular grid
              xStartLocalSlab += slabThickness;
              ++localSlabIdx;
            }
            slabIdx = localSlabIdx;
          } else {
            // TODO: better
            slabIdx = pointCloud[node0.start][0] / slabThickness;
            // TODO: shall we and how to handle case where quantization is bigger than thickness
          }
          int idxLaPoint = laPointIdx[slabIdx].size();
          laPointIdx[slabIdx].resize(idxLaPoint + node0.end - node0.start);
          for (auto idx = node0.start; idx < node0.end; idx++)
            laPointIdx[slabIdx][idxLaPoint++] = idx;
        }

        for (int i = 0; i < 8; i++) {
          if (!node0.childCounts[i]) {
            // child is empty: skip
            continue;
          }

          // if the bitstream is configured to represent unique points,
          // no point count is sent.
          if (gps.geom_unique_points_flag) {
            assert(node0.childCounts[i] == 1);
            continue;
          }

          encoder.encodePositionLeafNumPoints(node0.childCounts[i]);
        }

        // leaf nodes do not get split
        return;
      }

      if (!isLeafNode(effectiveChildSizeLog2)) {
        // push child nodes to fifo
        for (int i = 0; i < 2; ++i) {
          int childIndex = (nodeSliceIndex << 2) + (tubeIndex << 1) + i;
          bool occupiedChild = node0.childCounts[childIndex] > 0;
          if (!occupiedChild) {
            // child is empty: skip
            node0.predPointsStartIdx += node0.predCounts[childIndex];
          }
          else {
            // create new child
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

            // nodeSizeLog2 > 1: for each child:
            //  - determine elegibility for IDCM
            //  - directly code point positions if IDCM allowed and selected
            //  - otherwise, insert split children into fifo while updating neighbour state
            int childPointsStartIdx = node0.start;

            for (int j = 0; j < childIndex; ++j)
              childPointsStartIdx += node0.childCounts[j];

            child.start = childPointsStartIdx;
            childPointsStartIdx += node0.childCounts[childIndex];
            child.end = childPointsStartIdx;

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

            if (isInter) {
              child.mSOctreeNodeIdx = mSOctreeCurr.nodes[node0.mSOctreeNodeIdx].child[childIndex];
              assert((forTrisoup && flagNonPow2) || child.mSOctreeNodeIdx); // if not, missmatch between octrees
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
            //  idcmEnableMask = rotateRight(idcmEnableMask, 1);
            //}
#endif

            if (forTrisoup) {
              if (isLastDepth) {
                rste->leaves.emplace_back(child);
                rste->leaves.back().pos *= S;
                if (flagNonPow2) {
                  int maskS = (1 << S2) - 1;
                  for (int np = child.start; np < child.end; np++) {
                    pointCloud[np][0] = ((pointCloud[np][0] >> S2) * S) + (pointCloud[np][0] & maskS);
                    pointCloud[np][1] = ((pointCloud[np][1] >> S2) * S) + (pointCloud[np][1] & maskS);
                    pointCloud[np][2] = ((pointCloud[np][2] >> S2) * S) + (pointCloud[np][2] & maskS);
                  }

                  for (int np = child.predStart; np < child.predEnd; np++) {
                    compensatedPointCloud[np][0] = ((compensatedPointCloud[np][0] >> S2) * S) + (compensatedPointCloud[np][0] & maskS);
                    compensatedPointCloud[np][1] = ((compensatedPointCloud[np][1] >> S2) * S) + (compensatedPointCloud[np][1] & maskS);
                    compensatedPointCloud[np][2] = ((compensatedPointCloud[np][2] >> S2) * S) + (compensatedPointCloud[np][2] & maskS);
                  }
                }

                if (lastPos0 != INT_MIN && child.pos[0] != lastPos0)
                  rste->callTriSoupSlice(false); // TriSoup unpile slices (not final = false)
                lastPos0 = child.pos[0];
              }
            }

            numNodesNextLvl++;
          }
        }
      }
    });

    fifo.resize(0);
    fifo.swap(fifoNext);

    // calculate the number of points that would be decoded if decoding were
    // to stop at this point.
    if (gps.octree_point_count_list_present_flag) {
      int numPtsAtLvl = numNodesNextLvl + nextDmIdx - 1;
      gbh.footer.octree_lvl_num_points_minus1.push_back(numPtsAtLvl);
    }
  }

  if (!(gps.interPredictionEnabledFlag
    && gps.gof_geom_entropy_continuation_enabled_flag)
    && !(gps.trisoup_enabled_flag || gbh.entropy_continuation_flag)) {
    encoder.clearMap();
  }

  // the last element is the number of decoded points
  if (!gbh.footer.octree_lvl_num_points_minus1.empty())
    gbh.footer.octree_lvl_num_points_minus1.pop_back();

  // save the context state for re-use by a future slice if required
  //ctxtMem = encoder.getCtx(); // ctxtMem is now directly used

  if (forTrisoup) {
    // TriSoup final pass (true)
    rste->callTriSoupSlice(true);
    rste->finishSlice();
    return;
  } else if (nodesRemaining) {
    // return partial coding result
    //  - add missing levels to node positions
    //  - inverse quantise the point cloud
    // todo(df): this does not yet support inverse quantisation of node.pos
    auto nodeSizeLog2 = lvlNodeSizeLog2[maxDepth];
    for (auto& node : fifo) {
      node.pos <<= nodeSizeLog2;
      geometryScale(pointCloud, node, quantNodeSizeLog2);
    }
    *nodesRemaining = std::move(fifo);
    return;
  }

  ////
  // The following is to re-order the points according in the decoding
  // order since IDCM causes leaves to be coded earlier than they
  // otherwise would.
  if (!useLocalAttr) {
    // copy points with DM points first, the rest second
    nRecPoints = nextDmIdx;
    for (int i = 0; i < pointIdxToDmIdx.size(); i++) {
      int dstIdx = pointIdxToDmIdx[i];
      if (dstIdx == -1) {
        dstIdx = nRecPoints++;
      }
      else if (dstIdx == -2) {  // ignore duplicated points
        continue;
      }

      recPointCloud[dstIdx] = pointCloud[i];
      if (pointCloud.hasColors())
        recPointCloud.setColor(dstIdx, pointCloud.getColor(i));
      if (pointCloud.hasReflectances())
        recPointCloud.setReflectance(dstIdx, pointCloud.getReflectance(i));
    }
  } else {
  ////
  // The following is to render the last slab
    localPointCloud.appendPartition(pointCloud,laPointIdx[localSlabIdx]);
    auto nRecPointsLocal = localPointCloud.getPointCount();
    if (nRecPointsLocal) {
      rootEncoder.processNextSlabAttributes(localPointCloud, xStartLocalSlab, true);
      recPointCloud.setFromPartition(localPointCloud, 0, nRecPointsLocal, nRecPoints);
      nRecPoints += nRecPointsLocal;
    }
  }
  recPointCloud.resize(nRecPoints);
  swap(pointCloud, recPointCloud);
}

// instanciate for Trisoup
template void
encodeGeometryOctree<true>(
  const EncoderParams& encParams,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Encoder3& rootEncoder,
  RasterScanTrisoupEdgesEncoder* rste);

// instanciate for Octree
template void
encodeGeometryOctree<false>(
  const EncoderParams& encParams,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Encoder3& rootEncoder,
  RasterScanTrisoupEdgesEncoder* rste);

//============================================================================

void
encodeGeometryOctree(
  const EncoderParams& opt,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Encoder3& encoder)
{
  encodeGeometryOctree<false>(
    opt, gps, gbh, pointCloud, ctxtMem, arithmeticEncoders, nullptr, refFrame,
    sps, interPredParams, encoder);
}

//-------------------------------------------------------------------------

//============================================================================
}  // namespace pcc
