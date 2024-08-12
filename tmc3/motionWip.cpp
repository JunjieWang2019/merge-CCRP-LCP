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
#include "TMC3.h"

#include "motionWip.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <set>
#include <map>
#include <vector>
#include <unordered_map>

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "geometry_octree.h"
#include "PCCTMC3Encoder.h"

namespace pcc {

//============================================================================
struct MotionEntropy {
  AdaptiveBitModel splitPu;
  StaticBitModel mvSign;
  AdaptiveBitModel mvIsZero;
  AdaptiveBitModel mvIsOne;
  AdaptiveBitModel mvIsTwo;
  AdaptiveBitModel mvIsThree;
  StaticBitModel expGolombV0;
  AdaptiveBitModel expGolombV[6];
  AdaptiveBitModel _ctxLocalMV;

  MotionEntropy();
};

//----------------------------------------------------------------------------
class MotionEntropyEncoder : public MotionEntropy {
public:
  MotionEntropyEncoder(EntropyEncoder* arithmeticEncoder)
    : _arithmeticEncoder(arithmeticEncoder)
  {}

  // local
  void encodeSplitPu(int symbol);
  void encodeVector(const point_t& mv);

private:
  EntropyEncoder* _arithmeticEncoder;
};

//----------------------------------------------------------------------------
class MotionEntropyDecoder : public MotionEntropy {
public:
  MotionEntropyDecoder(EntropyDecoder* arithmeticDecoder)
    : _arithmeticDecoder(arithmeticDecoder)
  {}

  //local
  bool decodeSplitPu();
  void decodeVector(point_t* mv);

private:
  EntropyDecoder* _arithmeticDecoder;
};

//----------------------------------------------------------------------------

struct MotionEntropyEstimate {
  MotionEntropyEstimate(const MotionEntropy& codec, size_t wSize, int boundPrefix, int boundSuffix);

  double estimateVector(const point_t& mv) const;

  double estimateSplit(bool splitFlag) const {
    const double inv2powkFPP = 1.525878906250000e-05;  // thuis 2^-kFPP ; kFPP =16
    return double(hSplitPu[splitFlag])* inv2powkFPP;
  }

private:
  static constexpr int kFPP = 16;
  int64_t hMvIsZero[2];
  int64_t hMvIsOne[2];
  int64_t hMvIsTwo[2];
  int64_t hMvIsThree[2];
  int64_t hExpGolombV[6][2];
  int64_t hSplitPu[2];

  int64_t prepareEstimate(unsigned absval) const;
  void prepareEstimateVector();

  int64_t estimateComponent(unsigned absval) const {
    if (absval < LUT_MVestimate.size())
      return LUT_MVestimate[absval];
    auto res = dyn_MVestimate.emplace(std::make_pair(absval, 0LL));
    if (res.second)
      res.first->second = prepareEstimate(absval);
    return res.first->second;
  }
  int boundPrefix;
  int boundSuffix;
  std::vector<int64_t> LUT_MVestimate;
  mutable std::map<int,int64_t> dyn_MVestimate;
};

//============================================================================
MotionEntropy::MotionEntropy()
{
  //local
  splitPu.reset(true);
  mvIsZero.reset(true);
  mvIsOne.reset(true);
  mvIsTwo.reset(true);
  mvIsThree.reset(true);
  for (int i = 0; i < 6; i++)
    expGolombV[i].reset(true);
}

//----------------------------------------------------------------------------
MotionEntropyEstimate::MotionEntropyEstimate(const MotionEntropy& codec, size_t wSize, int boundPrefix, int boundSuffix)
: boundPrefix(boundPrefix)
, boundSuffix(boundSuffix)
, LUT_MVestimate(wSize+1)
{
  codec.splitPu.getEntropy<kFPP>(hSplitPu);
  codec.mvIsZero.getEntropy<kFPP>(hMvIsZero);
  codec.mvIsOne.getEntropy<kFPP>(hMvIsOne);
  codec.mvIsTwo.getEntropy<kFPP>(hMvIsTwo);
  codec.mvIsThree.getEntropy<kFPP>(hMvIsThree);
  for (int i = 0; i < 6; i++)
    codec.expGolombV[i].getEntropy<kFPP>(hExpGolombV[i]);

  prepareEstimateVector();
}

//----------------------------------------------------------------------------
inline void
MotionEntropyEncoder::encodeSplitPu(int symbol)
{
  _arithmeticEncoder->encode(symbol, splitPu);
}

//----------------------------------------------------------------------------

inline bool
MotionEntropyDecoder::decodeSplitPu()
{
  return _arithmeticDecoder->decode(splitPu);
}

//----------------------------------------------------------------------------

inline void
MotionEntropyEncoder::encodeVector(
  const point_t& mv)
{
  for (int comp = 0; comp < 3; comp++) {
    int v = mv[comp];
    if (v == 0) {
      _arithmeticEncoder->encode(1, mvIsZero);
    }
    else {
      _arithmeticEncoder->encode(0, mvIsZero);

      _arithmeticEncoder->encode(v < 0, mvSign);
      if (v < 0)
        v = -v;
      v--;
      _arithmeticEncoder->encode(v == 0, mvIsOne);

      if (!v) {
        continue;
      }
      v--;

      _arithmeticEncoder->encode(v == 0, mvIsTwo);
      if (!v) {
        continue;
      }
      v--;

      _arithmeticEncoder->encode(v == 0, mvIsThree);
      if (!v) {
        continue;
      }
      v--;

      // expGolomb on |v|-1 with truncation
      _arithmeticEncoder->encodeExpGolomb(uint32_t(v), 1, _ctxLocalMV);
    }
  }
}

//----------------------------------------------------------------------------

inline void
MotionEntropyDecoder::decodeVector(point_t* mv)
{
  for (int comp = 0; comp < 3; comp++) {
    if (_arithmeticDecoder->decode(mvIsZero)) {
      (*mv)[comp] = 0;
      continue;
    }
    bool sign = _arithmeticDecoder->decode(mvSign);

    if (_arithmeticDecoder->decode(mvIsOne)) {
      (*mv)[comp] = sign ? -1 : 1;
      continue;
    }
    if (_arithmeticDecoder->decode(mvIsTwo)) {
      (*mv)[comp] = sign ? -2 : 2;
      continue;
    }
    if (_arithmeticDecoder->decode(mvIsThree)) {
      (*mv)[comp] = sign ? -3 : 3;
      continue;
    }

    int v = 4 + _arithmeticDecoder->decodeExpGolomb(1, _ctxLocalMV);
    if (sign)
      v = -v;
    (*mv)[comp] = v;
  }
}



//----------------------------------------------------------------------------

int64_t
MotionEntropyEstimate::prepareEstimate(unsigned absval) const
{
  int v = absval;

  if (!v) {
    return hMvIsZero[1];
  }
  else {
    int64_t r = hMvIsZero[0] + (1ULL << kFPP);  // unpredictable sign
    v--;

    r += hMvIsOne[v == 0];
    if (!v) {
      return r;
    }
    v--;

    r += hMvIsTwo[v == 0];
    if (!v) {
      return r;
    }
    v--;

    r += hMvIsThree[v == 0];
    if (!v) {
      return r;
    }
    v--;

    int k = 1;  // initially, expgolomb order
    int ctx_number = 0;
    int num_bit_prefix = 0;
    while (1) {
      if (boundPrefix && v >= static_cast<unsigned int>(1 << k)) {
        // unary part
        r += hExpGolombV[ctx_number][1];
        v = v - (1 << k);
        k++;

        ctx_number++;
        if (ctx_number > 5)
          ctx_number = 5;
        num_bit_prefix++;
        continue;
      }

      // terminated zero of unary part + suffix
      if (num_bit_prefix < boundPrefix)
        r += hExpGolombV[ctx_number][0];
      if (num_bit_prefix == boundPrefix)
        k = boundSuffix;
      r += (int64_t(k) << kFPP);
      break;
    }
    return r;
  }
}

//----------------------------------------------------------------------------

void
MotionEntropyEstimate::prepareEstimateVector()
{
  for (int v0 = 0; v0 < LUT_MVestimate.size(); v0++) {
    LUT_MVestimate[v0] = prepareEstimate(v0);
  }
}


//----------------------------------------------------------------------------

double
MotionEntropyEstimate::estimateVector(
  const point_t& mv) const
{
  const double inv2powkFPP = 1.525878906250000e-05;  // thuis 2^-kFPP; kFPP =16
  return double(estimateComponent(std::abs(mv[0]))
    + estimateComponent(std::abs(mv[1]))
    + estimateComponent(std::abs(mv[2]))) * inv2powkFPP;
}

//============================================================================
static const int LUT_LOG2[64]{
  INT_MIN, 0,  16, 25, 32, 37, 41, 45, 48, 51, 53, 55, 57, 59, 61, 63,
  64,      65, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 79,
  80,      81, 81, 82, 83, 83, 84, 85, 85, 86, 86, 87, 87, 88, 88, 89,
  89,      90, 90, 91, 91, 92, 92, 93, 93, 93, 94, 94, 95, 95, 95, 96};

//----------------------------------------------------------------------------
inline int
plus1log2shifted4(int x)
{
  if (x < 62)
    return LUT_LOG2[x + 1];

  x++;
  int result = 0;
  while (x >= 64) {
    x >>= 1;
    result += 16;
  }

  return result + LUT_LOG2[x];
}


//----------------------------------------------------------------------------
int
roundIntegerHalfInf(const double x)
{
  return (x >= 0) ? int(x + 0.5) : -int(-x + 0.5);
}

//----------------------------------------- LOCAL MOTION -------------------


//============================================================================

int
deriveMotionMaxPrefixBits(int window_size)
{
  int max_MV =
    std::max(0, window_size - 1);
  if (max_MV >= 256)
    return 31;

  static int LUT_bound_prefix[256] = {
    0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7};

  return LUT_bound_prefix[max_MV];
}

int
deriveMotionMaxSuffixBits(int window_size)
{

  int max_MV =
    std::max(0, window_size  - 1);
  if (max_MV >= 256)
    return 31;

  static int LUT_bound_suffix[256] = {
    0, 1, 0, 1, 2, 2, 0, 1, 2, 2, 3, 3, 3, 3, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 1};

  return LUT_bound_suffix[max_MV];
}

//----------------------------------------------------------------------------

void
bounded_splitPU_MC(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  point_t nodeSizeLog2,
  point_t boundPos0,
  point_t boundPos1,
  PCCPointSet3* compensatedPointCloud)
{
  const int node_size = 1 << nodeSizeLog2[0];

  auto& puNode = mvField.puNodes[puNodeIdx];
  // --------------  non-split / terminal case  ----------------
  if (node_size <= param.motion_min_pu_size || !puNode._childsMask) {

    // use MV
    point_t MV = mvField.mvPool[puNode._mvIdx];

    auto node0pos0 = node0->pos << nodeSizeLog2[0];
    auto node0pos1 = (node0->pos + 1 << nodeSizeLog2[0]) - 1;

    for (int k = 0; k < 3; ++k) {
      node0pos0[k] = std::max(boundPos0[k], node0pos0[k]);
      node0pos1[k] = std::min(boundPos1[k], node0pos1[k]);
    }

    // motion compensated attribute projection
    mSOctree.apply_recolor_motion(MV, node0, *compensatedPointCloud);

    node0->isCompensated = true;
    return;
  }

  // --------------- split case ----------------------
}

//----------------------------------------------------------------------------

template <bool mcap>
void
encode_splitPU_MV_MC(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2)
{
  MotionEntropyEncoder motionEncoder(arithmeticEncoder);

  auto& puNode = mvField.puNodes[puNodeIdx];
  // --------------  non-split / terminal case  ----------------
  if (nodeSize <= param.motion_min_pu_size || !puNode._childsMask) {
    if (nodeSize > param.motion_min_pu_size) {
      motionEncoder.encodeSplitPu(0);
    }

    // encode MV
    point_t MV = mvField.mvPool[puNode._mvIdx];
    motionEncoder.encodeVector(MV);

    if (!mcap) {
      mSOctree.apply_motion(
        node0->pos * nodeSize,
        (node0->pos + 1) * nodeSize - 1,
        MV, node0, compensatedPointCloud, mSOctree.depth, flagNonPow2, S, S2);
    } else {
      mSOctree.apply_recolor_motion(MV, node0, *compensatedPointCloud);
    }
    node0->isCompensated = true;
    return;
  }

  // --------------- split case ----------------------
  motionEncoder.encodeSplitPu(1);
}

// instanciate for geometry
template void encode_splitPU_MV_MC<false>(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2);

// instanciate for mcap
template void encode_splitPU_MV_MC<true>(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2);

//----------------------------------------------------------------------------

template <bool mcap>
void
decode_splitPU_MV_MC(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2)
{
  // Note:
  // geometry mvField is stored only for attributes and local attributes processing
  // when attributes do not use dual motion
  MotionEntropyDecoder motionDecoder(arithmeticDecoder);

  auto& puNode = mvField.puNodes[puNodeIdx];
  puNode.set_pos0(node0->pos * nodeSize);
  puNode._puSizeLog2 = ilog2(uint32_t(nodeSize - 1)) + 1; // TODO: check if we need the true size at some point or clean

  // decode split flag
  bool split = false;
  if (nodeSize > param.motion_min_pu_size)
    split = motionDecoder.decodeSplitPu();

  if (!split) {  // not split
                 // decode MV
    point_t MV = 0;
    motionDecoder.decodeVector(&MV);
    puNode._childsMask = 0;
    puNode._mvIdx = mvField.mvPool.size();
    mvField.mvPool.emplace_back(MV);

    if (!mcap) {
      mSOctree.apply_motion(
        node0->pos * nodeSize,
        (node0->pos + 1) * nodeSize - 1,
        MV, node0, compensatedPointCloud, mSOctree.depth, flagNonPow2, S, S2);
    } else {
      mSOctree.apply_recolor_motion(MV, node0, *compensatedPointCloud);
    }

    node0->isCompensated = true;
    return;
  }

  // split; nothing to do
}

// instanciate for geometry
template void decode_splitPU_MV_MC<false>(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2);

// instanciate for mcap
template void decode_splitPU_MV_MC<true>(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2);

//============================================================================

struct MSOctreeStackElt {
  MSOctreeStackElt(int32_t nodeIdx, int32_t childIdx)
    : nodeIdx(nodeIdx), childIdx(childIdx) {}
  int32_t nodeIdx;
  int32_t childIdx;
};

MSOctree::MSOctree(
    PCCPointSet3* predPointCloud,
    point_t offsetOrigin,
    uint32_t leafSizeLog2
  )
  : pointCloud(predPointCloud)
  , offsetOrigin(offsetOrigin)
{
  PCCPointSet3& pointCloud(*this->pointCloud);
  if (!pointCloud.size())
    return;

  nodes.reserve(pointCloud.size());

  auto bbox = pointCloud.computeBoundingBox();

  assert(bbox.min[0] >= 0 && bbox.min[1] >= 0 && bbox.min[2] >= 0);

  maxDepth = numBits(std::max(std::max(bbox.max[0],bbox.max[1]), bbox.max[2]));
  depth = maxDepth - std::min(maxDepth, leafSizeLog2);

  // push the first node
  nodes.emplace_back();
  MSONode& node00 = nodes.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(pointCloud.getPointCount());
  node00.pos0 = point_t{0} + offsetOrigin;
  node00.parent = 0;
  node00.sizeMinus1 = int32_t((1 << maxDepth) - 1);

  // constructing depth first might be better for memory
  std::vector<MSOctreeStackElt> stack;
  stack.reserve(32);

  stack.push_back(MSOctreeStackElt(0,-1));

  while (!stack.empty()) {
    auto& curr = stack.back();

    const int currDepth = stack.size() - 1;
    const int childSizeLog2 = maxDepth - currDepth - 1;
    const int pointSortMask = 1 << childSizeLog2;
    const int childSizeMinus1 = pointSortMask - 1;

    if (curr.childIdx == -1) {
      MSONode& node0 = nodes[curr.nodeIdx]; // adress may be invalidated after emplace_back....
      std::array<int32_t, 8> childCounts = {};

      countingSort(
        PCCPointSet3::iterator(&pointCloud, node0.start),
        PCCPointSet3::iterator(&pointCloud, node0.end),
        childCounts, [=](const PCCPointSet3::Proxy& proxy) {
          const auto & point = *proxy;
          return !!(int(point[2]) & pointSortMask)
            | (!!(int(point[1]) & pointSortMask) << 1)
            | (!!(int(point[0]) & pointSortMask) << 2);
        });

      int32_t childNodeIdx = nodes.size();
      for (int i = 0; i < 8; ++i) {
        if (childCounts[i]) {
          node0.child[i] = childNodeIdx++;
        }
      }
      // create child nodes
      uint32_t childStart = node0.start;
      auto pos0 = node0.pos0;
      for (int i = 0; i < 8; ++i) {
        if (childCounts[i]) {
          nodes.emplace_back(); // node0 may be invalidated after that
          MSONode& child0 = nodes.back();
          uint32_t childEnd = childStart + childCounts[i];
          child0.start = childStart;
          child0.end = childEnd;
          childStart = childEnd;
          child0.pos0 = pos0 + (Vec3<int32_t>{i >> 2, (i >> 1) & 1, i & 1} << childSizeLog2);
          child0.sizeMinus1 = childSizeMinus1;
          child0.parent = curr.nodeIdx;
        }
      }

      if (currDepth == depth-1) {
        stack.pop_back();
        continue;
      }
      curr.childIdx = 0;
    }

    MSONode& node0 = nodes[curr.nodeIdx];

    while (curr.childIdx < 8 && !node0.child[curr.childIdx])
      curr.childIdx++;

    if (curr.childIdx < 8) {
      stack.push_back(MSOctreeStackElt(node0.child[curr.childIdx], -1));
      curr.childIdx++;
    }
    else {
      stack.pop_back();
    }
  }

  // now everything is built, we can offset the points
  for (int i=0; i < pointCloud.size(); ++i) {
    pointCloud[i] += offsetOrigin;
  }

  allocRingBuffers();
}

//----------------------------------------------------------------------------

int
MSOctree::nearestNeighbour_updateDMax(point_t pos, int32_t& d_max, bool approximate) const {
  if (approximate)
    return iApproximateNearestNeighbour_updateDMax(pos, d_max);
  else
    return iNearestNeighbour_updateDMax(pos, d_max);
}

//----------------------------------------------------------------------------

struct NNStackElt {
  NNStackElt() = default;
  NNStackElt(const NNStackElt&) = default;
  NNStackElt(int32_t nodeIdx, int16_t firstChildIdx)
    : firstChildIdx(firstChildIdx), childIdx(0), nodeIdx(nodeIdx) {}
  int32_t nodeIdx;
  int16_t childIdx;
  int16_t firstChildIdx;
};

inline
int
MSOctree::iNearestNeighbour_updateDMax(const point_t& pos, int32_t& d_max) const {
  std::array<NNStackElt,32> stack;
  int stack_last = -1;

  const int depthMax = depth;
  int32_t nodeIdx = 0;
  const MSONode* node = &nodes[nodeIdx];
  int currDepth = 0;
  const point_t posOf = pos - offsetOrigin;
  if (  (posOf[0] >= 0) & (posOf[0] < 1 << maxDepth)
      & (posOf[1] >= 0) & (posOf[1] < 1 << maxDepth)
      & (posOf[2] >= 0) & (posOf[2] < 1 << maxDepth)
  ) {
    int pointChildMask = 1 << maxDepth - 1;
    for (; currDepth < depthMax; ++currDepth) {
      const int childIdx
        = (!!((posOf[2]) & pointChildMask))
        | (!!((posOf[1]) & pointChildMask) << 1)
        | (!!((posOf[0]) & pointChildMask) << 2);

      if (!node->child[childIdx])
        break;

      stack[++stack_last] = NNStackElt(nodeIdx, childIdx);
      nodeIdx = node->child[childIdx];
      node = &nodes[nodeIdx];
      pointChildMask >>= 1;
    }
  }
  for (; currDepth < depthMax; ++currDepth) {
    int d_min = INT32_MAX;
    int i_min;
    for(int i = 0; i < 8; ++i)
      if (node->child[i]) {
        const MSONode& child = nodes[node->child[i]];

        const auto dPos0 = child.pos0 - pos;
        const auto dPos1 = child.sizeMinus1 + dPos0;
        int local_d_min
          = (dPos0[0] > 0 ? dPos0[0] : 0)
          + (dPos0[1] > 0 ? dPos0[1] : 0)
          + (dPos0[2] > 0 ? dPos0[2] : 0)
          - (dPos1[0] < 0 ? dPos1[0] : 0)
          - (dPos1[1] < 0 ? dPos1[1] : 0)
          - (dPos1[2] < 0 ? dPos1[2] : 0);

        if (local_d_min < d_min) {
          d_min = local_d_min;
          i_min = i;
        }
      }
    const int childNodeIdx = node->child[i_min];
    stack[++stack_last] = NNStackElt(nodeIdx, i_min);
    nodeIdx = childNodeIdx;
    node = &nodes[childNodeIdx];
  }

  int32_t local_d_max = INT32_MAX;

  int nearestIdx = node->start;

  for (int i = node->start; i < node->end; ++i) {
    auto dPoint = pos - (*pointCloud)[i];
    int32_t d
      = std::abs(dPoint[0])
      + std::abs(dPoint[1])
      + std::abs(dPoint[2]);
    if (d < local_d_max) {
      local_d_max = d;
      nearestIdx = i;
    }
  }

  if (local_d_max < d_max)
    d_max = local_d_max;

  if (!local_d_max)
    return nearestIdx;

  while (stack_last >= 0) {
    auto& sn = stack[stack_last];

    const MSONode& node = nodes[sn.nodeIdx];

    while (sn.childIdx < 8
        && (!node.child[sn.childIdx] || sn.childIdx == sn.firstChildIdx))
      ++sn.childIdx;

    if (sn.childIdx < 8) {
      const auto childNodeIdx = node.child[sn.childIdx];
      const MSONode& child = nodes[childNodeIdx];

      const auto dPos0 = child.pos0 - pos;
      const auto dPos1 = child.sizeMinus1 + dPos0;
      const auto d_min
        = (dPos0[0] > 0 ? dPos0[0] : 0)
        + (dPos0[1] > 0 ? dPos0[1] : 0)
        + (dPos0[2] > 0 ? dPos0[2] : 0)
        - (dPos1[0] < 0 ? dPos1[0] : 0)
        - (dPos1[1] < 0 ? dPos1[1] : 0)
        - (dPos1[2] < 0 ? dPos1[2] : 0);

      if (d_min < local_d_max) {
        if (stack_last < depthMax-1) {
          stack[++stack_last] = NNStackElt(childNodeIdx, -1);
        }
        else {
          for (int i = child.start; i < child.end; ++i) {
            auto dPoint = pos - (*pointCloud)[i];
            int32_t d
              = std::abs(dPoint[0])
              + std::abs(dPoint[1])
              + std::abs(dPoint[2]);
            if (d < local_d_max) {
              local_d_max = d;
              nearestIdx = i;
            }
          }
        }
      }
      ++sn.childIdx;
    } else {
      --stack_last;
    }
  }

  if (local_d_max < d_max)
    d_max = local_d_max;

  return nearestIdx;
}

//----------------------------------------------------------------------------

inline
int
MSOctree::iApproximateNearestNeighbour_updateDMax(const point_t& pos, int32_t& d_max) const {
  const int depthMax = depth;
  int32_t nodeIdx = 0;
  const MSONode* node = &nodes[nodeIdx];
  int currDepth = 0;
  const point_t posOf = pos - offsetOrigin;
  if (  (posOf[0] >= 0) & (posOf[0] < 1 << maxDepth)
      & (posOf[1] >= 0) & (posOf[1] < 1 << maxDepth)
      & (posOf[2] >= 0) & (posOf[2] < 1 << maxDepth)
  ) {
    int pointChildMask = 1 << maxDepth - 1;
    for (; currDepth < depthMax; ++currDepth) {
      const int childIdx
        = (!!((posOf[2]) & pointChildMask))
        | (!!((posOf[1]) & pointChildMask) << 1)
        | (!!((posOf[0]) & pointChildMask) << 2);

      if (!node->child[childIdx])
        break;

      nodeIdx = node->child[childIdx];
      node = &nodes[nodeIdx];
      pointChildMask >>= 1;
    }
  }
  for (; currDepth < depthMax; ++currDepth) {
    int d_min = INT32_MAX;
    int i_min;
    for(int i = 0; i < 8; ++i)
      if (node->child[i]) {
        const MSONode& child = nodes[node->child[i]];

        const auto dPos0 = child.pos0 - pos;
        const auto dPos1 = child.sizeMinus1 + dPos0;
        int local_d_min
          = (dPos0[0] > 0 ? dPos0[0] : 0)
          + (dPos0[1] > 0 ? dPos0[1] : 0)
          + (dPos0[2] > 0 ? dPos0[2] : 0)
          - (dPos1[0] < 0 ? dPos1[0] : 0)
          - (dPos1[1] < 0 ? dPos1[1] : 0)
          - (dPos1[2] < 0 ? dPos1[2] : 0);

        if (local_d_min < d_min) {
          d_min = local_d_min;
          i_min = i;
        }
      }
    const int childNodeIdx = node->child[i_min];
    nodeIdx = childNodeIdx;
    node = &nodes[childNodeIdx];
  }

  int32_t local_d_max = INT32_MAX;

  int nearestIdx = node->start;

  for (int i = node->start; i < node->end; ++i) {
    auto dPoint = pos - (*pointCloud)[i];
    int32_t d
      = std::abs(dPoint[0])
      + std::abs(dPoint[1])
      + std::abs(dPoint[2]);
    if (d < local_d_max) {
      local_d_max = d;
      nearestIdx = i;
    }
  }

  if (local_d_max < d_max)
    d_max = local_d_max;

  return nearestIdx;
}

//----------------------------------------------------------------------------

inline
int
MSOctree::iApproxNearestNeighbourAttr(const point_t& pos) const {
  const int depthMax = depth;
  int32_t nodeIdx = 0;
  const MSONode* node = &nodes[nodeIdx];
  int currDepth = 0;
  const point_t posOf = pos - offsetOrigin;
  if (  (posOf[0] >= 0) & (posOf[0] < 1 << maxDepth)
      & (posOf[1] >= 0) & (posOf[1] < 1 << maxDepth)
      & (posOf[2] >= 0) & (posOf[2] < 1 << maxDepth)
  ) {
    int pointChildMask = 1 << maxDepth - 1;
    for (; currDepth < depthMax; ++currDepth) {
      const int childIdx
        = (!!((posOf[2]) & pointChildMask))
        | (!!((posOf[1]) & pointChildMask) << 1)
        | (!!((posOf[0]) & pointChildMask) << 2);

      if (!node->child[childIdx])
        break;

      nodeIdx = node->child[childIdx];
      node = &nodes[nodeIdx];
      pointChildMask >>= 1;
    }
  }
  for (; currDepth < depthMax; ++currDepth) {
    int d_min = INT32_MAX;
    int i_min;
    for(int i = 0; i < 8; ++i)
      if (node->child[i]) {
        const MSONode& child = nodes[node->child[i]];

        const auto dPos0 = child.pos0 - pos;
        const auto dPos1 = child.sizeMinus1 + dPos0;
        int local_d_min
          = (dPos0[0] > 0 ? dPos0[0] : 0)
          + (dPos0[1] > 0 ? dPos0[1] : 0)
          + (dPos0[2] > 0 ? dPos0[2] : 0)
          - (dPos1[0] < 0 ? dPos1[0] : 0)
          - (dPos1[1] < 0 ? dPos1[1] : 0)
          - (dPos1[2] < 0 ? dPos1[2] : 0);

        if (local_d_min < d_min) {
          d_min = local_d_min;
          i_min = i;
        }
      }
    const int childNodeIdx = node->child[i_min];
    nodeIdx = childNodeIdx;
    node = &nodes[childNodeIdx];
  }

  int32_t d_max = INT32_MAX;

  int nearestIdx = node->start;

  for (int i = node->start; i < node->end; ++i) {
    auto dPoint = pos - (*pointCloud)[i];
    int32_t d
      = std::abs(dPoint[0])
      + std::abs(dPoint[1])
      + std::abs(dPoint[2]);
    if (d < d_max) {
      d_max = d;
      nearestIdx = i;
    }
  }

  return nearestIdx;
}

//----------------------------------------------------------------------------

double
MSOctree::find_motion(
  const EncodeMotionSearchParams& param,
  const ParameterSetMotion& mvPS,
  const MotionEntropyEstimate& motionEntropy,
  const PCCPointSet3& Block0,
  const point_t& xyz0,
  int local_size,
  MVField& mvField,
  uint32_t puNodeIdx) const // node Idx in mvField
{
  //if (!Window.size())
  //  return DBL_MAX;

  // ---------------------------- test no split --------------------
  // test V=0 and  dynamic planing at Block0-Window +/- window_size
  double cost_NoSplit = DBL_MAX;
  //int wSize = param.window_size;

  const point_t V00 = 0;

  //const int max_distance = 3 * wSize +10;//3 * wSize;
  point_t bestV_NoSplit = V00;
  point_t bestV_NoSplit_2nd = V00;

  int jumpBlock = 1 + (Block0.size() >> param.decimate);  // (kind of) random sampling of the original block to code

  int Dist = 0;

  // TODO: buffer vector difference or dmax to estimate dmax after motion is applied
  std::vector<int32_t> min_d0(1+Block0.size()/jumpBlock);
  std::vector<int32_t> min_d1(1+Block0.size()/jumpBlock);
  std::vector<int32_t> min_dK(1+Block0.size()/jumpBlock);
  std::vector<int32_t> min_dTmp(1+Block0.size()/jumpBlock);
  std::vector<int32_t> min_start;

  std::vector<int32_t> min_d0R;
  min_d0R.reserve(1+Block0.size()/jumpBlock);

  // initial vector
  point_t startMV = V00;
  point_t V0 = startMV;
  int NtestedPoints = 0;
  int NtestedPointsBlock0 = 0;
  startMV = {};
  for (int Nb = 0, idx = 0; Nb < Block0.size(); Nb += jumpBlock, ++idx) {
    int32_t min_d = INT32_MAX;
    auto p = Block0[Nb];
    int nearestPointIdx = nearestNeighbour_updateDMax(p + V0, min_d, param.approximate_nn);

    int dColor_forMinD = 0;
    if (Block0.hasColors()) {
      auto cW = pointCloud->getColor(nearestPointIdx);
      auto cB = Block0.getColor(Nb);
      dColor_forMinD += std::abs(cW[0] - cB[0]);
      dColor_forMinD += std::abs(cW[1] - cB[1]);
      dColor_forMinD += std::abs(cW[2] - cB[2]);
    }

    min_d0[idx] = min_d;

    NtestedPoints++;
    NtestedPointsBlock0++;
    startMV += (*pointCloud)[nearestPointIdx] - Block0[Nb];

    Dist += plus1log2shifted4(int(min_d + param.dgeom_color_factor * dColor_forMinD));  // 1/0.0625 = 16 times log
  }
  if (NtestedPoints)
    startMV /= NtestedPoints;

  double d = jumpBlock * Dist * 0.0625 + param.lambda * motionEntropy.estimateVector(V0);
  bestV_NoSplit = V0;
  bestV_NoSplit_2nd = V0;

  // set loop search parameters
  double best_d[2] = { d, d };
  std::set<int64_t> list_tested = { ((V0[0] + 32768LL) << 32) + ((V0[1] + 32768LL) << 16) + V0[2] + 32768LL };

  point_t VPrev = V00;
  V0 = startMV;
  min_dTmp = min_d0;
  min_d1 = min_d0;
  auto d_prev = std::numeric_limits<double>::max();
  for (int k = 0; k < param.K && d < d_prev; ++k) {
    d_prev = d;
    int64_t V1D = ((V0[0] + 32768LL) << 32) + ((V0[1] + 32768LL) << 16) + V0[2] + 32768LL;
    if(!list_tested.emplace(V1D).second)
      break;

    point_t meanV = {};
    NtestedPoints = 0;
    Dist = 0;
    for (int Nb = 0, idx = 0; Nb < Block0.size(); Nb += jumpBlock, ++idx) {
      auto p = Block0[Nb];
      auto offset = V0 - VPrev;
      int32_t min_d = min_dTmp[idx] + std::abs(offset[0]) + std::abs(offset[1]) + std::abs(offset[2]);
      int nearestPointIdx = nearestNeighbour_updateDMax(p + V0, min_d, param.approximate_nn);

      int dColor_forMinD = 0;
      if (Block0.hasColors()) {
        auto cW = pointCloud->getColor(nearestPointIdx);
        auto cB = Block0.getColor(Nb);
        dColor_forMinD += std::abs(cW[0] - cB[0]);
        dColor_forMinD += std::abs(cW[1] - cB[1]);
        dColor_forMinD += std::abs(cW[2] - cB[2]);
      }

      min_dTmp[idx] = min_d;

      NtestedPoints++;
      meanV += (*pointCloud)[nearestPointIdx] - Block0[Nb];

      Dist += plus1log2shifted4(int(min_d + param.dgeom_color_factor * dColor_forMinD));  // 1/0.0625 = 16 times log
    } // loop on points of block
    if (NtestedPoints)
      meanV /= NtestedPoints;

    d = jumpBlock * Dist * 0.0625 + param.lambda * motionEntropy.estimateVector(V0);

    if (d < best_d[0]) {
      best_d[1] = best_d[0];
      bestV_NoSplit_2nd = bestV_NoSplit;
      best_d[0] = d;
      bestV_NoSplit = V0;
      std::swap(min_d1, min_d0);
      min_d0 = min_dTmp;
    }
    else if (d < best_d[1]) {
      best_d[1] = d;
      bestV_NoSplit_2nd = V0;
      min_d1 = min_dTmp;
    }
    // start next step from meanV
    VPrev = V0;
    V0 = meanV;
  }

  const int searchPattern[3 * 18] = { 1,0,0,   0,1,0,  0,-1,0,  -1,0,0,  0,0,1, 0,0,-1,   1,1,0,     1,-1,0,  -1,1,0,  -1,-1,0,   0,1,1,    0,-1,1,  0,1,-1,  0,-1,-1,   -1,0,1,  1,0,-1,  -1,0,-1,  1,0,1 };
  //const int searchPattern[3 * 6] = {   1,0,0,   0,1,0,  0,-1,0,    -1,0,0,    0,0,1,   0,0,-1};
  //const int searchPattern[3 * 26] = { 1,0,0,   0,1,0,  0,-1,0,   -1,0,0,   0,0,1,   0,0,-1,  1,1,0,  1,-1,0,  -1,1,0,   -1,-1,0,  1,0,1,  0,1,1,  0,-1,1,  -1,0,1,  1,0,-1,  0,1,-1,  0,-1,-1,  -1,0,-1,  1,1,1,  1,1,-1,  1,-1,1,  1,-1,-1,  -1,1,1,  -1,1,-1,  -1,-1,1,  -1,-1,-1 };

  // loop MV search
  int see = 0;
  int Amotion = param.Amotion0;
  while (Amotion >= 1) {

    // pick one MV best seeds
    point_t Vs = see == 0 ? bestV_NoSplit : bestV_NoSplit_2nd;
    min_start = see == 0 ? min_d0 : min_d1;
    see = 1;

    // loop on searchPattern
    const int* pSearch = searchPattern;
    bool flagBetter = false;
    for (int t = 0; t < 18; t++, pSearch += 3) {

      point_t V = point_t(pSearch[0] * Amotion, pSearch[1] * Amotion, pSearch[2] * Amotion);
      V0 = Vs + V;

      int64_t V1D = ((V0[0] + 32768LL) << 32) + ((V0[1] + 32768LL) << 16) + V0[2] + 32768LL;
      if(!list_tested.emplace(V1D).second)
        continue;

      Dist = 0;
      for (int Nb = 0, idx = 0; Nb < Block0.size(); Nb += jumpBlock, ++idx) {
        auto p = Block0[Nb];
        auto offset = V;
        int32_t min_d = min_start[idx] + std::abs(offset[0]) + std::abs(offset[1]) + std::abs(offset[2]);
        int nearestPointIdx = nearestNeighbour_updateDMax(p + V0, min_d, param.approximate_nn);

        int dColor_forMinD = 0;
        if (Block0.hasColors()) {
          auto cW = pointCloud->getColor(nearestPointIdx);
          auto cB = Block0.getColor(Nb);
          dColor_forMinD += std::abs(cW[0] - cB[0]);
          dColor_forMinD += std::abs(cW[1] - cB[1]);
          dColor_forMinD += std::abs(cW[2] - cB[2]);
        }

        min_dK[idx] = min_d;

        Dist += plus1log2shifted4(int(min_d + param.dgeom_color_factor * dColor_forMinD));  // 1/0.0625 = 16 times log
      } // loop on points of block

      d = jumpBlock * Dist * 0.0625 + param.lambda * motionEntropy.estimateVector(V0);

      // keep 2 best MV
      if (d < best_d[0]) {
        best_d[1] = best_d[0];
        bestV_NoSplit_2nd = bestV_NoSplit;
        best_d[0] = d;
        bestV_NoSplit = V0;
        see = 0;
        flagBetter = true;
        std::swap(min_d1,min_d0);
        std::swap(min_d0,min_dK);

        break;
      }
      else if (d < best_d[1]) {
        best_d[1] = d;
        bestV_NoSplit_2nd = V0;
        std::swap(min_d1,min_dK);
      }

    }  // end loop on searchPattern

    // log reduction of search range
    if (!flagBetter && see == 1) {
      Amotion >>= 1;
      if (!bestV_NoSplit[0] && !bestV_NoSplit[1] && !bestV_NoSplit[2]) {
        bool flag = Amotion > 1;
        Amotion >>= 1;
        if (flag)
          Amotion = std::max(Amotion, 1);
      }
    }

  }  // end loop MV search

  cost_NoSplit = best_d[0];

  auto node = &mvField.puNodes[puNodeIdx];

  // ---------------------------- test split --------------------
  auto numPUsBeforeSplit = mvField.puNodes.size();
  auto numMVsBeforeSplit = mvField.mvPool.size();

  double cost_Split = DBL_MAX;

  if (local_size > mvPS.motion_min_pu_size && Block0.size() >= 8) {
    // condition on number of points for search acceleration
    int local_size1 = local_size >> 1;

    std::array<point_t, 8> list_xyz = {
      point_t(0, 0, 0) + xyz0,
      point_t(0, 0, local_size1) + xyz0,
      point_t(0, local_size1, 0) + xyz0,
      point_t(0, local_size1, local_size1) + xyz0,
      point_t(local_size1, 0, 0) + xyz0,
      point_t(local_size1, 0, local_size1) + xyz0,
      point_t(local_size1, local_size1, 0) + xyz0,
      point_t(local_size1, local_size1, local_size1) + xyz0
    };

    // loop on 8 child PU
    cost_Split = 0.;

    PCCPointSet3 Block1;
    Block1.reserve(Block0.size());

    // NOTE: points shall be already ordered
    //       and all belonging to the current node
    std::array<int32_t, 8> childCounts = {};
    // child idx
    int childIdx = 0;
    // child PU coordinates
    point_t xyz1 = list_xyz[childIdx];
    // block for child PU
    point_t xyz1High = xyz1 + local_size1;
    // TODO: might be simplified / accelerated
    // by dichotomic search or suited comparizon according to index
    // => childIdx outside and index in block0
    for (const auto& b : Block0) {
      while (childIdx < 8 && (
        b[2] < xyz1[2] || b[2] >= xyz1High[2]
        || b[1] < xyz1[1] || b[1] >= xyz1High[1]
        || b[0] < xyz1[0] || b[0] >= xyz1High[0]
      )) {
        ++childIdx;
        xyz1 = list_xyz[childIdx];
        xyz1High = xyz1 + local_size1;
      }
      ++childCounts[childIdx];
    }

    node->_firstChildIdx = numPUsBeforeSplit;
    node->_childsMask = 0;
    for (int childIdx = 0; childIdx < 8; childIdx++) {
      if(childCounts[childIdx]) {
        node->_childsMask += 1 << childIdx;
        // add node for current pu
        mvField.puNodes.emplace_back();
        auto& childNode = mvField.puNodes.back();
        // address may have changed with emplace_back()
        node = &mvField.puNodes[puNodeIdx];
        xyz1 = list_xyz[childIdx];
        childNode.set_pos0(xyz1);
        childNode._puSizeLog2 = node->_puSizeLog2 - 1;
      }
    }
    assert(node->_childsMask);

    int childStart = 0;
    uint32_t childNodeIdx = node->_firstChildIdx;
    for (int childIdx = 0; childIdx < 8; childIdx++) {
      cost_Split += 1.0;  // the cost due to not coding the occupancy with inter pred
      if(!childCounts[childIdx]) {  // empty PU
        continue;
      }
      Block1.resize(0);
      Block1.appendPartition(
        Block0, childStart, childStart + childCounts[childIdx]);
      childStart += childCounts[childIdx];

      xyz1 = list_xyz[childIdx];

      cost_Split += find_motion(
        param, mvPS, motionEntropy, Block1, xyz1, local_size1, mvField,
        childNodeIdx);

      ++childNodeIdx;
    }
  }

  // ---------------------------- choose split vs no split --------------------
  if (local_size > mvPS.motion_min_pu_size) {
    // cost no split flag
    cost_NoSplit += param.lambda * motionEntropy.estimateSplit(false);
  }
  // cost split flag
  cost_Split += param.lambda * motionEntropy.estimateSplit(true);

  if (local_size <= mvPS.motion_min_pu_size || cost_NoSplit <= cost_Split) {  // no split
    mvField.puNodes.resize(numPUsBeforeSplit);
    mvField.mvPool.resize(numMVsBeforeSplit + 1);
    // address may have changed
    node = &mvField.puNodes[puNodeIdx];
    //node._firstChildIdx = 0;
    node->_mvIdx = numMVsBeforeSplit;
    node->_childsMask = 0; // not split
    mvField.mvPool.back() = bestV_NoSplit;
    return cost_NoSplit;
  }
  else {
    return cost_Split;
  }
}

//----------------------------------------------------------------------------

bool
motionSearchForNode(
  const MSOctree& mSOctreeOrig,
  const MSOctree& mSOctree,
  const PCCOctree3Node* node0,
  const EncodeMotionSearchParams& param,
  const ParameterSetMotion& mvPS,
  int nodeSize,
  EntropyEncoder* arithmeticEncoder,
  MVField& mvField,
  uint32_t puNodeIdx, // node Idx in mvField
  bool flagNonPow2,
  int S,
  int S2)
{
  MotionEntropyEncoder motionEncoder(arithmeticEncoder);

  PCCPointSet3 Block0;
  Block0.appendPartition(*mSOctreeOrig.pointCloud, node0->start, node0->end);

  // entropy estimates
  int wSize = param.window_size;
  MotionEntropyEstimate mcEstimate(motionEncoder, wSize, param.max_prefix_bits, param.max_suffix_bits);

  // motion search
  Vec3<int32_t> pos = node0->pos * nodeSize;

  // TODO: scale point according to trisoup node size
  if (flagNonPow2) {
    int maskS = (1 << S2) - 1;
    for (int i = 0; i < Block0.size(); ++i) {
      Block0[i][0] = ((Block0[i][0] >> S2) * S) + (Block0[i][0] & maskS);
      Block0[i][1] = ((Block0[i][1] >> S2) * S) + (Block0[i][1] & maskS);
      Block0[i][2] = ((Block0[i][2] >> S2) * S) + (Block0[i][2] & maskS);
    }
  }

  auto& node = mvField.puNodes[puNodeIdx];
  node.set_pos0(pos);
  node._puSizeLog2 = ilog2(uint32_t(nodeSize - 1)) + 1; // TODO: check if we need the true size at some point or clean

  // MV search
  mSOctree.find_motion(
    param, mvPS, mcEstimate, Block0, pos, nodeSize, mvField, puNodeIdx);

  return true;
}

void
MSOctree::apply_motion(
  const point_t currNodePos0,
  const point_t currNodePos1,
  const point_t MVd,
  PCCOctree3Node* node0,
  PCCPointSet3* compensatedPointCloud,
  uint32_t depthMax,
  bool flagNonPow2,
  int S,
  int S2) const
{
  auto &fifo = a;
  assert(fifo.empty());
  fifo.clear();
  depthMax = std::min(depthMax, depth);
  const auto node0Pos0 = currNodePos0 + MVd;
  const auto node0Pos1 = currNodePos1 + MVd;
  const auto minNodeSizeMinus1 = (1 << maxDepth - depthMax) - 1;

  auto &local = b;
  assert(local.empty());
  local.clear();
  int addedPointCount = 0;

  fifo.push(0);
  while (!fifo.empty()) {
    const MSONode& node = nodes[fifo.front()];
    //
    const auto nodeSizeMinus1 = node.sizeMinus1;
    if (minNodeSizeMinus1 == nodeSizeMinus1)
      break;

    const auto nodePos0 = node.pos0;
    const auto nodePos1 = nodePos0 + nodeSizeMinus1;
    const auto dPos0 = nodePos0 - node0Pos0;
    const auto dPos1 = node0Pos1 - nodePos1;

    if ( (dPos0[0] | dPos0[1] | dPos0[2]
        | dPos1[0] | dPos1[1] | dPos1[2]) >= 0
    ) {
      local.push(fifo.front());
      addedPointCount += node.end - node.start;
    } else {
      int intersect = 0;
      for (int k=0; k < 3; ++k) {
        intersect |= std::min(nodePos1[k], node0Pos1[k]) - std::max(nodePos0[k], node0Pos0[k]);
      }
      if (intersect >= 0)
        for(int i = 0; i < 8; ++i)
          if (node.child[i])
            fifo.push(node.child[i]);
    }
    fifo.pop();
  }

  std::vector<int> indices;
  indices.reserve(pointCloud->size());
  indices.resize(addedPointCount);
  int i = 0;
  while (!local.empty()) {
    const MSONode& node = nodes[local.front()];
    for (int k=node.start; k<node.end; ++k)
      indices[i++] = k;
    local.pop();
  }

  while (!fifo.empty()) {
    const MSONode& node = nodes[fifo.front()];
    for (int i=node.start; i < node.end; ++i) {
      auto const & pt = (*pointCloud)[i];
      const auto dPos0 = pt - node0Pos0;
      const auto dPos1 = node0Pos1 - pt;
      if ( (dPos0[0] | dPos0[1] | dPos0[2]
          | dPos1[0] | dPos1[1] | dPos1[2]) >= 0)
        indices.push_back(i);
    }
    fifo.pop();
  }

  node0->predStart = compensatedPointCloud->size();
  compensatedPointCloud->appendPartition(*pointCloud, indices, false);
  node0->predEnd = compensatedPointCloud->size();

  for (int i = node0->predStart; i < node0->predEnd; ++i) {
    auto& predPoint = (*compensatedPointCloud)[i];
    predPoint[0] -= MVd[0];
    predPoint[1] -= MVd[1];
    predPoint[2] -= MVd[2];
  }
  // align points with octree nodes
  if (flagNonPow2) {
    // TODO:
    //  should we apply at octree node level to avoid division and multiplications ?
    const int factorS = (1 << S2) - S;
    for (int i = node0->predStart; i < node0->predEnd; ++i) {
      auto& predPoint = (*compensatedPointCloud)[i];

      int temp0 = predPoint[0] / S;
      int temp1 = predPoint[1] / S;
      int temp2 = predPoint[2] / S;
      predPoint[0] += temp0 * factorS;
      predPoint[1] += temp1 * factorS;
      predPoint[2] += temp2 * factorS;
    }
  }
}

void
MSOctree::apply_recolor_motion(
  point_t Mvd,
  PCCOctree3Node* node0,
  PCCPointSet3& pointCloud) const
{
  for (int i = node0->start; i < node0->end; ++i) {
    auto p = pointCloud[i];
    int nearestPointIdx = iApproxNearestNeighbourAttr(p + Mvd);
    pointCloud.setColor(i, this->pointCloud->getColor(nearestPointIdx));
 }
}

//----------------------------------------------------------------------------

int
MSOctree::nodeIdx(point_t nodePos0, uint32_t nodeSizeLog2) const
{
  uint32_t currSizeLog2 = maxDepth;
  uint32_t currNodeIdx = 0;

  while (currSizeLog2 > nodeSizeLog2) {
    const auto & node0 = nodes[currNodeIdx];
    const int childSizeLog2 = (currSizeLog2 - 1);
    int i
      = (((nodePos0[0] >> childSizeLog2)& 1) << 2)
      + (((nodePos0[1] >> childSizeLog2)& 1) << 1)
      + ((nodePos0[2] >> childSizeLog2)& 1);
    if (!node0.child[i])
      return -1;
    currNodeIdx = node0.child[i];
    currSizeLog2 = childSizeLog2;
  }
  return currNodeIdx;
}

//----------------------------------------------------------------------------
}  // namespace pcc
