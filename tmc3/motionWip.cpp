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
#include <vector>
#include <unordered_map>

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "geometry_octree.h"

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
  double hMvIsZero[2];
  double hMvIsOne[2];
  double hMvIsTwo[2];
  double hMvIsThree[2];
  double hExpGolombV[6][2];
  double hSplitPu[2];

  MotionEntropyEstimate(const MotionEntropy& codec, size_t wSize, int boundPrefix, int boundSuffix);

  double estimateVector(const point_t& mv) const;

private:
  double prepareEstimate(unsigned absval) const;
  void prepareEstimateVector();

  double estimateComponent(unsigned absval) const {
    if (absval < LUT_MVestimate.size())
      return LUT_MVestimate[absval];
    auto res = dyn_MVestimate.emplace(std::make_pair(absval, 0.));
    if (res.second)
      res.first->second = prepareEstimate(absval);
    return res.first->second;
  }
  int boundPrefix;
  int boundSuffix;
  std::vector<double> LUT_MVestimate;
  mutable std::map<int,double> dyn_MVestimate;
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
  codec.splitPu.getEntropy(hSplitPu);
  codec.mvIsZero.getEntropy(hMvIsZero);
  codec.mvIsOne.getEntropy(hMvIsOne);
  codec.mvIsTwo.getEntropy(hMvIsTwo);
  codec.mvIsThree.getEntropy(hMvIsThree);
  for (int i = 0; i < 6; i++)
    codec.expGolombV[i].getEntropy(hExpGolombV[i]);

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

double
MotionEntropyEstimate::prepareEstimate(unsigned absval) const
{
  int v = absval;

  if (!v) {
    return hMvIsZero[1];
  }
  else {
    double r = hMvIsZero[0] + 1.;  // unpredictable sign
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
      r += k;
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
  return estimateComponent(std::abs(mv[0])) + estimateComponent(std::abs(mv[1])) + estimateComponent(std::abs(mv[2]));
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
deriveMotionMaxPrefixBits(const GeometryParameterSet::Motion& param)
{
  int max_MV =
    std::max(0, param.motion_window_size - 1);
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
deriveMotionMaxSuffixBits(const GeometryParameterSet::Motion& param)
{

  int max_MV =
    std::max(0, param.motion_window_size  - 1);
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
extracPUsubtree(
  const GeometryParameterSet::Motion& param,
  PUtree* local_PU_tree,
  int block_size,
  int& pos_fs,
  int& pos_fp,
  int& pos_MV,
  PUtree* destination_tree)
{
  // non-split terminal case
  if (
    block_size <= param.motion_min_pu_size
    || !local_PU_tree->split_flags[pos_fs]) {
    if (block_size > param.motion_min_pu_size) {
      destination_tree->split_flags.push_back(0);
      pos_fs++;
    }

    destination_tree->MVs.push_back(local_PU_tree->MVs[pos_MV]);
    pos_MV++;
    return;
  }

  // split case
  destination_tree->split_flags.push_back(1);
  pos_fs++;

  // loop on 8 children
  for (int s = 0; s < 8; s++) {
    if (local_PU_tree->popul_flags[pos_fp]) {  // child populated
      destination_tree->popul_flags.push_back(1);
      pos_fp++;

      extracPUsubtree(
        param, local_PU_tree, block_size >> 1, pos_fs, pos_fp, pos_MV,
        destination_tree);
    } else {  // child not pouplated
      destination_tree->popul_flags.push_back(0);
      pos_fp++;
    }
  }
}


void
encode_splitPU_MV_MC(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  PUtree* local_PU_tree,
  const GeometryParameterSet::Motion& param,
  point_t nodeSizeLog2,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3* compensatedPointCloud,
  int numLPUPerLine,
  int log2MotionBlkSize,
  std::vector<MotionVector>& motionVectors)
{
  const int node_size = 1 << nodeSizeLog2[0];
  MotionEntropyEncoder motionEncoder(arithmeticEncoder);

  // --------------  non-split / terminal case  ----------------
  if (node_size <= param.motion_min_pu_size || !local_PU_tree->split_flags[0]) {
    if (node_size > param.motion_min_pu_size) {
      motionEncoder.encodeSplitPu(0);
    }

    // encode MV
    point_t MV = local_PU_tree->MVs[0];
    motionEncoder.encodeVector(MV);
    motionVectors.push_back({node0->pos * node_size, node_size, MV});

    point_t MVd = MV;

    mSOctree.apply_motion(MVd, node0, param, nodeSizeLog2[0], compensatedPointCloud, mSOctree.depth);
    node0->isCompensated = true;
    return;
  }

  // --------------- split case ----------------------
  motionEncoder.encodeSplitPu(1);
}

//----------------------------------------------------------------------------

void
decode_splitPU_MV_MC(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  const GeometryParameterSet::Motion& param,
  point_t nodeSizeLog2,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3* compensatedPointCloud,
  int numLPUPerLine,
  int log2MotionBlkSize,
  std::vector<MotionVector>& motionVectors)
{
  int node_size = 1 << nodeSizeLog2[0];
  MotionEntropyDecoder motionDecoder(arithmeticDecoder);

  // decode split flag
  bool split = false;
  if (node_size > param.motion_min_pu_size)
    split = motionDecoder.decodeSplitPu();

  if (!split) {  // not split
                 // decode MV
    point_t MV = 0;
    point_t MVd = 0.;
    motionDecoder.decodeVector(&MV);
    MVd = MV;

    motionVectors.push_back({node0->pos * node_size, node_size, MV});

    mSOctree.apply_motion(MVd, node0, param, nodeSizeLog2[0], compensatedPointCloud, mSOctree.depth);
    node0->isCompensated = true;
    return;
  }

  // split; nothing to do
}

//============================================================================

MSOctree::MSOctree(
    PCCPointSet3* predPointCloud,
    point_t offsetOrigin,
    uint32_t leafSizeLog2
  )
  : pointCloud(predPointCloud)
  , offsetOrigin(offsetOrigin)
{
  PCCPointSet3& pointCloud(*this->pointCloud);
  std::vector<MSONode> nodesNext;
  nodes.reserve(pointCloud.size());
  nodesNext.reserve(pointCloud.size());

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
  node00.sizeMinus1 = int32_t((1 << maxDepth) - 1);

  uint32_t nodesCurrNode = 0;
  for (int currDepth = 0; currDepth < depth; ++currDepth) {
    const uint32_t nodesCurrLvlEnd = nodes.size();
    const int childSizeLog2 = maxDepth - currDepth - 1;
    const int pointSortMask = 1 << childSizeLog2;
    const int childSizeMinus1 = pointSortMask - 1;

    for (; nodesCurrNode != nodesCurrLvlEnd; ++nodesCurrNode) {
      MSONode& node0 = nodes[nodesCurrNode];

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

      uint32_t childStart = node0.start;
      auto pos0 = node0.pos0;
      for (int i = 0; i < 8; ++i) {
        if (childCounts[i]) {
          node0.child[i] = nodesCurrLvlEnd + nodesNext.size();
          nodesNext.emplace_back();
          MSONode& child0 = nodesNext.back();
          uint32_t childEnd = childStart + childCounts[i];
          child0.start = childStart;
          child0.end = childEnd;
          childStart = childEnd;
          child0.pos0 = pos0 + (Vec3<int32_t>{i >> 2, (i >> 1) & 1, i & 1} << childSizeLog2);
          child0.sizeMinus1 = childSizeMinus1;
          child0.parent = nodesCurrNode;
        }
      }
    }
    nodes.insert(
      nodes.end(), nodesNext.begin(), nodesNext.end());
    nodesNext.clear();
  }

  // now everything is built, we can offset the points
  for (int i=0; i < pointCloud.size(); ++i) {
    pointCloud[i] += offsetOrigin;
  }
}

//----------------------------------------------------------------------------

int32_t
MSOctree::nearestNeighbour_estimateDMax(point_t pos, int32_t d_max, uint32_t depthMax) const {
  depthMax = std::min(depth, depthMax);

  const MSONode* node = &nodes[0];
  for (int currDepth = 0; currDepth < depthMax; ++currDepth) {
    const int childSizeLog2 = maxDepth - currDepth - 1;
    const int pointChildMask = 1 << childSizeLog2;
    const int childIdx
      = (!!((pos[2]+offsetOrigin[2]) & pointChildMask))
      | (!!((pos[1]+offsetOrigin[1]) & pointChildMask) << 1)
      | (!!((pos[0]+offsetOrigin[0]) & pointChildMask) << 2);

    if (!node->child[childIdx])
      break;

    node = &nodes[node->child[childIdx]];
  }

  const auto dPos0 = pos - node->pos0;
  const auto dPos1 = dPos0 - node->sizeMinus1;
  int32_t local_d_max
    = std::max(std::abs(dPos0[0]), std::abs(dPos1[0]))
    + std::max(std::abs(dPos0[1]), std::abs(dPos1[1]))
    + std::max(std::abs(dPos0[2]), std::abs(dPos1[2]));

  if (local_d_max < d_max)
    d_max = local_d_max;

  return d_max;
}

//----------------------------------------------------------------------------

std::tuple<int, int, int>
MSOctree::nearestNeighbour(point_t pos, int32_t d_max, uint32_t depthMax) const {
  fifo.push(0);
  int32_t d_min = d_max;
  int32_t local_d_max;
  int nearest_node_idx = -1;
  int nearest_point_idx = -1;
  const int32_t minSizeMinus1 = (1 << maxDepth - std::min(depth, depthMax)) - 1;

  // int currsizeMinus1 = nodes[fifo.front()].sizeMinus1;
  while (!fifo.empty() && nodes[fifo.front()].sizeMinus1 != minSizeMinus1) {
    const MSONode& node = nodes[fifo.front()];

    const auto dPos0 = node.pos0 - pos;
    const auto dPos1 = node.sizeMinus1 + dPos0;
    d_min
      = (dPos0[0] > 0 ? dPos0[0] : 0)
      + (dPos0[1] > 0 ? dPos0[1] : 0)
      + (dPos0[2] > 0 ? dPos0[2] : 0)
      - (dPos1[0] < 0 ? dPos1[0] : 0)
      - (dPos1[1] < 0 ? dPos1[1] : 0)
      - (dPos1[2] < 0 ? dPos1[2] : 0);

    if (d_min <= d_max) {
      /*if (d_min) {
        // point outside
        int d_child[8] = {};
        if (dPos0[0] > 0) {
          d_child[0]++;d_child[1]++;d_child[2]++;d_child[3]++;
        } else if (dPos1[0] < 0) {
          d_child[4]++;d_child[5]++;d_child[6]++;d_child[7]++;
        }
        if (dPos0[1] > 0) {
          d_child[0]++;d_child[1]++;d_child[4]++;d_child[5]++;
        } else if (dPos1[1] < 0) {
          d_child[2]++;d_child[3]++;d_child[6]++;d_child[7]++;
        }
        if (dPos0[2] > 0) {
          d_child[0]++;d_child[2]++;d_child[4]++;d_child[6]++;
        } else if (dPos1[2] < 0) {
          d_child[1]++;d_child[3]++;d_child[5]++;d_child[7]++;
        }
        int _max = 0;
        for(int i = 0; i < 8; ++i)
          if (node.child[i] && d_child[i] > _max)
            _max = d_child[i];
        for(int i = 0; i < 8; ++i)
          if (node.child[i] && d_child[i] >= _max-1)
            fifo.push(node.child[i]);
      }
      else*/
      for(int i = 0; i < 8; ++i)
        if (node.child[i])
          fifo.push(node.child[i]);

      local_d_max
        = std::max(std::abs(dPos0[0]), std::abs(dPos1[0]))
        + std::max(std::abs(dPos0[1]), std::abs(dPos1[1]))
        + std::max(std::abs(dPos0[2]), std::abs(dPos1[2]));

      d_max = local_d_max < d_max ? local_d_max : d_max;
    }

    fifo.pop();
  }
  local_d_max = d_max + 1; // to get point if d_max = min(d)
  while(!fifo.empty()) {
    const MSONode& node = nodes[fifo.front()];

    for (int i = node.start; i < node.end; ++i) {
      auto dPoint = pos - (*pointCloud)[i];
      int32_t d
        = std::abs(dPoint[0])
        + std::abs(dPoint[1])
        + std::abs(dPoint[2]);
      if (d < local_d_max) {
        local_d_max = d;
        nearest_node_idx = fifo.front();
        nearest_point_idx = i;
      }
    }
    fifo.pop();
  }
  if (local_d_max < d_max)
    d_max = local_d_max;
  return std::make_tuple(nearest_node_idx, nearest_point_idx, int(d_max));
}

//----------------------------------------------------------------------------

std::tuple<std::queue<uint32_t>, int>
MSOctree::nearestNodes(point_t node0Pos0, int32_t d_max, uint32_t node0SizeLog2) const {
  fifo.push(0);
  int32_t d_min = d_max;
  int32_t local_d_max;
  const int32_t minSizeMinus1 = (1 << std::max(node0SizeLog2, maxDepth - depth)) - 1;

  // int currsizeMinus1 = nodes[fifo.front()].sizeMinus1;
  while (!fifo.empty() && nodes[fifo.front()].sizeMinus1 > minSizeMinus1) {
    const MSONode& node = nodes[fifo.front()];

    const auto dPos0 = node.pos0 - node0Pos0;
    const auto dPos1 = node.sizeMinus1 + dPos0;
    d_min
      = (dPos0[0] > 0 ? dPos0[0] : 0)
      + (dPos0[1] > 0 ? dPos0[1] : 0)
      + (dPos0[2] > 0 ? dPos0[2] : 0)
      - (dPos1[0] < 0 ? dPos1[0] : 0)
      - (dPos1[1] < 0 ? dPos1[1] : 0)
      - (dPos1[2] < 0 ? dPos1[2] : 0);

    if (d_min <= d_max) {
      for(int i = 0; i < 8; ++i)
        if (node.child[i])
          fifo.push(node.child[i]);

      local_d_max
        = std::max(std::abs(dPos0[0]), std::abs(dPos1[0]))
        + std::max(std::abs(dPos0[1]), std::abs(dPos1[1]))
        + std::max(std::abs(dPos0[2]), std::abs(dPos1[2]));

      d_max = local_d_max < d_max ? local_d_max : d_max;
    }

    fifo.pop();
  }
  std::queue<std::tuple<int32_t,uint32_t>> intermediate;
  while(!fifo.empty()) {
    const MSONode& node = nodes[fifo.front()];

    const auto dPos0 = node.pos0 - node0Pos0;
    const auto dPos1 = node.sizeMinus1 + dPos0;
    d_min
      = (dPos0[0] > 0 ? dPos0[0] : 0)
      + (dPos0[1] > 0 ? dPos0[1] : 0)
      + (dPos0[2] > 0 ? dPos0[2] : 0)
      - (dPos1[0] < 0 ? dPos1[0] : 0)
      - (dPos1[1] < 0 ? dPos1[1] : 0)
      - (dPos1[2] < 0 ? dPos1[2] : 0);

    if (d_min <= d_max) {
      intermediate.push(std::make_tuple(int32_t(d_min), uint32_t(fifo.front())));

      local_d_max
        = std::max(std::abs(dPos0[0]), std::abs(dPos1[0]))
        + std::max(std::abs(dPos0[1]), std::abs(dPos1[1]))
        + std::max(std::abs(dPos0[2]), std::abs(dPos1[2]));

      d_max = local_d_max < d_max ? local_d_max : d_max;
    }

    fifo.pop();
  }
  std::queue<uint32_t> nearest_nodes_idx;
  while(!intermediate.empty()) {
    int d_min, node_idx;
    std::tie(d_min, node_idx) = intermediate.front();
    if (d_min <= d_max)
      nearest_nodes_idx.push(node_idx);
    intermediate.pop();
  }
  return std::make_tuple(nearest_nodes_idx, int(d_max));
}

//----------------------------------------------------------------------------

std::tuple<uint32_t, int>
MSOctree::nearestNode(point_t node0Pos0, int32_t d_max, uint32_t node0SizeLog2) const {
  fifo.push(0);
  int32_t d_min = d_max;
  uint32_t nearest_node_idx = 0;
  int32_t local_d_max;
  const int32_t minSizeMinus1 = (1 << std::max(node0SizeLog2, maxDepth - depth)) - 1;

  // int currsizeMinus1 = nodes[fifo.front()].sizeMinus1;
  while (!fifo.empty() && nodes[fifo.front()].sizeMinus1 > minSizeMinus1) {
    const MSONode& node = nodes[fifo.front()];

    const auto dPos0 = node.pos0 - node0Pos0;
    const auto dPos1 = node.sizeMinus1 + dPos0;
    d_min
      = (dPos0[0] > 0 ? dPos0[0] : 0)
      + (dPos0[1] > 0 ? dPos0[1] : 0)
      + (dPos0[2] > 0 ? dPos0[2] : 0)
      - (dPos1[0] < 0 ? dPos1[0] : 0)
      - (dPos1[1] < 0 ? dPos1[1] : 0)
      - (dPos1[2] < 0 ? dPos1[2] : 0);

    if (d_min <= d_max) {
      for(int i = 0; i < 8; ++i)
        if (node.child[i])
          fifo.push(node.child[i]);

      local_d_max
        = std::max(std::abs(dPos0[0]), std::abs(dPos1[0]))
        + std::max(std::abs(dPos0[1]), std::abs(dPos1[1]))
        + std::max(std::abs(dPos0[2]), std::abs(dPos1[2]));

      d_max = local_d_max < d_max ? local_d_max : d_max;
    }

    fifo.pop();
  }

  local_d_max = d_max + 1;

  while(!fifo.empty()) {
    const MSONode& node = nodes[fifo.front()];

    const auto dPos0 = node.pos0 - node0Pos0;
    const auto dPos1 = node.sizeMinus1 + dPos0;

    int32_t d
      = std::max(std::abs(dPos0[0]), std::abs(dPos1[0]))
      + std::max(std::abs(dPos0[1]), std::abs(dPos1[1]))
      + std::max(std::abs(dPos0[2]), std::abs(dPos1[2]));

    nearest_node_idx
      = d < local_d_max ? fifo.front() : nearest_node_idx;

    local_d_max = d < local_d_max ? d : local_d_max;

    fifo.pop();
  }
  d_max = local_d_max < d_max ? local_d_max : d_max;
  return std::make_tuple(nearest_node_idx, int(d_max));
}

//----------------------------------------------------------------------------

double
MSOctree::find_motion(
  const GeometryParameterSet::Motion& param,
  const MotionEntropyEstimate& motionEntropy,
  const MSOctree& mSOctreeOrig,
  uint32_t mSOctreeOrigNodeIdx,
  const PCCPointSet3& Block0,
  const point_t& xyz0,
  int local_size,
  PUtree* local_PU_tree) const
{
  //if (!Window.size())
  //  return DBL_MAX;

  // ---------------------------- test no split --------------------
  // test V=0 and  dynamic planing at Block0-Window +/- window_size
  double cost_NoSplit = DBL_MAX;
  //int wSize = param.motion_window_size;

  const point_t V00 = 0;

  //const int max_distance = 3 * wSize +10;//3 * wSize;
  point_t bestV_NoSplit = V00;
  point_t bestV_NoSplit_2nd = V00;

  int jumpBlock = 1 + (Block0.size() >> param.decimate);  // (kind of) random sampling of the original block to code

  int Dist = 0;
  double DistR = 0;
  int DColorBlock=0;

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
    int32_t d_max = INT32_MAX;
    auto p = Block0[Nb];
    if(idx > 0) {
      auto offset = Block0[Nb] - Block0[Nb - jumpBlock];
      d_max = min_d0[idx - 1] + std::abs(offset[0]) + std::abs(offset[1]) + std::abs(offset[2]);
    }
    else {
      d_max = nearestNeighbour_estimateDMax(p + V0, d_max);
    }
    int nearestNodeIdx, nearestPointIdx, min_d;
    std::tie(nearestNodeIdx, nearestPointIdx, min_d) = nearestNeighbour(p + V0, d_max, depth);
    assert(nearestNodeIdx >= 0 && nearestPointIdx >= 0);

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
    DistR = 0;
    DColorBlock=0;
    for (int Nb = 0, idx = 0; Nb < Block0.size(); Nb += jumpBlock, ++idx) {
      int32_t d_max = INT32_MAX;
      auto p = Block0[Nb];
      auto offset = V0 - VPrev;
      d_max = min_dTmp[idx] + std::abs(offset[0]) + std::abs(offset[1]) + std::abs(offset[2]);
      if(idx > 0) {
        auto offset = Block0[Nb] - Block0[Nb - jumpBlock];
        d_max = std::min(d_max, min_dTmp[idx - 1] + std::abs(offset[0]) + std::abs(offset[1]) + std::abs(offset[2]));
      }
      else {
        d_max = std::min(d_max, nearestNeighbour_estimateDMax(p + V0, d_max));
      }
      int nearestNodeIdx, nearestPointIdx, min_d;
      std::tie(nearestNodeIdx, nearestPointIdx, min_d) = nearestNeighbour(p + V0, d_max, depth);
      assert(nearestNodeIdx >= 0 && nearestPointIdx >= 0);

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
      DistR = 0;
      DColorBlock=0;
      for (int Nb = 0, idx = 0; Nb < Block0.size(); Nb += jumpBlock, ++idx) {
        int32_t d_max = INT32_MAX;
        auto p = Block0[Nb];
        auto offset = V;
        d_max = min_start[idx] + std::abs(offset[0]) + std::abs(offset[1]) + std::abs(offset[2]);
        if(idx > 0) {
          auto offset = Block0[Nb] - Block0[Nb - jumpBlock];
          d_max = std::min(d_max, min_dK[idx - 1] + std::abs(offset[0]) + std::abs(offset[1]) + std::abs(offset[2]));
        }
        else {
          d_max = std::min(d_max, nearestNeighbour_estimateDMax(p + V0, d_max));
        }
        int nearestNodeIdx, nearestPointIdx, min_d;
        std::tie(nearestNodeIdx, nearestPointIdx, min_d) = nearestNeighbour(p + V0, d_max, depth);
        assert(nearestNodeIdx >= 0 && nearestPointIdx >= 0);

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


  // ---------------------------- test split --------------------
  double cost_Split = DBL_MAX;
  PUtree* Split_PU_tree = new PUtree;  // local split tree

  if (local_size > param.motion_min_pu_size && Block0.size() >= 8) {
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
    std::vector<int32_t> blockIndices;
    blockIndices.reserve(Block0.size());

    for (int t = 0; t < 8; t++) {
      // child PU coordinates
      point_t xyz1 = list_xyz[t];
      // block for child PU
      point_t xyz1High = xyz1 + local_size1;
      blockIndices.resize(0);
      for (const auto& b : Block0) {
        // TODO: try to use formula in nearestNeighbour() to avoid
        // comparisons and branching operations
        if (b[0] >= xyz1[0] && b[0] < xyz1High[0]
            && b[1] >= xyz1[1] && b[1] < xyz1High[1]
            && b[2] >= xyz1[2] && b[2] < xyz1High[2])
          blockIndices.push_back(b.getIndex());
      }
      cost_Split += 1.0;  // the cost due to not coding the occupancy with inter pred

      if (!blockIndices.size()) {  // empty PU
        Split_PU_tree->popul_flags.push_back(0);
        continue;
      }
      Split_PU_tree->popul_flags.push_back(1);

      Block1.resize(0);
      Block1.appendPartition(Block0,blockIndices);

      uint32_t childNodeIdx = mSOctreeOrig.nodes[mSOctreeOrigNodeIdx].child[t];

      cost_Split += find_motion(param, motionEntropy, mSOctreeOrig, childNodeIdx, Block1, xyz1, local_size1, Split_PU_tree);
    }
  }

  // ---------------------------- choose split vs no split --------------------
  if (local_size > param.motion_min_pu_size) {
    cost_NoSplit +=
      param.lambda * motionEntropy.hSplitPu[0];  // cost no split flag
  }
  cost_Split += param.lambda * motionEntropy.hSplitPu[1];  // cost split flag

  if (local_size <= param.motion_min_pu_size|| cost_NoSplit <= cost_Split) {  // no split
    // push non split flag, only if size>size_min
    if (local_size > param.motion_min_pu_size) {
      local_PU_tree->split_flags.push_back(0);
    }
    // push MV
    local_PU_tree->MVs.push_back(bestV_NoSplit);

    delete Split_PU_tree;
    return cost_NoSplit;
  }
  else {
    // split
    local_PU_tree->split_flags.push_back(1);  // push split PU flag

    // append Split_PU_tree to  local_PU_tree
    for (const auto& f : Split_PU_tree->popul_flags)
      local_PU_tree->popul_flags.push_back(f);
    for (const auto& f : Split_PU_tree->split_flags)
      local_PU_tree->split_flags.push_back(f);
    for (const auto& v : Split_PU_tree->MVs)
      local_PU_tree->MVs.push_back(v);

    delete Split_PU_tree;
    return cost_Split;
  }
}

//----------------------------------------------------------------------------

bool
motionSearchForNode(
  const MSOctree& mSOctreeOrig,
  const MSOctree& mSOctree,
  const PCCOctree3Node* node0,
  const GeometryParameterSet::Motion& param,
  int nodeSizeLog2,
  EntropyEncoder* arithmeticEncoder,
  PUtree* local_PU_tree)
{
  MotionEntropyEncoder motionEncoder(arithmeticEncoder);

  PCCPointSet3 Block0;
  Block0.appendPartition(*mSOctreeOrig.pointCloud, node0->start, node0->end);

  // entropy estimates
  int wSize = param.motion_window_size;
  MotionEntropyEstimate mcEstimate(motionEncoder, wSize, param.motion_max_prefix_bits, param.motion_max_suffix_bits);

  // motion search
  Vec3<int32_t> pos = node0->pos << nodeSizeLog2;
  int mSOctreeOrigNodeIdx = node0->mSOctreeNodeIdx;

  // MV search
  mSOctree.find_motion(
    param, mcEstimate, mSOctreeOrig, mSOctreeOrigNodeIdx, Block0, pos, (1 << nodeSizeLog2), local_PU_tree);

  return true;
}

void
MSOctree::apply_motion(
  point_t MVd,
  PCCOctree3Node* node0,
  const GeometryParameterSet::Motion& param,
  int nodeSizeLog2,
  PCCPointSet3* compensatedPointCloud,
  uint32_t depthMax) const
{
  assert(fifo.empty());

  depthMax = std::min(depthMax, depth);
  const int32_t node0SizeMinus1 = (1 << nodeSizeLog2) - 1;
  const auto node0Pos0 = (node0->pos << nodeSizeLog2) + MVd;
  const auto node0Pos1 = node0Pos0 + node0SizeMinus1;
  const auto minNodeSizeMinus1 = (1 << maxDepth - depthMax) - 1;

  std::queue<int> local;
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
  compensatedPointCloud->appendPartition(*pointCloud, indices);
  node0->predEnd = compensatedPointCloud->size();

  for (int i = node0->predStart; i < node0->predEnd; ++i) {
    auto& predPoint = (*compensatedPointCloud)[i];
    predPoint[0] -= MVd[0];
    predPoint[1] -= MVd[1];
    predPoint[2] -= MVd[2];
  }

}

//----------------------------------------------------------------------------
}  // namespace pcc
