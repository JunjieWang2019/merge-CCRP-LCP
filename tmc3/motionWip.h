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

#pragma once
#include "TMC3.h"
#include <vector>

#include "PCCPointSet.h"
#include "entropy.h"
#include "hls.h"
#include "attr_tools.h"

#include <queue>
#include <tuple>

#include <queue>
#include <tuple>

#include <queue>
#include <tuple>

#include "ringbuf.h"

namespace pcc {

//============================================================================
static const unsigned int motionParamPrec = 16;
static const unsigned int motionParamScale = 1 << motionParamPrec;
static const unsigned int motionParamOffset = 1 << (motionParamPrec - 1);
int plus1log2shifted4(int x);
struct PCCOctree3Node;
struct MSOctree;
struct ParameterSetMotion;
struct EncodeMotionSearchParams;
//============================================================================

int roundIntegerHalfInf(const double x);

//----------------------------------------- LOCAL MOTION -------------------
struct MVField {
  static constexpr uint64_t kMaskX = 0x0FFFFF0000000000ULL;
  static constexpr uint64_t kMaskY = 0x000000FFFFF00000ULL;
  static constexpr uint64_t kMaskZ = 0x00000000000FFFFFULL;
  static constexpr int kOffsetX = 40;
  static constexpr int kOffsetY = 20;
  static constexpr int kOffsetZ = 0;

  static constexpr uint32_t kNotSetMVIdx = 0xFFFFFF;

  struct PUNode {
    PUNode()
      : _mvIdx(kNotSetMVIdx)
      , _childsMask(0)
    {}
    // should be 16 bytes aligned
    // TODO: use masks getter/setters to ensure bit order
    struct {
      union {
        uint32_t _firstChildIdx:24; // if _childsMask != 0
        uint32_t _mvIdx:24; // if _childsMask == 0, kNotSetMVIdx means not set
      };
      uint32_t _childsMask:8; // 0 means leaf node or not set
    };
    struct {
      uint32_t _reserved:27;
      uint32_t _puSizeLog2:5;
    };
    uint64_t _packedPos0; // z | y << 20 | x << 40

    point_t pos0() const
    { return point_t {
        int32_t((_packedPos0 & kMaskX) >> kOffsetX),
        int32_t((_packedPos0 & kMaskY) >> kOffsetY),
        int32_t((_packedPos0 & kMaskZ) >> kOffsetZ)};
    }

    void set_pos0(point_t _pos)
    {
      _packedPos0 =
        (int64_t(_pos[0]) << kOffsetX)
        + (int64_t(_pos[1]) << kOffsetY)
        + (int64_t(_pos[2]) << kOffsetZ);
    }
  };
  //uint32_t puSizeLog2;
  uint32_t numRoots;
  // (numRoots PUs in first layer)
  std::vector<PUNode> puNodes;
  // mvPool stores the motion vectors for child Nodes
  std::vector<point_t> mvPool;

  MVField() = default;
  MVField(const MVField&) = default;
  MVField(MVField&&) = default;
  MVField(const MVField& from, point_t begin, point_t end)
  {
    puNodes.reserve(from.puNodes.size());
    mvPool.reserve(from.mvPool.size());

    auto pos1 = end - 1;
    auto pos0 = begin;
    // For now taking entire node
    // TODO: generate new by only taking the intersection ? => not so easy
    for (int i = 0; i < from.numRoots; ++i) {
      const auto& nodeFrom = from.puNodes[i];
      auto nodePos0 = nodeFrom.pos0();
      auto nodePos1 = nodePos0 + ((1 << nodeFrom._puSizeLog2) - 1);
      int intersect = 0;
      for (int k=0; k < 3; ++k) {
        intersect |= std::min(nodePos1[k], pos1[k]) - std::max(nodePos0[k], pos0[k]);
      }
      if (intersect >= 0)
        puNodes.push_back(nodeFrom);
    }
    numRoots = puNodes.size();

    int N = puNodes.size();
    // TODO: depth first would be better
    for (int i = 0; i < N; ++i) {
      auto& node = puNodes[i];
      if (node._childsMask) {
        uint32_t childIdxFrom = node._firstChildIdx;
        node._firstChildIdx = N;
        for (int c = 0; c < 8; ++c) {
          if (node._childsMask & (1 << c)) {
            const auto& childNodeFrom =  from.puNodes[childIdxFrom++];
            puNodes.push_back(childNodeFrom);
            ++N;
          }
        }
      } else {
        uint32_t mvIdxFrom = node._mvIdx;
        node._mvIdx = mvPool.size();
        mvPool.push_back(from.mvPool[mvIdxFrom]);
      }
    }
  }
  MVField & operator =(const MVField&) = default;
  MVField & operator =(MVField&&) = default;
};

int deriveMotionMaxPrefixBits(int window_size);
int deriveMotionMaxSuffixBits(int window_size);

struct LPUwindow {
  Vec3<int> pos;
  Vec3<attr_t> color;
};


void bounded_splitPU_MC(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  point_t nodeSizeLog2,
  point_t boundPos0,
  point_t boundPos1,
  PCCPointSet3* compensatedPointCloud);

template <bool mcap>
void encode_splitPU_MV_MC(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2 = false,
  int S = -1,
  int S2 = -1);

// motion decoder

template <bool mcap>
void decode_splitPU_MV_MC(
  const MSOctree& mSOctree,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2 = false,
  int S = -1,
  int S2 = -1);

  //============================================================================

class MotionEntropyEstimate;
struct MSOctree {
  MSOctree& operator=(const MSOctree&) = default;
  MSOctree& operator=(MSOctree&&) = default;
  MSOctree() = default;
  MSOctree(const MSOctree&) = default;
  MSOctree(MSOctree&&) = default;
  MSOctree(
    PCCPointSet3* predPointCloud,
    point_t offsetOrigin,
    uint32_t leafSizeLog2 = 0
    );

  struct MSONode {
    uint32_t start;
    uint32_t end;
    std::array<uint32_t, 8> child = {}; // 0 means none
    int32_t sizeMinus1;
    point_t pos0;
    uint32_t parent;
    uint32_t reserved; // to align to 64 bytes (otherwise could be avoided)

    uint32_t numPoints() const { return end - start; }
  };

  point_t offsetOrigin; // offset applied to origin while construction the octree
  uint32_t maxDepth; // depth of full octree to get unitary sized nodes
  uint32_t depth; // depth of the motion search octree
  PCCPointSet3* pointCloud;
  std::vector<MSONode> nodes;

  void allocRingBuffers() {
    a = ringbuf<int>(nodes.size());
    b = ringbuf<int>(nodes.size());
  }

  int
  nearestNeighbour_updateDMax(point_t pos, int32_t& d_max, bool approximate = false) const;

  inline
  int
  iNearestNeighbour_updateDMax(const point_t& pos, int32_t& d_max) const;

  inline
  int
  iApproximateNearestNeighbour_updateDMax(const point_t& pos, int32_t& d_max) const;

  inline
  int
  iApproxNearestNeighbourAttr(const point_t& pos) const;

  double
  find_motion(
    const EncodeMotionSearchParams& param,
    const ParameterSetMotion& mvPS,
    const MotionEntropyEstimate& motionEntropy,
    const PCCPointSet3& Block0,
    const point_t& xyz0,
    int local_size,
    MVField& mvField,
    uint32_t puNodeIdx // node Idx in mvField
  ) const;

  void
  apply_motion(
    const point_t currNodePos0,
    const point_t currNodePos1,
    point_t Mvd,
    PCCOctree3Node* node0,
    PCCPointSet3* compensatedPointCloud,
    uint32_t depthMax = UINT32_MAX,
    bool flagNonPow2 = false,
    int S = -1,
    int S2 = -1
  ) const;

  int
  nodeIdx(point_t nodePos0, uint32_t nodeSizeLog2) const;

  void
  apply_recolor_motion(
    point_t Mvd,
    PCCOctree3Node* node0,
    PCCPointSet3& pointCloud
  ) const;

  mutable ringbuf<int> a; // for search
  mutable ringbuf<int> b; // for search
};

//----------------------------------------------------------------------------

bool
motionSearchForNode(
  const MSOctree& mSOctreeOrig,
  const MSOctree& mSOctree,
  const PCCOctree3Node* node0,
  const EncodeMotionSearchParams& msParams,
  const ParameterSetMotion& mvPS,
  int nodeSize,
  EntropyEncoder* arithmeticEncoder,
  MVField& mvField,
  uint32_t puNodeIdx, // node Idx in mvField
  bool flagNonPow2 = false,
  int S = -1,
  int S2 = -1
);

//============================================================================

struct InterPredParams {
  PCCPointSet3 referencePointCloud;
  MSOctree mSOctreeRef;
  PCCPointSet3 compensatedPointCloud;
  // Motion
  MVField mvField;
  // TMP hack
  mutable std::vector<int> attributes_mc;
};

//============================================================================

}  // namespace pcc
