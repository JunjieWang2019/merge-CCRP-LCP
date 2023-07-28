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

#ifndef PCCTMC3Common_h
#define PCCTMC3Common_h

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "constants.h"
#include "hls.h"

#include "nanoflann.hpp"

#include <cstdint>
#include <cstddef>
#include <memory>
#include <vector>

#include "attr_tools.h"

namespace pcc {

//============================================================================
// Hierachichal bounding boxes.
// Insert points (into the base layer), then generate the hierarchy via update.

template<int32_t BucketSizeLog2, int32_t LevelCount>
class BoxHierarchy {
public:
  void resize(const int32_t pointCount)
  {
    constexpr auto BucketSize = 1 << BucketSizeLog2;
    constexpr auto BucketSizeMinus1 = BucketSize - 1;
    int32_t count = pointCount;
    for (int i = 0; i < LevelCount; ++i) {
      count = (count + BucketSizeMinus1) >> BucketSizeLog2;
      _bBoxes[i].clear();
      _bBoxes[i].resize(count, Box3<int32_t>(INT32_MAX, INT32_MIN));
    }
  }

  void insert(const Vec3<int32_t>& point, const int32_t index)
  {
    const auto bindex = (index >> BucketSizeLog2);
    assert(bindex >= 0 && bindex < _bBoxes[0].size());
    _bBoxes[0][bindex].insert(point);
  }

  void update()
  {
    constexpr auto LevelCountMinus1 = LevelCount - 1;
    for (int i = 0; i < LevelCountMinus1; ++i) {
      for (int32_t j = 0, count = int32_t(_bBoxes[i].size()); j < count; ++j) {
        _bBoxes[i + 1][j >> BucketSizeLog2].merge(_bBoxes[i][j]);
      }
    }
  }

  const Box3<int32_t>& bBox(int32_t bindex, int32_t level) const
  {
    return _bBoxes[level][bindex];
  }

  int32_t bucketSizeLog2(int32_t level = 0) const
  {
    return BucketSizeLog2 * (1 + level);
  }

  int32_t bucketSize(int32_t level = 0) const
  {
    return 1 << bucketSizeLog2(level);
  }

private:
  std::vector<Box3<int32_t>> _bBoxes[LevelCount];
};

//============================================================================

class MortonIndexMap3d {
public:
  struct Range {
    int32_t start;
    int32_t end;
  };

  void resize(const int32_t cubeSizeLog2)
  {
    _cubeSizeLog2 = cubeSizeLog2;
    _cubeSize = 1 << cubeSizeLog2;
    _bufferSize = 1 << (3 * cubeSizeLog2);
    _mask = _bufferSize - 1;
    _buffer.reset(new Range[_bufferSize]);
  }

  void reserve(const uint32_t sz) { _updates.reserve(sz); }
  int cubeSize() const { return _cubeSize; }
  int cubeSizeLog2() const { return _cubeSizeLog2; }

  void init()
  {
    for (int32_t i = 0; i < _bufferSize; ++i) {
      _buffer[i] = {-1, -1};
    }
    _updates.resize(0);
  }

  void clearUpdates()
  {
    for (const auto index : _updates) {
      _buffer[index] = {-1, -1};
    }
    _updates.resize(0);
  }

  void set(const int64_t mortonCode, const int32_t index)
  {
    const int64_t mortonAddr = mortonCode & _mask;
    auto& unit = _buffer[mortonAddr];
    if (unit.start == -1) {
      unit.start = index;
    }
    unit.end = index + 1;
    _updates.push_back(mortonAddr);
  }

  Range get(const int64_t mortonCode) const
  {
    return _buffer[mortonCode & _mask];
  }

private:
  int32_t _cubeSize = 0;
  int32_t _cubeSizeLog2 = 0;
  int32_t _bufferSize = 0;
  int64_t _mask = 0;
  std::unique_ptr<Range[]> _buffer;

  // A list of indexes in _buffer that are dirty
  std::vector<int32_t> _updates;
};

//============================================================================

struct MortonCodeWithIndex {
  int64_t mortonCode;

  // The position used to generate the mortonCode
  Vec3<int32_t> position;

  int32_t index;

  bool operator<(const MortonCodeWithIndex& rhs) const
  {
    // NB: index used to maintain stable sort
    if (mortonCode == rhs.mortonCode)
      return index < rhs.index;
    return mortonCode < rhs.mortonCode;
  }
};

//---------------------------------------------------------------------------

struct PCCNeighborInfo {
  uint64_t weight;
  uint32_t predictorIndex;

  uint32_t pointIndex;
  bool interFrameRef;

  bool operator<(const PCCNeighborInfo& rhs) const
  {
    return weight < rhs.weight;
  }
};

//---------------------------------------------------------------------------

struct AttributeInterPredParams {
  PCCPointSet3 referencePointCloud;
  PCCPointSet3 compensatedPointCloud;
  std::vector<MotionVector> motionVectors;
  int frameDistance;
  bool enableAttrInterPred;
  int getPointCount() const { return referencePointCloud.getPointCount(); }
  void clear() { referencePointCloud.clear(); }
  bool hasLocalMotion() const { return compensatedPointCloud.getPointCount() > 0; }
};

//---------------------------------------------------------------------------

}  // namespace pcc

#endif /* PCCTMC3Common_h */
