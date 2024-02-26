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
#include <cstdint>

#include "FixedPoint.h"
#include "quantization.h"
#include "hls.h"
#include <vector>
#include "TMC3.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "PCCMisc.h"
#include "ModeCoder.h"

using VecAttr = std::vector<std::array<pcc::FixedPoint, 8>>;

namespace pcc {

void computeParentDc(
  const int64_t weights[],
  std::vector<int64_t>::const_iterator dc,
  std::vector<FixedPoint>& parentDc,
  bool integer_haar_enable_flag);

void translateLayer(
  std::vector<int64_t>& layerAttr,
  size_t layerDepth,
  size_t attrCount,
  size_t count_rf,
  size_t count_mc,
  int64_t* morton_rf,
  int64_t* morton_rf_transformed,
  int64_t* morton_mc,
  int* attr_mc,
  bool integer_haar_enable_flag,
  size_t layerSize = 0);

//============================================================================

struct UrahtNode {
  int64_t pos;
  int weight;
  Qps qp;
  bool decoded;
  int64_t offset;
  attr::Mode mode;
  attr::Mode _mode;
  attr::Mode mc;

  uint8_t occupancy;
  std::vector<UrahtNode>::iterator firstChild;
  std::vector<UrahtNode>::iterator lastChild;
};

enum NeighborType
{
  Self,
  Face,
  Edge,
  Other
};

struct ParentIndex {
  int index;
  NeighborType type;
  std::array<int, 8> children;
  std::array<NeighborType, 8> childrenType;

  ParentIndex() { reset(); }

  bool isvoid() const { return type == NeighborType::Other; }

  bool isvoid(int i) const { return childrenType[i] == NeighborType::Other; }

  int operator[](int i) const { return children[i]; }

  void set(int i, NeighborType t) { index = i; type = t;}
  void set(int j, int i, NeighborType t)
  {
    children[j] = i;
    childrenType[j] = t;
  }

  void operator=(NeighborType t) { type = t; }

  void setType(int i, NeighborType t) { childrenType[i] = t; }

  operator int() const { return index; }

  void reset()
  {
    type = NeighborType::Other;
    std::fill(childrenType.begin(), childrenType.end(), NeighborType::Other);
  }

  int getweight(const std::vector<int>& weights) const
  {
    switch (type) {
    case NeighborType::Self: return weights[0];
    case NeighborType::Face: return weights[1];
    case NeighborType::Edge: return weights[2];
    default: return 0;
    }
  }

  int getweight(int i, const std::vector<int>& weights) const
  {
    switch (childrenType[i]) {
    case NeighborType::Self: return 0;
    case NeighborType::Face: return weights[3];
    case NeighborType::Edge: return weights[4];
    default: return 0;
    }
  }
};

template<class T>
inline T
operator*(const T& a, const ParentIndex& b)
{
  return a * b.index;
}

template<class T>
inline T
operator*(const ParentIndex& a, const T& b)
{
  return a.index * b;
}

template<class T>
inline T
operator+(const T& a, const ParentIndex& b)
{
  return a + b.index;
}

template<class T>
inline T
operator+(const ParentIndex& a, const T& b)
{
  return a.index + b;
}

namespace attr {
  Mode getNeighborsMode(
    const bool& isEncoder,
    const int parentNeighIdx[19],
    const std::vector<UrahtNode>& weightsParent,
    int16_t& voteInterWeight,
    int16_t& voteIntraWeight,
    int16_t& voteInterLayerWeight,
    int16_t& voteIntraLayerWeight);

  Mode getInferredMode(
    int& ctxMode,
    bool enableIntraPrediction,
    bool enableInterPrediction,
    int childrenCount,
    Mode parent,
    Mode neighbors,
    const int numAttrs,
    const int64_t weights[8],
    std::vector<int64_t>::const_iterator dc);

  template<class Kernel>
  Mode choseMode(
    ModeEncoder& rdo,
    const VecAttr& transformBuf,
    const std::vector<Mode>& modes,
    const int64_t weights[],
    const int numAttrs,
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp);

}  // namespace AttrPrediction

//============================================================================
// Encapsulation of a RAHT transform stage.

class RahtKernel {
public:
  RahtKernel(int weightLeft, int weightRight, bool scaledWeights=false)
  {
    if (scaledWeights) {
      _a.val = weightLeft;
      _b.val = weightRight;
    } else {
      uint64_t w = weightLeft + weightRight;
      uint64_t isqrtW = irsqrt(w);
      _a.val =
        (isqrt(uint64_t(weightLeft) << (2 * _a.kFracBits)) * isqrtW) >> 40;
      _b.val =
        (isqrt(uint64_t(weightRight) << (2 * _b.kFracBits)) * isqrtW) >> 40;
    }
  }

  void fwdTransform(FixedPoint& lf, FixedPoint& hf)
  {
    // lf = right * b + left * a
    // hf = right * a - left * b
    auto tmp = lf.val * _b.val;
    lf.val = lf.val * _a.val + hf.val * _b.val;
    hf.val = hf.val * _a.val - tmp;
    //auto tmp = lf.val * _b.val;
    //lf.val *= _a.val;
    //lf.val += hf.val * _b.val;
    //hf.val *= _a.val;
    //hf.val -= tmp;

    lf.fixAfterMultiplication();
    hf.fixAfterMultiplication();
  }

  void invTransform(FixedPoint& left, FixedPoint& right)
  {
    // left = lf * a - hf * b
    // right = lf * b + hf * a
    auto tmp = right.val * _b.val;
    right.val = right.val * _a.val + left.val * _b.val;
    left.val = left.val *_a.val - tmp;
    //auto tmp = right.val * _b.val;
    //right.val *= _a.val;
    //right.val += left.val * _b.val;
    //left.val *= _a.val;
    //left.val -= tmp;

    right.fixAfterMultiplication();
    left.fixAfterMultiplication();
  }

  int64_t getW0() { return _a.val; }
  int64_t getW1() { return _b.val; }

private:
  FixedPoint _a, _b;
};

//============================================================================
// Encapsulation of an Integer Haar transform stage.

class HaarKernel {
public:
  HaarKernel(int weightLeft, int weightRight, bool scaledWeights=false) { }

  void fwdTransform(FixedPoint& lf, FixedPoint& hf)
  {
    hf.val -= lf.val;
    lf.val += ((hf.val >> (1 + hf.kFracBits)) << hf.kFracBits);
  }

  void invTransform(FixedPoint& left, FixedPoint& right)
  {
    left.val -= (((right.val >> (1 + right.kFracBits)) << right.kFracBits));
    right.val += left.val;
  }

  int64_t getW0() { return 1; }
  int64_t getW1() { return 1; }

};

//============================================================================
// In-place transform a set of sparse 2x2x2 blocks each using the same weights

template<class Kernel>
void
fwdTransformBlock222(
  const int numBufs, VecAttr::iterator buf, const int64_t weights[])
{
  static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};
  static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
  for (int i = 0, iw = 0; i < 12; i++, iw += 2) {
    int i0 = a[i];
    int i1 = b[i];
    int64_t w0 = weights[iw + 32];
    int64_t w1 = weights[iw + 33];

    if (!w0 && !w1)
      continue;

    // only one occupied, propagate to next level
    if (!w0 || !w1) {
      if (!w0) {
        for (int k = 0; k < numBufs; k++)
          std::swap(buf[k][i0], buf[k][i1]);
      }
      continue;
    }

    // actual transform
    Kernel kernel(w0, w1, true);
    for (int k = 0; k < numBufs; k++) {
      kernel.fwdTransform(buf[k][i0], buf[k][i1]);
    }
  }
}

template<class Kernel>
void
invTransformBlock222(
  const int numBufs, VecAttr::iterator buf, const int64_t weights[])
{
  static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};
  static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
  for (int i = 11, iw = 22; i >= 0; i--, iw -= 2) {
    int i0 = a[i];
    int i1 = b[i];
    int64_t w0 = weights[iw + 32];
    int64_t w1 = weights[iw + 33];

    if (!w0 && !w1)
      continue;

    // only one occupied, propagate to next level
    if (!w0 || !w1) {
      if (!w0) {
        for (int k = 0; k < numBufs; k++)
          std::swap(buf[k][i0], buf[k][i1]);
      }
      continue;
    }

    // actual transform
    Kernel kernel(w0, w1, true);
    for (int k = 0; k < numBufs; k++) {
      kernel.invTransform(buf[k][i0], buf[k][i1]);
    }
  }
}

} /* namespace pcc */
