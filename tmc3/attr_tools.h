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

#include "quantization.h"
#include "hls.h"
#include <vector>
#include "TMC3.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "PCCMisc.h"
#include "ModeCoder.h"

namespace pcc {

//============================================================================
struct UrahtNodeEncoder {
  int64_t pos;
  int weight;
  std::array<int16_t, 2> qp;
  bool decoded = false;
  attr::Mode mode = attr::Mode::size;
  attr::Mode _mode = attr::Mode::size;
  attr::Mode grand_mode = attr::Mode::size;
  uint8_t occupancy = 0;

  int firstChildIdx;
  uint8_t numChildren;
  int32_t sumAttr[3];
  int32_t sumAttrInter[3];
};

struct UrahtNodeDecoder {
  int64_t pos;
  int weight;
  std::array<int16_t, 2> qp;
  bool decoded = false;
  attr::Mode mode = attr::Mode::size;
  attr::Mode grand_mode = attr::Mode::size;
  uint8_t occupancy = 0;

  int firstChildIdx;
  uint8_t numChildren;
  int32_t sumAttrInter[3] = { 0,0,0 };
};

struct UrahtNodeDecoderHaar {
  int64_t pos;
  bool decoded = false;
  attr::Mode mode = attr::Mode::size;
  attr::Mode grand_mode = attr::Mode::size;
  uint8_t occupancy = 0;

  int firstChildIdx;
  uint8_t numChildren;
  int16_t sumAttrInter[3] = { 0,0,0 };

  // hack for template functions compatibility when non existing members are
  // referred conditionally (if constexpr could be used instead, if C++17
  // was allowed)
  static int weight;
  static std::array<int16_t, 2> qp;
};


//============================================================================

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

  Mode getNeighborsModeEncoder(
    const int parentNeighIdx[19],
    const std::vector<UrahtNodeEncoder>& weightsParent,
    int& voteInterWeight,
    int& voteIntraWeight,
    int& voteInterLayerWeight,
    int& voteIntraLayerWeight,
    std::vector<UrahtNodeEncoder>& weightsLf);

  int getInferredMode(
    bool enableIntraPred,
    bool enableInterPred,
    int childrenCount,
    Mode parent,
    Mode neighbors);

  template<class Kernel, int numAttrs, typename WeightsType>
  Mode choseMode(
    ModeEncoder& rdo,
    const int64_t* transformBuf,
    const std::vector<Mode>& modes,
    const WeightsType weights[],
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp,
    const bool inheritDC);

}  // namespace AttrPrediction

namespace RAHT {

//============================================================================

// Fixed Point precision for RAHT
static constexpr int kFPFracBits = 15;
static constexpr int64_t kFPOneHalf = 1 << (kFPFracBits - 1);
static constexpr int64_t kFPDecMask = ((1 << kFPFracBits) - 1);
static constexpr int64_t kFPIntMask = ~kFPDecMask;

//============================================================================
// Encapsulation of a RAHT transform stage.

class RahtKernel {
public:
  inline
  RahtKernel(int64_t weightLeft, int64_t weightRight, bool scaledWeights = false)
  {
    if (scaledWeights) {
      _a = weightLeft;
      _b = weightRight;
    }
    else {
      int64_t w = weightLeft + weightRight;
      int64_t isqrtW = fastIrsqrt(w);
      _a = fastIsqrt(weightLeft) * isqrtW >> 40;
      _b = fastIsqrt(weightRight) * isqrtW >> 40;
    }
  }

  void fwdTransform(int64_t& lf, int64_t& hf)
  {
    auto tmp = lf * _b;
    lf = fpReduce<kFPFracBits>(lf * _a + hf * _b);
    hf = fpReduce<kFPFracBits>(hf * _a - tmp);
  }

  void invTransform(int64_t& left, int64_t& right)
  {
    int64_t tmp = right * _b;
    right = fpReduce<kFPFracBits>(right * _a + left * _b);
    left = fpReduce<kFPFracBits>(left * _a - tmp);
  }

  int64_t getW0() { return _a; }
  int64_t getW1() { return _b; }

private:
  int64_t _a, _b;
};

//============================================================================
// Encapsulation of an Integer Haar transform stage.

class HaarKernel {
public:
  inline
  HaarKernel(int weightLeft, int weightRight, bool scaledWeights=false) { }

  void fwdTransform(int64_t& lf, int64_t& hf)
  {
    hf -= lf;
    lf += (hf >> 1 + kFPFracBits) << kFPFracBits;
  }

  void invTransform(int64_t& left, int64_t& right)
  {
    left -= (right >> 1 + kFPFracBits) << kFPFracBits;
    right += left;
  }

  int64_t getW0() { return 1; }
  int64_t getW1() { return 1; }

};

//============================================================================
// In-place transform a set of sparse 2x2x2 blocks each using the same weights

template<bool haarFlag>
struct FwdTransformBlock222 {
  template<class Kernel = RahtKernel>
  static inline void
  apply(
    const int numBufs, int64_t* buf, const int64_t weights[])
  {
    static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};
    static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
    const int64_t* w = &weights[32];
    for (int i = 0; i < 12; i++) {
      int64_t w0 = *w++;
      int64_t w1 = *w++;

      if (w0 || w1) {
        int i0 = a[i];
        int i1 = b[i];

        if (!w0) {
          for (int k = 0; k < numBufs; k++)
            std::swap(buf[8 * k + i0], buf[8 * k + i1]);
        }
        else if (w1) { // only one occupied, propagate to next level
          // actual transform
          Kernel kernel(w0, w1, true);
          for (int k = 0; k < numBufs; k++) {
            kernel.fwdTransform(buf[8 * k + i0], buf[8 * k + i1]);
          }
        }
      }
    }
  }

  template<int numBufs, class Kernel = RahtKernel>
  static void
  apply(int64_t* buf, const int64_t weights[])
  {
    apply<Kernel>(numBufs, buf, weights);
  }
};

template<>
struct FwdTransformBlock222<true> {
  template<class Kernel = HaarKernel>
  static inline void
  apply(const int numBufs, int64_t* buf, const bool weights[])
  {
    static const int a[4 + 4 + 4] = { 0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3 };
    static const int b[4 + 4 + 4] = { 1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7 };
    const bool* w = &weights[32];
    for (int i = 0; i < 12; i++) {
      bool w0 = *w++;
      bool w1 = *w++;

      if (w0 || w1) {
        int i0 = a[i];
        int i1 = b[i];

        if (!w0) {
          for (int k = 0; k < numBufs; k++)
            std::swap(buf[8 * k + i0], buf[8 * k + i1]);
        }
        else if (w1) { // only one occupied, propagate to next level
          for (int k = 0; k < numBufs; k++) {
            buf[8 * k + i1] -= buf[8 * k + i0];
            buf[8 * k + i0] += (buf[8 * k + i1] >> 1 + kFPFracBits) << kFPFracBits;
          }
        }
      }
    }
  }

  template<int numBufs, class Kernel = HaarKernel>
  static void
  apply(int64_t* buf, const bool weights[])
  {
    apply<Kernel>(numBufs, buf, weights);
  }
};



template<bool haarFlag>
struct InvTransformBlock222 {
  template<int numBufs, bool computekernel = false, class Kernel = RahtKernel>
  static void
  apply(int64_t* buf, const int64_t weights[])
  {
    static const int a[4 + 4 + 4] = { 0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3 };
    static const int b[4 + 4 + 4] = { 1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7 };
    const int64_t* w = &weights[33 + 22];
    for (int i = 11; i >= 0; i--) {

      int64_t w1 = *(w--);
      int64_t w0 = *(w--);

      if (w0 || w1) {
        int i0 = a[i];
        int i1 = b[i];

        if (!w0) {
          for (int k = 0; k < numBufs; k++)
            buf[8 * k + i1] = buf[8 * k + i0];
        }
        else if (w1) { // only one occupied, propagate to next level
          // actual transform
          Kernel kernel(w0, w1, !computekernel);
          for (int k = 0; k < numBufs; k++) {
            kernel.invTransform(buf[8 * k + i0], buf[8 * k + i1]);
          }
        }
      }
    }
  }
};

template<>
struct InvTransformBlock222<true> {
  template<int numBufs, bool computekernel = false, class Kernel = HaarKernel>
  static void
  apply(int64_t* buf, const bool weights[])
  {
    static const int a[4 + 4 + 4] = { 0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3 };
    static const int b[4 + 4 + 4] = { 1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7 };
    const bool* w = &weights[33 + 22];
    for (int i = 11; i >= 0; i--) {

      bool w1 = *(w--);
      bool w0 = *(w--);

      if (w0 || w1) {
        int i0 = a[i];
        int i1 = b[i];

        if (!w0) {
          for (int k = 0; k < numBufs; k++)
            buf[8 * k + i1] = buf[8 * k + i0];
        }
        else if (w1) { // only one occupied, propagate to next level
          for (int k = 0; k < numBufs; k++) {
            buf[8 * k + i0] -= (buf[8 * k + i1] >> 1 + kFPFracBits) << kFPFracBits;
            buf[8 * k + i1] += buf[8 * k + i0];
          }
        }
      }
    }
  }
};

} /* namespace RAHT */

} /* namespace pcc */
