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
#include "attr_tools.h"
#include "ply.h"
#include <vector>
#include "pointset_processing.h"
#include "ringbuf.h"

namespace pcc {

//============================================================================

void regionAdaptiveHierarchicalTransform(
  const RahtPredictionParams& rahtPredParams,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  const int attribCount,
  const int voxelCount,
  int64_t* mortonCode,
  attr_t* attributes,
  const attr_t* attributes_mc,
  int* coefficients,
  attr::ModeEncoder& encoder);

void regionAdaptiveHierarchicalInverseTransform(
  const RahtPredictionParams& rahtPredParams,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  const int attribCount,
  const int voxelCount,
  int64_t* mortonCode,
  attr_t* attributes,
  const attr_t* attributes_mc,
  int* coefficients,
  attr::ModeDecoder& decoder);

//============================================================================

namespace RAHT {

//============================================================================
// Find the neighbours of the node indicated by @t between @first and @last.
// The position weight of each found neighbour is stored in two arrays.

template<typename It>
int
findNeighbours(
  It first,
  It last,
  It it,
  int level,
  uint8_t occupancy,
  int parentNeighIdx[19])
{
  static const uint8_t neighMasks[19] = {255, 240, 204, 170, 192, 160, 136,
                                         3,   5,   15,  17,  51,  85,  10,
                                         34,  12,  68,  48,  80};

  // current position (discard extra precision)
  int64_t cur_pos = it->pos >> level;

  // the position of the parent, offset by (-1,-1,-1)
  int64_t base_pos = morton3dAdd(cur_pos, -1ll);

  // these neighbour offsets are relative to base_pos
  static const uint8_t neighOffset[19] = {0, 35, 21, 14, 49, 42, 28, 1,  2, 3,
                                          4, 5,  6,  10, 12, 17, 20, 33, 34};

  // special case for the direct parent (no need to search);
  parentNeighIdx[0] = std::distance(first, it);

  int parentNeighCount = parentNeighIdx[0] != -1;
  for (int i = 1; i < 19; i++) {
    // Only look for neighbours that have an effect
    if (!(occupancy & neighMasks[i])) {
      parentNeighIdx[i] = -1;
      continue;
    }

    // compute neighbour address to look for
    // the delta between it and the current position is
    int64_t neigh_pos = morton3dAdd(base_pos, neighOffset[i]);
    int64_t delta = neigh_pos - cur_pos;

    // find neighbour
    It start = first;
    It end = last;

    if (delta >= 0) {
      start = it;
      if ((delta + 1) < std::distance(it, last))
        end = std::next(it, delta + 1);
    }
    else {
      end = it;
      if ((-delta) < std::distance(first, it))
        start = std::prev(it, -delta);
    }
    It found = std::lower_bound(
      start, end, neigh_pos,
      [=](decltype(*it)& candidate, int64_t neigh_pos) {
        return (candidate.pos >> level) < neigh_pos; });
    parentNeighIdx[i] =
      found == end || (found->pos >> level) != neigh_pos
      ? -1 : std::distance(first, found);
    parentNeighCount += parentNeighIdx[i] != -1;
  }

  return parentNeighCount;
}

//============================================================================
// computeParentDc

template<bool haarFlag, int numAttrs>
void
computeParentDc(
  int64_t w,
  std::vector<int64_t>::const_iterator dc,
  int64_t parentDc[3]
)
{
  if (haarFlag) {
    for (int k = 0; k < numAttrs; k++) {
      parentDc[k] = dc[k];
    }
    return;
  }

  int shift = 5 * ((w > 1024) + (w > 1048576));
  int64_t rsqrtWeight = fastIrsqrt(w) >> (40 - shift - kFPFracBits);
  for (int k = 0; k < numAttrs; k++) {
    parentDc[k] = fpReduce<kFPFracBits>(
      (dc[k] >> shift) * rsqrtWeight);
  }
}

//============================================================================
// expand a set of eight weights into three levels

template<bool skipkernel, class Kernel>
void
mkWeightTree(int64_t weights[8 + 8 + 8 + 8 + 24])
{
  auto in = &weights[0];
  auto out = &weights[8];

  for (int i = 0; i < 4; i++) {
    out[0] = in[0] + in[1];
    out[4] = (in[0] && in[1]) * out[0];  // single node, no high frequencies
    in += 2;
    out++;
  }
  out += 4;
  for (int i = 0; i < 4; i++) {
    out[0] = in[0] + in[1];
    out[4] = (in[0] && in[1]) * out[0];  // single node, no high frequencies
    in += 2;
    out++;
  }
  out += 4;
  for (int i = 0; i < 4; i++) {
    out[0] = in[0] + in[1];
    out[4] = (in[0] && in[1]) * out[0];  // single node, no high frequencies
    in += 2;
    out++;
  }

  for (int i = 0; i < 24; i += 2) {
    weights[i + 32] = (!!weights[i]) << kFPFracBits;
    weights[i + 33] = (!!weights[i + 1]) << kFPFracBits;

    if (weights[i] && weights[i + 1]) {
      if (skipkernel) {
        weights[i + 32] = weights[i];
        weights[i + 33] = weights[i + 1];
      }
      else {
        Kernel w(weights[i], weights[i + 1]);
        weights[i + 32] = w.getW0();
        weights[i + 33] = w.getW1();
      }
    }
  }
}

//============================================================================
// Invoke mapFn(coefIdx) for each present coefficient in the transform

template<class T>
void
scanBlock(int64_t weights[], T mapFn)
{
  static const int8_t kRahtScanOrder[] = {0, 4, 2, 1, 6, 5, 3, 7};

  // there is always the DC coefficient (empty blocks are not transformed)
  mapFn(0);

  for (int i = 1; i < 8; i++) {
    if (weights[24 + kRahtScanOrder[i]])
      mapFn(kRahtScanOrder[i]);
  }
}

//============================================================================
// translateLayer

template<bool haarFlag, int numAttrs, typename UrahtNode>
void
translateLayer(
  std::vector<int64_t>& layerAttr,
  const attr_t* attr_mc,
  std::vector<UrahtNode>& weightsLf)
{
  // associate mean attribute of MC PC to each unique node
  layerAttr.resize(weightsLf.size() * numAttrs);
  auto layer = layerAttr.begin();
  for (int i = 0, j = 0;
      i < weightsLf.size();
      i++, layer += numAttrs) {

    for (int k = 0; k < numAttrs; k++)
      layer[k] = -1;

    int weight = weightsLf[i].weight;
    int jEnd = j + weight;

    auto attr = &attr_mc[numAttrs * j];

    std::array<int, numAttrs> sumAtt;
    std::fill_n(sumAtt.begin(), numAttrs, 0);

    for (; j < jEnd; ++j) {
      for (int k = 0; k < numAttrs; k++)
        sumAtt[k] += *attr++;
    }

    if (weight) {
      if (haarFlag)
        for (int k = 0; k < numAttrs; k++)
          layer[k] = int64_t(sumAtt[k]);
      else
        for (int k = 0; k < numAttrs; k++)
          layer[k] = int64_t(sumAtt[k]) << kFPFracBits;

      if (weight != 1) {
        for (int k = 0; k < numAttrs; k++) {
          layer[k] /= weight;
        }
      }

      if (haarFlag)
        for (int k = 0; k < numAttrs; k++)
          layer[k] <<= kFPFracBits;
    }
  }
}

//============================================================================
// Cross Chroma Component Prediction

struct PCCRAHTComputeCCCP {
  int8_t computeCrossChromaComponentPredictionCoeff(
    int m, int64_t coeffs[][3]);

  void reset()
  {
    sum = {0, 0};
    window.clear();
  }

private:
  struct Elt {
    int64_t k1k2;
    int64_t k1k1;
  };

  Elt sum {0, 0};
  ringbuf<Elt> window = ringbuf<Elt>(128);
};

//============================================================================
// Compute CCRP filter by floor division

static constexpr int kCCRPFiltPrecisionbits = 3;
static constexpr int kCCRPtemplatesize = 16;

template<int precbits=kCCRPFiltPrecisionbits, int buflen=kCCRPtemplatesize>
struct _CCRPFilter {

  struct Corr {
    int64_t yy;
    int64_t ycb;
    int64_t ycr;

    void operator -=(const Corr& other) {
      yy -= other.yy;
      ycb -= other.ycb;
      ycr -= other.ycr;
    }

    void operator +=(const Corr& other) {
      yy += other.yy;
      ycb += other.ycb;
      ycr += other.ycr;
    }
  };

  void update(const Corr& curCorr) {
    if (window.size() == buflen) {
      sum -= window.front();
      window.pop();
    }
    window.push(curCorr);
    sum += curCorr;

    _YCbFilt = getFilterTap(sum.yy+1, sum.ycb);
    _YCrFilt = getFilterTap(sum.yy+1, sum.ycr);
  }

  int getYCbFilt() const { return _YCbFilt; }
  int getYCrFilt() const { return _YCrFilt; }

private:
  Corr sum {0, 0, 0};

  int _YCbFilt = 0;
  int _YCrFilt = 0;

  ringbuf<Corr> window = ringbuf<Corr>(buflen);
public: // TMP
  static int getFilterTap (int64_t autocorr, int64_t crosscorr)
  {
    constexpr int unity = 1 << precbits;
    //binary search to replace divison. ( returns 128*crosscorr/autocorr)
    if (crosscorr == 0)
      return 0;

    const bool isneg = crosscorr < 0;
    crosscorr = std::abs(crosscorr);

    if (crosscorr == autocorr)
      return isneg ? -unity : unity;

    int tapint = 0, tapfrac = -1;

    //determine integer part by repeated subtraction
    while (crosscorr >= autocorr) {
      crosscorr -= autocorr;
      tapint += unity;
    }

    if (crosscorr == 0) {
      return isneg ? -tapint : tapint;
    }

    const int64_t numerator = crosscorr << precbits;
    const int64_t denominator = autocorr;

    for(int64_t cumsum = 0; cumsum <= numerator; cumsum += denominator) {
      tapfrac++;
    }

    return isneg ? -tapint - tapfrac : tapint + tapfrac;
  }
};

using CCRPFilter = _CCRPFilter<>;

//============================================================================

} /* namespace RAHT */

//============================================================================

} /* namespace pcc */
