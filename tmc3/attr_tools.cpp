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

#include "attr_tools.h"
#include "RAHT.h"

#include <cassert>
#include <cinttypes>
#include <climits>
#include <cstddef>
#include <utility>
#include <vector>
#include <stdio.h>

#include "PCCTMC3Common.h"
#include "PCCMisc.h"

namespace pcc {

void
computeParentDc(
  const int64_t weights[],
  std::vector<int64_t>::const_iterator dc,
  std::vector<FixedPoint>& parentDc)
{
  uint64_t w = weights[0];
  for (int i = 1; i < 8; i++)
    w += weights[i];

  FixedPoint rsqrtWeight;
  int shift = w > 1024 ? ilog2(w - 1) >> 1 : 0;
  rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
  for (int k = 0; k < parentDc.size(); k++) {
    parentDc[k].val = dc[k] >> shift;
    parentDc[k] *= rsqrtWeight;
  }
}

void
translateLayer(
  std::vector<int64_t>& layerAttr,
  size_t layerDepth,
  size_t attrCount,
  size_t count_rf,
  size_t count_mc,
  int64_t* morton_rf,
  int64_t* morton_rf_transformed,
  int64_t* morton_mc,
  int* attr_mc,
  size_t layerSize)
{
  size_t shift = layerDepth * 3;
  std::vector<MortonCodeWithIndex> packed;
  if (layerSize)
    packed.reserve(layerSize / attrCount);
  else
    packed.reserve(count_rf);

  int64_t prev = -1;
  size_t i = 0;
  MortonCodeWithIndex val;
  for (size_t n = 0; n < count_rf; n++) {
    int64_t curr = morton_rf[n] >> shift;
    if (curr != prev) {
      prev = curr;
      val.mortonCode = morton_rf_transformed[n] >> shift;
      val.index = i++;
      packed.push_back(val);
    }
  }
  std::sort(packed.begin(), packed.end());
  count_rf = packed.size();
  layerAttr.resize(count_rf * attrCount);

  i = 0;
  size_t j = 0;
  while (i < count_rf && j < count_mc) {
    prev = packed[i].mortonCode;

    while (j < count_mc && prev > (morton_mc[j] >> shift))
      j++;

    int64_t weight = 0;
    auto layer = std::next(layerAttr.begin(), attrCount * packed[i].index);
    for (size_t k = 0; k < attrCount; k++)
      layer[k] = 0;

    while (j < count_mc && prev == (morton_mc[j] >> shift)) {
      weight++;
      auto attr = &attr_mc[attrCount * j];
      for (size_t k = 0; k < attrCount; k++)
        layer[k] += static_cast<int64_t>(attr[k]) << FixedPoint::kFracBits;
      j++;
    }

    if (weight > 1) {
      for (size_t k = 0; k < attrCount; k++)
        layer[k] /= weight;
    } else if (!weight) {
      for (size_t k = 0; k < attrCount; k++)
        layer[k] = -1;
    }

    i++;
    while (i < count_rf && prev == packed[i].mortonCode) {
      std::copy(
        layer, layer + attrCount,
        std::next(layerAttr.begin(), attrCount * packed[i].index));
      i++;
    }
  }
}


namespace attr {
  Mode getNeighborsMode(
    const int parentNeighIdx[19],
    const std::vector<UrahtNode>& weightsParent)
  {
    int voteNull = 0;
    int voteIntra = 0;
    int voteInter = 0;

    for (int i = 1; i < 19; i++) {
      if (parentNeighIdx[i] == -1)
        continue;

      auto uncle = std::next(weightsParent.begin(), parentNeighIdx[i]);

      if (isNull(uncle->mode))
        voteNull += 1;
      else if (isIntra(uncle->mode))
        voteIntra += 1;
      else if (isInter(uncle->mode))
        voteInter += 1;

      if (uncle->decoded) {
        for (auto cousin = uncle->firstChild; cousin < uncle->lastChild;
             cousin++) {
          if (isNull(cousin->mode))
            voteNull += 3;
          else if (isIntra(cousin->mode))
            voteIntra += 3;
          else if (isInter(cousin->mode))
            voteInter += 3;
        }
      }
    }

    if (voteNull > voteIntra && voteNull > voteInter)
      return Mode::Null;
    if (voteIntra > voteInter)
      return Mode::Intra;
    return Mode::Inter;
  }

  Mode getInferredMode(
    int& ctxMode,
    bool enableIntraPrediction,
    bool enableInterPrediction,
    int childrenCount,
    Mode parent,
    Mode neighbors,
    const int numAttrs,
    const int64_t weights[8],
    std::vector<int64_t>::const_iterator dc)
  {
    // [3] Context: enabled flag
    if (enableIntraPrediction && enableInterPrediction)
      ctxMode = 2;
    else if (enableIntraPrediction)
      ctxMode = 1;
    else
      ctxMode = 0;

    // [3] Context: parent
    if (isInter(parent))
      ctxMode += 2 * 3;
    else if (isIntra(parent) || false)
      ctxMode += 1 * 3;

    // [4] Context: neighbors
    if (isInter(neighbors))
      ctxMode += 3 * 3 * 3;
    else if (isIntra(neighbors))
      ctxMode += 2 * 3 * 3;
    else if (isNull(neighbors))
      ctxMode += 1 * 3 * 3;

    // [3] Context: number of children
    if (childrenCount > 5)
      ctxMode += 2 * 4 * 3 * 3;
    else if (childrenCount > 3)
      ctxMode += 1 * 4 * 3 * 3;

    if (enableIntraPrediction)
      return Mode::Intra;
    return Mode::Null;
  }

  int estimateExpGolombBits(unsigned int symbol, int k)
  {
    int bits = 0;

    while (symbol >= (1u << k)) {
      bits++;
      symbol -= 1u << k;
      k++;
    }
    bits++;

    bits += k;
    return bits;
  }

  int rdoReconstruct(FixedPoint& val, FixedPoint& recons, Quantizer& q)
  {
    auto coeff = q.quantize(val.round() << kFixedPointAttributeShift);

    if (coeff)
      recons += divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
    else
      return 0;

    if (coeff < 0)
      coeff = -coeff;
    coeff--;

    if (coeff)
      return estimateExpGolombBits(coeff - 1, 1) + 2;
    else
      return 2;
  }

  Mode choseMode(
    ModeEncoder& rdo,
    const VecAttr& transformBuf,
    const std::vector<Mode>& modes,
    const int64_t weights[],
    const int numAttrs,
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp)
  {
    auto coefReal = transformBuf.begin();
    auto coefPred = coefReal + numAttrs;

    VecAttr reconsBuf(numAttrs * (modes.size() + 1));
    auto recReal = reconsBuf.begin();
    auto recPred = recReal + numAttrs;

    // Estimate rate
    std::vector<int> rate;
    rate.resize(modes.size(), 0);

    for (int k = 0; k < numAttrs; k++) {
      std::copy(coefReal[k].begin(), coefReal[k].end(), recReal[k].begin());
      for (int mode = 0; mode < modes.size(); mode++)
        recPred[numAttrs * mode + k][0] = coefReal[k][0];
    }

    auto w = weights + 8 + 8 + 8;
    static const std::vector<int> kRahtScanOrder = {4, 2, 1, 6, 5, 3, 7};
    for (auto i : kRahtScanOrder) {
      if (!w[i])
        continue;

      auto quantizers = qpset.quantizers(qpLayer, nodeQp[i]);
      for (int k = 0; k < numAttrs; k++) {
        auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

        auto real = coefReal[k][i];
        for (int mode = 0; mode < modes.size(); mode++) {
          auto attr = real;
          auto& rec = recPred[numAttrs * mode + k][i];
          if (mode) {
            rec = coefPred[numAttrs * (mode - 1) + k][i];
            attr -= rec;
          }
          rate[mode] += rdoReconstruct(attr, rec, q);
        }
      }
    }

    // Estimate distortion
    std::vector<double> error;
    error.resize(modes.size(), 0.0);

    invTransformBlock222(reconsBuf.size(), reconsBuf.begin(), weights);

    w = weights;
    int64_t sum_weight = 0;
    int childCount = 0;
    for (int i = 0; i < 8; i++) {
      if (!w[i])
        continue;

      sum_weight += w[i];
      childCount++;

      for (int k = 0; k < numAttrs; k++) {
        auto real = reconsBuf[k][i].round();
        for (int mode = 0; mode < modes.size(); mode++) {
          double diff = static_cast<double>(
            recPred[numAttrs * mode + k][i].round() - real);
          error[mode] += diff * diff;
        }
      }
    }

    assert(sum_weight > 1);
    assert(childCount > 1);

    /* A value in the interval [3,5] seens good */
    double lambda = rdo.getLambda(4);

    std::vector<double> costFun;
    costFun.resize(modes.size());
    for (int mode = 0; mode < modes.size(); mode++) {
      error[mode] /= sum_weight;
      costFun[mode] = error[mode]
        + lambda * (double(rate[mode]) + rdo.entropy[modes[mode]])
          / childCount;
    }
    rdo.update(
      error[Mode::Null],
      (double(rate[Mode::Null]) + rdo.entropy[Mode::Null]) / childCount);

    for (int i = modes.size() - 1; i > 0; i--) {
      bool selected = true;
      int j = 0;
      while (selected && j < i)
        selected = costFun[i] < costFun[j++];
      if (selected)
        return modes[i];
    }
    return Mode::Null;
  }
}  // namespace attr

//============================================================================
// Encapsulation of a RAHT transform stage.

RahtKernel::RahtKernel(int weightLeft, int weightRight)
{
  uint64_t w = weightLeft + weightRight;
  uint64_t isqrtW = irsqrt(w);
  _a.val = (isqrt(uint64_t(weightLeft) << (2 * _a.kFracBits)) * isqrtW) >> 40;
  _b.val = (isqrt(uint64_t(weightRight) << (2 * _b.kFracBits)) * isqrtW) >> 40;
}

void
RahtKernel::fwdTransform(
  FixedPoint& lf, FixedPoint& hf)
{
  // lf = right * b + left * a
  // hf = right * a - left * b
  auto tmp = lf.val * _b.val;
  lf.val *= _a.val;
  lf.val += hf.val * _b.val;
  hf.val *= _a.val;
  hf.val -= tmp;

  lf.fixAfterMultiplication();
  hf.fixAfterMultiplication();
}

void
RahtKernel::invTransform(
  FixedPoint& left, FixedPoint& right)
{
  // left = lf * a - hf * b
  // right = lf * b + hf * a
  auto tmp = right.val * _b.val;
  right.val *= _a.val;
  right.val += left.val * _b.val;
  left.val *= _a.val;
  left.val -= tmp;

  right.fixAfterMultiplication();
  left.fixAfterMultiplication();
}

//============================================================================
// In-place transform a set of sparse 2x2x2 blocks each using the same weights

void
fwdTransformBlock222(
  const int numBufs, VecAttr::iterator buf,
  const int64_t weights[])
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
    for (int k = 0; k < numBufs; k++) {
      auto& lf = buf[k][i0];
      auto& hf = buf[k][i1];

      auto tmp = lf.val * w1;
      lf.val = lf.val * w0 + hf.val * w1;
      hf.val = hf.val * w0 - tmp;

      lf.fixAfterMultiplication();
      hf.fixAfterMultiplication();
    }
  }
}

void
invTransformBlock222(
  const int numBufs, VecAttr::iterator buf,
  const int64_t weights[])
{
  static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};
  static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
  for (int i = 11, iw = 22; i >= 0; i--, iw -= 2) {
    int i0 = a[i];
    int i1 = b[i];
    int64_t w0 = weights[iw + 32];
    int64_t w1 = weights[iw + 33];

    if (!a && !b)
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
    for (int k = 0; k < numBufs; k++) {
      auto& left = buf[k][i0];
      auto& right = buf[k][i1];

      auto tmp = right.val * w1;
      right.val = right.val * w0 + left.val * w1;
      left.val = left.val * w0 - tmp;

      right.fixAfterMultiplication();
      left.fixAfterMultiplication();
    }
  }
}
}  // namespace pcc
