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
  std::vector<FixedPoint>& parentDc,
  bool integer_haar_enable_flag
  )
{
  if (integer_haar_enable_flag) {
    for (int k = 0; k < parentDc.size(); k++) {
      parentDc[k].val = dc[k];
    }
    return;
  }

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
  int64_t* morton_mc,
  int* attr_mc,
  bool integer_haar_enable_flag,
  size_t layerSize)
{
  size_t shift = layerDepth * 3;
  std::vector<int64_t> morton_layer;
  if (layerSize)
    morton_layer.reserve(layerSize / attrCount);
  else
    morton_layer.reserve(count_rf);

  int64_t prev = -1;
  size_t i = 0;
  for (size_t n = 0; n < count_rf; n++) {
    int64_t curr = morton_rf[n] >> shift;
    if (curr != prev) {
      prev = curr;
      morton_layer.push_back(curr);
    }
  }

  count_rf = morton_layer.size();
  layerAttr.resize(count_rf * attrCount);

  i = 0;
  size_t j = 0;
  while (i < count_rf && j < count_mc) {
    prev = morton_layer[i];

    while (j < count_mc && prev > (morton_mc[j] >> shift))
      j++;

    int64_t weight = 0;
    auto layer = std::next(layerAttr.begin(), attrCount * i);
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
      for (size_t k = 0; k < attrCount; k++) {
        layer[k] /= weight;
        if (integer_haar_enable_flag)
          layer[k] = (layer[k] >> FixedPoint::kFracBits) << FixedPoint::kFracBits;
      }
    } else if (!weight) {
      for (size_t k = 0; k < attrCount; k++)
        layer[k] = -1;
    }

    i++;
    while (i < count_rf && prev == morton_layer[i]) {
      std::copy(
        layer, layer + attrCount,
        std::next(layerAttr.begin(), attrCount * i));
      i++;
    }
  }
}


namespace attr {
  Mode getNeighborsMode(
    const bool& isEncoder,
    const int parentNeighIdx[19],
    const std::vector<UrahtNode>& weightsParent,
    int16_t& voteInterWeight,
    int16_t& voteIntraWeight,
    int16_t& voteInterLayerWeight,
    int16_t& voteIntraLayerWeight)
  {
    int voteNull = 0;
    int voteIntra = 0;
    int voteInter = 0;

    int voteNullLayer = 0;
    int voteIntraLayer = 0;
    int voteInterLayer = 0;

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

      if (isEncoder) {
        if (isNull(uncle->mode))
          voteNullLayer += 1;
        else if (isIntra(uncle->mode))
          voteIntraLayer += 1;
        else if (isInter(uncle->mode))
          voteInterLayer += 1;
      }

      if (uncle->decoded) {
        for (auto cousin = uncle->firstChild; cousin < uncle->lastChild;
             cousin++) {
          if (isNull(cousin->mode))
            voteNull += 3;
          else if (isIntra(cousin->mode))
            voteIntra += 3;
          else if (isInter(cousin->mode))
            voteInter += 3;

          if (isEncoder) {
            if (isNull(cousin->_mode))
              voteNullLayer += 3;
            else if (isIntra(cousin->_mode))
              voteIntraLayer += 3;
            else if (isInter(cousin->_mode))
              voteInterLayer += 3;
          }
        }
      }
    }

    voteInterWeight = voteInter * 2 + voteNull;
    voteIntraWeight = voteIntra * 2 + voteNull;
    if (isEncoder) {
      voteInterLayerWeight = voteInterLayer * 2 + voteNullLayer;
      voteIntraLayerWeight = voteIntraLayer * 2 + voteNullLayer;
    }

    if (1) {
      auto parent = std::next(weightsParent.begin(), parentNeighIdx[0]);

      if (isNull(parent->mode)) {
        voteIntraWeight += 1;
        voteInterWeight += 1;
      } else if (isIntra(parent->mode))
        voteIntraWeight += 2;
      else if (isInter(parent->mode))
        voteInterWeight += 2;

      if (isEncoder) {
        if (isNull(parent->mode)) {
          voteIntraLayerWeight += 1;
          voteInterLayerWeight += 1;
        } else if (isIntra(parent->mode))
          voteIntraLayerWeight += 2;
        else if (isInter(parent->mode))
          voteInterLayerWeight += 2;
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

  template<class Kernel>
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
    bool lossless = typeid(Kernel) == typeid(HaarKernel);

    VecAttr reconsBuf(numAttrs * (modes.size() + 1));
    auto recReal = reconsBuf.begin();
    auto recPred = recReal + numAttrs;

    // Estimate rate
    std::vector<int> rate;
    rate.resize(modes.size(), 0);
    std::vector<double> error;
    error.resize(modes.size(), 0.0);

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
          if(lossless){
            rate[mode] += rdoReconstruct(attr, rec, q);
          }
          else{
            rec = 0;
            rate[mode] += rdoReconstruct(attr, rec, q);
            double diff = static_cast<double>(attr.round() - rec.round());
            error[mode] += diff * diff;
          }
        }
      }
    }

    w = weights;
    int64_t sum_weight = 0;
    int childCount = 0;
    for (int i = 0; i < 8; i++) {
      if (!w[i])
        continue;

      sum_weight += w[i];
      childCount++;
    }

    assert(sum_weight > 1);
    assert(childCount > 1);

    /* A value in the interval [3,5] seens good */
    double lambda = lossless ? 1. : rdo.getLambda(4);

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

//============================================================================
// instanciate for Haar and Raht kernels

template Mode choseMode<HaarKernel>(
    ModeEncoder& rdo,
    const VecAttr& transformBuf,
    const std::vector<Mode>& modes,
    const int64_t weights[],
    const int numAttrs,
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp);

template Mode choseMode<RahtKernel>(
    ModeEncoder& rdo,
    const VecAttr& transformBuf,
    const std::vector<Mode>& modes,
    const int64_t weights[],
    const int numAttrs,
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp);

//============================================================================

}  // namespace attr
}  // namespace pcc
