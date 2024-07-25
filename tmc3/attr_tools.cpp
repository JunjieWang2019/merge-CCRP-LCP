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

using namespace RAHT;

namespace attr {
  Mode getNeighborsMode(
    const bool& isEncoder,
    const int parentNeighIdx[19],
    const std::vector<UrahtNode>& weightsParent,
    int& voteInterWeight,
    int& voteIntraWeight,
    int& voteInterLayerWeight,
    int& voteIntraLayerWeight)
  {
    int vote[4] = { 0,0,0,0 }; // Null, intra, inter, size;
    int voteLayer[4] = { 0,0,0,0 }; // Null, intra, inter, size;

    for (int i = 1; i < 19; i++) {
      if (parentNeighIdx[i] == -1)
        continue;

      auto uncle = std::next(weightsParent.begin(), parentNeighIdx[i]);

      vote[uncle->mode]++;

      if (isEncoder)
        voteLayer[uncle->mode]++;

      if (uncle->decoded) {
        auto cousin = uncle->firstChild;
        for (int t = 0; t < uncle->numChildren; t++, cousin++) {
          vote[cousin->mode] += 3;
          if (isEncoder)
            voteLayer[cousin->_mode] += 3;
        }
      }
    }

    auto parent = std::next(weightsParent.begin(), parentNeighIdx[0]);

    voteIntraWeight = vote[1] * 2 + vote[0];
    voteInterWeight = vote[2] * 2 + vote[0];
    voteIntraWeight += isNull(parent->mode) + 2 * isIntra(parent->mode);
    voteInterWeight += isNull(parent->mode) + 2 * isInter(parent->mode);

    if (isEncoder) {
      voteIntraLayerWeight = voteLayer[1] * 2 + voteLayer[0];
      voteInterLayerWeight = voteLayer[2] * 2 + voteLayer[0];
      voteIntraLayerWeight += isNull(parent->mode) + 2 * isIntra(parent->mode);
      voteInterLayerWeight += isNull(parent->mode) + 2 * isInter(parent->mode);
    }

    if (vote[0] > vote[1] && vote[0] > vote[2])
      return Mode::Null;
    if (vote[1] > vote[2])
      return Mode::Intra;
    return Mode::Inter;
  }

  int getInferredMode(
    bool enableIntraPred,
    bool enableInterPred,
    int childrenCount,
    Mode parent,
    Mode neighbors)
  {
    // [3] Context: enabled flag
    int ctxMode = enableIntraPred << enableInterPred;

    // [3] Context: parent
    ctxMode += 2 * 3 * isInter(parent);
    ctxMode += 1 * 3 * isIntra(parent);

    // [3] Context: neighbors
    ctxMode += 2 * 3 * 3 * isInter(neighbors);
    ctxMode += 1 * 3 * 3 * isIntra(neighbors);

    // [3] Context: number of children
    ctxMode += 1 * 3 * 3 * 3 * ((childrenCount > 3) + (childrenCount > 5));

    return (ctxMode << 1) + enableIntraPred;
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

  int rdoReconstruct(int64_t& val, int64_t& recons, Quantizer& q)
  {
    auto coeff = q.quantize(fpReduce<kFPFracBits>(val) << kFixedPointAttributeShift);

    if (coeff)
      recons += fpExpand<kFPFracBits>(
        divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift));
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
    const int64_t* transformBuf,
    const std::vector<Mode>& modes,
    const int64_t weights[],
    const int numAttrs,
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp,
    const bool inheritDC)
  {
    auto coefReal = transformBuf;
    auto coefPred = coefReal + 8 * numAttrs;
    bool lossless = typeid(Kernel) == typeid(HaarKernel);

    int64_t reconsBuf[8 * numAttrs * (Mode::size + 1)];
    auto recReal = reconsBuf;
    auto recPred = recReal + 8 * numAttrs;

    // Estimate rate
    std::vector<int> rate;
    rate.resize(modes.size(), 0);
    std::vector<double> error;
    error.resize(modes.size(), 0.0);

    std::copy_n(coefReal, 8 * numAttrs, recReal);
    for (int k = 0; k < numAttrs; k++) {
      for (int mode = !inheritDC; mode < modes.size(); mode++)
        recPred[8 * numAttrs * mode + k] = coefReal[8 * k];
    }

    auto w = weights + 8 + 8 + 8;
    static const std::vector<int> kRahtScanOrder = {4, 2, 1, 6, 5, 3, 7};
    static const std::vector<int> kRahtScanOrderRoot = {0, 4, 2, 1, 6, 5, 3, 7};
    for (auto i : inheritDC ? kRahtScanOrder : kRahtScanOrderRoot) {
      if (!w[i])
        continue;

      auto quantizers = qpset.quantizers(qpLayer, nodeQp[i]);
      for (int k = 0; k < numAttrs; k++) {
        auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

        auto real = coefReal[8 * k + i];
        for (int mode = 0; mode < modes.size(); mode++) {
          auto attr = real;
          auto& rec = recPred[8 * (numAttrs * mode + k) + i];
          if (mode) {
            rec = coefPred[8 * (numAttrs * (mode - 1) + k) + i];
            attr -= rec;
          }
          if(lossless){
            rate[mode] += rdoReconstruct(attr, rec, q);
          }
          else{
            rec = 0;
            rate[mode] += rdoReconstruct(attr, rec, q);
            double diff = static_cast<double>(
              fpReduce<kFPFracBits>(attr) - fpReduce<kFPFracBits>(rec));
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

    rdo.update(
      error[Mode::Null] / sum_weight,
      fpToDouble<32>(
        ((int64_t(rate[Mode::Null]) << 32) + rdo.entropy[modes[Mode::Null]])
          / childCount));
    /* A value in the interval [3,5] seens good */
    double lambda = lossless ? 1. : rdo.getLambda(4);

    std::vector<double> costFun;
    costFun.resize(modes.size());
    for (int mode = 0; mode < modes.size(); mode++) {
      error[mode] /= sum_weight;
      costFun[mode] = error[mode]
        + lambda * fpToDouble<32>(
          ((int64_t(rate[mode]) << 32) + rdo.entropy[modes[mode]])
            / childCount);
    }

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
    const int64_t* transformBuf,
    const std::vector<Mode>& modes,
    const int64_t weights[],
    const int numAttrs,
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp,
    const bool inheritDC);

template Mode choseMode<RahtKernel>(
    ModeEncoder& rdo,
    const int64_t* transformBuf,
    const std::vector<Mode>& modes,
    const int64_t weights[],
    const int numAttrs,
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp,
    const bool inheritDC);

//============================================================================

}  // namespace attr
}  // namespace pcc
