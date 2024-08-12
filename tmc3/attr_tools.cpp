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
  Mode getNeighborsModeEncoder(
    const int parentNeighIdx[19],
    const std::vector<UrahtNodeEncoder>& weightsParent,
    int& voteInterWeight,
    int& voteIntraWeight,
    int& voteInterLayerWeight,
    int& voteIntraLayerWeight,
    std::vector<UrahtNodeEncoder>& weightsLf)
  {
    int vote[4] = { 0,0,0,0 }; // Null, intra, inter, size;
    int voteLayer[4] = { 0,0,0,0 }; // Null, intra, inter, size;

    for (int i = 1; i < 19; i++) {
      if (parentNeighIdx[i] == -1)
        continue;

      auto uncle = std::next(weightsParent.begin(), parentNeighIdx[i]);
      vote[uncle->mode]++;
      voteLayer[uncle->mode]++;

      if (uncle->decoded) {
        auto cousin = &weightsLf[uncle->firstChildIdx];
        for (int t = 0; t < uncle->numChildren; t++, cousin++) {
          vote[cousin->mode] += 3;
          voteLayer[cousin->_mode] += 3;
        }
      }
    }

    auto parent = std::next(weightsParent.begin(), parentNeighIdx[0]);

    voteIntraWeight = vote[1] * 2 + vote[0];
    voteInterWeight = vote[2] * 2 + vote[0];
    voteIntraWeight += isNull(parent->mode) + 2 * isIntra(parent->mode);
    voteInterWeight += isNull(parent->mode) + 2 * isInter(parent->mode);

    voteIntraLayerWeight = voteLayer[1] * 2 + voteLayer[0];
    voteInterLayerWeight = voteLayer[2] * 2 + voteLayer[0];
    voteIntraLayerWeight += isNull(parent->mode) + 2 * isIntra(parent->mode);
    voteInterLayerWeight += isNull(parent->mode) + 2 * isInter(parent->mode);

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
    int bits = -k;
    while (symbol >= (1u << k)) {
      symbol -= 1u << k;
      k++;
    }
    return bits + 2 * k + 1;
  }

  inline
  int rdoReconstruct(int64_t coeff)
  {
    static const int LUT_RDO[20] = {
      0, 2, 4, 4, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10};

    if (coeff < 20)
      return LUT_RDO[coeff];
    else
      return estimateExpGolombBits(coeff - 2, 1) + 2;
  }

  template<class Kernel, int numAttrs, typename WeightsType>
  Mode choseMode(
    ModeEncoder& rdo,
    const int64_t* transformBuf,
    const std::vector<Mode>& modes,
    const WeightsType weights[],
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp,
    const bool inheritDC)
  {
    auto coefReal = transformBuf;
    auto coefPred = coefReal + 8 * numAttrs;
    constexpr bool lossless = std::is_same<Kernel, HaarKernel>::value;

    int64_t reconsBuf[8 * numAttrs * (Mode::size + 1)];
    auto recReal = reconsBuf;
    auto recPred = recReal + 8 * numAttrs;

    // Estimate rate
    const int NMode = modes.size();
    std::vector<int> rate(NMode, 0);
    std::vector<int64_t> error(NMode, 0);

    std::copy_n(coefReal, 8 * numAttrs, recReal);
    for (int mode = !inheritDC; mode < NMode; mode++)
      for (int k = 0; k < numAttrs; k++)
        recPred[8 * numAttrs * mode + k] = coefReal[8 * k];

    auto w = weights + 8 + 8 + 8;
    for (int i = 1; i < 8; i++ ) { //only AC coeff
      if (w[i]) {
        auto quantizers = qpset.quantizers(qpLayer, nodeQp[i]);

        for (int k = 0; k < numAttrs; k++) {
          const auto real = coefReal[8 * k + i];
          for (int mode = 0; mode < NMode; mode++) {
            auto attr = real;
            auto& rec = recPred[8 * (numAttrs * mode + k) + i];
            if (mode) {
              rec = coefPred[8 * (numAttrs * (mode - 1) + k) + i];
              attr -= rec;
            }

            if (lossless) {
              auto coeff = fpReduce<kFPFracBits>(attr);
              rate[mode] += rdoReconstruct(std::abs(coeff));
            }
            else {
              auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];
              auto coeff = q.quantize(
                fpReduce<kFPFracBits>(attr) << kFixedPointAttributeShift);
              rate[mode] += rdoReconstruct(std::abs(coeff));
              rec = divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
              int64_t diff = fpReduce<kFPFracBits>(attr) - rec;
              error[mode] += diff * diff;
            }
          }
        }
      }
    }

    int64_t sum_weight = 0;
    int childCount = 0;
    for (int i = 0; i < 8; i++) {
      sum_weight += weights[i];
      childCount += !!weights[i];
    }

    assert(sum_weight > 1);
    assert(childCount > 1);
    const double inv2pow32 = 2.328306436538696e-10;
    const double inv_sum_weight = inv2pow32 / double(sum_weight);

    static const double invChildOver2pow32[9] = {
      0,                     2.328306436538696e-10, 1.164153218269348e-10,
      7.761021455128987e-11, 5.820766091346741e-11, 4.656612873077393e-11,
      3.880510727564494e-11, 3.326152052198137e-11, 2.910383045673370e-11 };
    const double inv_childCount = invChildOver2pow32[childCount];

    std::vector<double> rate_d(NMode);
    std::vector<double> error_d(NMode);
    for (int mode = 0; mode < NMode; mode++) {
      rate_d[mode] = inv_childCount * double(
        (int64_t(rate[mode]) << 32) + rdo.entropy[modes[mode]]);
      error_d[mode] = double(error[mode]) * inv_sum_weight;
    }


    rdo.update(error_d[Mode::Null], rate_d[Mode::Null]);
    // A value in the interval [3,5] seens good
    double lambda = lossless ? 1. : rdo.getLambda(4);

    // sleect best mode
    int selectedModeIdx = 0;
    double minCost = error_d[0] + lambda * rate_d[0];
    for (int mode = 1; mode <= NMode - 1; mode++) {
      double costFun = error_d[mode] + lambda * rate_d[mode];
      if (costFun < minCost) {
        minCost = costFun;
        selectedModeIdx = mode;
      }
    }

    return modes[selectedModeIdx];
  }

//============================================================================
// instanciate for Haar and Raht kernels with 3 attributes

template Mode choseMode<HaarKernel, 3, bool>(
    ModeEncoder& rdo,
    const int64_t* transformBuf,
    const std::vector<Mode>& modes,
    const bool weights[],
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp,
    const bool inheritDC);

template Mode choseMode<RahtKernel, 3, int64_t>(
    ModeEncoder& rdo,
    const int64_t* transformBuf,
    const std::vector<Mode>& modes,
    const int64_t weights[],
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp,
    const bool inheritDC);

#if defined _MSC_VER && defined _DEBUG
// These instanciations are not used but required by MSVC in debug mode
template Mode choseMode<HaarKernel, 3, int64_t>(
    ModeEncoder& rdo,
    const int64_t* transformBuf,
    const std::vector<Mode>& modes,
    const int64_t weights[],
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp,
    const bool inheritDC);

template Mode choseMode<RahtKernel, 3, bool>(
    ModeEncoder& rdo,
    const int64_t* transformBuf,
    const std::vector<Mode>& modes,
    const bool weights[],
    const QpSet& qpset,
    const int qpLayer,
    const Qps* nodeQp,
    const bool inheritDC);
#endif

//============================================================================

}  // namespace attr
}  // namespace pcc
