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
#include "attr_tools.h"

using pcc::attr::Mode;
using namespace pcc::RAHT;

//============================================================================

// Return the value of `level / 3`,
// where `level` must satisfy: `0 <= level < 512`
template<class T>
inline T div3level(T level)
{
  assert(level >= 0 && level < 512);
  return level * 171 >> 9;
}

// Return the value of `level % 3`, and compute the value `layer = level / 3`
// where `level` must satisfy: `0 <= level < 512`
template<class T>
inline T rem3level(T level, T& layer)
{
  layer = div3level(level);
  return level - layer * 3;
}

namespace pcc {

//============================================================================

struct PCCRAHTACCoefficientEntropyEstimate {
  PCCRAHTACCoefficientEntropyEstimate()
  {
    init();
  }

  PCCRAHTACCoefficientEntropyEstimate(
    const PCCRAHTACCoefficientEntropyEstimate& other) = default;

  PCCRAHTACCoefficientEntropyEstimate&
    operator=(const PCCRAHTACCoefficientEntropyEstimate&) = default;

  void resStatUpdate(int32_t values, int k);
  void init();
  void updateCostBits(int32_t values, int k);
  int64_t costBits() { return sumCostBits; }
  void resetCostBits() { sumCostBits = 0; }

private:
  // Encoder side residual cost calculation
  static constexpr int log2scaleRes = 20;
  static constexpr unsigned scaleRes = 1 << log2scaleRes;
  static constexpr unsigned windowLog2 = 6;
  int probResGt0[3];  //prob of residuals larger than 0: 1 for each component
  int probResGt1[3];  //prob of residuals larger than 1: 1 for each component
  int64_t sumCostBits;
};

//============================================================================

void
PCCRAHTACCoefficientEntropyEstimate::init()
{
  for (int k = 0; k < 3; k++)
    probResGt0[k] = probResGt1[k] = (scaleRes >> 1);
  sumCostBits = 0;
}

//---------------------------------------------------------------------------

void
PCCRAHTACCoefficientEntropyEstimate::updateCostBits(int32_t value, int k)
{
  int64_t bits = 0;
  bits += value ? -fpLog2<log2scaleRes, 32>(probResGt0[k])
    : -fpLog2<log2scaleRes, 32>(scaleRes - probResGt0[k]);  //Gt0
  int mag = abs(value);
  if (mag) {
    bits += mag > 1 ? -fpLog2<log2scaleRes, 32>(probResGt1[k])
      : -fpLog2<log2scaleRes, 32>(scaleRes - probResGt1[k]);  //Gt1
    bits += 1ULL << 32;  //sign bit.
    if (mag > 1)
      bits += 2 * fpLog2<0, 32>(mag - 1.0) + 1.0;  //EG0 approximation.
  }
  sumCostBits += bits;
}

//----------------------------------------------------------------------------

void
PCCRAHTACCoefficientEntropyEstimate::resStatUpdate(int32_t value, int k)
{
  probResGt0[k] += value ? scaleRes - probResGt0[k] >> windowLog2
    : -(probResGt0[k] >> windowLog2);
  if (value)
    probResGt1[k] += abs(value) > 1 ? scaleRes - probResGt1[k] >> windowLog2
      : -(probResGt1[k] >> windowLog2);
}

//============================================================================

int8_t
PCCRAHTComputeCCCP::computeCrossChromaComponentPredictionCoeff(
  int m, int64_t coeffs[][3])
{
  Elt sum12 {0, 0};

  for (size_t coeffIdx = 0; coeffIdx < m; ++coeffIdx) {
    const auto& attr = coeffs[coeffIdx];
    sum12.k1k2 += attr[1] * attr[2];
    sum12.k1k1 += attr[1] * attr[1];
  }

  if (window.size() == 128) {
    const auto& removed = window.front();
    sum.k1k2 -= removed.k1k2;
    sum.k1k1 -= removed.k1k1;
    window.pop_front();
  }

  sum.k1k2 += sum12.k1k2;
  sum.k1k1 += sum12.k1k1;
  window.push(sum12);

  int scale = 0;

  if (sum.k1k2 && sum.k1k1) {
    // sign(sum.k1k2) * sign(sum.k1k1)
    scale = divApproxRoundHalfInf(sum.k1k2, sum.k1k1, 4);
  }

  // NB: coding range is limited to +-4
  return PCCClip(scale, -16, 16);
}

//============================================================================
// remove any non-unique leaves from a level in the uraht tree
template<bool haarFlag, int numAttrs>
int
reduceUnique(
  int numNodes,
  std::vector<UrahtNodeLight>* weightsIn,
  std::vector<UrahtNodeLight>* weightsOut,
  std::vector<int64_t>* attrsIn,
  std::vector<int64_t>* attrsOut)
{
  // process a single level of the tree
  int64_t posPrev = -1;
  auto weightsInWrIt = weightsIn->begin();
  auto weightsInRdIt = weightsIn->cbegin();
  auto attrsInWrIt = attrsIn->begin();
  auto attrsInRdIt = attrsIn->begin();
  for (int i = 0; i < numNodes; i++) {
    const auto& node = *weightsInRdIt++;

    // copy across unique nodes
    if (node.pos != posPrev) {
      posPrev = node.pos;
      *weightsInWrIt++ = node;
      // Encoder only
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;

      continue;
    }

    // duplicate node
    (weightsInWrIt - 1)->weight += node.weight;
    weightsOut->push_back(node);

    // Encoder only
    for (int k = 0; k < numAttrs; k++) {
      if (haarFlag) {
        attrsOut->push_back(
          *attrsInRdIt++ - *(attrsInWrIt - numAttrs + k));
        *(attrsInWrIt - numAttrs + k) +=
          attrsOut->back() >> 1;
      } else {
        *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;
        attrsOut->push_back(*attrsInRdIt++);
      }
    }
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================
// Split a level of values into sum and difference pairs.
template<bool haarFlag, int numAttrs>
int
reduceLevel(
  int level,
  int numNodes,
  std::vector<UrahtNodeLight>* weightsIn,
  std::vector<UrahtNodeLight>* weightsOut,
  std::vector<int64_t>* attrsIn,
  std::vector<int64_t>* attrsOut)
{
  // process a single level of the tree
  int64_t posPrev = -1;
  auto weightsInWrIt = weightsIn->begin();
  auto weightsInRdIt = weightsIn->cbegin();
  auto attrsInWrIt = attrsIn->begin();
  auto attrsInRdIt = attrsIn->begin();
  for (int i = 0; i < numNodes; i++) {
    auto& node = *weightsInRdIt++;
    bool newPair = (posPrev ^ node.pos) >> level != 0;
    posPrev = node.pos;
    if (newPair) {
      *weightsInWrIt++ = node;
      // Encoder only
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;
    } else {
      auto& left = *(weightsInWrIt - 1);
      left.weight += node.weight;
      left.qp[0] = (left.qp[0] + node.qp[0]) >> 1;
      left.qp[1] = (left.qp[1] + node.qp[1]) >> 1;
      weightsOut->push_back(node);
      // Encoder only
      for (int k = 0; k < numAttrs; k++) {
        if (haarFlag) {
          attrsOut->push_back(*attrsInRdIt++ - *(attrsInWrIt - numAttrs + k));
          *(attrsInWrIt - numAttrs + k) += attrsOut->back() >> 1;
        } else {
          *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;
          attrsOut->push_back(*attrsInRdIt++);
        }
      }
    }
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================
// Merge sum and difference values to form a tree.
template<bool haarFlag, int numAttrs>
void
expandLevel(
  int level,
  int numNodes,
  std::vector<UrahtNodeLight>* weightsIn,   // expand by numNodes before expand
  std::vector<UrahtNodeLight>* weightsOut,  // shrink after expand
  std::vector<int64_t>* attrsIn,
  std::vector<int64_t>* attrsOut)
{
  if (numNodes == 0)
    return;

  // process a single level of the tree
  auto weightsInWrIt = weightsIn->rbegin();
  auto weightsInRdIt = std::next(weightsIn->crbegin(), numNodes);
  auto weightsOutRdIt = weightsOut->crbegin();
  auto attrsInWrIt = attrsIn->rbegin();
  auto attrsInRdIt = std::next(attrsIn->crbegin(), numNodes * numAttrs);
  auto attrsOutRdIt = attrsOut->crbegin();
  for (int i = 0; i < numNodes;) {
    bool isPair = (weightsOutRdIt->pos ^ weightsInRdIt->pos) >> level == 0;
    if (!isPair) {
      weightsInWrIt->mode = Mode::Null;
      *weightsInWrIt++ = *weightsInRdIt++;
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;
      continue;
    }

    // going to process a pair
    i++;

    // Out node is inserted before In node.
    const auto& nodeDelta = *weightsInWrIt++ = *weightsOutRdIt++;
    auto curAttrIt = attrsInWrIt;
    for (int k = 0; k < numAttrs; k++)
      *attrsInWrIt++ = *attrsOutRdIt++;

    // move In node to correct position, subtracting delta
    *weightsInWrIt = *weightsInRdIt++;
    weightsInWrIt->mode = Mode::Null;
    (weightsInWrIt++)->weight -= nodeDelta.weight;
    for (int k = 0; k < numAttrs; k++) {
      *attrsInWrIt = *attrsInRdIt++;
      if (haarFlag) {
        *attrsInWrIt -= *curAttrIt >> 1;
        *curAttrIt++ += *attrsInWrIt++;
      }
      else {
        *attrsInWrIt++ -= *curAttrIt++;
      }
    }
  }
}

//============================================================================

template<typename It, typename It2>
void
findNeighboursChildren(
  It first,
  It2 firstChild,
  int level,
  uint8_t occupancy,
  const int parentNeighIdx[19],
  int childNeighIdx[12][8])
{
  memset(childNeighIdx, -1, 96 * sizeof(int));

  static const uint8_t occuMasks[12] = {3,  5,  15, 17, 51, 85,
                                        10, 34, 12, 68, 48, 80};
  static const uint8_t occuShift[12] = {6, 5, 4, 3, 2, 1, 3, 1, 2, 1, 2, 3};

  for (int i = 0; i < 9; i++) {
    if (parentNeighIdx[7 + i] == -1)
      continue;

    auto neiIt = first + parentNeighIdx[7 + i];
    uint8_t mask =
      (neiIt->occupancy >> occuShift[i]) & occupancy & occuMasks[i];
    if (!mask)
      continue;
    auto it = neiIt->firstChild;
    for (int t = 0; t < neiIt->numChildren; t++, it++) {
      int nodeIdx = ((it->pos >> level) & 0x7) - occuShift[i];
      if ((nodeIdx >= 0) && ((mask >> nodeIdx) & 1)) {
        childNeighIdx[i][nodeIdx] = std::distance(firstChild, it);
      }
    }
  }

  for (int i = 9; i < 12; i++) {
    if (parentNeighIdx[7 + i] == -1)
      continue;

    auto neiIt = first + parentNeighIdx[7 + i];
    uint8_t mask =
      (neiIt->occupancy << occuShift[i]) & occupancy & occuMasks[i];
    if (!mask)
      continue;

    auto it = neiIt->firstChild;
    for (int t = 0; t < neiIt->numChildren; t++, it++) {
      int nodeIdx = ((it->pos >> level) & 0x7) + occuShift[i];
      if ((nodeIdx < 8) && ((mask >> nodeIdx) & 1)) {
        childNeighIdx[i][nodeIdx] = std::distance(firstChild, it);
      }
    }
  }
}

//============================================================================
// Generate the spatial prediction of a block.

template<bool haarFlag, int numAttrs, typename It, typename It2 >
void
intraDcPred(
  bool isEncoder,
  const int parentNeighIdx[19],
  const int childNeighIdx[12][8],
  int occupancy,
  It first,
  It firstChild,
  It2 intraLayerFirstChild,
  It2 interLayerFirstChild,
  int64_t* predBuf,
  int64_t* intraLayerPredBuf,
  int64_t* interLayerPredBuf,
  const RahtPredictionParams &rahtPredParams,
  const bool& enableLayerCoding)
{
  static const uint8_t predMasks[19] = {255, 240, 204, 170, 192, 160, 136,
                                        3,   5,   15,  17,  51,  85,  10,
                                        34,  12,  68,  48,  80};

  const auto& predWeightParent = rahtPredParams.predWeightParent;
  const auto& predWeightChild = rahtPredParams.predWeightChild;
  int weightSum[8] = {0, 0, 0, 0, 0, 0, 0, 0};

  std::fill_n(predBuf, 8 * numAttrs, 0);

  if (isEncoder && enableLayerCoding) {
    std::fill_n(intraLayerPredBuf, 8 * numAttrs, 0);
    std::fill_n(interLayerPredBuf, 8 * numAttrs, 0);
  }

  int64_t neighValue[3];
  int64_t childNeighValue[3];
  int64_t intraLayerChildNeighValue[3];
  int64_t interLayerChildNeighValue[3];
  int64_t limitLow = 0;
  int64_t limitHigh = 0;

  const auto parentOnlyCheckMaxIdx =
    rahtPredParams.subnode_prediction_enabled_flag ? 7 : 19;
  for (int i = 0; i < parentOnlyCheckMaxIdx; i++) {
    if (parentNeighIdx[i] == -1)
      continue;

    auto neighValueIt = std::next(first, numAttrs * parentNeighIdx[i]);
    for (int k = 0; k < numAttrs; k++)
      neighValue[k] = *neighValueIt++;

    // skip neighbours that are outside of threshold limits
    if (i) {
      if (10 * neighValue[0] <= limitLow || 10 * neighValue[0] >= limitHigh)
        continue;
    } else {
      constexpr int ratioThreshold1 = 2;
      constexpr int ratioThreshold2 = 25;
      limitLow = ratioThreshold1 * neighValue[0];
      limitHigh = ratioThreshold2 * neighValue[0];
    }

    // apply weighted neighbour value to masked positions
    for (int k = 0; k < numAttrs; k++)
      neighValue[k] *= predWeightParent[i];

    int mask = predMasks[i] & occupancy;
    for (int j = 0; mask; j++, mask >>= 1) {
      if (mask & 1) {
        weightSum[j] += predWeightParent[i];
        for (int k = 0; k < numAttrs; k++) {
          predBuf[8 * k + j] += neighValue[k];
          if (isEncoder && enableLayerCoding) {
            intraLayerPredBuf[8 * k + j] += neighValue[k];
            interLayerPredBuf[8 * k + j] += neighValue[k];
          }
        }
      }
    }
  }
  if (rahtPredParams.subnode_prediction_enabled_flag) {
    for (int i = 0; i < 12; i++) {
      if (parentNeighIdx[7 + i] == -1)
        continue;

      auto neighValueIt = std::next(first, numAttrs * parentNeighIdx[7 + i]);
      for (int k = 0; k < numAttrs; k++)
        neighValue[k] = *neighValueIt++;

      // skip neighbours that are outside of threshold limits
      if (10 * neighValue[0] <= limitLow || 10 * neighValue[0] >= limitHigh)
        continue;

      // apply weighted neighbour value to masked positions
      for (int k = 0; k < numAttrs; k++)
        neighValue[k] *= predWeightParent[7 + i];

      int mask = predMasks[7 + i] & occupancy;
      for (int j = 0; mask; j++, mask >>= 1) {
        if (mask & 1) {
          if (childNeighIdx[i][j] != -1) {
            weightSum[j] += predWeightChild[i];
            auto childNeighValueIt =
              std::next(firstChild, numAttrs * childNeighIdx[i][j]);
            for (int k = 0; k < numAttrs; k++)
              childNeighValue[k] = (*childNeighValueIt++)
                * predWeightChild[i];

            for (int k = 0; k < numAttrs; k++)
              predBuf[8 * k + j] += childNeighValue[k];

            if (isEncoder && enableLayerCoding) {
              auto intraChildNeighValueIt =
                std::next(intraLayerFirstChild, numAttrs * childNeighIdx[i][j]);
              for (int k = 0; k < numAttrs; k++)
                intraLayerChildNeighValue[k] =
                (*intraChildNeighValueIt++) * predWeightChild[i];
              for (int k = 0; k < numAttrs; k++)
                intraLayerPredBuf[8 * k + j] += intraLayerChildNeighValue[k];

              auto interChildNeighValueIt =
                std::next(interLayerFirstChild, numAttrs * childNeighIdx[i][j]);
              for (int k = 0; k < numAttrs; k++)
                interLayerChildNeighValue[k] =
                (*interChildNeighValueIt++) * predWeightChild[i];
              for (int k = 0; k < numAttrs; k++) {
                interLayerPredBuf[8 * k + j] += interLayerChildNeighValue[k];
              }
            }
          } else {
            weightSum[j] += predWeightParent[7 + i];
            for (int k = 0; k < numAttrs; k++) {
              predBuf[8 * k + j] += neighValue[k];
              if (isEncoder && enableLayerCoding) {
                intraLayerPredBuf[8 * k + j] += neighValue[k];
                interLayerPredBuf[8 * k + j] += neighValue[k];
              }
            }
          }
        }
      }
    }
  }

  // normalise
  for (int i = 0; i < 8; i++, occupancy >>= 1) {
    if (occupancy & 1) {
      int w = weightSum[i];
      if (w > 1) {
        ApproxNormalize div(w);
        for (int k = 0; k < numAttrs; k++) {
          div(predBuf[8 * k + i]);
          if (isEncoder && enableLayerCoding) {
            div(intraLayerPredBuf[8 * k + i]);
            div(interLayerPredBuf[8 * k + i]);
          }
        }
      }
      if (haarFlag) {
        for (int k = 0; k < numAttrs; k++) {
          predBuf[8 * k + i] = fpExpand<kFPFracBits>(
            fpReduce<kFPFracBits>(predBuf[8 * k + i]));
          if (isEncoder && enableLayerCoding) {
            intraLayerPredBuf[8 * k + i] = fpExpand<kFPFracBits>(
              fpReduce<kFPFracBits>(intraLayerPredBuf[8 * k + i]));
            interLayerPredBuf[8 * k + i] = fpExpand<kFPFracBits>(
              fpReduce<kFPFracBits>(interLayerPredBuf[8 * k + i]));
          }
        }
      }
    }
  }
}

//============================================================================
// Tests if two positions are siblings at the given tree level

inline static bool
isSibling(int64_t pos0, int64_t pos1, int level)
{
  return ((pos0 ^ pos1) >> level) == 0;
}

//============================================================================
int getRate(int trainZeros)
{
  static const int LUTbins[11] = { 1,2,3, 5,5, 7,7, 9,9 ,11 ,11 };
  int Rate = LUTbins[trainZeros > 10 ? 10 : trainZeros];
  if (trainZeros > 10) {
    int temp = trainZeros - 11;
    // prefix k =2
    temp += 1;
    int a = 0;
    while (temp) {
      a++;
      temp >>= 1;

    }
    Rate += 2 * a - 1;
    // suffix  k=2
    Rate += 2;
  }
  return Rate;
}

//============================================================================
// Core transform process (for encoder)

template<bool haarFlag, int numAttrs, class GetMode>
inline void
uraht_process_encoder(
  const RahtPredictionParams& rahtPredParams,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int numPoints,
  int64_t* positions,
  attr_t* attributes,
  const attr_t* attributes_mc,
  int32_t* coeffBufIt,
  attr::ModeEncoder& coder,
  GetMode getMode)
{
  // coefficients are stored in three planar arrays.  coeffBufItK is a set
  // of iterators to each array.
  int32_t* coeffBufItK[3] = {
    coeffBufIt,
    coeffBufIt + numPoints,
    coeffBufIt + numPoints * 2,
  };

  if (numPoints == 1) {
    auto quantizers = qpset.quantizers(0, pointQpOffsets[0]);
    for (int k = 0; k < numAttrs; k++) {
      auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

      auto coeff = attributes[k];
      assert(coeff <= INT_MAX && coeff >= INT_MIN);
      *coeffBufItK[k]++ = coeff =
        q.quantize(coeff << kFixedPointAttributeShift);
      attributes[k] =
        divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
    }
    return;
  }

  std::vector<UrahtNodeLight> weightsLf, weightsHf;
  std::vector<int64_t> attrsLf, attrsHf;

  bool enableACInterPred =
    rahtPredParams.enable_inter_prediction && attributes_mc;
  bool enableACRDOInterPred =
    rahtPredParams.raht_enable_inter_intra_layer_RDO
    && enableACInterPred && rahtPredParams.prediction_enabled_flag;

  coder.setInterEnabled(
    rahtPredParams.prediction_enabled_flag && enableACInterPred);

  const bool CbCrEnabled =
    !coder.isInterEnabled()
    && rahtPredParams.cross_chroma_component_prediction_flag;

  const bool CCRPEnabled =
    rahtPredParams.cross_component_residual_prediction_flag;

  const int maxlevelCCPenabled = (rahtPredParams.numlayer_CCRP_enabled - 1) * 3;

  int RDOCodingDepth = abh.attr_layer_code_mode.size();
  std::vector<int64_t> interTree;

  weightsLf.reserve(numPoints);
  attrsLf.reserve(numPoints * numAttrs);

  int regionQpShift = 4;

  // copy positions into internal form
  for (int i = 0; i < numPoints; i++) {
    UrahtNodeLight node;
    node.pos = positions[i];
    node.weight = 1;
    node.qp = {
      int16_t(pointQpOffsets[i][0] << regionQpShift),
      int16_t(pointQpOffsets[i][1] << regionQpShift)};
    weightsLf.emplace_back(node);
    for (int k = 0; k < numAttrs; k++) {
      attrsLf.push_back(attributes[i * numAttrs + k]);
    }
  }

  // ascend tree
  weightsHf.reserve(numPoints);
  attrsHf.reserve(numPoints * numAttrs);
  std::vector<int> levelHfPos;

  for (int level = 0, numNodes = weightsLf.size(); numNodes > 1; level++) {
    levelHfPos.push_back(weightsHf.size());
    if (level == 0) {
      // process any duplicate points
      numNodes = reduceUnique<haarFlag, numAttrs>(
        numNodes, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    } else {
      // normal level reduction
      numNodes = reduceLevel<haarFlag, numAttrs>(
        level, numNodes, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    }
  }

  assert(weightsLf[0].weight == numPoints);
  weightsLf[0].mode = Mode::Null;

  // reconstruction buffers
  std::vector<int32_t> attrRec, attrRecParent;
  std::vector<int32_t> interLayerAttrRec, intraLayerAttrRec;
  attrRec.resize(numPoints * numAttrs);
  attrRecParent.resize(numPoints * numAttrs);
  if (enableACRDOInterPred) {
    interLayerAttrRec.resize(numPoints * numAttrs);
    intraLayerAttrRec.resize(numPoints * numAttrs);
  }

  std::vector<int64_t> attrRecUs, attrRecParentUs, intraLayerAttrRecUs, interLayerAttrRecUs;
  attrRecUs.resize(numPoints * numAttrs);
  attrRecParentUs.resize(numPoints * numAttrs);
  if (enableACRDOInterPred) {
    intraLayerAttrRecUs.resize(numPoints* numAttrs);
    interLayerAttrRecUs.resize(numPoints* numAttrs);
  }

  std::vector<UrahtNode> weightsParent(1);
  weightsParent[0] = weightsLf[0];

  std::vector<int8_t> numParentNeigh, numGrandParentNeigh;
  numParentNeigh.resize(numPoints);
  numGrandParentNeigh.resize(numPoints);

  // indexes of the neighbouring parents
  int parentNeighIdx[19];
  int childNeighIdx[12][8];

  // Prediction buffers
  int64_t SampleDomainBuff[3 * numAttrs * 8];
  int64_t transformBuf[3 * numAttrs * 8];

  const int numBuffers = (2 + coder.isInterEnabled()) * numAttrs;
  constexpr size_t sizeofBuf = sizeof(int64_t) * 8;
  const size_t sizeofBufs = sizeofBuf * numBuffers;

  int64_t* attrPred;
  int64_t* attrPredTransform;
  int64_t* attrReal;
  int64_t* attrRealTransform;
  int64_t* attrPredIntra;
  int64_t* attrPredInter;
  int64_t* attrBestPredIt;
  int64_t* attrPredIntraTransformIt;
  int64_t* attrPredInterTransformIt;
  int64_t* attrBestPredTransformIt;

  std::vector<Mode> modes;

  int64_t SampleLayerBuf[2 * numAttrs * 8];
  int64_t transformLayerBuf[2 * numAttrs * 8];

  const int numLayerBuffers = 2 * coder.isInterEnabled() * numAttrs;
  const size_t sizeofLayerBuf = sizeofBuf * numLayerBuffers;

  int64_t attrInterPredLayer[1 * numAttrs * 8];
  int64_t attrInterPredLayerTransform[1 * numAttrs * 8];
  int64_t* attrPredIntraLayer = SampleLayerBuf + 8 * numAttrs;
  int64_t* attrPredInterLayer = SampleLayerBuf;
  int64_t* attrOrgPredInterLayer = attrInterPredLayer;
  int64_t* attrPredIntraLayerTransformIt = transformLayerBuf + 8 * numAttrs;
  int64_t* attrPredInterLayerTransformIt = transformLayerBuf;
  int64_t* attrOrgPredInterLayerTransformIt = attrInterPredLayerTransform;

  attrReal = SampleDomainBuff;
  attrRealTransform = transformBuf;
  attrPred = SampleDomainBuff + 8 * numAttrs;
  attrPredTransform = transformBuf + 8 * numAttrs;

  modes.push_back(Mode::Null);

  if (coder.isInterEnabled()) {
    attrPredInter = attrPred;
    attrPredInterTransformIt = attrPredTransform;
    modes.push_back(Mode::Inter);
    attrPred += 8 * numAttrs;
    attrPredTransform += 8 * numAttrs;
  }

  attrPredIntra = attrPred;
  attrPredIntraTransformIt = attrPredTransform;
  modes.push_back(Mode::Intra);

  // quant layer selection
  auto qpLayer = 0;

  // descend tree
  weightsLf.resize(1);
  attrsLf.resize(numAttrs);
  int trainZeros = 0;
  // NB: rootLayer = ceil((levelHfPos.size() - 1)/3.0)
  int rootLayer = div3level(levelHfPos.size() + 1);
  int CccpCoeff = 0;
  PCCRAHTComputeCCCP curlevelCccp;
  int intraLayerTrainZeros = 0;
  int interLayerTrainZeros = 0;
  PCCRAHTACCoefficientEntropyEstimate intraLayerEstimate;
  PCCRAHTACCoefficientEntropyEstimate interLayerEstimate;
  PCCRAHTACCoefficientEntropyEstimate curEstimate;
  int depth = 0;
  std::vector<int> intraLayerACCoeffcients, interLayerACCoeffcients;
  if (enableACRDOInterPred) {
    intraLayerACCoeffcients.resize(numPoints * numAttrs);
    interLayerACCoeffcients.resize(numPoints * numAttrs);
  }
  // number of coded coefficients (=nodes) in the layer
  // n.b. in first layer there is also the DC coefficient (=root node) to count
  int sumNodes = 1;

  int preLayerCodeMode = 0;

  for (int level = levelHfPos.size() - 1, isFirst = 1; level > 0; /*nop*/) {
    int numNodes = weightsHf.size() - levelHfPos[level];
    sumNodes += numNodes;
    weightsLf.resize(weightsLf.size() + numNodes);
    attrsLf.resize(attrsLf.size() + numNodes * numAttrs);
    expandLevel<haarFlag, numAttrs>(
      level, numNodes, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    weightsHf.resize(levelHfPos[level]);
    attrsHf.resize(levelHfPos[level] * numAttrs);

    // expansion of level is complete, processing is now on the next level
    level--;

    // every three levels, perform transform
    int layer;
    if (rem3level(level, layer))
      continue;

    CccpCoeff = 0;
    curlevelCccp.reset();

    int distanceToRoot = rootLayer - layer;
    int layerD = RDOCodingDepth - distanceToRoot;
    int predCtxLevel = 0;
    if (rahtPredParams.enable_inter_prediction) {
      predCtxLevel = layerD - rahtPredParams.mode_level;
      if (predCtxLevel >= NUMBER_OF_LEVELS_MODE)
        predCtxLevel = NUMBER_OF_LEVELS_MODE - 1;
    } else if (rahtPredParams.prediction_enabled_flag) {
      predCtxLevel = layerD - rahtPredParams.intra_mode_level;
      if (predCtxLevel >= NUMBER_OF_LEVELS_MODE)
        predCtxLevel = NUMBER_OF_LEVELS_MODE - 1;
    }

    bool upperInferMode = coder.isInterEnabled()
      && distanceToRoot < rahtPredParams.upper_mode_level
      && distanceToRoot < RDOCodingDepth - rahtPredParams.mode_level + 1;

    const bool enableAveragePredictionLevel =
      rahtPredParams.enable_average_prediction
      && distanceToRoot >= (rootLayer - rahtPredParams.mode_level
        - rahtPredParams.upper_mode_level_for_average_prediction + 1)
      && distanceToRoot < (rootLayer - rahtPredParams.mode_level
        + rahtPredParams.lower_mode_level_for_average_prediction + 1)
      && !upperInferMode && coder.isInterEnabled();

    bool realInferInLowerLevel =
      rahtPredParams.enable_average_prediction
      ? (distanceToRoot >= (RDOCodingDepth - rahtPredParams.mode_level
          + rahtPredParams.lower_mode_level_for_average_prediction + 1)
        && !upperInferMode && coder.isInterEnabled())
      : predCtxLevel < 0 && !upperInferMode && coder.isInterEnabled();

    // Motion compensation
    if (coder.isInterEnabled()) {
      translateLayer<haarFlag, numAttrs>(
        interTree, attributes_mc, weightsLf);
    }

    // initial scan position of the coefficient buffer
    //  -> first level = all coeffs
    //  -> otherwise = ac coeffs only
    bool inheritDc = !isFirst;
    bool enableIntraPredInLvl =
      inheritDc && rahtPredParams.prediction_enabled_flag;

    auto it = weightsLf.begin();
    for (auto i = 0; i < weightsParent.size(); i++) {
      weightsParent[i].occupancy = 0;
      weightsParent[i].decoded = 0;
      weightsParent[i].firstChild = it++;
      int numChildren = 1;

      if (preLayerCodeMode)
        weightsParent[i].mode = (pcc::attr::Mode)(3 - preLayerCodeMode); // preLayerCodeMode == 2 -> 1 INTRA, preLayerCodeMode == 1 -> 2 INTER

      while (it != weightsLf.end()
          && !((it->pos ^ weightsParent[i].pos) >> (level + 3))) {
        it++;
        numChildren++;
      }
      weightsParent[i].numChildren = numChildren;

      if (isFirst)
        weightsParent[i].mode = Mode::Null;
    }

    isFirst = 0;

    // select quantiser according to transform layer
    qpLayer = std::min(qpLayer + 1, int(qpset.layers.size()) - 1);

    // prepare reconstruction buffers
    //  previous reconstruction -> attrRecParent
    std::swap(attrRec, attrRecParent);
    std::swap(attrRecUs, attrRecParentUs);
    std::swap(numParentNeigh, numGrandParentNeigh);
    auto numGrandParentNeighIt = numGrandParentNeigh.cbegin();

    if (attributes_mc) {
      assert(interTree.size() <= attrRec.size());
      std::copy(interTree.begin(), interTree.end(), attrRec.begin());
    }

    int32_t* intraLayerCoeffBufItK[3] = {
      intraLayerACCoeffcients.data(),
      intraLayerACCoeffcients.data() + sumNodes,
      intraLayerACCoeffcients.data() + sumNodes * 2,
    };
    int32_t* intraLayerCoeffBufItBeginK[3] = {
      intraLayerCoeffBufItK[0],
      intraLayerCoeffBufItK[1],
      intraLayerCoeffBufItK[2],
    };
    int32_t* interLayerCoeffBufItK[3] = {
      interLayerACCoeffcients.data(),
      interLayerACCoeffcients.data() + sumNodes,
      interLayerACCoeffcients.data() + sumNodes * 2,
    };
    int32_t* interLayerCoeffBufItBeginK[3] = {
      interLayerCoeffBufItK[0],
      interLayerCoeffBufItK[1],
      interLayerCoeffBufItK[2],
    };

    int32_t* coeffBufItBeginK[3] = {
      coeffBufItK[0],
      coeffBufItK[1],
      coeffBufItK[2],
    };

    bool curLevelEnableLayerModeCoding = false;
    bool curLevelEnableACInterPred = false;
    enableACRDOInterPred =
      rahtPredParams.raht_enable_inter_intra_layer_RDO
      && enableACInterPred && rahtPredParams.prediction_enabled_flag
      && depth < RDOCodingDepth;

    if (enableACRDOInterPred) {
      curLevelEnableACInterPred = rahtPredParams.prediction_enabled_flag;
      curLevelEnableLayerModeCoding = curLevelEnableACInterPred;
    }

    if (enableACRDOInterPred) {
      coder.restoreStates();
      coder.resetModeBits();
    }

    double distinterLayer = 0;
    double distintraLayer = 0;
    double distMCPredLayer = 0;
    double dlambda = 1.0;

    //CCRP Parameters
    const bool CCRPFlag =
      !haarFlag && CCRPEnabled && level <= maxlevelCCPenabled && !CbCrEnabled;

    CCRPFilter ccrpFilter;
    CCRPFilter ccrpFilterIntra;
    CCRPFilter ccrpFilterInter;

    int i = 0;
    for (auto weightsParentIt = weightsParent.begin();
        weightsParentIt < weightsParent.end();
        weightsParentIt++) {

      memset(SampleDomainBuff, 0, sizeofBufs);
      memset(transformBuf, 0, sizeofBufs);

      if (enableACRDOInterPred) {
        memset(SampleLayerBuf, 0, sizeofLayerBuf);
        memset(transformLayerBuf, 0, sizeofLayerBuf);
      }

      int64_t transformIntraLayerBuf[8 * numAttrs] = {};
      int64_t transformInterLayerBuf[8 * numAttrs] = {};

      int64_t weights[8 + 8 + 8 + 8 + 24] = {};
      int64_t sqrtweightsbuf[8] = {};
      int64_t normsqrtsbuf[8] = {};
      bool skipinverse = !haarFlag;
      bool skipinverseinterlayer = skipinverse;
      bool skipinverseintralayer = skipinverse;


      int64_t PredDC[3] = { 0,0,0 };
      int64_t PredDCIntraLayer[3] = { 0,0,0 };
      int64_t PredDCInterLayer[3] = { 0,0,0 };

      Qps nodeQp[8] = {};
      uint8_t occupancy = 0;

      // generate weights, occupancy mask, and fwd transform buffers
      // for all siblings of the current node.
      int nodeCnt = weightsParentIt->numChildren;
      for (int t = 0, j0 = i; t < nodeCnt; t++, j0++) {

        int nodeIdx = (weightsLf[j0].pos >> level) & 0x7;
        weights[nodeIdx] = weightsLf[j0].weight;
        nodeQp[nodeIdx][0] = weightsLf[j0].qp[0] >> regionQpShift;
        nodeQp[nodeIdx][1] = weightsLf[j0].qp[1] >> regionQpShift;

        occupancy |= 1 << nodeIdx;

        //store grandmode
        if (depth >= 1) {
          weightsLf[j0].grand_mode = weightsParentIt->mode;
        } else {
          weightsLf[j0].grand_mode = Mode::Null;
        }

        for (int k = 0; k < numAttrs; k++)
          attrReal[8 * k + nodeIdx] = fpExpand<kFPFracBits>(attrsLf[j0 * numAttrs + k]);
      }
      int64_t sumweights = weightsParentIt->weight;
      mkWeightTree<false, RahtKernel>(weights);

      if (!inheritDc) {
        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])
            continue;
          numParentNeigh[j++] = 19;
        }
      }

      auto attrRecParentUsIt = std::next(
        attrRecParentUs.cbegin(),
        std::distance(weightsParent.begin(), weightsParentIt) * numAttrs);

      // Inter-level prediction:
      //  - Find the parent neighbours of the current node
      //  - Generate prediction for all attributes into transformIntraBuf
      //  - Subtract transformed coefficients from forward transform
      //  - The transformIntraBuf is then used for reconstruction
      bool enableIntraPred = false;
      bool enableInterPred = false;
      //< layer mode coding
      bool enableIntraLayerPred = false;
      bool enableInterLayerPred = false;
      if (enableACRDOInterPred) {
        enableIntraPred =
          rahtPredParams.enable_inter_prediction
          ? enableIntraPredInLvl && (nodeCnt > 1) && (distanceToRoot > 2)
          : enableIntraPredInLvl;
        enableInterPred = coder.isInterEnabled() && (nodeCnt > 1);

        enableIntraLayerPred = enableIntraPredInLvl && (nodeCnt > 1);
        enableInterLayerPred = curLevelEnableACInterPred && (nodeCnt > 1);
      } else {
        enableIntraPred =
          rahtPredParams.enable_inter_prediction
          ? enableIntraPredInLvl && (nodeCnt > 1) && (distanceToRoot > 2)
          : enableIntraPredInLvl;
        enableInterPred = coder.isInterEnabled() && (nodeCnt > 1);
      }

      // inter prediction
      weightsParentIt->occupancy = occupancy;
      Mode neighborsMode = Mode::size;
      if (enableInterPred || enableInterLayerPred) {
        enableInterPred = true;

        auto pred = attrPredInter;
        auto inter = std::next(
          interTree.begin(),
          std::distance(weightsLf.begin(), weightsParentIt->firstChild)
            * numAttrs);

        for (int nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if ((occupancy >> nodeIdx) & 0x1)
            for (int k = 0; k < numAttrs; k++)
              pred[8 * k + nodeIdx] = *(inter++);
        }
      }

      int voteInterWeight = 1, voteIntraWeight = 1;
      int voteInterLayerWeight = 1, voteIntraLayerWeight = 1;

      if (enableIntraPred || enableIntraLayerPred) {
        int parentNeighCount = 0;
        if (rahtPredParams.enable_inter_prediction
            || (!(rahtPredParams.prediction_skip1_flag && nodeCnt == 1)
              && !(*numGrandParentNeighIt < rahtPredParams.prediction_threshold0))) {

          parentNeighCount = findNeighbours(
            weightsParent.begin(), weightsParent.end(), weightsParentIt,
            level + 3, occupancy, parentNeighIdx);

          if (rahtPredParams.subnode_prediction_enabled_flag)
            findNeighboursChildren(
              weightsParent.begin(), weightsLf.begin(), level, occupancy,
              parentNeighIdx, childNeighIdx);

          if (rahtPredParams.prediction_enabled_flag)
            neighborsMode = attr::getNeighborsMode(
              true, parentNeighIdx, weightsParent, voteInterWeight,
              voteIntraWeight, voteInterLayerWeight, voteIntraLayerWeight);
        }
        if (!rahtPredParams.enable_inter_prediction
            && rahtPredParams.prediction_skip1_flag && nodeCnt == 1) {
          enableIntraPred = false;
          enableIntraLayerPred = false;
          parentNeighCount = 19;
        } else if (!rahtPredParams.enable_inter_prediction
            && *numGrandParentNeighIt < rahtPredParams.prediction_threshold0) {
          enableIntraPred = false;
          enableIntraLayerPred = false;
        } else if (parentNeighCount < rahtPredParams.prediction_threshold1) {
            enableIntraPred = false;
            enableIntraLayerPred = false;
        } else {
          intraDcPred<haarFlag, numAttrs>(
            true, parentNeighIdx, childNeighIdx, occupancy,
            attrRecParent.begin(), attrRec.begin(), intraLayerAttrRec.begin(),
            interLayerAttrRec.begin(), attrPredIntra, attrPredIntraLayer,
            attrPredInterLayer, rahtPredParams, enableACRDOInterPred);
        }
        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])
            continue;
          numParentNeigh[j++] = parentNeighCount;
        }
      }

      if (inheritDc) {
        numGrandParentNeighIt++;
        weightsParentIt->decoded = true;
      }

      if (!enableIntraPred) {
        neighborsMode = enableInterPred ? Mode::Inter : Mode::Null;
      }

      if (enableACRDOInterPred) {
        bool interpredlayeravailable =
          enableInterLayerPred && enableInterPred
          || enableIntraLayerPred;
        bool intrapredlayeravailable = enableIntraLayerPred;
        skipinverseinterlayer = skipinverseinterlayer && interpredlayeravailable;
        skipinverseintralayer = skipinverseintralayer && intrapredlayeravailable;
      }

      if (haarFlag) {
        memcpy(transformBuf, SampleDomainBuff, sizeofBufs);
        fwdTransformBlock222<HaarKernel>(numBuffers, transformBuf, weights);
        if (enableACRDOInterPred) {
          memcpy(transformLayerBuf, SampleLayerBuf, sizeofLayerBuf);
          fwdTransformBlock222<HaarKernel>(
            numLayerBuffers, transformLayerBuf, weights);
        }
      } else {
        // normalise coefficients
        for (int childIdx = 0; childIdx < 8; childIdx++) {
          if (weights[childIdx] <= 1) {
            sqrtweightsbuf[childIdx] = 32768;
            continue;
          }

          // Summed attribute values
          uint64_t w = weights[childIdx];
          int shift = 5 * ((w > 1024) + (w > 1048576));
          int64_t rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;
          for (int k = 0; k < numAttrs; k++) {
            SampleDomainBuff[8 * k + childIdx] = fpReduce<kFPFracBits>(
              (SampleDomainBuff[8 * k + childIdx] >> shift) * rsqrtWeight);
          }

          // Predicted attribute values
          int64_t sqrtWeight = fastIsqrt(weights[childIdx]);
          sqrtweightsbuf[childIdx] = sqrtWeight;

          for (auto buf = SampleDomainBuff + 8 * numAttrs;
              buf < SampleDomainBuff + 8 * numBuffers;
              buf += 8 * numAttrs) {
            for (int k = 0; k < numAttrs; k++)
              buf[8 * k + childIdx] = fpReduce<kFPFracBits>(
                buf[8 * k + childIdx] * sqrtWeight);
          }
          if (enableACRDOInterPred) {
            for (auto buf = SampleLayerBuf;
                buf < SampleLayerBuf + 8 * numLayerBuffers;
                buf += 8 * numAttrs) {
              for (int k = 0; k < numAttrs; k++)
                buf[8 * k + childIdx] = fpReduce<kFPFracBits>(
                  buf[8 * k + childIdx] * sqrtWeight);
            }
          }
        }
        memcpy(transformBuf, SampleDomainBuff, sizeofBufs);
        fwdTransformBlock222<RahtKernel>(numBuffers, transformBuf, weights);
        if (enableACRDOInterPred) {
          memcpy(transformLayerBuf, SampleLayerBuf, sizeofLayerBuf);
          fwdTransformBlock222<RahtKernel>(
            numLayerBuffers, transformLayerBuf, weights);
        }
      } //else normalize

      if (enableACRDOInterPred) {
        std::copy_n(transformBuf, 8 * numAttrs, transformIntraLayerBuf);
        std::copy_n(transformBuf, 8 * numAttrs, transformInterLayerBuf);
      }

      const bool enableAveragePrediction =
        enableAveragePredictionLevel
        && (enableIntraPred || enableIntraLayerPred)
        && enableInterPred;
      const bool enableAverageLayerPrediction =
        curLevelEnableLayerModeCoding && enableInterLayerPred
        && enableIntraLayerPred && enableInterPred;
      int64_t weightIntra, weightInter;
      int64_t weightIntraLayer, weightInterLayer;

      if (enableAveragePrediction) {
        bool flag = !isInter(weightsParentIt->mode) && !isIntra(weightsParentIt->mode);
        voteInterWeight += 12 * isInter(weightsParentIt->mode) + 6 * flag;
        voteIntraWeight += 12 * isIntra(weightsParentIt->mode) + 6 * flag;

        if (depth > 1) {
          bool flag_grand = !isInter(weightsParentIt->grand_mode)
            && !isIntra(weightsParentIt->grand_mode);
          voteInterWeight +=
            28 * isInter(weightsParentIt->grand_mode) + 14 * flag_grand;
          voteIntraWeight +=
            28 * isIntra(weightsParentIt->grand_mode) + 14 * flag_grand;
        }

        weightIntra = divApprox(voteIntraWeight << kFPFracBits, voteInterWeight + voteIntraWeight, 0);
        weightInter = (1 << kFPFracBits) - weightIntra;



        if (enableAverageLayerPrediction) {
          std::copy_n(attrPredInter, 8 * numAttrs, attrOrgPredInterLayer);
          std::copy_n(attrPredInterTransformIt, 8 * numAttrs, attrOrgPredInterLayerTransformIt);

          voteInterLayerWeight += 12 * isInter(weightsParentIt->mode) + 6 * flag;
          voteIntraLayerWeight += 12 * isIntra(weightsParentIt->mode) + 6 * flag;

          if (depth > 1) {
            bool flag_grand = !isInter(weightsParentIt->grand_mode)
              && !isIntra(weightsParentIt->grand_mode);
            voteInterLayerWeight +=
              28 * isInter(weightsParentIt->grand_mode) + 14 * flag_grand;
            voteIntraLayerWeight +=
              28 * isIntra(weightsParentIt->grand_mode) + 14 * flag_grand;
          }

          weightIntraLayer = divApprox(
            voteIntraLayerWeight << kFPFracBits,
            voteInterLayerWeight + voteIntraLayerWeight, 0);
          weightInterLayer = (1 << kFPFracBits) - weightIntraLayer;
        }

        for (int nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          for (int k = 0; k < numAttrs; k++) {
            if (enableIntraPred && enableInterPred) {
              if (predCtxLevel < 0) {
                attrPredIntra[8 * k + nodeIdx] = fpReduce<kFPFracBits>(
                  attrPredIntra[8 * k + nodeIdx] * weightIntra
                  + attrPredInter[8 * k + nodeIdx] * weightInter);
                attrPredIntraTransformIt[8 * k + nodeIdx] = fpReduce<kFPFracBits>(
                  attrPredIntraTransformIt[8 * k + nodeIdx] * weightIntra
                  + attrPredInterTransformIt[8 * k + nodeIdx] * weightInter);
                if (haarFlag) {
                  attrPredIntra[8 * k + nodeIdx] &= kFPIntMask;
                  attrPredIntraTransformIt[8 * k + nodeIdx] &= kFPIntMask;
                }
              }
              else {
                attrPredInter[8 * k + nodeIdx] = fpReduce<kFPFracBits>(
                  attrPredIntra[8 * k + nodeIdx] * weightIntra
                  + attrPredInter[8 * k + nodeIdx] * weightInter);
                attrPredInterTransformIt[8 * k + nodeIdx] = fpReduce<kFPFracBits>(
                  attrPredIntraTransformIt[8 * k + nodeIdx] * weightIntra
                  + attrPredInterTransformIt[8 * k + nodeIdx] * weightInter);
                if (haarFlag) {
                  attrPredInter[8 * k + nodeIdx] &= kFPIntMask;
                  attrPredInterTransformIt[8 * k + nodeIdx] &= kFPIntMask;
                }
              }
            }
            if (enableAverageLayerPrediction) {
              attrPredInterLayer[8 * k + nodeIdx] = fpReduce<kFPFracBits>(
                attrPredInterLayer[8 * k + nodeIdx] * weightIntraLayer
                + attrOrgPredInterLayer[8 * k + nodeIdx] * weightInterLayer);
              attrPredInterLayerTransformIt[8 * k + nodeIdx] = fpReduce<kFPFracBits>(
                attrPredInterLayerTransformIt[8 * k + nodeIdx] * weightIntraLayer
                + attrOrgPredInterLayerTransformIt[8 * k + nodeIdx] * weightInterLayer);
              if (haarFlag) {
                attrPredInterLayer[8 * k + nodeIdx] &= kFPIntMask;
                attrPredInterLayerTransformIt[8 * k + nodeIdx] &= kFPIntMask;
              }
            }
          }
        }
      }
      else {  //!enableAveragePrediction
        if (enableInterLayerPred && enableInterPred) {
          std::copy_n(attrPredInter, 8 * numAttrs, attrPredInterLayer);
          std::copy_n(attrPredInterTransformIt, 8 * numAttrs, attrPredInterLayerTransformIt);
        }
      }

      Mode predMode = Mode::Null;
      predMode =
        rahtPredParams.prediction_enabled_flag
        ? getMode(
          coder, nodeCnt, predCtxLevel, enableIntraPred,
          enableInterPred, weightsParentIt->mode, neighborsMode, numAttrs,
          weights, attrRecParentUsIt, transformBuf, modes, qpLayer, nodeQp,
          upperInferMode, inheritDc)
        : Mode::Null;

      auto weightsChild = weightsParentIt->firstChild;
      for (int t = 0; t < weightsParentIt->numChildren; t++, weightsChild++) {
        if (int(predMode) >= Mode::Inter)
          weightsChild->mode = Mode::Inter;
        else
          weightsChild->mode = predMode;
        if (curLevelEnableLayerModeCoding)
        {
          if (enableInterLayerPred && enableInterPred)
            weightsChild->_mode = Mode::Inter;
          else if (enableIntraLayerPred)
            weightsChild->_mode = Mode::Intra;
          else
            weightsChild->_mode = Mode::Null;
        }
      }

      if (attr::isNull(predMode)) {
        skipinverse = false;
        attrBestPredIt = SampleDomainBuff + 8 * numAttrs;
        attrBestPredTransformIt = transformBuf + 8 * numAttrs;
        std::fill_n(attrBestPredIt, 8 * numAttrs, 0);
        std::fill_n(attrBestPredTransformIt, 8 * numAttrs, 0);
        for (int k = 0; k < numAttrs; k++)
          PredDC[k] = 0;
      } else if (attr::isIntra(predMode)) {
        attrBestPredIt = attrPredIntra;
        attrBestPredTransformIt = attrPredIntraTransformIt;
      } else {
        attrBestPredIt = attrPredInter;
        attrBestPredTransformIt = attrPredInterTransformIt;
      }
      //compute DC of prediction signals
      if (!attr::isNull(predMode) || enableACRDOInterPred) {
        if (!haarFlag) {
          int64_t rsqrtweightsum = fastIrsqrt(sumweights);

          for (int childIdx = 0; childIdx < 8; childIdx++) {
            if (weights[childIdx] == 0)
              continue;

            int64_t normSqrtW = sqrtweightsbuf[childIdx] * rsqrtweightsum >> 40;
            normsqrtsbuf[childIdx] = normSqrtW;
            for (int k = 0; k < numAttrs; k++){
              if (!attr::isNull(predMode)) {
                PredDC[k] += fpReduce<kFPFracBits>(
                  normSqrtW * attrBestPredIt[8 * k + childIdx]);
              }
              if (enableACRDOInterPred) {
                PredDCIntraLayer[k] += fpReduce<kFPFracBits>(
                  normSqrtW * attrPredIntraLayer[8 * k + childIdx]);
                PredDCInterLayer[k] += fpReduce<kFPFracBits>(
                  normSqrtW * attrPredInterLayer[8 * k + childIdx]);
              }
            }
          }
        }
      }

      // per-coefficient operations:
      //  - subtract transform domain prediction (encoder)
      //  - write out/read in quantised coefficients
      //  - inverse quantise + add transform domain prediction
      int64_t BestRecBuf[8 * numAttrs] = { 0 };
      int64_t BestRecBufIntraLayer[8 * numAttrs] = { 0 };
      int64_t BestRecBufInterLayer[8 * numAttrs] = { 0 };

      int64_t CoeffRecBuf[8][numAttrs] = { 0 };
      int nodelvlSum = 0;
      int64_t transformRecBuf[numAttrs] = { 0 };

      CCRPFilter::Corr curCorr = { 0, 0, 0 };
      CCRPFilter::Corr curCorrIntra  = { 0, 0, 0 };
      CCRPFilter::Corr curCorrInter  = { 0, 0, 0 };

      scanBlock(weights, [&](int idx) {
        // skip the DC coefficient unless at the root of the tree
        if (inheritDc && !idx)
          return;

        // subtract transformed prediction (skipping DC)
        if (!attr::isNull(predMode)) {
          for (int k = 0; k < numAttrs; k++) {
            transformBuf[8 * k + idx] -= attrBestPredTransformIt[8 * k + idx];
          }
        }

        if (enableACRDOInterPred) {
          for (int k = 0; k < numAttrs; k++) {
            if (!realInferInLowerLevel)
              transformIntraLayerBuf[8 * k + idx] -= attrPredIntraLayerTransformIt[8 * k + idx];

            if (!upperInferMode)
              transformInterLayerBuf[8 * k + idx] -= attrPredInterLayerTransformIt[8 * k + idx];
          }
        }

        // decision for RDOQ
        int64_t sumCoeff = 0;
        const int LUTlog[16] = {0,   256, 406, 512, 594, 662, 719,  768,
                                812, 850, 886, 918, 947, 975, 1000, 1024};
        bool flagRDOQ = false;
        int64_t Qcoeff;
        int64_t intraLayerSumCoeff = 0;
        bool intraLayerFlagRDOQ = false;

        int64_t interLayerSumCoeff = 0;
        bool interLayerFlagRDOQ = false;

        if (!haarFlag) {
          int64_t recDist2 = 0;
          int64_t Dist2 = 0;
          int Ratecoeff = 0;
          int64_t lambda0;
          int64_t lambda;
          int64_t coeff = 0;

          int64_t intraLayerRecDist2 = 0;
          int64_t intraLayerDist2 = 0;
          int intraLayerRatecoeff = 0;

          int64_t interLayerRecDist2 = 0;
          int64_t interLayerDist2 = 0;
          int interLayerRatecoeff = 0;

          for (int k = 0; k < numAttrs; k++) {
            //auto q = Quantizer(qpLayer[std::min(k, int(quantizers.size()) - 1)] + nodeQp[idx]);
            auto quantizers = qpset.quantizers(qpLayer, nodeQp[idx]);
            auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];
            coeff = fpReduce<kFPFracBits>(transformBuf[8 * k + idx]);
            Dist2 += coeff * coeff;

            if (!haarFlag && rahtPredParams.cross_chroma_component_prediction_flag && !coder.isInterEnabled()) {
              if (k != 2) {
                Qcoeff = q.quantize(coeff << kFixedPointAttributeShift);
                transformRecBuf[k] = fpExpand<kFPFracBits>(
                  divExp2RoundHalfUp(q.scale(Qcoeff), kFixedPointAttributeShift));
              } else {
                transformRecBuf[k] = transformBuf[8 * k + idx] - ((CccpCoeff * transformRecBuf[1]) >> 4);
                coeff = fpReduce<kFPFracBits>(transformRecBuf[k]);
                Qcoeff = q.quantize((coeff) << kFixedPointAttributeShift);
              }
            } else
              Qcoeff = q.quantize(coeff << kFixedPointAttributeShift);

            auto recCoeff =
              divExp2RoundHalfUp(q.scale(Qcoeff), kFixedPointAttributeShift);
            recDist2 += (coeff - recCoeff) * (coeff - recCoeff);

            sumCoeff += std::abs(Qcoeff);
            //Ratecoeff += !!Qcoeff; // sign
            Ratecoeff +=
              std::abs(Qcoeff) < 15 ? LUTlog[std::abs(Qcoeff)] : LUTlog[15];
            if (!k)
              lambda0 = q.scale(1);
            lambda = lambda0 * lambda0 * (numAttrs == 1 ? 25 : 35);
            if (enableACRDOInterPred) {
              if (!upperInferMode) {
                auto interLayerCoeff = fpReduce<kFPFracBits>(
                  transformInterLayerBuf[8 * k + idx]);
                interLayerDist2 += interLayerCoeff * interLayerCoeff;
                auto interLayerQcoeff =
                  q.quantize(interLayerCoeff << kFixedPointAttributeShift);

                auto recInterCoeff =
                  divExp2RoundHalfUp(q.scale(interLayerQcoeff), kFixedPointAttributeShift);
                interLayerRecDist2 +=
                  (interLayerCoeff - recInterCoeff) * (interLayerCoeff - recInterCoeff);

                interLayerSumCoeff += std::abs(interLayerQcoeff);
                interLayerRatecoeff +=
                  std::abs(interLayerQcoeff) < 15
                  ? LUTlog[std::abs(interLayerQcoeff)] : LUTlog[15];
              }

              if (!realInferInLowerLevel) {
                auto intraLayerCoeff = fpReduce<kFPFracBits>(
                  transformIntraLayerBuf[8 * k + idx]);
                intraLayerDist2 += intraLayerCoeff * intraLayerCoeff;
                auto intraLayerQcoeff =
                  q.quantize(intraLayerCoeff << kFixedPointAttributeShift);

                auto recIntraCoeff =
                  divExp2RoundHalfUp(q.scale(intraLayerQcoeff), kFixedPointAttributeShift);
                intraLayerRecDist2 +=
                  (intraLayerCoeff - recIntraCoeff) * (intraLayerCoeff - recIntraCoeff);

                intraLayerSumCoeff += std::abs(intraLayerQcoeff);
                intraLayerRatecoeff +=
                  std::abs(intraLayerQcoeff) < 15
                  ? LUTlog[std::abs(intraLayerQcoeff)] : LUTlog[15];
              }
            }
          }
          dlambda = (double)lambda;
          if (sumCoeff < 3) {
            int Rate = getRate(trainZeros);
            Rate += (Ratecoeff + 128) >> 8;
            flagRDOQ = (Dist2 << 26) < (lambda * Rate + (recDist2 << 26));
          }

          if (enableACRDOInterPred) {
            if (!realInferInLowerLevel && intraLayerSumCoeff < 3) {
              int Rate = getRate(intraLayerTrainZeros);
              Rate += (intraLayerRatecoeff + 128) >> 8;
              intraLayerFlagRDOQ =
                (intraLayerDist2 << 26) < (lambda * Rate + (intraLayerRecDist2 << 26));
            }

            if (!upperInferMode && interLayerSumCoeff < 3) {
              int Rate = getRate(interLayerTrainZeros);
              Rate += (interLayerRatecoeff + 128) >> 8;
              interLayerFlagRDOQ =
                (interLayerDist2 << 26) < (lambda * Rate + (interLayerRecDist2 << 26));
            }
          }
        }

        // Track RL for RDOQ
        if (flagRDOQ || sumCoeff == 0)
          trainZeros++;
        else
          trainZeros = 0;

        if (enableACRDOInterPred) {
          if (!realInferInLowerLevel) {
            if (intraLayerFlagRDOQ || intraLayerSumCoeff == 0)
              intraLayerTrainZeros++;
            else
              intraLayerTrainZeros = 0;
          }

          if (!upperInferMode) {
            if (interLayerFlagRDOQ || interLayerSumCoeff == 0)
              interLayerTrainZeros++;
            else
              interLayerTrainZeros = 0;
          }
        }

        // The RAHT transform
        auto quantizers = qpset.quantizers(qpLayer, nodeQp[idx]);

        int64_t quantizedLuma = 0, quantizedChroma = 0;
        int64_t quantizedLumaIntra = 0, quantizedChromaIntra = 0;
        int64_t quantizedLumaInter = 0, quantizedChromaInter = 0;

        for (int k = 0; k < numAttrs; k++) {
          if (flagRDOQ) { // apply RDOQ
            transformBuf[8 * k + idx] = 0;
            transformRecBuf[k] = 0;
          }

          if (enableACRDOInterPred) {
            if (!realInferInLowerLevel) {
              if (intraLayerFlagRDOQ)
                transformIntraLayerBuf[8 * k + idx] = 0;
            }

            if (!upperInferMode) {
              if (interLayerFlagRDOQ)
                transformInterLayerBuf[8 * k + idx] = 0;
            }
          }

          auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

          int64_t CCRPPred = 0, CCRPPredIntra = 0, CCRPPredInter = 0;

          if(k && CCRPFlag){
            int64_t CCRPFilt =
              (k == 1) ? ccrpFilter.getYCbFilt() : ccrpFilter.getYCrFilt();

            CCRPPred = quantizedLuma * CCRPFilt >> kCCRPFiltPrecisionbits;

            transformBuf[8 * k + idx] =
              transformBuf[8 * k + idx] - fpExpand<kFPFracBits>(CCRPPred);

            if(enableACRDOInterPred){
              int64_t CCRPFiltIntra =
                (k == 1)
                ? ccrpFilterIntra.getYCbFilt()
                : ccrpFilterIntra.getYCrFilt();

              CCRPPredIntra =
                quantizedLumaIntra * CCRPFiltIntra >> kCCRPFiltPrecisionbits;

              transformIntraLayerBuf[8 * k + idx] =
                transformIntraLayerBuf[8 * k + idx]
                - fpExpand<kFPFracBits>(CCRPPredIntra);

              int64_t CCRPFiltInter =
                (k == 1)
                ? ccrpFilterInter.getYCbFilt()
                : ccrpFilterInter.getYCrFilt();

              CCRPPredInter =
                quantizedLumaInter * CCRPFiltInter >> kCCRPFiltPrecisionbits;

              transformInterLayerBuf[8 * k + idx] =
                transformInterLayerBuf[8 * k + idx]
                - fpExpand<kFPFracBits>(CCRPPredInter);
            }
          }

          auto coeff = fpReduce<kFPFracBits>(transformBuf[8 * k + idx]);
          assert(coeff <= INT_MAX && coeff >= INT_MIN);
          coeff =
            q.quantize(coeff << kFixedPointAttributeShift);

          if (!haarFlag && rahtPredParams.cross_chroma_component_prediction_flag && !coder.isInterEnabled()) {
            if (k != 2) {
              BestRecBuf[8 * k + idx] = transformRecBuf[k];
            } else {
              coeff = fpReduce<kFPFracBits>(transformRecBuf[k]);
              coeff = q.quantize(coeff << kFixedPointAttributeShift);
              transformRecBuf[k] = fpExpand<kFPFracBits>(
                divExp2RoundHalfUp(
                  q.scale(coeff), kFixedPointAttributeShift));
              transformRecBuf[k]+= CccpCoeff * transformRecBuf[1] >> 4;
              BestRecBuf[8 * k + idx] = transformRecBuf[k];
            }
            CoeffRecBuf[nodelvlSum][k] = fpReduce<kFPFracBits>(
              transformRecBuf[k]);
          } else {
            int64_t quantizedValue =
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
            BestRecBuf[8 * k + idx] = fpExpand<kFPFracBits>(quantizedValue);

            if(CCRPFlag){
              if (k == 0) {
                quantizedLuma = quantizedValue;
                curCorr.yy += quantizedLuma * quantizedLuma;
              } else {
                quantizedChroma = quantizedValue;
                quantizedChroma += CCRPPred;
                BestRecBuf[8 * k + idx] = fpExpand<kFPFracBits>(quantizedChroma);
                skipinverse = skipinverse && (quantizedChroma == 0);
                if (k == 1)
                  curCorr.ycb += quantizedLuma * quantizedChroma;
                else
                  curCorr.ycr += quantizedLuma * quantizedChroma;
              }
            }
          }

          if (enableACRDOInterPred)
            curEstimate.updateCostBits(coeff, k);

          *coeffBufItK[k]++ = coeff;
          skipinverse = skipinverse && (coeff == 0);

          if (enableACRDOInterPred)
            curEstimate.resStatUpdate(coeff, k);

          if (enableACRDOInterPred) {
            if (!realInferInLowerLevel) {
              auto intraLayerCoeff = fpReduce<kFPFracBits>(
                transformIntraLayerBuf[8 * k + idx]);
              assert(intraLayerCoeff <= INT_MAX && intraLayerCoeff >= INT_MIN);
              intraLayerCoeff =
                q.quantize(intraLayerCoeff << kFixedPointAttributeShift);

              intraLayerEstimate.updateCostBits(intraLayerCoeff, k);

              *intraLayerCoeffBufItK[k]++ = intraLayerCoeff;
              skipinverseintralayer = skipinverseintralayer && (intraLayerCoeff == 0);
              int64_t quantizedValueIntra =
                divExp2RoundHalfUp(q.scale(intraLayerCoeff), kFixedPointAttributeShift);
              BestRecBufIntraLayer[8 * k + idx] = fpExpand<kFPFracBits>(quantizedValueIntra);
              intraLayerEstimate.resStatUpdate(intraLayerCoeff, k);
              if(CCRPFlag){
                if (k == 0) {
                  quantizedLumaIntra = quantizedValueIntra;
                  curCorrIntra.yy += quantizedLumaIntra * quantizedLumaIntra;
                } else {
                  quantizedChromaIntra = quantizedValueIntra;
                  quantizedChromaIntra += CCRPPredIntra;
                  BestRecBufIntraLayer[8 * k + idx] = fpExpand<kFPFracBits>(quantizedChromaIntra);
                  skipinverseintralayer = skipinverseintralayer && (quantizedChromaIntra == 0);
                  if (k == 1)
                    curCorrIntra.ycb += quantizedLumaIntra * quantizedChromaIntra;
                  else
                    curCorrIntra.ycr += quantizedLumaIntra * quantizedChromaIntra;
                }
              }
            }

            if (!upperInferMode) {
              auto interLayerCoeff = fpReduce<kFPFracBits>(
                transformInterLayerBuf[8 * k + idx]);
              assert(interLayerCoeff <= INT_MAX && interLayerCoeff >= INT_MIN);
              interLayerCoeff =
                q.quantize(interLayerCoeff << kFixedPointAttributeShift);

              interLayerEstimate.updateCostBits(interLayerCoeff, k);

              skipinverseinterlayer = skipinverseinterlayer && (interLayerCoeff == 0);
              *interLayerCoeffBufItK[k]++ = interLayerCoeff;
              // still in RAHT domain
              int64_t quantizedValueInter =
                divExp2RoundHalfUp(q.scale(interLayerCoeff), kFixedPointAttributeShift);
              BestRecBufInterLayer[8 * k + idx] = fpExpand<kFPFracBits>(quantizedValueInter);

              interLayerEstimate.resStatUpdate(interLayerCoeff, k);

              if(CCRPFlag){
                if (k == 0) {
                  quantizedLumaInter = quantizedValueInter;
                  curCorrInter.yy += quantizedLumaInter * quantizedLumaInter;
                } else {
                  quantizedChromaInter = quantizedValueInter;
                  quantizedChromaInter += CCRPPredInter;
                  BestRecBufInterLayer[8 * k + idx] = fpExpand<kFPFracBits>(quantizedChromaInter);
                  skipinverseinterlayer = skipinverseinterlayer && (quantizedChromaInter == 0);
                  if (k == 1)
                    curCorrInter.ycb += quantizedLumaInter * quantizedChromaInter;
                  else
                    curCorrInter.ycr += quantizedLumaInter * quantizedChromaInter;
                }
              }
            }

            int64_t iresidueinterLayer = 0;
            if (!upperInferMode)
              // quantization reconstruction error
              iresidueinterLayer = fpReduce<kFPFracBits>(
                transformInterLayerBuf[8 * k + idx]
                - BestRecBufInterLayer[8 * k + idx]);

            int64_t iresidueintraLayer = 0;
            if (!realInferInLowerLevel)
              // quantization reconstruction error
              iresidueintraLayer = fpReduce<kFPFracBits>(
                transformIntraLayerBuf[8 * k + idx]
                - BestRecBufIntraLayer[8 * k + idx]);

            int64_t iresidueMCPred = fpReduce<kFPFracBits>(
              transformBuf[8 * k + idx] - BestRecBuf[8 * k + idx]);

            int64_t idistinterLayer;
            if (!upperInferMode) {
              idistinterLayer = iresidueinterLayer * iresidueinterLayer;
            }

            int64_t idistintraLayer;
            if (!realInferInLowerLevel) {
              idistintraLayer = iresidueintraLayer * iresidueintraLayer;
            }

            int64_t idistMCPredLayer = iresidueMCPred * iresidueMCPred;

            if (!upperInferMode)
              distinterLayer += (double)idistinterLayer;

            if (!realInferInLowerLevel)
              distintraLayer += (double)idistintraLayer;

            distMCPredLayer += (double)idistMCPredLayer;
          }

          if (haarFlag) {
            // Transform Domain Pred for Lossless case
            BestRecBuf[8 * k + idx] += attrBestPredTransformIt[8 * k + idx];
            if (enableACRDOInterPred) {
              if (!realInferInLowerLevel)
                BestRecBufIntraLayer[8 * k + idx] += attrPredIntraLayerTransformIt[8 * k + idx];

              if (!upperInferMode)
                BestRecBufInterLayer[8 * k + idx] += attrPredInterLayerTransformIt[8 * k + idx];
            }
          }
        }
        nodelvlSum++;
      });

      if (CCRPFlag) {
        ccrpFilter.update(curCorr);
        if (enableACRDOInterPred) {
          ccrpFilterIntra.update(curCorrIntra);
          ccrpFilterInter.update(curCorrInter);
        }
      }

       // compute last component coefficient
      if (numAttrs == 3 && nodeCnt > 1 && !haarFlag
          && rahtPredParams.cross_chroma_component_prediction_flag
          && !coder.isInterEnabled()) {
        CccpCoeff = curlevelCccp.computeCrossChromaComponentPredictionCoeff(nodelvlSum, CoeffRecBuf);
      }
      // replace DC coefficient with parent if inheritable
      if (inheritDc) {
        for (int k = 0; k < numAttrs; k++) {
          if (haarFlag) {
            BestRecBuf[8 * k + 0] = attrRecParentUsIt[k];
          } else {
            BestRecBuf[8 * k + 0] = attrRecParentUsIt[k] - PredDC[k];
          }

          if (enableACRDOInterPred ) {
            if (haarFlag) {
              if (!realInferInLowerLevel)
                BestRecBufIntraLayer[8 * k + 0] = attrRecParentUsIt[k];

              if (!upperInferMode)
                BestRecBufInterLayer[8 * k + 0] = attrRecParentUsIt[k];
            } else {
              if (!realInferInLowerLevel)
                BestRecBufIntraLayer[8 * k + 0] = attrRecParentUsIt[k] - PredDCIntraLayer[k];

              if (!upperInferMode)
                BestRecBufInterLayer[8 * k + 0] = attrRecParentUsIt[k] - PredDCInterLayer[k];
            }
          }
        }
      }

      if (haarFlag) {
        invTransformBlock222<numAttrs, false, HaarKernel>(BestRecBuf, weights);
        if (enableACRDOInterPred) {
          if (!realInferInLowerLevel)
            invTransformBlock222<numAttrs, false, HaarKernel>(BestRecBufIntraLayer, weights);
          if (!upperInferMode)
            invTransformBlock222<numAttrs, false, HaarKernel>(BestRecBufInterLayer, weights);
        }
      } else {
        if (skipinverse) {
          int64_t DCerror[3];
          for (int k = 0; k < numAttrs; k++) {
            DCerror[k] = BestRecBuf[8 * k + 0];
            BestRecBuf[8 * k + 0] = 0;
          }
          for (int cidx = 0; cidx < 8; cidx++) {
            if (!weights[cidx])
              continue;

            for (int k = 0; k < numAttrs; k++) {
              int64_t Correctionterm = fpReduce<kFPFracBits>(
                normsqrtsbuf[cidx] * DCerror[k]);
              BestRecBuf[8 * k + cidx] = Correctionterm;
            }
          }
        } else {
          invTransformBlock222<numAttrs, false, RahtKernel>(BestRecBuf, weights);
        }

        if (enableACRDOInterPred) {
          if (!realInferInLowerLevel) {
            //intralayer
            if (skipinverseintralayer) {
              int64_t DCerror[3];
              for (int k = 0; k < numAttrs; k++) {
                DCerror[k] = BestRecBufIntraLayer[8 * k + 0];
                BestRecBufIntraLayer[8 * k + 0] = 0;
              }
              for (int cidx = 0; cidx < 8; cidx++) {
                if (!weights[cidx])
                  continue;

                for (int k = 0; k < numAttrs; k++) {
                  int64_t Correctionterm = fpReduce<kFPFracBits>(
                    normsqrtsbuf[cidx] * DCerror[k]);
                  BestRecBufIntraLayer[8 * k + cidx] = Correctionterm;
                }
              }
            } else {
              invTransformBlock222<numAttrs, false, RahtKernel>(BestRecBufIntraLayer, weights);
            }
          }

          if (!upperInferMode) {
            //interlayer
            if (skipinverseinterlayer) {
              int64_t DCerror[3];
              for (int k = 0; k < numAttrs; k++) {
                DCerror[k] = BestRecBufInterLayer[8 * k + 0];
                BestRecBufInterLayer[8 * k + 0] = 0;
              }
              for (int cidx = 0; cidx < 8; cidx++) {
                if (!weights[cidx])
                  continue;

                for (int k = 0; k < numAttrs; k++) {
                  int64_t Correctionterm = fpReduce<kFPFracBits>(
                    normsqrtsbuf[cidx] * DCerror[k]);
                  BestRecBufInterLayer[8 * k + cidx] = Correctionterm;
                }
              }
            } else {
              invTransformBlock222<numAttrs, false, RahtKernel>(BestRecBufInterLayer, weights);
            }
          }
        }
      }

      for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
        if (!weights[nodeIdx])
          continue;

        for (int k = 0; k < numAttrs; k++) {
          if (haarFlag) {
            attrRecUs[j * numAttrs + k] = BestRecBuf[8 * k + nodeIdx];
          } else {
            BestRecBuf[8 * k + nodeIdx] += attrBestPredIt[8 * k + nodeIdx]; //Sample Domain Reconstruction
            attrRecUs[j * numAttrs + k] = BestRecBuf[8 * k + nodeIdx];
          }

          if (enableACRDOInterPred ) {
            if (haarFlag) {
              if (!realInferInLowerLevel)
                intraLayerAttrRecUs[j * numAttrs + k] = BestRecBufIntraLayer[8 * k + nodeIdx];

              if (!upperInferMode)
                interLayerAttrRecUs[j * numAttrs + k] = BestRecBufInterLayer[8 * k + nodeIdx];
            } else {
              if (!realInferInLowerLevel) {
                BestRecBufIntraLayer[8 * k + nodeIdx] += attrPredIntraLayer[8 * k + nodeIdx];
                intraLayerAttrRecUs[j * numAttrs + k] = BestRecBufIntraLayer[8 * k + nodeIdx];
              }

              if (!upperInferMode) {
                BestRecBufInterLayer[8 * k + nodeIdx] += attrPredInterLayer[8 * k + nodeIdx];
                interLayerAttrRecUs[j * numAttrs + k] = BestRecBufInterLayer[8 * k + nodeIdx];
              }
            }
          }
        }

        // scale values for next level
        if (!haarFlag) {
          if (weights[nodeIdx] > 1) {
            uint64_t w = weights[nodeIdx];
            int shift = 5 * ((w > 1024) + (w > 1048576));
            uint64_t rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;
            for (int k = 0; k < numAttrs; k++) {
              BestRecBuf[8 * k + nodeIdx] = fpReduce<kFPFracBits>(
                (BestRecBuf[8 * k + nodeIdx] >> shift) * rsqrtWeight);
              if (enableACRDOInterPred) {
                if (!realInferInLowerLevel)
                  BestRecBufIntraLayer[8 * k + nodeIdx] = fpReduce<kFPFracBits>(
                    (BestRecBufIntraLayer[8 * k + nodeIdx] >> shift) * rsqrtWeight);
                if (!upperInferMode)
                  BestRecBufInterLayer[8 * k + nodeIdx] = fpReduce<kFPFracBits>(
                    (BestRecBufInterLayer[8 * k + nodeIdx] >> shift) * rsqrtWeight);
              }
            }
          }
        }

        for (int k = 0; k < numAttrs; k++) {
          attrRec[j * numAttrs + k] = BestRecBuf[8 * k + nodeIdx];
          if(enableACRDOInterPred) {
            if (!realInferInLowerLevel)
              intraLayerAttrRec[j * numAttrs + k] = BestRecBufIntraLayer[8 * k + nodeIdx];

            if (!upperInferMode)
              interLayerAttrRec[j * numAttrs + k] = BestRecBufInterLayer[8 * k + nodeIdx];
          }
        }
        j++;
      }
      i += nodeCnt;
    }


    if (enableACRDOInterPred) {
      double curCost = fpToDouble<32>(
        curEstimate.costBits() + coder.getModeBits());

      double intraLayerCost;
      if (!realInferInLowerLevel)
        intraLayerCost = fpToDouble<32>(intraLayerEstimate.costBits());

      double interLayerCost;
      if (!upperInferMode)
        interLayerCost = fpToDouble<32>(interLayerEstimate.costBits());

      int64_t ifactor = 1 << 24;
      double dfactor = (double)(ifactor);
      double rdcostMCPred = distMCPredLayer * dfactor + dlambda * curCost;

      double rdcostintraLayer;
      if (inheritDc && !realInferInLowerLevel)
        rdcostintraLayer = distintraLayer * dfactor + dlambda * intraLayerCost;
      else
        rdcostintraLayer = std::numeric_limits<double>::infinity();

      double rdcostinterLayer;
      if (!upperInferMode)
        rdcostinterLayer = distinterLayer * dfactor + dlambda * interLayerCost;
      else
        rdcostinterLayer = std::numeric_limits<double>::infinity();

      if (rdcostinterLayer < rdcostintraLayer && rdcostinterLayer < rdcostMCPred) {
        for (int k = 0; k < numAttrs; ++k)
          std::copy_n(interLayerCoeffBufItBeginK[k], sumNodes, coeffBufItBeginK[k]);
        std::swap(interLayerAttrRec, attrRec);
        std::swap(interLayerAttrRecUs, attrRecUs);
        curEstimate = interLayerEstimate;
        intraLayerEstimate = interLayerEstimate;
        abh.attr_layer_code_mode[depth] = 1;
        trainZeros = interLayerTrainZeros;
        intraLayerTrainZeros = interLayerTrainZeros;
        preLayerCodeMode = 1;
        coder.reloadPrevStates();
      }
      else if (rdcostintraLayer < rdcostinterLayer && rdcostintraLayer < rdcostMCPred) {
        for (int k = 0; k < numAttrs; ++k)
          std::copy_n(intraLayerCoeffBufItBeginK[k], sumNodes, coeffBufItBeginK[k]);
        std::swap(intraLayerAttrRec, attrRec);
        std::swap(intraLayerAttrRecUs, attrRecUs);
        curEstimate = intraLayerEstimate;
        interLayerEstimate = intraLayerEstimate;
        abh.attr_layer_code_mode[depth] = 2;
        trainZeros = intraLayerTrainZeros;
        interLayerTrainZeros = intraLayerTrainZeros;
        preLayerCodeMode = 2;
        coder.reloadPrevStates();
      }
      else {
        intraLayerEstimate = curEstimate;
        interLayerEstimate = curEstimate;
        abh.attr_layer_code_mode[depth] = 0;
        intraLayerTrainZeros = trainZeros;
        interLayerTrainZeros = trainZeros;
        preLayerCodeMode = 0;
      }
      curEstimate.resetCostBits();
      intraLayerEstimate.resetCostBits();
      interLayerEstimate.resetCostBits();
    }

    // preserve current weights/positions for later search
    //weightsParent = weightsLf;
    weightsParent.resize(weightsLf.size());
    for (int t = 0; t < weightsLf.size(); t++) {
      weightsParent[t] = weightsLf[t];
    }

    sumNodes = 0;
    ++depth;
  }

  // process duplicate points at level 0
  std::swap(attrRec, attrRecParent);
  auto attrRecParentIt = attrRecParent.cbegin();
  auto attrsHfIt = attrsHf.cbegin();

  for (int i = 0, out = 0, iEnd = weightsLf.size(); i < iEnd; i++) {
    int weight = weightsLf[i].weight;

    // unique points have weight = 1
    if (weight == 1) {
      for (int k = 0; k < numAttrs; k++)
        attrRec[out++] = *attrRecParentIt++;
      continue;
    }

    // duplicates
    Qps nodeQp = {
      weightsLf[i].qp[0] >> regionQpShift,
      weightsLf[i].qp[1] >> regionQpShift};

    int64_t attrSum[3];
    int64_t attrRecDc[3];
    int64_t sqrtWeight = fastIsqrt(uint64_t(weight));

    int64_t sumCoeff = 0;
    for (int k = 0; k < numAttrs; k++) {
      attrSum[k] = attrsLf[i * numAttrs + k];
      attrRecDc[k] = *attrRecParentIt++;
      if (!haarFlag) {
        attrRecDc[k] = fpReduce<kFPFracBits>(
          attrRecDc[k] * sqrtWeight);
      }
    }

    int64_t rsqrtWeight;
    for (int w = weight - 1; w > 0; w--) {
      RahtKernel kernel(w, 1);
      HaarKernel haarkernel(w, 1);
      int shift = 5 * ((w > 1024) + (w > 1048576));
      rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;

      auto quantizers = qpset.quantizers(qpLayer, nodeQp);
      for (int k = 0; k < numAttrs; k++) {
        auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

        int64_t transformBuf[2];

        // invert the initial reduction (sum)
        // NB: read from (w-1) since left side came from attrsLf.
        transformBuf[1] = fpExpand<kFPFracBits>(
          attrsHfIt[(w - 1) * numAttrs + k]);
        if (haarFlag) {
          attrSum[k] -= transformBuf[1] >> 1;
          transformBuf[1] += attrSum[k];
          transformBuf[0] = attrSum[k];
        } else {
          attrSum[k] -= transformBuf[1];
          transformBuf[0] = attrSum[k];

          // NB: weight of transformBuf[1] is by construction 1.
          transformBuf[0] = fpReduce<kFPFracBits>(
            (transformBuf[0] >> shift) * rsqrtWeight);
        }

        if (haarFlag) {
          haarkernel.fwdTransform(transformBuf[0], transformBuf[1]);
        } else {
          kernel.fwdTransform(transformBuf[0], transformBuf[1]);
        }

        auto coeff = fpReduce<kFPFracBits>(transformBuf[1]);
        assert(coeff <= INT_MAX && coeff >= INT_MIN);
        *coeffBufItK[k]++ = coeff =
          q.quantize(coeff << kFixedPointAttributeShift);
        transformBuf[1] = fpExpand<kFPFracBits>(
          divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift));

        sumCoeff += std::abs(q.quantize(coeff << kFixedPointAttributeShift));

        // inherit the DC value
        transformBuf[0] = attrRecDc[k];

        if (haarFlag) {
          haarkernel.invTransform(transformBuf[0], transformBuf[1]);
        } else {
          kernel.invTransform(transformBuf[0], transformBuf[1]);
        }

        attrRecDc[k] = transformBuf[0];
        attrRec[out + w * numAttrs + k] = transformBuf[1];
        if (w == 1)
          attrRec[out + k] = transformBuf[0];
      }

      // Track RL for RDOQ
      if (sumCoeff == 0)
        trainZeros++;
      else
        trainZeros = 0;
    }

    attrsHfIt += (weight - 1) * numAttrs;
    out += weight * numAttrs;
  }

  // write-back reconstructed attributes
  assert(attrRec.size() == numAttrs * numPoints);
  auto attrOut = attributes;
  for (auto attr : attrRec) {
    auto v = attr + kFPOneHalf >> kFPFracBits;
    *attrOut++ = PCCClip(v, 0, std::numeric_limits<attr_t>::max());
  }
}

//============================================================================
/*
 * RAHT Fixed Point
 *
 * Inputs:
 * quantStepSizeLuma = Quantization step
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 *
 * Outputs:
 * weights = list of 'voxelCount' weights associated with each transform coefficient
 * coefficients = quantized transformed attributes array, in column-major order
 * binaryLayer = binary layer where each coefficient was generated
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void
regionAdaptiveHierarchicalTransform(
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
  attr::ModeEncoder& encoder)
{
  auto getMode =
    [&qpset, &rahtPredParams](
      attr::ModeEncoder& encoder,
      int nodeCnt, int predCtxLevel,
      bool enableIntraPred, bool enableInterPred, Mode parentMode,
      Mode neighborsMode, int numAttrs, int64_t weights[],
      std::vector<int64_t>::const_iterator attrRecParent,
      int64_t* transformBuf, std::vector<Mode> modes, const int qpLayer,
      const Qps* nodeQp, bool upperInferMode, bool inheritDC) {
      if (nodeCnt > 1) {
        if (upperInferMode) {
          if (enableInterPred)
            return Mode::Inter;
          else if (enableIntraPred)
            return Mode::Intra;
          else
            return Mode::Null;
        }

        if (predCtxLevel < 0)
          return enableIntraPred ? Mode::Intra : Mode::Null;

        int predCtxMode = attr::getInferredMode(enableIntraPred,
          enableInterPred, nodeCnt, parentMode, neighborsMode);

        if (encoder.isInterEnabled()) {
          if (!enableIntraPred && enableInterPred) {
            modes.resize(2);
          }
          else if (!enableIntraPred && !enableInterPred)
            return Mode::Null;
        } else {
          if (!enableIntraPred)
            return Mode::Null;
        }

        encoder.getEntropy(predCtxMode, predCtxLevel);
        Mode predMode;
        if(rahtPredParams.integer_haar_enable_flag) {
          predMode = attr::choseMode<HaarKernel>(
            encoder, transformBuf, modes, weights, numAttrs, qpset, qpLayer,
            nodeQp, inheritDC);
        } else {
          predMode = attr::choseMode<RahtKernel>(
            encoder, transformBuf, modes, weights, numAttrs, qpset, qpLayer,
            nodeQp, inheritDC);
        }
        encoder.encode(predCtxMode, predCtxLevel, predMode);
        encoder.updateModeBits(predMode);
        return predMode;
      }
      return Mode::Null;
    };

  switch (attribCount) {
  case 3:
    if (!rahtPredParams.integer_haar_enable_flag)
      uraht_process_encoder<false,3>(
        rahtPredParams, abh,qpset, pointQpOffsets, voxelCount, mortonCode,
        attributes, attributes_mc, coefficients, encoder, getMode);
    else
      uraht_process_encoder<true, 3>(
        rahtPredParams, abh, qpset, pointQpOffsets, voxelCount, mortonCode,
        attributes, attributes_mc, coefficients, encoder, getMode);
    break;
  default:
    throw std::runtime_error("attribCount != 3 not tested yet");
  }
}

//============================================================================

}  // namespace pcc
