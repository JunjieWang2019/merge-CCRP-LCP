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
#include <cfloat>
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
  if (value) {
    sumCostBits += fpEntropyProbaLUT<log2scaleRes, 32>(probResGt0[k]);  //Gt0
    probResGt0[k] +=  scaleRes - probResGt0[k] >> windowLog2;
  }
  else {
    sumCostBits += fpEntropyProbaLUT<log2scaleRes, 32>(scaleRes - probResGt0[k]);  //Gt0
    probResGt0[k] -= probResGt0[k] >> windowLog2;
  }


  uint32_t mag = abs(value);
  if (mag) {
    sumCostBits += 1ULL << 32;  //sign bit.
    if (mag > 1) {
      sumCostBits += fpEntropyProbaLUT<log2scaleRes, 32>(probResGt1[k]);  //Gt1
      probResGt1[k] += scaleRes - probResGt1[k] >> windowLog2 ;

      mag--;
      sumCostBits += 2 * (int64_t(pcc::ilog2(mag)) << 32) + 1;  //EG0 approximation.

    }
    else {
      sumCostBits += fpEntropyProbaLUT<log2scaleRes, 32>(scaleRes - probResGt1[k]);  //Gt1
      probResGt1[k] -=  probResGt1[k] >> windowLog2;
    }
  }
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
reduceUniqueEncoder(
  int numNodes,
  std::vector<UrahtNodeEncoder>* weightsIn,
  std::vector<UrahtNodeEncoder>* weightsOut)
{
  // process a single level of the tree
  int64_t posPrev = -1;
  auto weightsInWrIt = weightsIn->begin();
  auto weightsInRdIt = weightsIn->cbegin();
  for (int i = 0; i < numNodes; i++) {
    const auto& node = *weightsInRdIt++;

    // copy across unique nodes
    if (node.pos != posPrev) {
      posPrev = node.pos;
      *weightsInWrIt++ = node;
      continue;
    }

    // duplicate node
    (weightsInWrIt - 1)->weight += node.weight;
    weightsOut->push_back(node);
    if (haarFlag) {
      for (int k = 0; k < numAttrs; k++) {
        auto temp = node.sumAttr[k] - (weightsInWrIt - 1)->sumAttr[k];
        (weightsInWrIt - 1)->sumAttr[k] += temp >> 1;
        weightsOut->back().sumAttr[k] = temp;
      }
    }
    else {
      for (int k = 0; k < numAttrs; k++)
        (weightsInWrIt - 1)->sumAttr[k] += node.sumAttr[k];
    }
    for (int k = 0; k < numAttrs; k++)
      (weightsInWrIt - 1)->sumAttrInter[k] += node.sumAttrInter[k];
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================
struct HaarNode {
  int64_t pos;
  int32_t attr[3];
  int32_t attrInter[3];
};

template<bool haarFlag, int numAttrs>
int
reduceDepthEncoder(
  int level,
  int numNodes,
  std::vector<UrahtNodeEncoder>* weightsIn,
  std::vector<UrahtNodeEncoder>* weightsOut)
{
  // process a single level of the tree
  int64_t posPrev = -1;
  auto weightsInRdIt = weightsIn->begin();
  for (int i = 0; i < weightsIn->size(); ) {

    // this is a new node
    UrahtNodeEncoder last = weightsInRdIt[i];
    posPrev = last.pos;
    last.firstChildIdx = i;

    // look for same node
    int i2 = i + 1;
    for (; i2 < weightsIn->size(); i2++)
      if ((posPrev ^ weightsInRdIt[i2].pos) >> level)
        break;

    // process same nodes
    last.numChildren = i2 - i;

    for (int j = i + 1; j < i2; j++) {
      const auto node = weightsInRdIt[j];
      last.weight += node.weight;
      // TODO: fix local qp to be same in encoder and decoder
      last.qp[0] = (last.qp[0] + node.qp[0]) >> 1;
      last.qp[1] = (last.qp[1] + node.qp[1]) >> 1;

      if (!haarFlag) {
        for (int k = 0; k < numAttrs; k++)
          last.sumAttr[k] += node.sumAttr[k];
        for (int k = 0; k < numAttrs; k++)
          last.sumAttrInter[k] += node.sumAttrInter[k];
      }
    }

    //attribute processign for Haar per direction in the interval [i, i2[
    if (haarFlag) {
      HaarNode haarNode[4];
      // first direction (at most 8 nodes)
      int numNode = 0;
      int64_t posPrevH = -1;
      for (int j = i; j < i2; j++) {
        const auto node = weightsInRdIt[j];
        bool newPair = (posPrevH ^ node.pos) >> (level - 2) != 0;
        posPrevH = node.pos;

        if (newPair) {
          haarNode[numNode].pos = node.pos;
          for (int k = 0; k < numAttrs; k++) {
            haarNode[numNode].attr[k] = node.sumAttr[k];
            haarNode[numNode].attrInter[k] = node.sumAttrInter[k];
          }
          numNode++;
        }
        else {
          auto& lastH = haarNode[numNode - 1];
          for (int k = 0; k < numAttrs; k++) {
            auto temp = node.sumAttr[k] - lastH.attr[k];
            lastH.attr[k] += temp >> 1;
            temp = node.sumAttrInter[k] - lastH.attrInter[k];
            lastH.attrInter[k] += temp >> 1;
          }
        }
      }

      // second direction (at most 4 nodes)
      int numNode2 = 0;
      posPrevH = -1;
      for (int j = 0; j < numNode; j++) {
        const auto node = haarNode[j];
        bool newPair = (posPrevH ^ node.pos) >> (level - 1) != 0;
        posPrevH = node.pos;

        if (newPair) {
          haarNode[numNode2].pos = node.pos;
          for (int k = 0; k < numAttrs; k++) {
            haarNode[numNode2].attr[k] = node.attr[k];
            haarNode[numNode2].attrInter[k] = node.attrInter[k];
          }
          numNode2++;
        }
        else {
          auto& lastH = haarNode[numNode2 - 1];
          for (int k = 0; k < numAttrs; k++) {
            auto temp = node.attr[k] - lastH.attr[k];
            lastH.attr[k] += temp >> 1;
            temp = node.attrInter[k] - lastH.attrInter[k];
            lastH.attrInter[k] += temp >> 1;
          }
        }
      }

      // third direction (at most 2 nodes).
      auto& lastH = haarNode[0];
      for (int k = 0; k < numAttrs; k++) {
        last.sumAttr[k] = lastH.attr[k];
        last.sumAttrInter[k] = lastH.attrInter[k];
      }

      if (numNode2 == 2) {
        lastH = haarNode[1];
        for (int k = 0; k < numAttrs; k++) {
          auto temp = lastH.attr[k] - last.sumAttr[k];
          last.sumAttr[k] += temp >> 1;
          temp = lastH.attrInter[k] - last.sumAttrInter[k];
          last.sumAttrInter[k] += temp >> 1;
        }
      }
    } // end Haar attributes

    weightsOut->push_back(last);
    i = i2;
  }

  // number of nodes in next level
  return weightsOut->size();
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
int getRate(int trainZeros)
{
  static const int LUTbins[11] = { 1, 2, 3, 5, 5, 7, 7, 9, 9, 11 ,11 };
  int Rate = LUTbins[trainZeros > 10 ? 10 : trainZeros];
  if (trainZeros > 10) {
    int temp = trainZeros - 11 + 1; // prefix k =2
    while (temp) {
      Rate += 2;
      temp >>= 1;
    }
    Rate += -1 + 2; // suffix  k=2
  }
  return Rate;
}

//============================================================================
// Core transform process (for encoder)

template<bool haarFlag, int numAttrs>
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
  attr::ModeEncoder& coder)
{
  // coefficients are stored in three planar arrays.  coeffBufItK is a set
  // of iterators to each array.
  int32_t* coeffBufItK[3] = {
    coeffBufIt, coeffBufIt + numPoints, coeffBufIt + numPoints * 2,
  };

  if (numPoints == 1) {
    auto quantizers = qpset.quantizers(0, pointQpOffsets[0]);
    for (int k = 0; k < numAttrs; k++) {
      auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

      auto coeff = attributes[k];
      assert(coeff <= INT_MAX && coeff >= INT_MIN);
      *coeffBufItK[k]++ = coeff = q.quantize(coeff << kFixedPointAttributeShift);
      attributes[k] = divExp2RoundHalfUp(
        q.scale(coeff), kFixedPointAttributeShift);
    }
    return;
  }


  // --------- ascend tree per depth  -----------------
  // create leaf nodes
  int regionQpShift = 4;
  std::vector<UrahtNodeEncoder>  weightsHf;
  std::vector < std::vector<UrahtNodeEncoder>> weightsLfStack;

  weightsLfStack.emplace_back();
  weightsLfStack.back().reserve(numPoints);
  auto weightsLfRef = &weightsLfStack.back();
  auto attr = attributes;
  auto attrPredictor = attributes_mc;

  for (int i = 0; i < numPoints; i++) {
    UrahtNodeEncoder node;
    node.pos = positions[i];
    node.weight = 1;
    node.qp = {
      int16_t(pointQpOffsets[i][0] << regionQpShift),
      int16_t(pointQpOffsets[i][1] << regionQpShift) };
    for (int k = 0; k < numAttrs; k++) {
      node.sumAttr[k] = (*attr++);
      node.sumAttrInter[k] = attributes_mc ? (*attrPredictor++) : 0;
    }
    weightsLfRef->emplace_back(node);
  }

  // -----------  bottom up per depth  --------------------
  int numNodes = weightsLfRef->size();
  // for duplicates, skipable if it is known there is no duplicate
  numNodes = reduceUniqueEncoder<haarFlag, numAttrs>(numNodes, weightsLfRef, &weightsHf);
  const bool flagNoDuplicate = weightsHf.size() == 0;
  int numDepth = 0;
  for (int levelD = 3; numNodes > 1; levelD += 3) {
    // one depth reduction
    weightsLfStack.emplace_back();
    weightsLfStack.back().reserve(numNodes / 3);
    weightsLfRef = &weightsLfStack.back();

    auto weightsLfRefold = &weightsLfStack[weightsLfStack.size() - 2];
    numNodes = reduceDepthEncoder<haarFlag, numAttrs>(levelD, numNodes, weightsLfRefold, weightsLfRef);
    numDepth++;
  }

  // --------- initialize stuff ----------------
  // root node
  auto& rootNode = weightsLfStack.back()[0];
  rootNode.mode = Mode::Null;
  rootNode.firstChildIdx = 0;
  assert(rootNode.weight == numPoints);

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

  std::vector<int8_t> numParentNeigh, numGrandParentNeigh;
  if (!rahtPredParams.enable_inter_prediction) {
    numParentNeigh.resize(numPoints);
    numGrandParentNeigh.resize(numPoints);
  }

  // indexes of the neighbouring parents
  int parentNeighIdx[19];
  int childNeighIdx[12][8];

  // Prediction buffers
  int64_t SampleDomainBuff[3 * numAttrs * 8];
  int64_t transformBuf[3 * numAttrs * 8];

  const int numBuffers = (2 + coder.isInterEnabled()) * numAttrs;
  constexpr size_t sizeofBuf = sizeof(int64_t) * 8;
  const size_t sizeofBufs = sizeofBuf * numBuffers;

  int64_t* attrPred = SampleDomainBuff + 8 * numAttrs;;
  int64_t* attrPredTransform = transformBuf + 8 * numAttrs;;
  int64_t* attrReal = SampleDomainBuff;
  int64_t* attrPredIntra;
  int64_t* attrPredInter;
  int64_t* attrBestPredIt;
  int64_t* attrPredIntraTransformIt;
  int64_t* attrPredInterTransformIt;
  int64_t* attrBestPredTransformIt;

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

  std::vector<Mode> modes;
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
  int trainZeros = 0;
  int CccpCoeff = 0;
  PCCRAHTComputeCCCP curlevelCccp;
  int intraLayerTrainZeros = 0;
  int interLayerTrainZeros = 0;
  PCCRAHTACCoefficientEntropyEstimate intraLayerEstimate;
  PCCRAHTACCoefficientEntropyEstimate interLayerEstimate;
  PCCRAHTACCoefficientEntropyEstimate curEstimate;
  std::vector<int> intraLayerACCoeffcients, interLayerACCoeffcients;
  if (enableACRDOInterPred) {
    intraLayerACCoeffcients.resize(numPoints * numAttrs);
    interLayerACCoeffcients.resize(numPoints * numAttrs);
  }

  // -------------- descend tree, loop on depth --------------
  int rootLayer = numDepth;
  int depth = 0;
  int sumNodes = 1; // number of coded coefficients (=nodes) in the layer; NB: in first layer there is also the DC coefficient (=root node) to count
  int preLayerCodeMode = 0;
  for (int levelD = numDepth, isFirst = 1; levelD > 0; /*nop*/) {
    // references
    std::vector<UrahtNodeEncoder>& weightsParent = weightsLfStack[levelD];
    std::vector<UrahtNodeEncoder>& weightsLf = weightsLfStack[levelD - 1];
    sumNodes += weightsLf.size() - weightsParent.size();

    levelD--;
    int level = 3 * levelD;

    CccpCoeff = 0;
    curlevelCccp.reset();

    int distanceToRoot = rootLayer - levelD;
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

    // initial scan position of the coefficient buffer
    //  -> first level = all coeffs
    //  -> otherwise = ac coeffs only
    bool inheritDc = !isFirst;
    bool enableIntraPredInLvl =
      inheritDc && rahtPredParams.prediction_enabled_flag;

    //!!!! needed for RDO ???
    for (auto i = 0; i < weightsParent.size(); i++) {
      if (preLayerCodeMode)
        // preLayerCodeMode == 2 -> 1 INTRA, preLayerCodeMode == 1 -> 2 INTER
        weightsParent[i].mode = (pcc::attr::Mode)(3 - preLayerCodeMode);
    }

    if (isFirst)
      weightsParent[0].mode = Mode::Null; // root node has null mode
    isFirst = 0;

    // select quantiser according to transform layer
    qpLayer = std::min(qpLayer + 1, int(qpset.layers.size()) - 1);

    // prepare reconstruction buffers
    //  previous reconstruction -> attrRecParent
    std::swap(attrRec, attrRecParent);
    std::swap(attrRecUs, attrRecParentUs);
    std::swap(numParentNeigh, numGrandParentNeigh);
    auto numGrandParentNeighIt = numGrandParentNeigh.cbegin();

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
    enableACRDOInterPred = rahtPredParams.raht_enable_inter_intra_layer_RDO
      && enableACInterPred && rahtPredParams.prediction_enabled_flag
      && depth < RDOCodingDepth;

    if (enableACRDOInterPred) {
      curLevelEnableACInterPred = rahtPredParams.prediction_enabled_flag;
      curLevelEnableLayerModeCoding = curLevelEnableACInterPred;

      coder.restoreStates();
      coder.resetModeBits();
    }

    int64_t distinterLayer = 0;
    int64_t distintraLayer = 0;
    int64_t distMCPredLayer = 0;
    double dlambda = 1.;

    //CCRP Parameters
    const bool CCRPFlag = !haarFlag && CCRPEnabled && level <= maxlevelCCPenabled && !CbCrEnabled;

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

      using WeightsType =
        typename std::conditional<haarFlag, bool, int64_t>::type;
      WeightsType weights[8 + 8 + 8 + 8 + 24] = {};

      int64_t interPredictor[8 * 3] = { 0 };
      int childTable[8] = { };

      int64_t sqrtweightsbuf[8] = { 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768 }; //value for weigth is 1
      int64_t normsqrtsbuf[8] = {};
      bool skipinverse = !haarFlag;
      bool skipinverseinterlayer = skipinverse;
      bool skipinverseintralayer = skipinverse;

      Qps nodeQp[8] = {};
      uint8_t occupancy = 0;

      // generate weights, occupancy mask, and fwd transform buffers
      // for all siblings of the current node.
      int nodeCnt = weightsParentIt->numChildren;
      for (int t = 0, j0 = i; t < nodeCnt; t++, j0++) {
        int nodeIdx = (weightsLf[j0].pos >> level) & 0x7;
        childTable[t] = nodeIdx;
        weights[nodeIdx] = weightsLf[j0].weight;
        nodeQp[nodeIdx][0] = weightsLf[j0].qp[0] >> regionQpShift;
        nodeQp[nodeIdx][1] = weightsLf[j0].qp[1] >> regionQpShift;

        // inter predictor
        if (attributes_mc) {
          auto pred = &interPredictor[nodeIdx];
          if (haarFlag) {
            for (int k = 0; k < numAttrs; k++)
              pred[8 * k] = int64_t(weightsLf[j0].sumAttrInter[k]) << kFPFracBits;
          }
          else {
            //int64_t divisor = (int64_t(1) << 30) / weights[nodeIdx];
            int64_t w = weights[nodeIdx];
            int shift = 5 * (1 + (w > 1024) + (w > 1048576));
            int64_t rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;
            int64_t divisor = fpReduce<kFPFracBits>(((int64_t(1) << 30) >> shift) * rsqrtWeight);
            divisor = fpReduce<kFPFracBits>((divisor >> shift) * rsqrtWeight);

            for (int k = 0; k < numAttrs; k++)
              pred[8 * k] = (int64_t(weightsLf[j0].sumAttrInter[k]) * divisor) >> (30 - kFPFracBits);
          }
        }

        occupancy |= 1 << nodeIdx;

        //store grandmode
        weightsLf[j0].grand_mode = depth >= 1 ? weightsParentIt->mode : Mode::Null;

        for (int k = 0; k < numAttrs; k++) {
          attrReal[8 * k + nodeIdx] = fpExpand<kFPFracBits>(int64_t(weightsLf[j0].sumAttr[k]));
        }
      }
      weightsParentIt->occupancy = occupancy;
      int64_t sumweights = weightsParentIt->weight;
      using Kernel =
        typename std::conditional<haarFlag, HaarKernel, RahtKernel>::type;
      mkWeightTree<haarFlag>::template apply<false, Kernel>(weights);

      if (!inheritDc && !rahtPredParams.enable_inter_prediction) {
        for (int j = i, n = 0; n < nodeCnt; n++)
          numParentNeigh[j++] = 19;
      }

      auto attrRecParentUsIt = std::next(
        attrRecParentUs.cbegin(),
        std::distance(weightsParent.begin(), weightsParentIt) * numAttrs);

      // Inter-level prediction:
      //  - Find the parent neighbours of the current node
      //  - Generate prediction for all attributes into transformIntraBuf
      //  - Subtract transformed coefficients from forward transform
      //  - The transformIntraBuf is then used for reconstruction
      bool enableIntraPred = rahtPredParams.enable_inter_prediction
        ? enableIntraPredInLvl && (nodeCnt > 1) && (distanceToRoot > 2)
        : enableIntraPredInLvl; ;
      bool enableInterPred = coder.isInterEnabled() && (nodeCnt > 1);
      //< layer mode coding
      bool enableIntraLayerPred = false;
      bool enableInterLayerPred = false;
      if (enableACRDOInterPred) {
        enableIntraLayerPred = enableIntraPredInLvl && (nodeCnt > 1);
        enableInterLayerPred = curLevelEnableACInterPred && (nodeCnt > 1);
      }

      // inter prediction
      Mode neighborsMode = Mode::size;
      if (enableInterPred || enableInterLayerPred) {
        enableInterPred = true;
        memcpy(attrPredInter, interPredictor, 8 * numAttrs * sizeof(int64_t));
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
              parentNeighIdx, childNeighIdx, weightsLf);

          if (rahtPredParams.prediction_enabled_flag)
            neighborsMode = attr::getNeighborsModeEncoder(parentNeighIdx, weightsParent, voteInterWeight,
              voteIntraWeight, voteInterLayerWeight, voteIntraLayerWeight, weightsLf);
        }
        if (!rahtPredParams.enable_inter_prediction && rahtPredParams.prediction_skip1_flag && nodeCnt == 1) {
          enableIntraPred = false;
          enableIntraLayerPred = false;
          parentNeighCount = 19;
        }
        else if (!rahtPredParams.enable_inter_prediction && *numGrandParentNeighIt < rahtPredParams.prediction_threshold0) {
          enableIntraPred = false;
          enableIntraLayerPred = false;
        }
        else if (parentNeighCount < rahtPredParams.prediction_threshold1) {
          enableIntraPred = false;
          enableIntraLayerPred = false;
        }
        else {
          intraDcPred<haarFlag, numAttrs>(
            true, parentNeighIdx, childNeighIdx, occupancy,
            attrRecParent.begin(), attrRec.begin(), intraLayerAttrRec.begin(),
            interLayerAttrRec.begin(), attrPredIntra, attrPredIntraLayer,
            attrPredInterLayer, rahtPredParams, enableACRDOInterPred);
        }

        if (!rahtPredParams.enable_inter_prediction) {
          for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
            if (!weights[nodeIdx])
              continue;
            numParentNeigh[j++] = parentNeighCount;
          }
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
        FwdTransformBlock222<haarFlag>
        ::template apply(numBuffers, transformBuf, weights);
        if (enableACRDOInterPred) {
          memcpy(transformLayerBuf, SampleLayerBuf, sizeofLayerBuf);
          FwdTransformBlock222<haarFlag>
          ::template apply(numLayerBuffers, transformLayerBuf, weights);
        }
      }
      else {
        // normalise coefficients
        for (int n = 0; n < nodeCnt; n++) {
          const int childIdx = childTable[n];
          int64_t w = weights[childIdx];
          if (w > 1) {
            // Summed attribute values
            int shift = 5 * ((w > 1024) + (w > 1048576));
            int64_t rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;
            for (int k = 0; k < numAttrs; k++) {
              SampleDomainBuff[8 * k + childIdx] = fpReduce<kFPFracBits>((SampleDomainBuff[8 * k + childIdx] >> shift) * rsqrtWeight);
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
        }
        memcpy(transformBuf, SampleDomainBuff, sizeofBufs);
        FwdTransformBlock222<haarFlag>
        ::template apply(numBuffers, transformBuf, weights);
        if (enableACRDOInterPred) {
          memcpy(transformLayerBuf, SampleLayerBuf, sizeofLayerBuf);
          FwdTransformBlock222<haarFlag>
          ::template apply(numLayerBuffers, transformLayerBuf, weights);
        }
      } //else normalize

      if (enableACRDOInterPred) {
        std::copy_n(transformBuf, 8 * numAttrs, transformIntraLayerBuf);
        std::copy_n(transformBuf, 8 * numAttrs, transformInterLayerBuf);
      }

      const bool enableAveragePrediction = enableAveragePredictionLevel
        && (enableIntraPred || enableIntraLayerPred) && enableInterPred;

      if (enableAveragePrediction) {
        bool flag = !isInter(weightsParentIt->mode) && !isIntra(weightsParentIt->mode);
        voteInterWeight += 12 * isInter(weightsParentIt->mode) + 6 * flag;
        voteIntraWeight += 12 * isIntra(weightsParentIt->mode) + 6 * flag;

        if (depth > 1) {
          bool flag_grand = !isInter(weightsParentIt->grand_mode) && !isIntra(weightsParentIt->grand_mode);
          voteInterWeight += 28 * isInter(weightsParentIt->grand_mode) + 14 * flag_grand;
          voteIntraWeight += 28 * isIntra(weightsParentIt->grand_mode) + 14 * flag_grand;
        }

        int64_t weightIntra = divApprox(voteIntraWeight << kFPFracBits, voteInterWeight + voteIntraWeight, 0);
        int64_t weightInter = (1 << kFPFracBits) - weightIntra;


        const bool enableAverageLayerPrediction = curLevelEnableLayerModeCoding && enableInterLayerPred
          && enableIntraLayerPred && enableInterPred;

        if (enableAverageLayerPrediction) {
          std::copy_n(attrPredInter, 8 * numAttrs, attrOrgPredInterLayer);
          std::copy_n(attrPredInterTransformIt, 8 * numAttrs, attrOrgPredInterLayerTransformIt);
        }


        if (enableIntraPred && enableInterPred) {
          auto aP = predCtxLevel < 0 ? attrPredIntra : attrPredInter;
          auto aPT = predCtxLevel < 0 ? attrPredIntraTransformIt : attrPredInterTransformIt;
          for (int t = 0; t < 8 * numAttrs; t++) {
            aP[t] = fpReduce<kFPFracBits>(attrPredIntra[t] * weightIntra + attrPredInter[t] * weightInter);
            aPT[t] = fpReduce<kFPFracBits>(attrPredIntraTransformIt[t] * weightIntra + attrPredInterTransformIt[t] * weightInter);
            if (haarFlag) {
              aP[t] &= kFPIntMask;
              aPT[t] &= kFPIntMask;
            }
          }
        }

        if (enableAverageLayerPrediction) {
          voteInterLayerWeight += 12 * isInter(weightsParentIt->mode) + 6 * flag;
          voteIntraLayerWeight += 12 * isIntra(weightsParentIt->mode) + 6 * flag;

          if (depth > 1) {
            bool flag_grand = !isInter(weightsParentIt->grand_mode) && !isIntra(weightsParentIt->grand_mode);
            voteInterLayerWeight += 28 * isInter(weightsParentIt->grand_mode) + 14 * flag_grand;
            voteIntraLayerWeight += 28 * isIntra(weightsParentIt->grand_mode) + 14 * flag_grand;
          }

          int64_t weightIntraLayer = divApprox(voteIntraLayerWeight << kFPFracBits, voteInterLayerWeight + voteIntraLayerWeight, 0);
          int64_t weightInterLayer = (1 << kFPFracBits) - weightIntraLayer;

          for (int t = 0; t < 8 * numAttrs; t++) {
            attrPredInterLayer[t] = fpReduce<kFPFracBits>(attrPredInterLayer[t] * weightIntraLayer + attrOrgPredInterLayer[t] * weightInterLayer);
            attrPredInterLayerTransformIt[t] = fpReduce<kFPFracBits>(attrPredInterLayerTransformIt[t] * weightIntraLayer + attrOrgPredInterLayerTransformIt[t] * weightInterLayer);
            if (haarFlag) {
              attrPredInterLayer[t] &= kFPIntMask;
              attrPredInterLayerTransformIt[t] &= kFPIntMask;
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

      // ---------- determine best prediction mode -------------
      Mode predMode = Mode::Null;
      if (rahtPredParams.prediction_enabled_flag && nodeCnt > 1) {
        if (upperInferMode) {
          if (enableInterPred)
            predMode = Mode::Inter;
          else if (enableIntraPred)
            predMode = Mode::Intra;
        }
        else {
          if (predCtxLevel < 0)
            predMode = enableIntraPred ? Mode::Intra : Mode::Null;
          else {
            bool modeIsNull = !enableIntraPred && !enableInterPred;
            if (!modeIsNull) {
              std::vector<Mode> modes2 = modes;
              if (enableInterPred && !enableIntraPred)
                modes2.resize(2);

              int predCtxMode = attr::getInferredMode(enableIntraPred,
                enableInterPred, nodeCnt, weightsParentIt->mode, neighborsMode);
              coder.getEntropy(predCtxMode, predCtxLevel);
              if (haarFlag) {
                predMode = attr::choseMode<HaarKernel, numAttrs, WeightsType>(
                  coder, transformBuf, modes2, weights, qpset, qpLayer, nodeQp, inheritDc);
              }
              else {
                predMode = attr::choseMode<RahtKernel, numAttrs, WeightsType>(
                  coder, transformBuf, modes2, weights, qpset, qpLayer, nodeQp, inheritDc);
              }
              coder.encode(predCtxMode, predCtxLevel, predMode);
              coder.updateModeBits(predMode);
            }
          }
        }

      }


      // store pred mode in child nodes, to determine best mode at next depth
      auto predMode2 = int(predMode) >= Mode::Inter ? Mode::Inter : predMode;
      auto _predMode2 = Mode::size;
      if (curLevelEnableLayerModeCoding)
        _predMode2 = enableInterLayerPred && enableInterPred ? Mode::Inter : (enableIntraLayerPred ? Mode::Intra : Mode::Null);
      auto ChildPt = &weightsLf[weightsParentIt->firstChildIdx];
      for (int t = 0; t < nodeCnt; t++, ChildPt++) {
        ChildPt->mode = predMode2;
        ChildPt->_mode = _predMode2;
      }

      if (attr::isNull(predMode)) {
        skipinverse = false;
        attrBestPredIt = SampleDomainBuff + 8 * numAttrs;
        attrBestPredTransformIt = transformBuf + 8 * numAttrs;
        std::fill_n(attrBestPredIt, 8 * numAttrs, 0);
        std::fill_n(attrBestPredTransformIt, 8 * numAttrs, 0);
      }
      else if (attr::isIntra(predMode)) {
        attrBestPredIt = attrPredIntra;
        attrBestPredTransformIt = attrPredIntraTransformIt;
      }
      else {
        attrBestPredIt = attrPredInter;
        attrBestPredTransformIt = attrPredInterTransformIt;
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
      CCRPFilter::Corr curCorrIntra = { 0, 0, 0 };
      CCRPFilter::Corr curCorrInter = { 0, 0, 0 };


      // subtract transformed prediction (should skipping DC, but ok)
      if (!attr::isNull(predMode)) {
        for (int t = 0; t < numAttrs * 8; t++)
          transformBuf[t] -= attrBestPredTransformIt[t];
      }
      if (enableACRDOInterPred) {
        if (!realInferInLowerLevel)
          for (int t = 0; t < numAttrs * 8; t++)
            transformIntraLayerBuf[t] -= attrPredIntraLayerTransformIt[t];
        if (!upperInferMode)
          for (int t = 0; t < numAttrs * 8; t++)
            transformInterLayerBuf[t] -= attrPredInterLayerTransformIt[t];
      }

      // ----------scan blocks -------
      for (int idxB = 0; idxB < 8; idxB++) {
        if ((idxB == 0 || weights[24 + idxB]) // there is always the DC coefficient (empty blocks are not transformed)
          && !(inheritDc && !idxB)) {  // skip the DC coefficient unless at the root of the tree

          // decision for RDOQ
          int64_t sumCoeff = 0;
          const int LUTlog[16] = { 0,   256, 406, 512, 594, 662, 719,  768, 812, 850, 886, 918, 947, 975, 1000, 1024 };
          bool flagRDOQ = false;
          int64_t Qcoeff;
          int64_t intraLayerSumCoeff = 0;
          bool intraLayerFlagRDOQ = false;
          int64_t interLayerSumCoeff = 0;
          bool interLayerFlagRDOQ = false;

          auto quantizers = qpset.quantizers(qpLayer, nodeQp[idxB]);
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
              auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

              const int k8idx = 8 * k + idxB;
              coeff = fpReduce<kFPFracBits>(transformBuf[k8idx]);
              Dist2 += coeff * coeff;

              if (rahtPredParams.cross_chroma_component_prediction_flag && !coder.isInterEnabled()) {
                if (k != 2) {
                  Qcoeff = q.quantize(coeff << kFixedPointAttributeShift);
                  transformRecBuf[k] = fpExpand<kFPFracBits>(divExp2RoundHalfUp(q.scale(Qcoeff), kFixedPointAttributeShift));
                }
                else {
                  transformRecBuf[k] = transformBuf[k8idx] - ((CccpCoeff * transformRecBuf[1]) >> 4);
                  coeff = fpReduce<kFPFracBits>(transformRecBuf[k]);
                  Qcoeff = q.quantize(coeff << kFixedPointAttributeShift);
                }
              }
              else
                Qcoeff = q.quantize(coeff << kFixedPointAttributeShift);

              auto recCoeff = divExp2RoundHalfUp(q.scale(Qcoeff), kFixedPointAttributeShift);
              recDist2 += (coeff - recCoeff) * (coeff - recCoeff);

              sumCoeff += std::abs(Qcoeff);
              Ratecoeff += std::abs(Qcoeff) < 15 ? LUTlog[std::abs(Qcoeff)] : LUTlog[15];
              if (!k)
                lambda0 = q.scale(1);
              lambda = lambda0 * lambda0 * (numAttrs == 1 ? 25 : 35);

              if (enableACRDOInterPred) {
                if (!upperInferMode) {
                  auto interLayerCoeff = fpReduce<kFPFracBits>(transformInterLayerBuf[k8idx]);
                  interLayerDist2 += interLayerCoeff * interLayerCoeff;
                  auto interLayerQcoeff = q.quantize(interLayerCoeff << kFixedPointAttributeShift);

                  auto recInterCoeff = divExp2RoundHalfUp(q.scale(interLayerQcoeff), kFixedPointAttributeShift);
                  interLayerRecDist2 += (interLayerCoeff - recInterCoeff) * (interLayerCoeff - recInterCoeff);

                  interLayerSumCoeff += std::abs(interLayerQcoeff);
                  interLayerRatecoeff += std::abs(interLayerQcoeff) < 15 ? LUTlog[std::abs(interLayerQcoeff)] : LUTlog[15];
                }

                if (!realInferInLowerLevel) {
                  auto intraLayerCoeff = fpReduce<kFPFracBits>(transformIntraLayerBuf[k8idx]);
                  intraLayerDist2 += intraLayerCoeff * intraLayerCoeff;
                  auto intraLayerQcoeff = q.quantize(intraLayerCoeff << kFixedPointAttributeShift);

                  auto recIntraCoeff = divExp2RoundHalfUp(q.scale(intraLayerQcoeff), kFixedPointAttributeShift);
                  intraLayerRecDist2 += (intraLayerCoeff - recIntraCoeff) * (intraLayerCoeff - recIntraCoeff);

                  intraLayerSumCoeff += std::abs(intraLayerQcoeff);
                  intraLayerRatecoeff += std::abs(intraLayerQcoeff) < 15 ? LUTlog[std::abs(intraLayerQcoeff)] : LUTlog[15];
                }
              }
            }

            dlambda = double(lambda);
            if (sumCoeff < 3) {
              int Rate = getRate(trainZeros);
              Rate += (Ratecoeff + 128) >> 8;
              flagRDOQ = (Dist2 << 26) < (lambda * Rate + (recDist2 << 26));
            }

            if (enableACRDOInterPred) {
              if (!realInferInLowerLevel && intraLayerSumCoeff < 3) {
                int Rate = getRate(intraLayerTrainZeros);
                Rate += (intraLayerRatecoeff + 128) >> 8;
                intraLayerFlagRDOQ = (intraLayerDist2 << 26) < (lambda * Rate + (intraLayerRecDist2 << 26));
              }

              if (!upperInferMode && interLayerSumCoeff < 3) {
                int Rate = getRate(interLayerTrainZeros);
                Rate += (interLayerRatecoeff + 128) >> 8;
                interLayerFlagRDOQ = (interLayerDist2 << 26) < (lambda * Rate + (interLayerRecDist2 << 26));
              }
            }
          }

          // Track RL for RDOQ
          trainZeros++;
          trainZeros *= (flagRDOQ || sumCoeff == 0);

          if (enableACRDOInterPred) {
            if (!realInferInLowerLevel) {
              intraLayerTrainZeros++;
              intraLayerTrainZeros *= (intraLayerFlagRDOQ || intraLayerSumCoeff == 0);
            }

            if (!upperInferMode) {
              interLayerTrainZeros++;
              interLayerTrainZeros *= (interLayerFlagRDOQ || interLayerSumCoeff == 0);
            }
          }

          // The RAHT transform
          int64_t quantizedLuma = 0, quantizedChroma = 0;
          int64_t quantizedLumaIntra = 0, quantizedChromaIntra = 0;
          int64_t quantizedLumaInter = 0, quantizedChromaInter = 0;

          int64_t flagRDOnot = !flagRDOQ;
          int64_t flagRDOintra = !intraLayerFlagRDOQ;
          int64_t flagRDOinter = !interLayerFlagRDOQ;

          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + idxB;
            // apply RDOQ
            transformBuf[k8idx] *= flagRDOnot;
            transformRecBuf[k] *= flagRDOnot;
            transformIntraLayerBuf[k8idx] *= flagRDOintra;
            transformInterLayerBuf[k8idx] *= flagRDOinter;

            int64_t CCRPPred = 0, CCRPPredIntra = 0, CCRPPredInter = 0;
            if (k && CCRPFlag) {
              if (k == 1) {
                CCRPPred = quantizedLuma * ccrpFilter.getYCbFilt() >> kCCRPFiltPrecisionbits;
                transformBuf[k8idx] -= fpExpand<kFPFracBits>(CCRPPred);
                if (enableACRDOInterPred) {
                  CCRPPredIntra = quantizedLumaIntra * ccrpFilterIntra.getYCbFilt() >> kCCRPFiltPrecisionbits;
                  transformIntraLayerBuf[k8idx] -= fpExpand<kFPFracBits>(CCRPPredIntra);
                  CCRPPredInter = quantizedLumaInter * ccrpFilterInter.getYCbFilt() >> kCCRPFiltPrecisionbits;
                  transformInterLayerBuf[k8idx] -= fpExpand<kFPFracBits>(CCRPPredInter);
                }
              }
              if (k == 2) {
                CCRPPred = quantizedLuma * ccrpFilter.getYCrFilt() >> kCCRPFiltPrecisionbits;
                transformBuf[k8idx] -= fpExpand<kFPFracBits>(CCRPPred);
                if (enableACRDOInterPred) {
                  CCRPPredIntra = quantizedLumaIntra * ccrpFilterIntra.getYCrFilt() >> kCCRPFiltPrecisionbits;
                  transformIntraLayerBuf[k8idx] -= fpExpand<kFPFracBits>(CCRPPredIntra);
                  CCRPPredInter = quantizedLumaInter * ccrpFilterInter.getYCrFilt() >> kCCRPFiltPrecisionbits;
                  transformInterLayerBuf[k8idx] -= fpExpand<kFPFracBits>(CCRPPredInter);
                }
              }
            }

            auto coeff = fpReduce<kFPFracBits>(transformBuf[k8idx]);
            assert(coeff <= INT_MAX && coeff >= INT_MIN);
            auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];
            coeff = q.quantize(coeff << kFixedPointAttributeShift); // does nothing for Haar

            if (!haarFlag && rahtPredParams.cross_chroma_component_prediction_flag && !coder.isInterEnabled()) {
              if (k == 2) {
                coeff = fpReduce<kFPFracBits>(transformRecBuf[k]);
                coeff = q.quantize(coeff << kFixedPointAttributeShift);
                transformRecBuf[k] = fpExpand<kFPFracBits>(divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift));
                transformRecBuf[k] += CccpCoeff * transformRecBuf[1] >> 4;
                BestRecBuf[k8idx] = transformRecBuf[k];
              }
              BestRecBuf[k8idx] = transformRecBuf[k];
              CoeffRecBuf[nodelvlSum][k] = fpReduce<kFPFracBits>(transformRecBuf[k]);
            }
            else {
              int64_t quantizedValue = divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
              BestRecBuf[k8idx] = fpExpand<kFPFracBits>(quantizedValue);

              if (CCRPFlag) {
                if (k == 0) {
                  quantizedLuma = quantizedValue;
                  curCorr.yy += quantizedLuma * quantizedLuma;
                }
                else {
                  quantizedChroma = quantizedValue + CCRPPred;
                  BestRecBuf[k8idx] = fpExpand<kFPFracBits>(quantizedChroma);
                  skipinverse = skipinverse && (quantizedChroma == 0);
                  if (k == 1)
                    curCorr.ycb += quantizedLuma * quantizedChroma;
                  else
                    curCorr.ycr += quantizedLuma * quantizedChroma;
                }
              }
            }

            *coeffBufItK[k]++ = coeff;
            skipinverse = skipinverse && (coeff == 0);

            if (enableACRDOInterPred) {
              curEstimate.updateCostBits(coeff, k);

              if (!realInferInLowerLevel) {
                auto intraLayerCoeff = fpReduce<kFPFracBits>(transformIntraLayerBuf[k8idx]);
                assert(intraLayerCoeff <= INT_MAX && intraLayerCoeff >= INT_MIN);
                intraLayerCoeff = q.quantize(intraLayerCoeff << kFixedPointAttributeShift); // does nothing for Haar

                intraLayerEstimate.updateCostBits(intraLayerCoeff, k);

                *intraLayerCoeffBufItK[k]++ = intraLayerCoeff;
                skipinverseintralayer = skipinverseintralayer && (intraLayerCoeff == 0);

                int64_t quantizedValueIntra = divExp2RoundHalfUp(q.scale(intraLayerCoeff), kFixedPointAttributeShift);
                BestRecBufIntraLayer[k8idx] = fpExpand<kFPFracBits>(quantizedValueIntra);

                if (CCRPFlag) {
                  if (k == 0) {
                    quantizedLumaIntra = quantizedValueIntra;
                    curCorrIntra.yy += quantizedLumaIntra * quantizedLumaIntra;
                  }
                  else {
                    quantizedChromaIntra = quantizedValueIntra + CCRPPredIntra;
                    BestRecBufIntraLayer[k8idx] = fpExpand<kFPFracBits>(quantizedChromaIntra);
                    skipinverseintralayer = skipinverseintralayer && (quantizedChromaIntra == 0);
                    if (k == 1)
                      curCorrIntra.ycb += quantizedLumaIntra * quantizedChromaIntra;
                    else
                      curCorrIntra.ycr += quantizedLumaIntra * quantizedChromaIntra;
                  }
                }
              }

              if (!upperInferMode) {
                auto interLayerCoeff = fpReduce<kFPFracBits>(transformInterLayerBuf[k8idx]);
                assert(interLayerCoeff <= INT_MAX && interLayerCoeff >= INT_MIN);
                interLayerCoeff = q.quantize(interLayerCoeff << kFixedPointAttributeShift); // does nothing for Haar

                interLayerEstimate.updateCostBits(interLayerCoeff, k);

                skipinverseinterlayer = skipinverseinterlayer && (interLayerCoeff == 0);
                *interLayerCoeffBufItK[k]++ = interLayerCoeff;
                // still in RAHT domain
                int64_t quantizedValueInter =
                  divExp2RoundHalfUp(q.scale(interLayerCoeff), kFixedPointAttributeShift);
                BestRecBufInterLayer[k8idx] = fpExpand<kFPFracBits>(quantizedValueInter);

                if (CCRPFlag) {
                  if (k == 0) {
                    quantizedLumaInter = quantizedValueInter;
                    curCorrInter.yy += quantizedLumaInter * quantizedLumaInter;
                  }
                  else {
                    quantizedChromaInter = quantizedValueInter + CCRPPredInter;
                    BestRecBufInterLayer[k8idx] = fpExpand<kFPFracBits>(quantizedChromaInter);
                    skipinverseinterlayer = skipinverseinterlayer && (quantizedChromaInter == 0);
                    if (k == 1)
                      curCorrInter.ycb += quantizedLumaInter * quantizedChromaInter;
                    else
                      curCorrInter.ycr += quantizedLumaInter * quantizedChromaInter;
                  }
                }
              }

              // quantization reconstruction error
              if (!upperInferMode) {
                int64_t iresidueinterLayer = fpReduce<kFPFracBits>(transformInterLayerBuf[k8idx] - BestRecBufInterLayer[k8idx]);
                distinterLayer += iresidueinterLayer * iresidueinterLayer;
              }
              if (!realInferInLowerLevel) {
                int64_t iresidueintraLayer = fpReduce<kFPFracBits>(transformIntraLayerBuf[k8idx] - BestRecBufIntraLayer[k8idx]);
                distintraLayer += iresidueintraLayer * iresidueintraLayer;
              }
              int64_t iresidueMCPred = fpReduce<kFPFracBits>(transformBuf[k8idx] - BestRecBuf[k8idx]);
              distMCPredLayer += iresidueMCPred * iresidueMCPred;

            }
          } //end loop on attributes

          nodelvlSum++;
        }
      }// end of scan block

      // Transform Domain Pred for Lossless case
      if (haarFlag) {
        for (int t = 0; t < numAttrs * 8; t++)
          BestRecBuf[t] += attrBestPredTransformIt[t];
        if (enableACRDOInterPred) {
          if (!realInferInLowerLevel)
            for (int t = 0; t < numAttrs * 8; t++)
              BestRecBufIntraLayer[t] += attrPredIntraLayerTransformIt[t];
          if (!upperInferMode)
            for (int t = 0; t < numAttrs * 8; t++)
              BestRecBufInterLayer[t] += attrPredInterLayerTransformIt[t];
        }
      }


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


      //compute DC of prediction signals
      int64_t PredDC[3] = { 0,0,0 };
      int64_t PredDCIntraLayer[3] = { 0,0,0 };
      int64_t PredDCInterLayer[3] = { 0,0,0 };
      if (!haarFlag && (!attr::isNull(predMode) || enableACRDOInterPred)) {
        int64_t rsqrtweightsum = fastIrsqrt(sumweights);
        for (int n = 0; n < nodeCnt; n++) {
          const int childIdx = childTable[n];

          int64_t normSqrtW = sqrtweightsbuf[childIdx] * rsqrtweightsum >> 40;
          normsqrtsbuf[childIdx] = normSqrtW;
          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + childIdx;
            if (!attr::isNull(predMode))
              PredDC[k] += fpReduce<kFPFracBits>(normSqrtW * attrBestPredIt[k8idx]);

            if (enableACRDOInterPred) {
              PredDCIntraLayer[k] += fpReduce<kFPFracBits>(normSqrtW * attrPredIntraLayer[k8idx]);
              PredDCInterLayer[k] += fpReduce<kFPFracBits>(normSqrtW * attrPredInterLayer[k8idx]);
            }
          }
        }
      }

      // replace DC coefficient with parent if inheritable
      if (inheritDc) {
        for (int k = 0; k < numAttrs; k++) {
          const int k8 = 8 * k;
          BestRecBuf[k8] = attrRecParentUsIt[k];
          if (!haarFlag)
            BestRecBuf[k8] -= PredDC[k];

          if (enableACRDOInterPred) {
            if (!realInferInLowerLevel) {
              BestRecBufIntraLayer[k8] = attrRecParentUsIt[k];
              if (!haarFlag)
                BestRecBufIntraLayer[k8] -= PredDCIntraLayer[k];
            }
            if (!upperInferMode) {
              BestRecBufInterLayer[k8] = attrRecParentUsIt[k];
              if (!haarFlag)
                BestRecBufInterLayer[k8] -= PredDCInterLayer[k];
            }
          }
        }
      }

      if (haarFlag) {
        InvTransformBlock222<haarFlag>
        ::template apply<numAttrs>(BestRecBuf, weights);
        if (enableACRDOInterPred) {
          if (!realInferInLowerLevel)
            InvTransformBlock222<haarFlag>
            ::template apply<numAttrs>(BestRecBufIntraLayer, weights);
          if (!upperInferMode)
            InvTransformBlock222<haarFlag>
            ::template apply<numAttrs>(BestRecBufInterLayer, weights);
        }
      } else {
        if (skipinverse) {
          int64_t DCerror[numAttrs];
          for (int k = 0; k < numAttrs; k++) {
            DCerror[k] = BestRecBuf[8 * k];
          }
          for (int n = 0; n < nodeCnt; n++) {
            const int childIdx = childTable[n];
            for (int k = 0; k < numAttrs; k++)
              BestRecBuf[8 * k + childIdx] = fpReduce<kFPFracBits>(normsqrtsbuf[childIdx] * DCerror[k]);
          }
        }
        else {
          InvTransformBlock222<haarFlag>
          ::template apply<numAttrs>(BestRecBuf, weights);
        }

        if (enableACRDOInterPred) {
          if (!realInferInLowerLevel) {
            //intralayer
            if (skipinverseintralayer) {
              int64_t DCerror[numAttrs];
              for (int k = 0; k < numAttrs; k++) {
                DCerror[k] = BestRecBufIntraLayer[8 * k];
              }
              for (int n = 0; n < nodeCnt; n++) {
                const int childIdx = childTable[n];
                for (int k = 0; k < numAttrs; k++)
                  BestRecBufIntraLayer[8 * k + childIdx] = fpReduce<kFPFracBits>(normsqrtsbuf[childIdx] * DCerror[k]);
              }
            } else {
              InvTransformBlock222<haarFlag>
              ::template apply<numAttrs>(BestRecBufIntraLayer, weights);
            }
          }

          if (!upperInferMode) {
            //interlayer
            if (skipinverseinterlayer) {
              int64_t DCerror[numAttrs];
              for (int k = 0; k < numAttrs; k++) {
                DCerror[k] = BestRecBufInterLayer[8 * k];
              }
              for (int n = 0; n < nodeCnt; n++) {
                const int childIdx = childTable[n];
                for (int k = 0; k < numAttrs; k++)
                  BestRecBufInterLayer[8 * k + childIdx] = fpReduce<kFPFracBits>(normsqrtsbuf[childIdx] * DCerror[k]);
              }
            } else {
              InvTransformBlock222<haarFlag>
              ::template apply<numAttrs>(BestRecBufInterLayer, weights);
            }
          }
        }
      }

      for (int j = i,  n = 0; n < nodeCnt; j++, n++) {
        const int nodeIdx = childTable[n];

        for (int k = 0; k < numAttrs; k++) {
          const int k8idx = 8 * k + nodeIdx;
          const int recidx = j * numAttrs + k;

          if (!haarFlag)
            BestRecBuf[k8idx] += attrBestPredIt[k8idx]; //Sample Domain Reconstruction
          attrRecUs[recidx] = BestRecBuf[k8idx];

          if (enableACRDOInterPred ) {
            if (!realInferInLowerLevel) {
              if (!haarFlag)
                BestRecBufIntraLayer[k8idx] += attrPredIntraLayer[k8idx];
              intraLayerAttrRecUs[recidx] = BestRecBufIntraLayer[k8idx];
            }
            if (!upperInferMode) {
              if (!haarFlag)
                BestRecBufInterLayer[k8idx] += attrPredInterLayer[k8idx];
              interLayerAttrRecUs[recidx] = BestRecBufInterLayer[k8idx];
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
              const int k8idx = 8 * k + nodeIdx;
              BestRecBuf[k8idx] = fpReduce<kFPFracBits>((BestRecBuf[k8idx] >> shift) * rsqrtWeight);
              if (enableACRDOInterPred) {
                if (!realInferInLowerLevel)
                  BestRecBufIntraLayer[k8idx] = fpReduce<kFPFracBits>((BestRecBufIntraLayer[k8idx] >> shift) * rsqrtWeight);
                if (!upperInferMode)
                  BestRecBufInterLayer[k8idx] = fpReduce<kFPFracBits>((BestRecBufInterLayer[k8idx] >> shift) * rsqrtWeight);
              }
            }
          }
        }

        for (int k = 0; k < numAttrs; k++) {
          const int k8idx = 8 * k + nodeIdx;
          const int recidx = j * numAttrs + k;
          attrRec[recidx] = BestRecBuf[k8idx];
          if(enableACRDOInterPred) {
            if (!realInferInLowerLevel)
              intraLayerAttrRec[recidx] = BestRecBufIntraLayer[k8idx];

            if (!upperInferMode)
              interLayerAttrRec[recidx] = BestRecBufInterLayer[k8idx];
          }
        }
      }
      i += nodeCnt;
    }


    if (enableACRDOInterPred) {
      const double inv2pow32 = 2.328306436538696e-10;
      double curCost = double(curEstimate.costBits() + coder.getModeBits()) * inv2pow32;

      const double dfactor = 16777216.;  //(1<<24)
      double rdcostMCPred = double(distMCPredLayer) * dfactor + dlambda * curCost;

      double rdcostintraLayer = rdcostintraLayer = DBL_MAX;
      if (inheritDc && !realInferInLowerLevel) {
        double intraLayerCost = double(intraLayerEstimate.costBits()) * inv2pow32;
        rdcostintraLayer = double(distintraLayer) * dfactor + dlambda * intraLayerCost;
      }

      double rdcostinterLayer = DBL_MAX;
      if (!upperInferMode) {
        double interLayerCost = double(interLayerEstimate.costBits()) * inv2pow32;
        rdcostinterLayer = double(distinterLayer) * dfactor + dlambda * interLayerCost;
      }

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

    sumNodes = 0;
    ++depth;
    weightsParent.clear();
  }

  // -------------- process duplicate points at level 0 --------------
  if (flagNoDuplicate) { // write-back reconstructed attributes
    auto attrOut = attributes;
    for (auto attr : attrRec) {
      auto v = attr + kFPOneHalf >> kFPFracBits;
      *attrOut++ = PCCClip(v, 0, std::numeric_limits<attr_t>::max());
    }
    return;
  }

  // case there are duplicates
  std::swap(attrRec, attrRecParent);
  auto attrRecParentIt = attrRecParent.cbegin();
  auto attrsHfIt = weightsHf.cbegin();

  std::vector<UrahtNodeEncoder>& weightsLf = weightsLfStack[0];
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
      attrSum[k] = int64_t(weightsLf[i].sumAttr[k]);
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
        transformBuf[1] = fpExpand<kFPFracBits>(int64_t(attrsHfIt[w - 1].sumAttr[k]));
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
        *coeffBufItK[k]++ = coeff = q.quantize(coeff << kFixedPointAttributeShift); // quantize does nothing for Haar
        transformBuf[1] = fpExpand<kFPFracBits>(divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift));

        sumCoeff += std::abs(q.quantize(coeff << kFixedPointAttributeShift)); // quantize does nothing for Haar

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

    attrsHfIt += (weight - 1);
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
  switch (attribCount) {
  case 3:
    if (!rahtPredParams.integer_haar_enable_flag)
      uraht_process_encoder<false,3>(
        rahtPredParams, abh,qpset, pointQpOffsets, voxelCount, mortonCode,
        attributes, attributes_mc, coefficients, encoder);
    else
      uraht_process_encoder<true, 3>(
        rahtPredParams, abh, qpset, pointQpOffsets, voxelCount, mortonCode,
        attributes, attributes_mc, coefficients, encoder);
    break;
  default:
    throw std::runtime_error("attribCount != 3 not tested yet");
  }
}

//============================================================================

}  // namespace pcc
