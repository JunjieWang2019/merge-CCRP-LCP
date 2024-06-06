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

using pcc::attr::Mode;
using namespace pcc::RAHT;

namespace pcc {

//============================================================================
// remove any non-unique leaves from a level in the uraht tree

int
reduceUniqueDecoder(
  int numNodes,
  std::vector<UrahtNodeDecoder>* weightsIn,
  std::vector<UrahtNodeDecoder>* weightsOut)
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
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================
// Split a level of values into sum and difference pairs.

int
reduceLevelDecoder(
  int level,
  int numNodes,
  std::vector<UrahtNodeLight>* weightsIn,
  std::vector<UrahtNodeLight>* weightsOut)
{
  // process a single level of the tree
  int64_t posPrev = -1;
  auto weightsInWrIt = weightsIn->begin();
  auto weightsInRdIt = weightsIn->cbegin();
  for (int i = 0; i < numNodes; i++) {
    auto& node = *weightsInRdIt++;
    bool newPair = (posPrev ^ node.pos) >> level;
    posPrev = node.pos;
    if (newPair) {
      *weightsInWrIt++ = node;
    }
    else {
      auto& left = *(weightsInWrIt - 1);
      left.weight += node.weight;
      left.qp[0] = (left.qp[0] + node.qp[0]) >> 1;
      left.qp[1] = (left.qp[1] + node.qp[1]) >> 1;
      weightsOut->push_back(node);
    }
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

int
reduceDepthDecoder(
  int level,
  int numNodes,
  std::vector<UrahtNodeDecoder>* weightsIn,
  std::vector<UrahtNodeDecoder>* weightsOut)
{
  // process a single level of the tree
  int64_t posPrev = -1;
  auto weightsInRdIt = weightsIn->begin();
  for (int i = 0; i < weightsIn->size(); i++) {
    auto& node = *weightsInRdIt++;
    bool newNode = (posPrev ^ node.pos) >> level;
    posPrev = node.pos;
    if (newNode) {
      weightsOut->push_back(node);
      weightsOut->back().numChildren = 1;
      weightsOut->back().firstChildIdx = i;
    }
    else {
      auto& last = weightsOut->back();
      last.weight += node.weight;
      last.numChildren++;
      // TODO: fix local qp to be same in encoder and decoder
      last.qp[0] = (last.qp[0] + node.qp[0]) >> 1;
      last.qp[1] = (last.qp[1] + node.qp[1]) >> 1;
    }
  }

  // number of nodes in next level
  return weightsOut->size();
}

//============================================================================

template<typename T>
Mode getNeighborsModeDecoder(
  const int parentNeighIdx[19],
  const std::vector<T>& weightsParent,
  int& voteInterWeight,
  int& voteIntraWeight,
  std::vector<UrahtNodeDecoder>& weightsLf,
  pcc::attr::Mode modeParents)
{
  int vote[4] = { 0, 0, 0, 0 }; // Null, intra, inter, size;

  for (int i = 1; i < 19; i++) {
    if (parentNeighIdx[i] == -1)
      continue;

    auto uncle = std::next(weightsParent.begin(), parentNeighIdx[i]);

    auto uncleMode =
      modeParents == attr::Mode::size ? uncle->mode : modeParents;
    vote[uncleMode]++;


    if (uncle->decoded) {
      auto cousin = &weightsLf[uncle->firstChildIdx];
      for (int t = 0; t < uncle->numChildren; t++, cousin++) {
        vote[cousin->mode] += 3;
      }
    }
  }

  auto parent = std::next(weightsParent.begin(), parentNeighIdx[0]);
  auto parentMode =
    modeParents == attr::Mode::size ? parent->mode : modeParents;

  voteIntraWeight = vote[1] * 2 + vote[0];
  voteInterWeight = vote[2] * 2 + vote[0];
  voteIntraWeight += isNull(parentMode) + 2 * isIntra(parentMode);
  voteInterWeight += isNull(parentMode) + 2 * isInter(parentMode);

  if (vote[0] > vote[1] && vote[0] > vote[2])
    return Mode::Null;
  if (vote[1] > vote[2])
    return Mode::Intra;
  return Mode::Inter;
}

//============================================================================

template<typename It, typename It2>
void
findNeighboursChildrenDecoder(
  It first,
  It2 firstChild,
  int level,
  uint8_t occupancy,
  const int parentNeighIdx[19],
  int childNeighIdx[12][8],
  std::vector<UrahtNodeDecoder>& weightsLf)
{
  memset(childNeighIdx, -1, 96 * sizeof(int));
  static const uint8_t occuMasks[12] = {3,  5,  15, 17, 51, 85,
                                        10, 34, 12, 68, 48, 80};
  static const uint8_t occuShift[12] = {6, 5, 4, 3, 2, 1, 3, 1, 2, 1, 2, 3};

  const int* pN = &parentNeighIdx[7];
  for (int i = 0; i < 9; i++, pN++) {
    if (*pN == -1)
      continue;

    auto neiIt = first + *pN;
    uint8_t mask =
      (neiIt->occupancy >> occuShift[i]) & occupancy & occuMasks[i];
    if (mask) {
      auto it = &weightsLf[neiIt->firstChildIdx];
      for (int t = 0; t < neiIt->numChildren; t++, it++) {
        int nodeIdx = ((it->pos >> level) & 0x7) - occuShift[i];
        if ((nodeIdx >= 0) && ((mask >> nodeIdx) & 1)) {
          childNeighIdx[i][nodeIdx] = t + neiIt->firstChildIdx;
        }
      }
    }
  }

  for (int i = 9; i < 12; i++, pN++) {
    if (*pN == -1)
      continue;

    auto neiIt = first + *pN;
    uint8_t mask =
      (neiIt->occupancy << occuShift[i]) & occupancy & occuMasks[i];
    if (mask) {
      auto it = &weightsLf[neiIt->firstChildIdx];;
      for (int t = 0; t < neiIt->numChildren; t++, it++) {
        int nodeIdx = ((it->pos >> level) & 0x7) + occuShift[i];
        if ((nodeIdx < 8) && ((mask >> nodeIdx) & 1)) {
          childNeighIdx[i][nodeIdx] = t + neiIt->firstChildIdx;
        }
      }
    }
  }
}

//============================================================================
// translateLayer

template<bool haarFlag, int numAttrs>
bool
translateLayerDecoder(
  std::vector<int64_t>& layerAttr,
  size_t level,
  size_t count_mc,
  int64_t* morton_mc,
  int* attr_mc,
  std::vector<UrahtNodeDecoder>& weightsLf)
{
  bool flagMCmatchCurrent = true;
  // associate mean attribute of MC PC to each unique node
  auto layer = layerAttr.begin();
  layerAttr.resize(weightsLf.size() * numAttrs);
  for (int i = 0, j = 0;
      i < weightsLf.size() && j < count_mc;
      i++, layer += numAttrs) {
    int64_t pos = weightsLf[i].pos >> level;

    while (j < count_mc && pos > (morton_mc[j] >> level))
      j++;

    for (size_t k = 0; k < numAttrs; k++)
      layer[k] = -1;

    int64_t weight = 0;
    auto attr = &attr_mc[numAttrs * j];
    int64_t sumAtt[3] = { 0,0,0 };
    while (j < count_mc && pos == (morton_mc[j] >> level)) {
      weight++;
      for (size_t k = 0; k < numAttrs; k++)
        sumAtt[k] += static_cast<int64_t>(*attr++) << kFPFracBits;
      j++;
    }

    if (weight) {
      for (int k = 0; k < numAttrs; k++)
        layer[k] = sumAtt[k];

      if (weight != 1) {
        for (int k = 0; k < numAttrs; k++) {
          layer[k] /= weight;
          if (haarFlag)
            layer[k] = (layer[k] >> kFPFracBits) << kFPFracBits;
        }
      }
    }

    flagMCmatchCurrent = flagMCmatchCurrent && (weight > 0);
  }
  return flagMCmatchCurrent;
}

//============================================================================
// Generate the spatial prediction of a block.

template<bool haarFlag, int numAttrs, typename It>
void
intraDcPredDecoder(
  const int parentNeighIdx[19],
  const int childNeighIdx[12][8],
  int occupancy,
  It first,
  It firstChild,
  int64_t* predBuf,
  const RahtPredictionParams& rahtPredParams)
{
  static const uint8_t predMasks[19] = {255, 240, 204, 170, 192, 160, 136,
                                        3,   5,   15,  17,  51,  85,  10,
                                        34,  12,  68,  48,  80};

  const auto& predWeightParent = rahtPredParams.predWeightParent;
  const auto& predWeightChild = rahtPredParams.predWeightChild;
  int weightSum[8] = {0, 0, 0, 0, 0, 0, 0, 0};

  memset(predBuf, 0, 8 * numAttrs * sizeof(predBuf));

  int64_t neighValue[3];
  int64_t childNeighValue[3];
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
    }
    else {
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
              childNeighValue[k] = *childNeighValueIt++ * predWeightChild[i];

            for (int k = 0; k < numAttrs; k++)
              predBuf[8 * k + j] += childNeighValue[k];
          }
          else {
            weightSum[j] += predWeightParent[7 + i];
            for (int k = 0; k < numAttrs; k++)
              predBuf[8 * k + j] += neighValue[k];
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
        for (int k = 0; k < numAttrs; k++)
          div(predBuf[8*k+i]);
      }
      if (haarFlag) {
        for (int k = 0; k < numAttrs; k++)
          predBuf[8*k+i] = fpExpand<kFPFracBits>(
            fpReduce<kFPFracBits>(predBuf[8*k+i]));
      }
    }
  }
}

//============================================================================
// Core transform process (for decoder)

template<bool haarFlag, int numAttrs>
inline void
uraht_process_decoder(
  const RahtPredictionParams& rahtPredParams,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  const int numPoints,
  int64_t* positions,
  int* attributes,
  const int numPoints_mc,
  int64_t* positions_mc,
  int* attributes_mc,
  int32_t* coeffBufIt,
  attr::ModeDecoder& coder)
{
  // coefficients are stored in three planar arrays.
  // coeffBufItK is a set of iterators to each array.
  int32_t* coeffBufItK[3] =
    {coeffBufIt, coeffBufIt + numPoints, coeffBufIt + numPoints * 2};

  // early termination only one point
  if (numPoints == 1) {
    auto quantizers = qpset.quantizers(0, pointQpOffsets[0]);
    for (int k = 0; k < numAttrs; k++) {
      auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];
      int64_t coeff = *coeffBufItK[k]++;
      attributes[k] = divExp2RoundHalfUp(
        q.scale(coeff), kFixedPointAttributeShift);
    }
    return;
  }


  // --------- ascend tree per depth  -----------------
  // create leaf nodes
  int regionQpShift = 4;
  std::vector<UrahtNodeDecoder>  weightsHf;
  std::vector < std::vector<UrahtNodeDecoder>> weightsLfStack;

  weightsLfStack.emplace_back();
  weightsLfStack.back().reserve(numPoints);
  auto weightsLfRef = &weightsLfStack.back();

  for (int i = 0; i < numPoints; i++) {
    UrahtNodeDecoder node;
    node.pos = positions[i];
    node.weight = 1;
    node.qp = {
      int16_t(pointQpOffsets[i][0] << regionQpShift),
      int16_t(pointQpOffsets[i][1] << regionQpShift)};
    weightsLfRef->emplace_back(node);
  }

  // -----------  bottom up per depth  --------------------
  int numNodes = weightsLfRef->size();
  // for duplicates, skipable if it is known there is no duplicate
  numNodes = reduceUniqueDecoder(numNodes, weightsLfRef, &weightsHf);
  const bool flagNoDuplicate = weightsHf.size() == 0;
  int numDepth = 0;
  for (int levelD = 3; numNodes > 1; levelD += 3) {
    // one depth reduction
    weightsLfStack.emplace_back();
    weightsLfStack.back().reserve(numNodes / 3);
    weightsLfRef = &weightsLfStack.back();

    auto weightsLfRefold = &weightsLfStack[weightsLfStack.size() - 2];
    numNodes = reduceDepthDecoder(
      levelD, numNodes, weightsLfRefold, weightsLfRef);
    numDepth++;
  }

  // --------- initialize stuff ----------------
  // root node
  auto& rootNode = weightsLfStack.back()[0];
  rootNode.mode = Mode::Null;
  rootNode.firstChildIdx = 0;
  assert(rootNode.weight == numPoints);

  const int RDOCodingDepth = abh.attr_layer_code_mode.size();
  const bool enableACInterPred =
    rahtPredParams.enable_inter_prediction && (numPoints_mc > 0);

  coder.setInterEnabled(
    rahtPredParams.prediction_enabled_flag && enableACInterPred);

  // reconstruction buffers
  std::vector<int32_t> attrRec, attrRecParent;
  attrRec.resize(numPoints * numAttrs);
  attrRecParent.resize(numPoints * numAttrs);

  std::vector<int64_t> attrRecUs, attrRecParentUs;
  attrRecUs.resize(numPoints * numAttrs);
  attrRecParentUs.resize(numPoints * numAttrs);

  std::vector<int8_t> numParentNeigh, numGrandParentNeigh;
  numParentNeigh.resize(numPoints);
  numGrandParentNeigh.resize(numPoints);

  // indexes of the neighbouring parents
  int parentNeighIdx[19];
  int childNeighIdx[12][8];

  // Prediction buffers
  int64_t SampleDomainBuff[2 * numAttrs * 8];
  int64_t transformBuf[2 * numAttrs * 8];

  int64_t* attrPredIntra = SampleDomainBuff;
  int64_t* attrPredInter;
  int64_t* attrBestPred;
  int64_t* attrPredIntraTransform = transformBuf;
  int64_t* attrPredInterTransform;
  int64_t* attrBestPredTransform;

  // modes, inter
  std::vector<Mode> modes;
  modes.push_back(Mode::Null);
  std::vector<int64_t> interTree;
  if (coder.isInterEnabled()) {
    attrPredInter = attrPredIntra;
    attrPredInterTransform = attrPredIntraTransform;
    modes.push_back(Mode::Inter);
    attrPredIntra += 8 * numAttrs;
    attrPredIntraTransform += 8 * numAttrs;
    interTree.reserve(numPoints * numAttrs);
  }
  modes.push_back(Mode::Intra);

  // quant layer selection
  auto qpLayer = 0;

  // cross channel prediction
  int CccpCoeff = 0;
  PCCRAHTComputeCCCP curlevelCccp;


  // -------------- descend tree, loop on depth --------------
  int rootLayer = numDepth;
  int depth = 0;
  int preLayerCodeMode = 0;
  pcc::attr::Mode modeParents = pcc::attr::Mode::size;
  for (int levelD = numDepth, isFirst = 1; levelD > 0; /*nop*/) {
    // references
    std::vector<UrahtNodeDecoder>& weightsParent = weightsLfStack[levelD];
    std::vector<UrahtNodeDecoder>& weightsLf = weightsLfStack[levelD - 1];

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
    }
    else if (rahtPredParams.prediction_enabled_flag) {
      predCtxLevel = layerD - rahtPredParams.intra_mode_level;
      if (predCtxLevel >= NUMBER_OF_LEVELS_MODE)
        predCtxLevel = NUMBER_OF_LEVELS_MODE - 1;
    }

    bool upperInferMode =
      coder.isInterEnabled()
      && distanceToRoot < rahtPredParams.upper_mode_level
      && distanceToRoot < RDOCodingDepth - rahtPredParams.mode_level + 1;

    const bool enableAveragePredictionLevel =
      rahtPredParams.enable_average_prediction
      && distanceToRoot >= (
        rootLayer - rahtPredParams.mode_level
        - rahtPredParams.upper_mode_level_for_average_prediction + 1)
      && distanceToRoot < (
        rootLayer - rahtPredParams.mode_level
        + rahtPredParams.lower_mode_level_for_average_prediction + 1)
      && !upperInferMode && coder.isInterEnabled();

    // Motion compensation
    bool flagMCmatchCurrent = true;
    if (coder.isInterEnabled()) {
      flagMCmatchCurrent = translateLayerDecoder<haarFlag, numAttrs>(
        interTree, level, numPoints_mc, positions_mc, attributes_mc,
        weightsLf);
    }

    // initial scan position of the coefficient buffer
    //  -> first level = all coeffs
    //  -> otherwise = ac coeffs only
    bool inheritDc = !isFirst;
    bool enableIntraPredInLvl =
      inheritDc && rahtPredParams.prediction_enabled_flag;

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

    bool curLevelEnableLayerModeCoding = false;
    bool curLevelEnableACInterPred = false;
    bool enableACRDOInterPred =
      rahtPredParams.raht_enable_inter_intra_layer_RDO && enableACInterPred
      && rahtPredParams.prediction_enabled_flag && depth < RDOCodingDepth;

    if (enableACRDOInterPred && rahtPredParams.prediction_enabled_flag) {
      curLevelEnableLayerModeCoding = abh.attr_layer_code_mode[depth];
      curLevelEnableACInterPred = abh.attr_layer_code_mode[depth] == 1;
    }

    // -------------- loop on nodes of the depth --------------
    int i = 0;
    auto attrRecParentUsIt = attrRecParentUs.begin();
    for (auto weightsParentIt = weightsParent.begin();
        weightsParentIt < weightsParent.end();
        weightsParentIt++, attrRecParentUsIt += numAttrs) {

      if (modeParents != attr::Mode::size)
        // replaces parent mode for the whole depth;
        // done also in function getNeighborsModeDecoder
        weightsParentIt->mode = modeParents;

      int64_t weights[8 + 8 + 8 + 8 + 24] = {};
      bool skipinverse = !haarFlag;

      // generate weights, occupancy mask, and fwd transform
      // for all siblings of the current node.
      Qps nodeQp[8] = {};
      uint8_t occupancy = 0;
      const int nodeCnt = weightsParentIt->numChildren;
      for (int t = 0, j0 = i; t < nodeCnt; t++, j0++) {
        int nodeIdx = (weightsLf[j0].pos >> level) & 0x7;
        weights[nodeIdx] = weightsLf[j0].weight;
        nodeQp[nodeIdx][0] = weightsLf[j0].qp[0] >> regionQpShift;
        nodeQp[nodeIdx][1] = weightsLf[j0].qp[1] >> regionQpShift;
        occupancy |= 1 << nodeIdx;
      }
      weightsParentIt->occupancy = occupancy;
      int64_t sumweights = weightsParentIt->weight;

      // Inter-level prediction:
      memset(SampleDomainBuff, 0, sizeof(SampleDomainBuff));
      memset(transformBuf, 0, sizeof(transformBuf));

      bool enableIntraPred = false;
      bool enableInterPred = false;
      bool enableIntraLayerPred = false;
      bool enableInterLayerPred = false;

      if (enableACRDOInterPred && curLevelEnableLayerModeCoding) {
        enableIntraLayerPred = enableIntraPredInLvl && nodeCnt > 1;
        enableInterLayerPred = curLevelEnableACInterPred && nodeCnt > 1;
      } else {
        enableInterPred = coder.isInterEnabled() && nodeCnt > 1;
        enableIntraPred =
          rahtPredParams.enable_inter_prediction
          ? enableIntraPredInLvl && nodeCnt > 1 && distanceToRoot > 2
          : enableIntraPredInLvl;

      }

      // inter prediction
      if (enableInterPred || enableInterLayerPred) {
        if (!flagMCmatchCurrent) {
          auto inter = std::next(
            interTree.begin(), weightsParentIt->firstChildIdx * numAttrs);
          int availablePrediction = 0;
          for (int nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
            if (!((occupancy >> nodeIdx) & 0x1))
              continue;

            if (!(*inter < 0))
              availablePrediction |= 0x1 << nodeIdx;
            inter += numAttrs;
          }
          // always true if projection on decoded geo.
          enableInterPred = availablePrediction == occupancy;
        }
        else
          enableInterPred = true;
      }


      // ------- compute enable intra/iter and neighbor counts ----------
      int parentNeighCount = 19;
      if (enableIntraPred || enableIntraLayerPred) {
        bool flag1 =
          !rahtPredParams.enable_inter_prediction
          && rahtPredParams.prediction_skip1_flag && nodeCnt == 1;
        bool flag2 =
          !rahtPredParams.enable_inter_prediction
          && *numGrandParentNeighIt < rahtPredParams.prediction_threshold0;

        parentNeighCount = flag1 ? 19 : 0;
        if (!flag1 && !flag2) {
          parentNeighCount = findNeighbours(
            weightsParent.begin(), weightsParent.end(), weightsParentIt,
            level + 3, occupancy, parentNeighIdx);
        }
        else {
          enableIntraPred = false;
          enableIntraLayerPred = false;
        }

        if (parentNeighCount < rahtPredParams.prediction_threshold1) {
          enableIntraPred = false;
          enableIntraLayerPred = false;
        }
      }

      // store number of neighbors
      for (int t = 0; t < nodeCnt; t++) {
        numParentNeigh[i + t] = parentNeighCount;
      }

      if (inheritDc) {
        numGrandParentNeighIt++;
        weightsParentIt->decoded = true;
      }

      // ---------- determine best prediction mode -------------
      Mode predMode = Mode::Null;
      Mode neighborsMode = Mode::size;
      int voteInterWeight = 1, voteIntraWeight = 1;
      if (curLevelEnableLayerModeCoding) {
        if (enableInterLayerPred && enableInterPred)
          predMode = Mode::Inter;
        else if (enableIntraLayerPred)
          predMode = Mode::Intra;
      }
      else {
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
                neighborsMode = getNeighborsModeDecoder(
                  parentNeighIdx, weightsParent, voteInterWeight,
                  voteIntraWeight, weightsLf, modeParents);
                int predCtxMode = attr::getInferredMode(
                  enableIntraPred, enableInterPred, nodeCnt,
                  weightsParentIt->mode, neighborsMode);
                predMode = coder.decode(predCtxMode, predCtxLevel);
              }
            }
          }
        }
      }

      // store pred mode in child nodes, to determine best mode at next depth
      auto predMode2 = int(predMode) >= Mode::Inter ? Mode::Inter : predMode;
      auto ChildPt = &weightsLf[weightsParentIt->firstChildIdx];
      for (int t = 0; t < weightsParentIt->numChildren; t++, ChildPt++) {
        ChildPt->mode = predMode2;
      }

      const bool enableAveragePrediction =
        enableAveragePredictionLevel
        && (enableIntraPred || enableIntraLayerPred) && enableInterPred;
      const bool enableAverageLayerPrediction =
        curLevelEnableLayerModeCoding
        && enableInterLayerPred && enableIntraLayerPred && enableInterPred;
      const bool isAveragePred =
        enableAveragePrediction
        && (attr::isIntra(predMode) && predCtxLevel < 0
          || attr::isInter(predMode) && predCtxLevel >= 0
          || enableAverageLayerPrediction);


      // ---------- prepare best predictor attrBestPred attrBestPredTransform -------------
      if (attr::isNull(predMode)) { // null mode
        skipinverse = false;
        if (haarFlag) {
          attrBestPredTransform = transformBuf;
          memset(attrBestPredTransform, 0, numAttrs * 8 * sizeof(int64_t));
        }
        else {
          attrBestPred = SampleDomainBuff;
          memset(attrBestPred, 0, 8 * numAttrs * sizeof(int64_t));
        }
      }
      else if (attr::isIntra(predMode) || isAveragePred) { // intra or average mode
        if (rahtPredParams.subnode_prediction_enabled_flag)
          findNeighboursChildrenDecoder(
            weightsParent.begin(), weightsLf.begin(), level, occupancy,
            parentNeighIdx, childNeighIdx, weightsLf);

        intraDcPredDecoder<haarFlag, numAttrs>(
          parentNeighIdx, childNeighIdx, occupancy, attrRecParent.begin(),
          attrRec.begin(), attrPredIntra, rahtPredParams);
        attrBestPred = attrPredIntra;
        attrBestPredTransform = attrPredIntraTransform;
      }
      else { // inter mode
        attrBestPred = attrPredInter;
        attrBestPredTransform = attrPredInterTransform;
      }

      // prepare inter and average mode
      if (attr::isInter(predMode) || isAveragePred) {
        auto pred = attrPredInter;
        auto inter = std::next(
          interTree.begin(), weightsParentIt->firstChildIdx * numAttrs);

        if (!flagMCmatchCurrent) {
          bool notCalculatedParentDc = true;
          int64_t parentDc[3];
          for (int nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
            if (!((occupancy >> nodeIdx) & 0x1))
              continue;

            if (*inter < 0) { // always wrong if projection on decoded geo.
              inter += numAttrs;
              if (notCalculatedParentDc) {
                computeParentDc<haarFlag, numAttrs>(
                  sumweights, attrRecParentUsIt, parentDc);
                notCalculatedParentDc = false;
              }
              for (int k = 0; k < numAttrs; k++)
                pred[8 * k + nodeIdx] = parentDc[k];
            }
            else {
              for (int k = 0; k < numAttrs; k++) {
                pred[8 * k + nodeIdx] = *(inter++);
              }
            }
          }
        }
        else { // MC matches current
          for (int nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
            if ((occupancy >> nodeIdx) & 0x1)
              for (int k = 0; k < numAttrs; k++)
                pred[8 * k + nodeIdx] = *(inter++);
          }
        }

      }


      // -------  denormalize best predictor   ----------
      int64_t PredDC[3] = {0, 0, 0};
      int64_t sqrtweightsbuf[8] = {};
      int64_t normsqrtsbuf[8] = {};
      mkWeightTree<true, RahtKernel>(weights);

      if (!attr::isNull(predMode)) { // only if not NULL
        if (haarFlag) {
          memcpy(
            attrBestPredTransform, attrBestPred,
            8 * numAttrs * sizeof(int64_t));
          fwdTransformBlock222<numAttrs, HaarKernel>(
            attrBestPredTransform, weights);
          if (isAveragePred) {
            memcpy(
              attrPredInterTransform, attrPredInter,
              8 * numAttrs * sizeof(int64_t));
            fwdTransformBlock222<numAttrs, HaarKernel>(
              attrPredInterTransform, weights);
          }
        }
        else {
          // normalise predicted attribute values
          for (int childIdx = 0; childIdx < 8; childIdx++) {
            if (weights[childIdx] <= 1) {
              sqrtweightsbuf[childIdx] = 32768;
              continue;
            }

            int64_t sqrtWeight = fastIsqrt(uint64_t(weights[childIdx]));
            sqrtweightsbuf[childIdx] = sqrtWeight;
            for (int k = 0; k < numAttrs; k++) {
              attrBestPred[8 * k + childIdx] = fpReduce<kFPFracBits>(
                attrBestPred[8 * k + childIdx] * sqrtWeight);
              if (isAveragePred)
                attrPredInter[8 * k + childIdx] = fpReduce<kFPFracBits>(
                  attrPredInter[8 * k + childIdx] * sqrtWeight);
            }
          }
        }

        if (isAveragePred) {
          if (neighborsMode == Mode::size)
            neighborsMode = getNeighborsModeDecoder(
              parentNeighIdx, weightsParent, voteInterWeight, voteIntraWeight,
              weightsLf, modeParents);

          bool flag =
            !isInter(weightsParentIt->mode)
            && !isIntra(weightsParentIt->mode);

          voteInterWeight += 12 * isInter(weightsParentIt->mode) + 6 * flag;
          voteIntraWeight += 12 * isIntra(weightsParentIt->mode) + 6 * flag;

          int64_t weightIntra = divApprox(
            voteIntraWeight << kFPFracBits,
            voteInterWeight + voteIntraWeight, 0);
          int64_t weightInter = (1 << kFPFracBits) - weightIntra;

          for (int nodeIdx = 0; nodeIdx < 8; nodeIdx++)
            for (int k = 0; k < numAttrs; k++) {
              if (haarFlag) {
                attrBestPredTransform[8 * k + nodeIdx] =
                  fpReduce<kFPFracBits>(
                    attrPredInterTransform[8 * k + nodeIdx] * weightInter
                    + attrBestPredTransform[8 * k + nodeIdx] * weightIntra);
                attrBestPredTransform[8 * k + nodeIdx] &= kFPIntMask;
              }
              else {
                attrBestPred[8 * k + nodeIdx] = fpReduce<kFPFracBits>(
                  attrPredInter[8 * k + nodeIdx] * weightInter
                  + attrBestPred[8 * k + nodeIdx] * weightIntra);
              }
            }
        }

        // Compute DC of best pred
        if (!haarFlag) {
          int64_t rsqrtweightsum = fastIrsqrt(sumweights);
          for (int childIdx = 0; childIdx < 8; childIdx++) {
            if (weights[childIdx]) {
              int64_t normalizedsqrtweight =
                sqrtweightsbuf[childIdx] * rsqrtweightsum >> 40;
              normsqrtsbuf[childIdx] = normalizedsqrtweight;
              for (int k = 0; k < numAttrs; k++)
                PredDC[k] += fpReduce<kFPFracBits>(
                  normalizedsqrtweight * attrBestPred[8 * k + childIdx]);
            }
          }
        }
      }


      // per-coefficient operations:
      //  - read in quantised coefficients
      //  - inverse quantise + add transform domain prediction
      int64_t BestRecBuf[3][8] = {0};
      int64_t CoeffRecBuf[8][3] = {0};
      int nodelvlSum = 0;
      int64_t transformRecBuf[3] = {0, 0, 0};
      bool enableCCCP =
        !haarFlag && rahtPredParams.cross_chroma_component_prediction_flag
        && !coder.isInterEnabled();

      scanBlock(weights, [&](int idx) {
        // skip the DC coefficient unless at the root of the tree
        if (inheritDc && !idx)
          return;

        // The RAHT transform
        auto quantizers = qpset.quantizers(qpLayer, nodeQp[idx]);
        for (int k = 0; k < numAttrs; k++) {
          auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

          int64_t coeff = *coeffBufItK[k]++;
          skipinverse = skipinverse && (coeff == 0);
          CoeffRecBuf[nodelvlSum][k] = divExp2RoundHalfUp(
            q.scale(coeff), kFixedPointAttributeShift);
          transformRecBuf[k] = fpExpand<kFPFracBits>(
            CoeffRecBuf[nodelvlSum][k]);
          BestRecBuf[k][idx] = transformRecBuf[k];

          if (haarFlag) {
            // Transform Domain Pred for Lossless case
            BestRecBuf[k][idx] += attrBestPredTransform[8 * k + idx];
          }
        }

        // cross chroma
        if (enableCCCP && numAttrs == 3) {
          BestRecBuf[2][idx] += (CccpCoeff * transformRecBuf[1]) >> 4;
          CoeffRecBuf[nodelvlSum][2] = fpReduce<kFPFracBits>(
            BestRecBuf[2][idx]);
        }

        nodelvlSum++;
        });

      // compute last component coefficient
      if (numAttrs == 3 && nodeCnt > 1 && enableCCCP) {
        CccpCoeff = curlevelCccp.computeCrossChromaComponentPredictionCoeff(
          nodelvlSum, CoeffRecBuf);
      }
      // replace DC coefficient with parent if inheritable
      if (inheritDc) {
        for (int k = 0; k < numAttrs; k++) {
          if (haarFlag)
            BestRecBuf[k][0] = attrRecParentUsIt[k];
          else
            BestRecBuf[k][0] = attrRecParentUsIt[k] - PredDC[k];
        }
      }

      if (haarFlag) {
        invTransformBlock222<numAttrs, false, HaarKernel>(
          BestRecBuf, weights);
      }
      else {
        if (skipinverse) {
          int64_t DCerror[3];
          for (int k = 0; k < numAttrs; k++) {
            DCerror[k] = BestRecBuf[k][0];
            BestRecBuf[k][0] = 0;
          }
          for (int cidx = 0; cidx < 8; cidx++) {
            if (weights[cidx]) {
              for (int k = 0; k < numAttrs; k++)
                BestRecBuf[k][cidx] = fpReduce<kFPFracBits>(
                  normsqrtsbuf[cidx] * DCerror[k]);
            }
          }
        }
        else {
          invTransformBlock222<numAttrs, true, RahtKernel>(
            BestRecBuf, weights);
        }
      }

      for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
        if (!weights[nodeIdx])
          continue;

        for (int k = 0; k < numAttrs; k++) {
          if (!haarFlag)
            //Sample Domain Reconstruction
            BestRecBuf[k][nodeIdx] += attrBestPred[8 * k + nodeIdx];
          attrRecUs[j * numAttrs + k] = BestRecBuf[k][nodeIdx];
        }

        // scale values for next level
        if (!haarFlag) {
          if (weights[nodeIdx] > 1) {
            uint64_t w = weights[nodeIdx];
            int shift = 5 * ((w > 1024) + (w > 1048576));
            uint64_t rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;
            for (int k = 0; k < numAttrs; k++)
              BestRecBuf[k][nodeIdx] = fpReduce<kFPFracBits>(
                (BestRecBuf[k][nodeIdx] >> shift) * rsqrtWeight);
          }
        }

        for (int k = 0; k < numAttrs; k++)
          attrRec[j * numAttrs + k] = BestRecBuf[k][nodeIdx];

        j++;
      }

      i += nodeCnt;
    } // end loop on nodes of depth

    // parent mode for next depth
    if (enableACRDOInterPred) {
      preLayerCodeMode = 0;
      if (abh.attr_layer_code_mode[depth])
        preLayerCodeMode = 1 + (abh.attr_layer_code_mode[depth] != 1);

      // preLayerCodeMode == 2 -> 1 INTRA, preLayerCodeMode == 1 -> 2 INTER
      modeParents = (pcc::attr::Mode)(3 - preLayerCodeMode);
      ++depth;
    }
  } // end loop on depth


  // -------------- process duplicate points at level 0 --------------
  if (flagNoDuplicate) { // write-back reconstructed attributes
    auto attrOut = attributes;
    for (auto attr : attrRec)
      *attrOut++ = attr + kFPOneHalf >> kFPFracBits;
    return;
  }

  // case there are duplicates
  std::swap(attrRec, attrRecParent);
  auto attrRecParentIt = attrRecParent.cbegin();

  std::vector<int64_t> attrsHf;
  attrsHf.resize(weightsHf.size()* numAttrs);
  auto attrsHfIt = attrsHf.cbegin();

  std::vector<UrahtNodeDecoder>& weightsLf = weightsLfStack[0];
  for (int i = 0, out = 0, iEnd = weightsLf.size(); i < iEnd; i++) {
    // unique points have weight = 1
    int weight = weightsLf[i].weight;
    if (weight == 1) {
      for (int k = 0; k < numAttrs; k++)
        attrRec[out++] = *attrRecParentIt++;
      continue;
    }

    // duplicates
    Qps nodeQp = {
      weightsLf[i].qp[0] >> regionQpShift,
      weightsLf[i].qp[1] >> regionQpShift};

    int64_t attrRecDc[3];
    int64_t sqrtWeight = fastIsqrt(uint64_t(weight));

    for (int k = 0; k < numAttrs; k++) {
      attrRecDc[k] = *attrRecParentIt++;
      if (!haarFlag) {
        attrRecDc[k] = fpReduce<kFPFracBits>(
          attrRecDc[k] * sqrtWeight);
      }
    }

    for (int w = weight - 1; w > 0; w--) {
      RahtKernel kernel(w, 1);
      HaarKernel haarkernel(w, 1);

      auto quantizers = qpset.quantizers(qpLayer, nodeQp);
      for (int k = 0; k < numAttrs; k++) {
        auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

        int64_t transformBuf[2];
        int64_t coeff = *coeffBufItK[k]++;
        transformBuf[1] = fpExpand<kFPFracBits>(
          divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift));
        transformBuf[0] = attrRecDc[k]; // inherit the DC value

        if (haarFlag)
          haarkernel.invTransform(transformBuf[0], transformBuf[1]);
        else
          kernel.invTransform(transformBuf[0], transformBuf[1]);

        attrRecDc[k] = transformBuf[0];
        attrRec[out + w * numAttrs + k] = transformBuf[1];
        if (w == 1)
          attrRec[out + k] = transformBuf[0];
      }
    }

    attrsHfIt += (weight - 1) * numAttrs;
    out += weight * numAttrs;
  }

  // write-back reconstructed attributes
  assert(attrRec.size() == numAttrs * numPoints);
  auto attrOut = attributes;
  for (auto attr : attrRec)
    *attrOut++ = attr + kFPOneHalf >> kFPFracBits;
}

//============================================================================
/*
 * inverse RAHT Fixed Point
 *
 * Inputs:
 * quantStepSizeLuma = Quantization step
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 * coefficients = quantized transformed attributes array, in column-major order
 *
 * Outputs:
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void
regionAdaptiveHierarchicalInverseTransform(
  const RahtPredictionParams& rahtPredParams,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  const int attribCount,
  const int voxelCount,
  int64_t* mortonCode,
  int* attributes,
  const int voxelCount_mc,
  int64_t* mortonCode_mc,
  int* attributes_mc,
  int* coefficients,
  attr::ModeDecoder& decoder)
{
  switch (attribCount) {
  case 3:
    if (!rahtPredParams.integer_haar_enable_flag)
      uraht_process_decoder<false, 3>(
        rahtPredParams, abh,qpset, pointQpOffsets, voxelCount, mortonCode,
        attributes, voxelCount_mc, mortonCode_mc, attributes_mc, coefficients,
        decoder);
    else
      uraht_process_decoder<true, 3>(
        rahtPredParams, abh, qpset, pointQpOffsets, voxelCount, mortonCode,
        attributes, voxelCount_mc, mortonCode_mc, attributes_mc, coefficients,
        decoder);
    break;
  default:
    throw std::runtime_error("attribCount != 3 not tested yet");
  }
}

//============================================================================

}  // namespace pcc
