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

#if defined _MSC_VER && defined _DEBUG
// These variables are not used but must be instanciated for MSVC in debug mode
int UrahtNodeDecoderHaar::weight;
std::array<int16_t, 2> UrahtNodeDecoderHaar::qp;
#endif

//============================================================================
// remove any non-unique leaves from a level in the uraht tree

template<int numAttrs, typename Node>
int
reduceUniqueDecoder(
  const bool isInter,
  int numNodes,
  std::vector<Node>* weightsIn,
  std::vector<Node>* weightsOut)
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
    if (isInter)
      for (int k = 0; k < numAttrs; k++)
        (weightsInWrIt - 1)->sumAttrInter[k] += node.sumAttrInter[k];

    weightsOut->push_back(node);
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================

template <int numAttrs, typename Node>
struct reduceDepthDecoder_ {
  template <typename It>
  static inline void process(
    const bool isInter, const int level, const int i, const int i2,
    It& weightsInRdIt, Node& last)
  {
    throw std::runtime_error("Not implemented");
  }
};


template <int numAttrs>
struct reduceDepthDecoder_<numAttrs, UrahtNodeDecoder> {
  template <typename It>
  static inline void process(
    const bool isInter, const int level, const int i, const int i2,
    It& weightsInRdIt, UrahtNodeDecoder& last)
  {
    for (int j = i + 1; j < i2; j++) {
      const auto node = weightsInRdIt[j];
      last.weight += node.weight;
      // TODO: fix local qp to be same in encoder and decoder
      last.qp[0] = (last.qp[0] + node.qp[0]) >> 1;
      last.qp[1] = (last.qp[1] + node.qp[1]) >> 1;

      if (isInter)
        for (int k = 0; k < numAttrs; k++)
          last.sumAttrInter[k] += node.sumAttrInter[k];
    }
  }
};

template <int numAttrs>
struct reduceDepthDecoder_<numAttrs, UrahtNodeDecoderHaar> {
  template <typename It>
  static inline void process(
    const bool isInter, const int level, const int i, const int i2,
    It& weightsInRdIt, UrahtNodeDecoderHaar& last)
  {
    if (!isInter)
      return;

    //attribute processign for Haar per direction in the interval [i, i2[
    struct {
      int64_t pos;
      int32_t attrInter[3];
    } haarNode[4];

    // first direction (at most 8 nodes)
    int numNode = 0;
    int64_t posPrevH = -1;
    for (int j = i; j < i2; j++) {
      const auto node = weightsInRdIt[j];
      bool newPair = (posPrevH ^ node.pos) >> (level - 2) != 0;
      posPrevH = node.pos;

      if (newPair) {
        haarNode[numNode].pos = node.pos;
        for (int k = 0; k < numAttrs; k++)
          haarNode[numNode].attrInter[k] = node.sumAttrInter[k];
        numNode++;
      } else {
        auto& lastH = haarNode[numNode - 1];
        for (int k = 0; k < numAttrs; k++) {
          auto temp = node.sumAttrInter[k] - lastH.attrInter[k];
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
        for (int k = 0; k < numAttrs; k++)
          haarNode[numNode2].attrInter[k] = node.attrInter[k];
        numNode2++;
      } else {
        auto& lastH = haarNode[numNode2 - 1];
        for (int k = 0; k < numAttrs; k++) {
          auto temp = node.attrInter[k] - lastH.attrInter[k];
          lastH.attrInter[k] += temp >> 1;
        }
      }
    }

    // third direction (at most 2 nodes).
    auto& lastH = haarNode[0];
    for (int k = 0; k < numAttrs; k++)
      last.sumAttrInter[k] = lastH.attrInter[k];

    if (numNode2 == 2) {
      lastH = haarNode[1];
      for (int k = 0; k < numAttrs; k++) {
        auto temp = lastH.attrInter[k] - last.sumAttrInter[k];
        last.sumAttrInter[k] += temp >> 1;
      }
    }
  }
};

template<typename TreeNode, bool haarFlag, int numAttrs>
int
reduceDepthDecoder(
  int level,
  int numNodes,
  std::vector<TreeNode>* weightsIn,
  std::vector<TreeNode>* weightsOut,
  const bool isInter)
{
  // process a single level of the tree
  int64_t posPrev = -1;
  auto weightsInRdIt = weightsIn->begin();
  for (int i = 0; i < weightsIn->size(); ) {

    // this is a new node
    TreeNode last = weightsInRdIt[i];
    posPrev = last.pos;
    last.firstChildIdx = i;

    // look for same node
    int i2 = i + 1;
    for (; i2 < weightsIn->size(); i2++)
      if ((posPrev ^ weightsInRdIt[i2].pos) >> level)
        break;

    // process same nodes
    last.numChildren = i2 - i;

    reduceDepthDecoder_<numAttrs, TreeNode>::process(
      isInter, level, i, i2, weightsInRdIt, last);

    weightsOut->push_back(last);
    i = i2;
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
  std::vector<T>& weightsLf,
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

  memset(predBuf, 0, 8 * numAttrs * sizeof(*predBuf));

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

template<bool haarFlag, int numAttrs, typename Node>
inline void
uraht_process_decoder(
  const RahtPredictionParams& rahtPredParams,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  const int numPoints,
  int64_t* positions,
  attr_t* attributes,
  const attr_t* attributes_mc,
  int32_t* coeffBufIt,
  attr::ModeDecoder& coder)
{
  // coefficients are stored in three planar arrays.
  // coeffBufItK is a set of iterators to each array.
  int32_t* coeffBufItK[3] =
    {coeffBufIt, coeffBufIt + numPoints, coeffBufIt + numPoints * 2};

  // early termination only one point
  if (numPoints == 1) {
    if (haarFlag) {
      for (int k = 0; k < numAttrs; k++) {
        attributes[k] = *coeffBufItK[k]++;
      }
    } else {
      auto quantizers = qpset.quantizers(0, pointQpOffsets[0]);
      for (int k = 0; k < numAttrs; k++) {
        auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];
        int64_t coeff = *coeffBufItK[k]++;
        attributes[k] = divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
      }
    }
    return;
  }

  const bool isInter = attributes_mc;

  // --------- ascend tree per depth  -----------------
  // create leaf nodes
  int regionQpShift = 4;
  std::vector<Node> weightsHf;
  std::vector<std::vector<Node>> weightsLfStack;

  weightsLfStack.emplace_back();
  weightsLfStack.back().reserve(numPoints);
  auto weightsLfRef = &weightsLfStack.back();
  auto attrPredictor = attributes_mc;

  for (int i = 0; i < numPoints; i++) {
    Node node;
    node.pos = positions[i];
    if (!haarFlag) {
      node.weight = 1;
      node.qp = {
        int16_t(pointQpOffsets[i][0] << regionQpShift),
        int16_t(pointQpOffsets[i][1] << regionQpShift)};
    }
    if (isInter)
      for (int k = 0; k < numAttrs; k++)
        node.sumAttrInter[k] = *attrPredictor++;
    weightsLfRef->emplace_back(node);
  }

  // -----------  bottom up per depth  --------------------
  int numNodes = weightsLfRef->size();
#if 0 // duplicates won't work anymore with haar Nodes
      // Do we still need that with GeS-TM, if yes,
      // TODO: fix duplicates
  // for duplicates, skipable if it is known there is no duplicate
  numNodes =
    reduceUniqueDecoder<numAttrs>(isInter, numNodes, weightsLfRef, &weightsHf);
#endif

  const bool flagNoDuplicate = weightsHf.size() == 0;
  int numDepth = 0;
  for (int levelD = 3; numNodes > 1; levelD += 3) {
    // one depth reduction
    weightsLfStack.emplace_back();
    weightsLfStack.back().reserve(numNodes / 3);
    weightsLfRef = &weightsLfStack.back();

    auto weightsLfRefold = &weightsLfStack[weightsLfStack.size() - 2];
    numNodes = reduceDepthDecoder<Node, haarFlag, numAttrs>(
      levelD, numNodes, weightsLfRefold, weightsLfRef, isInter);
    numDepth++;
  }

  // --------- initialize stuff ----------------
  // root node
  auto& rootNode = weightsLfStack.back()[0];
  rootNode.mode = Mode::Null;
  rootNode.firstChildIdx = 0;
  if (!haarFlag)
    assert(rootNode.weight == numPoints);

  const int RDOCodingDepth = abh.attr_layer_code_mode.size();
  const bool enableACInterPred =
    rahtPredParams.enable_inter_prediction && isInter;

  coder.setInterEnabled(
    rahtPredParams.prediction_enabled_flag && enableACInterPred);

  const bool CbCrEnabled =
    !coder.isInterEnabled()
    && rahtPredParams.cross_chroma_component_prediction_flag;

  const bool CCRPEnabled =
    rahtPredParams.cross_component_residual_prediction_flag;

  const int maxlevelCCPenabled = (rahtPredParams.numlayer_CCRP_enabled - 1) * 3;

  // reconstruction buffers
  std::vector<int32_t> attrRec, attrRecParent;
  attrRec.resize(numPoints * numAttrs);
  attrRecParent.resize(numPoints * numAttrs);

  std::vector<int64_t> attrRecUs, attrRecParentUs;
  if (!haarFlag) {
    attrRecUs.resize(numPoints * numAttrs);
    attrRecParentUs.resize(numPoints * numAttrs);
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
  int64_t SampleDomainBuff[2 * numAttrs * 8];

  const size_t sizeofBuf = sizeof(int64_t) * numAttrs * 8 * (coder.isInterEnabled() ? 2 : 1);

  int64_t* attrPredIntra = SampleDomainBuff;
  int64_t* attrPredInter;
  int64_t* attrBestPred;

  // modes, inter
  std::vector<Mode> modes;
  modes.push_back(Mode::Null);
  if (coder.isInterEnabled()) {
    attrPredInter = attrPredIntra;
    modes.push_back(Mode::Inter);
    attrPredIntra += 8 * numAttrs;
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
    std::vector<Node>& weightsParent = weightsLfStack[levelD];
    std::vector<Node>& weightsLf = weightsLfStack[levelD - 1];

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
    if (!haarFlag)
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
    //CCRP parameters
    const bool CCRPFlag =
      !haarFlag && CCRPEnabled && level <= maxlevelCCPenabled && !CbCrEnabled;
    CCRPFilter ccrpFilter;

    // -------------- loop on nodes of the depth --------------
    int i = 0;
    auto attrRecParentUsIt = conditional_value<haarFlag>
      ::from(attrRecParent.begin(), attrRecParentUs.begin());
    for (auto weightsParentIt = weightsParent.begin();
        weightsParentIt < weightsParent.end();
        weightsParentIt++, attrRecParentUsIt += numAttrs) {

      if (modeParents != attr::Mode::size)
        // replaces parent mode for the whole depth;
        // done also in function getNeighborsModeDecoder
        weightsParentIt->mode = modeParents;

      typename std::conditional<haarFlag, bool, int64_t>::type
        weights[8 + 8 + 8 + 8 + 24] = {};
      int64_t interPredictor[8 * 3] = { 0 };
      int childTable[8] = { };
      bool skipinverse = !haarFlag;

      // generate weights, occupancy mask, and fwd transform
      // for all siblings of the current node.
      Qps nodeQp[8] = {};
      uint8_t occupancy = 0;
      const int nodeCnt = weightsParentIt->numChildren;
      for (int t = 0, j0 = i; t < nodeCnt; t++, j0++) {
        int nodeIdx = (weightsLf[j0].pos >> level) & 0x7;
        childTable[t] = nodeIdx;
        if (haarFlag) {
          weights[nodeIdx] = true;
        } else {
          weights[nodeIdx] = weightsLf[j0].weight;
          nodeQp[nodeIdx][0] = weightsLf[j0].qp[0] >> regionQpShift;
          nodeQp[nodeIdx][1] = weightsLf[j0].qp[1] >> regionQpShift;
        }

        // inter predictor
        if (isInter) {
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
      }
      weightsParentIt->occupancy = occupancy;
      int64_t sumweights = haarFlag ? 0 : weightsParentIt->weight;

      // Inter-level prediction:
      memset(SampleDomainBuff, 0, sizeofBuf);

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
      enableInterPred = enableInterPred || enableInterLayerPred; // inter prediction

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
      if (!rahtPredParams.enable_inter_prediction) {
        for (int t = 0; t < nodeCnt; t++)
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
                if (!enableIntraPred) {
                  neighborsMode = Mode::Inter;
                } else {
                  neighborsMode = getNeighborsModeDecoder(
                    parentNeighIdx, weightsParent, voteInterWeight,
                    voteIntraWeight, weightsLf, modeParents);
                }
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
      for (int t = 0; t < nodeCnt; t++, ChildPt++) {
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


      // ---------- prepare best predictor attrBestPred  -------------
      if (attr::isNull(predMode)) { // null mode
        skipinverse = false;
        attrBestPred = SampleDomainBuff;
        memset(attrBestPred, 0, numAttrs * 8 * sizeof(int64_t));
      }
      else if (attr::isIntra(predMode) || isAveragePred) { // intra or average mode
        if (rahtPredParams.subnode_prediction_enabled_flag)
          findNeighboursChildren(weightsParent.begin(), weightsLf.begin(), level, occupancy,
            parentNeighIdx, childNeighIdx, weightsLf);

        intraDcPredDecoder<haarFlag, numAttrs>(
          parentNeighIdx, childNeighIdx, occupancy, attrRecParent.begin(),
          attrRec.begin(), attrPredIntra, rahtPredParams);
        attrBestPred = attrPredIntra;
      }
      else { // inter mode
        attrBestPred = attrPredInter;
      }

      // prepare inter and average mode
      if (attr::isInter(predMode) || isAveragePred) {
        memcpy(attrPredInter, interPredictor, 8 * numAttrs * sizeof(int64_t));
      }


      // -------  denormalize best predictor   ----------
      int64_t PredDC[3] = {0, 0, 0};
      int64_t sqrtweightsbuf[8] = {};
      int64_t normsqrtsbuf[8] = {};
      using Kernel =
        typename std::conditional<haarFlag, HaarKernel, RahtKernel>::type;
      mkWeightTree<haarFlag>::template apply<true, Kernel>(weights);

      if (!attr::isNull(predMode)) { // only if not NULL
        if (haarFlag) {
          FwdTransformBlock222<haarFlag>
          ::template apply<numAttrs>(attrBestPred, weights);
          if (isAveragePred) {
            FwdTransformBlock222<haarFlag>
            ::template apply<numAttrs>(attrPredInter, weights);
          }
        }
        else {
          // normalise predicted attribute values
          for (int n = 0; n < nodeCnt; n++) {
            const int childIdx = childTable[n];
            int64_t w = weights[childIdx];
            if (w <= 1) {
              sqrtweightsbuf[childIdx] = 32768;
              continue;
            }

            int64_t sqrtWeight = fastIsqrt(uint64_t(w));
            sqrtweightsbuf[childIdx] = sqrtWeight;
            for (int k = 0; k < numAttrs; k++) {
              const int k8idx = 8 * k + childIdx;
              attrBestPred[k8idx] = fpReduce<kFPFracBits>(attrBestPred[k8idx] * sqrtWeight);
              if (isAveragePred)
                attrPredInter[k8idx] = fpReduce<kFPFracBits>(attrPredInter[k8idx] * sqrtWeight);
            }
          }
        }

        if (isAveragePred) {
          if (neighborsMode == Mode::size)
            neighborsMode = getNeighborsModeDecoder(
              parentNeighIdx, weightsParent, voteInterWeight, voteIntraWeight,
              weightsLf, modeParents);

          bool flag = !attr::isInter(weightsParentIt->mode) && !attr::isIntra(weightsParentIt->mode);
          voteInterWeight += 12 * attr::isInter(weightsParentIt->mode) + 6 * flag;
          voteIntraWeight += 12 * attr::isIntra(weightsParentIt->mode) + 6 * flag;

          if (depth > 1) {
            bool flag_grand = !attr::isInter(weightsParentIt->grand_mode) && !attr::isIntra(weightsParentIt->grand_mode);
            voteInterWeight += 28 * attr::isInter(weightsParentIt->grand_mode) + 14 * flag_grand;
            voteIntraWeight += 28 * attr::isIntra(weightsParentIt->grand_mode) + 14 * flag_grand;
          }

          int64_t weightIntra = divApprox(voteIntraWeight << kFPFracBits, voteInterWeight + voteIntraWeight, 0);
          int64_t weightInter = (1 << kFPFracBits) - weightIntra;

          for (int t = 0; t < 8 * numAttrs; t++) {
            attrBestPred[t] = fpReduce<kFPFracBits>(attrPredInter[t] * weightInter + attrBestPred[t] * weightIntra);
            if (haarFlag)
              attrBestPred[t] &= kFPIntMask;
          }
        }

        // Compute DC of best pred
        if (!haarFlag) {
          int64_t rsqrtweightsum = fastIrsqrt(sumweights);
          for (int n = 0; n < nodeCnt; n++) {
            int childIdx = childTable[n];
            int64_t normalizedsqrtweight = sqrtweightsbuf[childIdx] * rsqrtweightsum >> 40;
            normsqrtsbuf[childIdx] = normalizedsqrtweight;
            for (int k = 0; k < numAttrs; k++)
              PredDC[k] += fpReduce<kFPFracBits>(normalizedsqrtweight * attrBestPred[8 * k + childIdx]);

          }
        }
      }


      // per-coefficient operations:
      //  - read in quantised coefficients
      //  - inverse quantise + add transform domain prediction
      int64_t BestRecBuf[8 * numAttrs] = {0};
      int64_t CoeffRecBuf[8][numAttrs] = {0};
      int nodelvlSum = 0;
      int64_t transformRecBuf[numAttrs] = {0};
      bool enableCCCP =
        !haarFlag && rahtPredParams.cross_chroma_component_prediction_flag
        && !coder.isInterEnabled();

      CCRPFilter::Corr curCorr = { 0, 0, 0 };

      // ----------scan blocks -------
      for (int idxB = 0; idxB < 8; idxB++) {
        if ((idxB == 0 || weights[24 + idxB]) // there is always the DC coefficient (empty blocks are not transformed)
             && !(inheritDc && !idxB)) {  // skip the DC coefficient unless at the root of the tree
          // The RAHT transform
          auto quantizers = qpset.quantizers(qpLayer, nodeQp[idxB]);
          int64_t quantizedLuma = 0, quantizedChroma = 0;

          for (int k = 0; k < numAttrs; k++) {
            auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

            int64_t coeff = *coeffBufItK[k]++;
            skipinverse = skipinverse && (coeff == 0);
            CoeffRecBuf[nodelvlSum][k] = divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift); //div and scale do nothing for Haar
            int64_t quantizedValue = CoeffRecBuf[nodelvlSum][k];
            transformRecBuf[k] = fpExpand<kFPFracBits>(CoeffRecBuf[nodelvlSum][k]);

            const int k8idx = 8 * k + idxB;
            BestRecBuf[k8idx] = transformRecBuf[k];

            if (haarFlag) {
              // Transform Domain Pred for Lossless case
              BestRecBuf[k8idx] += attrBestPred[k8idx];
            }
            if (CCRPFlag) {
              if (k == 0) {
                quantizedLuma = quantizedValue;
                curCorr.yy += quantizedLuma * quantizedLuma;
              }
              if (k == 1) {
                int64_t CCRPPred = quantizedLuma * ccrpFilter.getYCbFilt() >> kCCRPFiltPrecisionbits;
                quantizedChroma = quantizedValue + CCRPPred;
                BestRecBuf[k8idx] = fpExpand<kFPFracBits>(quantizedChroma);
                skipinverse = skipinverse && (quantizedChroma == 0);
                curCorr.ycb += quantizedLuma * quantizedChroma;
              }
              if (k == 2) {
                int64_t CCRPPred = quantizedLuma * ccrpFilter.getYCrFilt() >> kCCRPFiltPrecisionbits;
                quantizedChroma = quantizedValue + CCRPPred;
                BestRecBuf[k8idx] = fpExpand<kFPFracBits>(quantizedChroma);
                skipinverse = skipinverse && (quantizedChroma == 0);
                curCorr.ycr += quantizedLuma * quantizedChroma;
              }
            }
          }

          // cross chroma
          if (enableCCCP && numAttrs == 3) {
            BestRecBuf[16 + idxB] += (CccpCoeff * transformRecBuf[1]) >> 4;
            CoeffRecBuf[nodelvlSum][2] = fpReduce<kFPFracBits>(BestRecBuf[16 + idxB]);
          }
          nodelvlSum++;
        }
      } // end scan

      if (CCRPFlag) {
        ccrpFilter.update(curCorr);
      }

      // compute last component coefficient
      if (numAttrs == 3 && nodeCnt > 1 && enableCCCP) {
        CccpCoeff = curlevelCccp.computeCrossChromaComponentPredictionCoeff(
          nodelvlSum, CoeffRecBuf);
      }
      // replace DC coefficient with parent if inheritable
      if (inheritDc) {
        for (int k = 0; k < numAttrs; k++) {
          if (haarFlag)
            BestRecBuf[8 * k] = attrRecParentUsIt[k];
          else
            BestRecBuf[8 * k] = attrRecParentUsIt[k] - PredDC[k];
        }
      }

      if (haarFlag) {
        InvTransformBlock222<haarFlag>
        ::template apply<numAttrs>(BestRecBuf, weights);
      }
      else {
        if (skipinverse) {
          int64_t DCerror[numAttrs];
          for (int k = 0; k < numAttrs; k++) {
            DCerror[k] = BestRecBuf[8 * k];
          }
          for (int n = 0; n < nodeCnt; n++) {
            int childIdx = childTable[n];
            for (int k = 0; k < numAttrs; k++)
              BestRecBuf[8 * k + childIdx] = fpReduce<kFPFracBits>(normsqrtsbuf[childIdx] * DCerror[k]);
          }
        }
        else {
          InvTransformBlock222<haarFlag>
          ::template apply<numAttrs, true>(BestRecBuf, weights);
        }
      }

      for (int j = i, n = 0; n < nodeCnt; j++, n++) {
        const int childIdx = childTable[n];

        // scale values for next level
        if (!haarFlag) {
          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + childIdx;
            //Sample Domain Reconstruction
            BestRecBuf[k8idx] += attrBestPred[k8idx];
            attrRecUs[j * numAttrs + k] = BestRecBuf[k8idx];
          }

          uint64_t w = weights[childIdx];
          if (w> 1) {
            int shift = 5 * ((w > 1024) + (w > 1048576));
            uint64_t rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;
            for (int k = 0; k < numAttrs; k++) {
              const int k8idx = 8 * k + childIdx;
              BestRecBuf[k8idx] = fpReduce<kFPFracBits>((BestRecBuf[k8idx] >> shift) * rsqrtWeight);
            }
          }
        }

        for (int k = 0; k < numAttrs; k++) {
          const int k8idx = 8 * k + childIdx;
          attrRec[j * numAttrs + k] = BestRecBuf[k8idx];
        }
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
    }
    else
      modeParents = attr::Mode::size;

    ++depth;
    weightsParent.clear();
  } // end loop on depth


  // -------------- process duplicate points at level 0 --------------
  if (flagNoDuplicate) { // write-back reconstructed attributes
    auto attrOut = attributes;
    for (auto attr : attrRec) {
      auto v = attr + kFPOneHalf >> kFPFracBits;
      *attrOut++ = PCCClip(v, 0, std::numeric_limits<attr_t>::max());
    }
    return;
  }

#if 0 // This code won't work with haar Nodes
      // Do we still need that with GeS-TM, if yes,
      // TODO: fix duplicates
      // One way to fix it would be to keep the duplicates count in a sevarate
      // buffer instead of using weight (which is not required for the rest of
      // haar)
  // case there are duplicates
  std::swap(attrRec, attrRecParent);
  auto attrRecParentIt = attrRecParent.cbegin();

  std::vector<int64_t> attrsHf;
  attrsHf.resize(weightsHf.size()* numAttrs);
  auto attrsHfIt = attrsHf.cbegin();

  std::vector<Node>& weightsLf = weightsLfStack[0];
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

    if (haarFlag) {
      for (int w = weight - 1; w > 0; w--) {
        HaarKernel haarkernel(w, 1);

        for (int k = 0; k < numAttrs; k++) {

          int64_t transformBuf[2];
          transformBuf[1] = *coeffBufItK[k]++;
          transformBuf[0] = attrRecDc[k]; // inherit the DC value

          haarkernel.invTransform(transformBuf[0], transformBuf[1]);

          attrRecDc[k] = transformBuf[0];
          attrRec[out + w * numAttrs + k] = transformBuf[1];
          if (w == 1)
            attrRec[out + k] = transformBuf[0];
        }
      }
    } else {
      for (int w = weight - 1; w > 0; w--) {
        RahtKernel kernel(w, 1);

        auto quantizers = qpset.quantizers(qpLayer, nodeQp);
        for (int k = 0; k < numAttrs; k++) {
          auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

          int64_t transformBuf[2];
          int64_t coeff = *coeffBufItK[k]++;
          transformBuf[1] = fpExpand<kFPFracBits>(divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift));
          transformBuf[0] = attrRecDc[k]; // inherit the DC value

          kernel.invTransform(transformBuf[0], transformBuf[1]);

          attrRecDc[k] = transformBuf[0];
          attrRec[out + w * numAttrs + k] = transformBuf[1];
          if (w == 1)
            attrRec[out + k] = transformBuf[0];
        }
      }

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
#endif
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
  attr_t* attributes,
  const attr_t* attributes_mc,
  int* coefficients,
  attr::ModeDecoder& decoder)
{
  switch (attribCount) {
  case 3:
    if (!rahtPredParams.integer_haar_enable_flag)
      uraht_process_decoder<false, 3, UrahtNodeDecoder>(
        rahtPredParams, abh,qpset, pointQpOffsets, voxelCount, mortonCode,
        attributes, attributes_mc, coefficients,
        decoder);
    else
      uraht_process_decoder<true, 3, UrahtNodeDecoderHaar>(
        rahtPredParams, abh, qpset, pointQpOffsets, voxelCount, mortonCode,
        attributes, attributes_mc, coefficients,
        decoder);
    break;
  default:
    throw std::runtime_error("attribCount != 3 not tested yet");
  }
}

//============================================================================

}  // namespace pcc
