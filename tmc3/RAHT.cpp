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

namespace pcc {

//============================================================================
// remove any non-unique leaves from a level in the uraht tree

int
reduceUnique(
  int numNodes,
  int numAttrs,
  std::vector<UrahtNode>* weightsIn,
  std::vector<UrahtNode>* weightsOut,
  std::vector<int>* attrsIn,
  std::vector<int>* attrsOut)
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
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;
      continue;
    }

    // duplicate node
    (weightsInWrIt - 1)->weight += node.weight;
    weightsOut->push_back(node);
    for (int k = 0; k < numAttrs; k++) {
      *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;
      attrsOut->push_back(*attrsInRdIt++);
    }
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================
// Split a level of values into sum and difference pairs.

int
reduceLevel(
  int level,
  int numNodes,
  int numAttrs,
  std::vector<UrahtNode>* weightsIn,
  std::vector<UrahtNode>* weightsOut,
  std::vector<int>* attrsIn,
  std::vector<int>* attrsOut)
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
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;
    } else {
      auto& left = *(weightsInWrIt - 1);
      left.weight += node.weight;
      left.qp[0] = (left.qp[0] + node.qp[0]) >> 1;
      left.qp[1] = (left.qp[1] + node.qp[1]) >> 1;
      weightsOut->push_back(node);

      for (int k = 0; k < numAttrs; k++) {
        *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;
        attrsOut->push_back(*attrsInRdIt++);
      }
    }
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================
// Merge sum and difference values to form a tree.

void
expandLevel(
  int level,
  int numNodes,
  int numAttrs,
  std::vector<UrahtNode>* weightsIn,   // expand by numNodes before expand
  std::vector<UrahtNode>* weightsOut,  // shrink after expand
  std::vector<int>* attrsIn,
  std::vector<int>* attrsOut)
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
      *attrsInWrIt++ -= *curAttrIt++;
    }
  }
}

//============================================================================
// Search for neighbour with @value in the ordered list [first, last).
//
// If distance is positive, search [from, from+distance].
// If distance is negative, search [from-distance, from].

template<typename It, typename T, typename T2, typename Cmp>
It
findNeighbour(It first, It last, It from, T value, T2 distance, Cmp compare)
{
  It start = first;
  It end = last;

  if (distance >= 0) {
    start = from;
    if ((distance + 1) < std::distance(from, last))
      end = std::next(from, distance + 1);
  } else {
    end = from;
    if ((-distance) < std::distance(first, from))
      start = std::prev(from, -distance);
  }

  auto found = std::lower_bound(start, end, value, compare);
  if (found == end)
    return last;
  return found;
}

//============================================================================
// Find the neighbours of the node indicated by @t between @first and @last.
// The position weight of each found neighbour is stored in two arrays.


#define MASK_Self 255
#define MASK_PX 240
#define MASK_NX 15
#define MASK_PY 204
#define MASK_NY 51
#define MASK_PZ 170
#define MASK_NZ 85
#define MASK_PXPY 192
#define MASK_PXNY 48
#define MASK_NXPY 12
#define MASK_NXNY 3
#define MASK_PYPZ 136
#define MASK_PYNZ 68
#define MASK_NYPZ 34
#define MASK_NYNZ 17
#define MASK_PZPX 160
#define MASK_PZNX 10
#define MASK_NZPX 80
#define MASK_NZNX 5
#define MASK_Other 0
#define MASKS \
  { \
    MASK_Self, \
    MASK_PX,   MASK_PY,   MASK_PZ, \
    MASK_NX,   MASK_NY,   MASK_NZ, \
    MASK_PXPY, MASK_PZPX, MASK_PYPZ, \
    MASK_NXNY, MASK_NZNX, MASK_NYNZ, \
    MASK_PZNX, MASK_NYPZ, MASK_NXPY, \
    MASK_PYNZ, MASK_PXNY, MASK_NZPX \
  };

template<typename It>
void
findNeighbours(
  It first,
  It last,
  It it,
  It firstChild,
  It lastChild,
  int level,
  uint8_t occupancy,
  std::array<ParentIndex, 19>& parentNeighIdx,
  bool subnodePredictionEnabled)
{
  static const uint8_t neighMasks[19] = MASKS;

  // current position (discard extra precision)
  int64_t cur_pos = it->pos >> level;

  // the position of the parent, offset by (-1,-1,-1)
  int64_t base_pos = morton3dAdd(cur_pos, -1ll);

  // these neighbour offsets are relative to base_pos
  static const uint8_t neighOffset[19] = {
    0, 35, 21, 14, 3, 5, 6, 49, 42, 28, 1, 2, 4, 10, 12, 17, 20, 33, 34};

  static const NeighborType type[19] = {
    Self, Face, Face, Face, Face, Face, Face, Edge, Edge, Edge,
    Edge, Edge, Edge, Edge, Edge, Edge, Edge, Edge, Edge};

  for (auto& neighbor : parentNeighIdx)
    neighbor.reset();

  // special case for the direct parent (no need to search);
  parentNeighIdx[0].set(std::distance(first, it), NeighborType::Self);

  for (int i = 1; i < 19; i++) {
    // Only look for neighbours that have an effect
    if (!(occupancy & neighMasks[i])) {
      continue;
    }

    // compute neighbour address to look for
    // the delta between it and the current position is
    int64_t neigh_pos = morton3dAdd(base_pos, neighOffset[i]);
    int64_t delta = neigh_pos - cur_pos;

    // find neighbour
    auto found = findNeighbour(
      first, last, it, neigh_pos, delta,
      [=](decltype(*it)& candidate, int64_t neigh_pos) {
        return (candidate.pos >> level) < neigh_pos;
      });

    if (found == last) {
      continue;
    }

    if ((found->pos >> level) != neigh_pos) {
      continue;
    }

    parentNeighIdx[i].set(std::distance(first, found), type[i]);
  }

  if (subnodePredictionEnabled) {
    static const int idxFace[6][4][4] = {
      {{0, 4, 6, 5}, {1, 5, 4, 7}, {2, 6, 4, 7}, {3, 7, 5, 6}},
      {{0, 2, 3, 6}, {1, 3, 2, 7}, {4, 6, 2, 7}, {5, 7, 3, 6}},
      {{0, 1, 3, 5}, {2, 3, 1, 7}, {4, 5, 1, 7}, {6, 7, 3, 5}},
      {{4, 0, 1, 2}, {5, 1, 0, 3}, {6, 2, 0, 3}, {7, 3, 1, 2}},
      {{2, 0, 1, 4}, {3, 1, 0, 5}, {6, 4, 0, 5}, {7, 5, 1, 4}},
      {{1, 0, 2, 4}, {3, 2, 0, 6}, {5, 4, 0, 6}, {7, 6, 2, 4}}};
    static const int idxEdge[12][2][2] = {
      {{0, 6}, {1, 7}}, {{0, 5}, {2, 7}}, {{0, 3}, {4, 7}}, {{6, 0}, {7, 1}},
      {{5, 0}, {7, 2}}, {{3, 0}, {7, 4}}, {{4, 1}, {6, 3}}, {{2, 1}, {6, 5}},
      {{4, 2}, {5, 3}}, {{1, 2}, {5, 6}}, {{2, 4}, {3, 5}}, {{1, 4}, {3, 6}}};

    int curLevel = level - 3;

    for (int i = 1; i < 7; i++) {
      auto& neigh = parentNeighIdx[i];
      if (neigh.isvoid())
        continue;

      auto neiIt = first + neigh;

      if (!neiIt->decoded)
        continue;

      int j = 0;
      auto IDX = idxFace[i - 1];
      auto idx = IDX[0];
      auto it = neiIt->firstChild;
      while (it != neiIt->lastChild) {
        int nodeIdx = (it->pos >> curLevel) & 0x7;

        if (nodeIdx > idx[0]) {
          j++;
          if (j < 4)
            idx = IDX[j];
          else
            break;
          continue;
        }
        if (nodeIdx < idx[0]) {
          it++;
          continue;
        }

        auto distance = std::distance(firstChild, it);
        neigh.set(idx[1], distance, NeighborType::Face);
        if (neigh.isvoid(idx[2]))
          neigh.set(idx[2], distance, NeighborType::Edge);
        if (neigh.isvoid(idx[3]))
          neigh.set(idx[3], distance, NeighborType::Edge);

        it++;
        j++;
        if (j < 4)
          idx = IDX[j];
        else
          break;
      }
    }

    for (int i = 7; i < 19; i++) {
      auto& neigh = parentNeighIdx[i];
      if (neigh.isvoid())
        continue;

      auto neiIt = first + neigh;

      if (!neiIt->decoded)
        continue;

      int j = 0;
      auto IDX = idxEdge[i - 7];
      auto idx = IDX[0];
      auto it = neiIt->firstChild;
      while (it != neiIt->lastChild) {
        int nodeIdx = (it->pos >> curLevel) & 0x7;

        if (nodeIdx > idx[0]) {
          j++;
          if (j < 2)
            idx = IDX[j];
          else
            break;
          continue;
        }
        if (nodeIdx < idx[0]) {
          it++;
          continue;
        }

        auto distance = std::distance(firstChild, it);
        neigh.set(idx[1], distance, NeighborType::Edge);

        it++;
        j++;
        if (j < 2)
          idx = IDX[j];
        else
          break;
      }
    }
  }
}

//============================================================================
// Generate the spatial prediction of a block.

template<typename It>
void
intraDcPred(
  int numAttrs,
  const std::array<ParentIndex, 19>& parentNeighIdx,
  int occupancy,
  It first,
  It firstChild,
  VecAttr::iterator predBuf,
  std::vector<UrahtNode>::iterator firstNode,
  const RahtPredictionParams& rahtPredParams)
{
  static const uint8_t predMasks[19] = MASKS;

  static const int kDivisors[64] = {
    32768, 16384, 10923, 8192, 6554, 5461, 4681, 4096, 3641, 3277, 2979,
    2731,  2521,  2341,  2185, 2048, 1928, 1820, 1725, 1638, 1560, 1489,
    1425,  1365,  1311,  1260, 1214, 1170, 1130, 1092, 1057, 1024, 993,
    964,   936,   910,   886,  862,  840,  819,  799,  780,  762,  745,
    728,   712,   697,   683,  669,  655,  643,  630,  618,  607,  596,
    585,   575,   565,   555,  546,  537,  529,  520,  512};

  int weightSum[8] = {-1, -1, -1, -1, -1, -1, -1, -1};

  int64_t neighValue[3];
  int64_t childNeighValue[3];
  int64_t limitLow = 0;
  int64_t limitHigh = 0;

  for (int i = 0; i < 19; i++) {
    auto& neigbor = parentNeighIdx[i];
    int predWeightParent = neigbor.getweight(rahtPredParams.prediction_weights);
    if (!predWeightParent)
      continue;

    auto neighValueIt = std::next(first, numAttrs * neigbor);
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

    auto node = std::next(firstNode, neigbor);

    // apply weighted neighbour value to masked positions
    for (int k = 0; k < numAttrs; k++)
      neighValue[k] *= predWeightParent;

    int mask = predMasks[i] & occupancy;
    for (int j = 0; mask; j++, mask >>= 1) {
      if (mask & 1) {
        int predWeightChild = 0;
        if (rahtPredParams.subnode_prediction_enabled_flag || node->decoded)
          predWeightChild =
            neigbor.getweight(j, rahtPredParams.prediction_weights);
        if (predWeightChild) {
          weightSum[j] += predWeightChild;
          auto childNeighValueIt =
            std::next(firstChild, numAttrs * parentNeighIdx[i][j]);
          for (int k = 0; k < numAttrs; k++)
            childNeighValue[k] = (*childNeighValueIt++) * predWeightChild;

          for (int k = 0; k < numAttrs; k++)
            predBuf[k][j].val += childNeighValue[k];
        } else {
          weightSum[j] += predWeightParent;
          for (int k = 0; k < numAttrs; k++)
            predBuf[k][j].val += neighValue[k];
        }
      }
    }
  }

  // normalise
  FixedPoint div;
  for (int i = 0; i < 8; i++, occupancy >>= 1) {
    if (occupancy & 1) {
      div.val = kDivisors[weightSum[i]];
      for (int k = 0; k < numAttrs; k++)
        predBuf[k][i] *= div;
    }
  }
}

//============================================================================
// expand a set of eight weights into three levels

template<class Kernel>
void
mkWeightTree(int64_t weights[8 + 8 + 8 + 8 + 24])
{
  auto in = &weights[0];
  auto out = &weights[8];

  for (int i = 0; i < 4; i++) {
    out[0] = out[4] = in[0] + in[1];
    if (!in[0] || !in[1])
      out[4] = 0;  // single node, no high frequencies
    in += 2;
    out++;
  }
  out += 4;
  for (int i = 0; i < 4; i++) {
    out[0] = out[4] = in[0] + in[1];
    if (!in[0] || !in[1])
      out[4] = 0;  // single node, no high frequencies
    in += 2;
    out++;
  }
  out += 4;
  for (int i = 0; i < 4; i++) {
    out[0] = out[4] = in[0] + in[1];
    if (!in[0] || !in[1])
      out[4] = 0;  // single node, no high frequencies
    in += 2;
    out++;
  }

  for (int i = 0; i < 24; i += 2) {
    if (!weights[i]) {
      if (weights[i + 1])
        weights[i + 33] = 0x1ll << FixedPoint::kFracBits;
    } else if (!weights[i + 1]) {
      weights[i + 32] = 0x1ll << FixedPoint::kFracBits;
    }
    else {
      Kernel w(weights[i], weights[i + 1]);
      weights[i + 32] = w.getW0();
      weights[i + 33] = w.getW1();
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
    if (!weights[24 + kRahtScanOrder[i]])
      continue;

    mapFn(kRahtScanOrder[i]);
  }
}

//============================================================================
// Tests if two positions are siblings at the given tree level

static bool
isSibling(int64_t pos0, int64_t pos1, int level)
{
  return ((pos0 ^ pos1) >> level) == 0;
}

//============================================================================
// Core transform process (for encoder/decoder)

template<bool isEncoder>
void
uraht_process(
  const RahtPredictionParams& rahtPredParams,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int numPoints,
  int numAttrs,
  int64_t* positions,
  int* attributes,
  int32_t* coeffBufIt)
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

      if (isEncoder) {
        auto coeff = attributes[k];
        assert(coeff <= INT_MAX && coeff >= INT_MIN);
        *coeffBufItK[k]++ = coeff =
          q.quantize(coeff << kFixedPointAttributeShift);
        attributes[k] =
          divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
      } else {
        int64_t coeff = *coeffBufItK[k]++;
        attributes[k] =
          divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
      }
    }
    return;
  }

  std::vector<UrahtNode> weightsLf, weightsHf;
  std::vector<int> attrsLf, attrsHf;

  weightsLf.reserve(numPoints);
  attrsLf.reserve(numPoints * numAttrs);

  int regionQpShift = 4;

  // copy positions into internal form
  // todo(df): lift to api
  for (int i = 0; i < numPoints; i++) {
    weightsLf.emplace_back(UrahtNode{
      positions[i],
      1,
      {pointQpOffsets[i][0] << regionQpShift,
       pointQpOffsets[i][1] << regionQpShift}});
    for (int k = 0; k < numAttrs; k++) {
      attrsLf.push_back(attributes[i * numAttrs + k]);
    }
  }

  weightsHf.reserve(numPoints);
  attrsHf.reserve(numPoints * numAttrs);

  // ascend tree
  std::vector<int> levelHfPos;

  for (int level = 0, numNodes = weightsLf.size(); numNodes > 1; level++) {
    levelHfPos.push_back(weightsHf.size());
    if (level == 0) {
      // process any duplicate points
      numNodes = reduceUnique(
        numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    } else {
      // normal level reduction
      numNodes = reduceLevel(
        level, numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    }
  }

  assert(weightsLf[0].weight == numPoints);
  weightsLf[0].mode = Mode::Null;
  weightsLf[0].offset = 0;

  // reconstruction buffers
  std::vector<int64_t> attrRec, attrRecParent;
  attrRec.resize(numPoints * numAttrs);
  attrRecParent.resize(numPoints * numAttrs);

  std::vector<int64_t> attrRecUs, attrRecParentUs;
  attrRecUs.resize(numPoints * numAttrs);
  attrRecParentUs.resize(numPoints * numAttrs);

  std::vector<UrahtNode> weightsParent;
  weightsParent.resize(1, weightsLf[0]);

  std::vector<int> numParentNeigh, numGrandParentNeigh;
  numParentNeigh.resize(numPoints);
  numGrandParentNeigh.resize(numPoints);

  // indexes of the neighbouring parents
  static std::array<ParentIndex, 19> parentNeighIdx;

  // Prediction buffers
  VecAttr transformBuf;
  VecAttr attrRecIntra;
  VecAttr::iterator attrPred;
  VecAttr::iterator attrReal;
  VecAttr::iterator attrPredIntra;
  std::vector<FixedPoint> parentDc(numAttrs);

  transformBuf.resize(2 * numAttrs);

  attrReal = transformBuf.begin();
  attrPred = std::next(attrReal, numAttrs);

  attrPredIntra = attrPred;

  // quant layer selection
  auto qpLayer = 0;

  // descend tree
  weightsLf.resize(1);
  attrsLf.resize(numAttrs);
  int trainZeros = 0;
  for (int level = levelHfPos.size() - 1, isFirst = 1; level > 0; /*nop*/) {
    int numNodes = weightsHf.size() - levelHfPos[level];
    weightsLf.resize(weightsLf.size() + numNodes);
    attrsLf.resize(attrsLf.size() + numNodes * numAttrs);
    expandLevel(
      level, numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    weightsHf.resize(levelHfPos[level]);
    attrsHf.resize(levelHfPos[level] * numAttrs);

    // expansion of level is complete, processing is now on the next level
    level--;

    // every three levels, perform transform
    if (level % 3)
      continue;

    // initial scan position of the coefficient buffer
    //  -> first level = all coeffs
    //  -> otherwise = ac coeffs only
    bool inheritDc = !isFirst;
    bool enableIntraPredictionInLvl =
      inheritDc && rahtPredParams.prediction_enabled_flag;

    if (enableIntraPredictionInLvl) {
      for (auto& ele : weightsParent)
        ele.occupancy = 0;

      const int parentCount = weightsParent.size();
      auto it = weightsLf.begin();
      for (auto i = 0; i < parentCount; i++) {
        weightsParent[i].decoded = 0;
        weightsParent[i].firstChild = it++;

        while (it != weightsLf.end()
               && !((it->pos ^ weightsParent[i].pos) >> (level + 3)))
          it++;
        weightsParent[i].lastChild = it;
        if (isFirst)
          weightsParent[i].mode = Mode::Null;
      }
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

    int i = 0;
    for (auto weightsParentIt = weightsParent.begin();
         weightsParentIt < weightsParent.end(); weightsParentIt++) {

      for (auto& buf : transformBuf) {
        std::fill(buf.begin(), buf.end(), FixedPoint(0));
      }

      int64_t weights[8 + 8 + 8 + 8 + 24] = {};
      Qps nodeQp[8] = {};
      uint8_t occupancy = 0;

      // generate weights, occupancy mask, and fwd transform buffers
      // for all siblings of the current node.
      int nodeCnt = 0;
      while (!isSibling(weightsLf[i].pos, weightsParentIt->pos, level + 3))
        i++;
      for (int j = i; j < weightsLf.size(); j++) {
        int nextNode =
          j > i && !isSibling(weightsLf[j].pos, weightsLf[i].pos, level + 3);
        if (nextNode)
          break;

        int nodeIdx = (weightsLf[j].pos >> level) & 0x7;
        weights[nodeIdx] = weightsLf[j].weight;
        nodeQp[nodeIdx][0] = weightsLf[j].qp[0] >> regionQpShift;
        nodeQp[nodeIdx][1] = weightsLf[j].qp[1] >> regionQpShift;

        occupancy |= 1 << nodeIdx;

        if (rahtPredParams.prediction_skip1_flag)
          nodeCnt++;

        if (isEncoder) {
          for (int k = 0; k < numAttrs; k++)
            attrReal[k][nodeIdx] = attrsLf[j * numAttrs + k];
        }
      }

      mkWeightTree<RahtKernel>(weights);

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
      bool enableIntraPrediction = enableIntraPredictionInLvl;

      if (enableIntraPrediction) {
        weightsParentIt->occupancy = occupancy;

        int parentNeighCount = 0;
        if (rahtPredParams.prediction_skip1_flag && nodeCnt == 1) {
          enableIntraPrediction = false;
          parentNeighCount = 19;
        } else if (
          *numGrandParentNeighIt < rahtPredParams.prediction_threshold0) {
          enableIntraPrediction = false;
        } else {
          findNeighbours(
            weightsParent.begin(), weightsParent.end(), weightsParentIt,
            weightsLf.begin(), weightsLf.begin() + i, level + 3, occupancy,
            parentNeighIdx, rahtPredParams.subnode_prediction_enabled_flag);

          parentNeighCount = std::count_if(
            parentNeighIdx.begin(), parentNeighIdx.end(),
            [](const int idx) { return idx >= 0; });
          if (parentNeighCount < rahtPredParams.prediction_threshold1) {
            enableIntraPrediction = false;
          } else if (isEncoder) {
            intraDcPred(
              numAttrs, parentNeighIdx, occupancy, attrRecParent.begin(),
              attrRec.begin(), attrPredIntra, weightsParent.begin(),
              rahtPredParams);
          }
        }
        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])
            continue;
          numParentNeigh[j++] = parentNeighCount;
        }
      }

      int parentWeight = 0;
      if (inheritDc) {
        numGrandParentNeighIt++;
        parentWeight = weightsParentIt->weight;
        weightsParentIt->decoded = true;
      }

      if (isEncoder) {
        // normalise coefficients
        for (int childIdx = 0; childIdx < 8; childIdx++) {
          if (weights[childIdx] <= 1)
            continue;

          // Summed attribute values
          FixedPoint rsqrtWeight;
          uint64_t w = weights[childIdx];
          int shift = w > 1024 ? ilog2(w - 1) >> 1 : 0;
          rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
          for (int k = 0; k < numAttrs; k++) {
            transformBuf[k][childIdx].val >>= shift;
            transformBuf[k][childIdx] *= rsqrtWeight;
          }

          // Predicted attribute values
          FixedPoint sqrtWeight;
          sqrtWeight.val =
            isqrt(uint64_t(weights[childIdx]) << (2 * FixedPoint::kFracBits));
          attrPred = transformBuf.begin() + numAttrs;
          while (attrPred < transformBuf.end()) {
            for (int k = 0; k < numAttrs; k++)
              attrPred[k][childIdx] *= sqrtWeight;
            attrPred += numAttrs;
          }
        }

        fwdTransformBlock222(
          transformBuf.size(), transformBuf.begin(), weights);
      }

      Mode predMode = nodeCnt > 1 && enableIntraPrediction ? Mode::Intra : Mode::Null;

      for (auto weightsChild = weightsParentIt->firstChild;
           weightsChild < weightsParentIt->lastChild; weightsChild++) {
        weightsChild->mode = predMode;
      }

      if (isEncoder) {
        if (attr::isNull(predMode)) {
          attrPred = std::next(transformBuf.begin(), numAttrs);
          for (int k = 0; k < numAttrs; k++)
            std::fill(attrPred[k].begin(), attrPred[k].end(), FixedPoint(0));
        } else {
          attrPred = attrPredIntra;
        }
      }

      if (!isEncoder) {
        // prediction
        if (attr::isNull(predMode)) {
          attrPred = std::next(transformBuf.begin(), numAttrs);
          for (int k = 0; k < numAttrs; k++)
            std::fill(attrPred[k].begin(), attrPred[k].end(), FixedPoint(0));
        } else {
          intraDcPred(
            numAttrs, parentNeighIdx, occupancy, attrRecParent.begin(),
            attrRec.begin(), attrPredIntra, weightsParent.begin(),
            rahtPredParams);
          attrPred = attrPredIntra;
        }

        if (!attr::isNull(predMode)) {
          // normalise predicted attribute values
          for (int childIdx = 0; childIdx < 8; childIdx++) {
            if (weights[childIdx] <= 1)
              continue;

            FixedPoint sqrtWeight;
            sqrtWeight.val = isqrt(
              uint64_t(weights[childIdx]) << (2 * FixedPoint::kFracBits));
            for (int k = 0; k < numAttrs; k++)
              attrPred[k][childIdx] *= sqrtWeight;
          }

          fwdTransformBlock222(numAttrs, attrPred, weights);
        }
      }

      // per-coefficient operations:
      //  - subtract transform domain prediction (encoder)
      //  - write out/read in quantised coefficients
      //  - inverse quantise + add transform domain prediction
      scanBlock(weights, [&](int idx) {
        // skip the DC coefficient unless at the root of the tree
        if (inheritDc && !idx)
          return;

        // subtract transformed prediction (skipping DC)
        if (isEncoder && !attr::isNull(predMode)) {
          for (int k = 0; k < numAttrs; k++) {
            transformBuf[k][idx] -= attrPred[k][idx];
          }
        }

        // decision for RDOQ
        int64_t sumCoeff = 0;
        bool flagRDOQ = false;
        if (isEncoder) {
          int64_t Dist2 = 0;
          int Ratecoeff = 0;
          int64_t lambda0;

          for (int k = 0; k < numAttrs; k++) {
            //auto q = Quantizer(qpLayer[std::min(k, int(quantizers.size()) - 1)] + nodeQp[idx]);
            auto quantizers = qpset.quantizers(qpLayer, nodeQp[idx]);
            auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];
            auto coeff = transformBuf[k][idx].round();
            Dist2 += coeff * coeff;
            auto Qcoeff = q.quantize(coeff << kFixedPointAttributeShift);
            sumCoeff += std::abs(Qcoeff);
            //Ratecoeff += !!Qcoeff; // sign
            Ratecoeff += int(std::log2(1 + std::abs(Qcoeff)) * 256.);
            if (!k)
              lambda0 = q.scale(1);
          }

          if (sumCoeff < 3) {
            int LUTbins[11] = { 1,2,3, 5,5, 7,7, 9,9 ,11 ,11 };
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
            //Rate = Rate / std::max(1, trainZeros);
            Rate += (Ratecoeff + 128) >> 8;

            int64_t lambda = lambda0 * lambda0 * (numAttrs == 1 ? 25 : 35);
            flagRDOQ = (Dist2 << 26) < lambda * Rate;
          }

        }

        // Track RL for RDOQ
        if (flagRDOQ || sumCoeff == 0)
          trainZeros++;
        else
          trainZeros = 0;


        // The RAHT transform
        auto quantizers = qpset.quantizers(qpLayer, nodeQp[idx]);
        for (int k = 0; k < numAttrs; k++) {
          if (flagRDOQ) // apply RDOQ
            transformBuf[k][idx].val = 0;

          auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

          if (isEncoder) {
            auto coeff = transformBuf[k][idx].round();
            assert(coeff <= INT_MAX && coeff >= INT_MIN);
            *coeffBufItK[k]++ = coeff =
              q.quantize(coeff << kFixedPointAttributeShift);
            attrPred[k][idx] +=
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
          } else {
            int64_t coeff = *coeffBufItK[k]++;
            attrPred[k][idx] +=
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
          }
        }
      });

      // replace DC coefficient with parent if inheritable
      if (inheritDc) {
        for (int k = 0; k < numAttrs; k++) {
          attrPred[k][0].val = attrRecParentUsIt[k];
        }
      }

      invTransformBlock222(numAttrs, attrPred, weights);

      for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
        if (!weights[nodeIdx])
          continue;

        for (int k = 0; k < numAttrs; k++)
          attrRecUs[j * numAttrs + k] = attrPred[k][nodeIdx].val;

        // scale values for next level
        if (weights[nodeIdx] > 1) {
          FixedPoint rsqrtWeight;
          uint64_t w = weights[nodeIdx];
          int shift = w > 1024 ? ilog2(w - 1) >> 1 : 0;
          rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
          for (int k = 0; k < numAttrs; k++) {
            attrPred[k][nodeIdx].val >>= shift;
            attrPred[k][nodeIdx] *= rsqrtWeight;
          }
        }

        for (int k = 0; k < numAttrs; k++) {
          attrRec[j * numAttrs + k] = attrPred[k][nodeIdx].val;
        }
        j++;
      }
    }

    // preserve current weights/positions for later search
    weightsParent = weightsLf;
  }

  // process duplicate points at level 0
  std::swap(attrRec, attrRecParent);
  auto attrRecParentIt = attrRecParent.cbegin();
  auto attrsHfIt = attrsHf.cbegin();

  for (int i = 0, out = 0, iEnd = weightsLf.size(); i < iEnd; i++) {
    int weight = weightsLf[i].weight;
    Qps nodeQp = {
      weightsLf[i].qp[0] >> regionQpShift,
      weightsLf[i].qp[1] >> regionQpShift};

    // unique points have weight = 1
    if (weight == 1) {
      for (int k = 0; k < numAttrs; k++)
        attrRec[out++] = *attrRecParentIt++;
      continue;
    }

    // duplicates
    FixedPoint attrSum[3];
    FixedPoint attrRecDc[3];
    FixedPoint sqrtWeight;
    sqrtWeight.val = isqrt(uint64_t(weight) << (2 * FixedPoint::kFracBits));

    int64_t sumCoeff = 0;
    for (int k = 0; k < numAttrs; k++) {
      if (isEncoder)
        attrSum[k] = attrsLf[i * numAttrs + k];
      attrRecDc[k].val = *attrRecParentIt++;
      attrRecDc[k] *= sqrtWeight;
    }

    FixedPoint rsqrtWeight;
    for (int w = weight - 1; w > 0; w--) {
      RahtKernel kernel(w, 1);
      int shift = w > 1024 ? ilog2(uint32_t(w - 1)) >> 1 : 0;
      if (isEncoder)
        rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);

      auto quantizers = qpset.quantizers(qpLayer, nodeQp);
      for (int k = 0; k < numAttrs; k++) {
        auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

        FixedPoint transformBuf[2];
        if (isEncoder) {
          // invert the initial reduction (sum)
          // NB: read from (w-1) since left side came from attrsLf.
          transformBuf[1] = attrsHfIt[(w - 1) * numAttrs + k];
          attrSum[k] -= transformBuf[1];
          transformBuf[0] = attrSum[k];

          // NB: weight of transformBuf[1] is by construction 1.
          transformBuf[0].val >>= shift;
          transformBuf[0] *= rsqrtWeight;

          kernel.fwdTransform(transformBuf[0], transformBuf[1]);

          auto coeff = transformBuf[1].round();
          assert(coeff <= INT_MAX && coeff >= INT_MIN);
          *coeffBufItK[k]++ = coeff =
            q.quantize(coeff << kFixedPointAttributeShift);
          transformBuf[1] =
            divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);

          sumCoeff += std::abs(q.quantize(coeff << kFixedPointAttributeShift));
        } else {
          int64_t coeff = *coeffBufItK[k]++;
          transformBuf[1] =
            divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
        }

        // inherit the DC value
        transformBuf[0] = attrRecDc[k];

        kernel.invTransform(transformBuf[0], transformBuf[1]);

        attrRecDc[k] = transformBuf[0];
        attrRec[out + w * numAttrs + k] = transformBuf[1].val;
        if (w == 1)
          attrRec[out + k] = transformBuf[0].val;
      }

      // Track RL for RDOQ
      if (isEncoder) {
        if (sumCoeff == 0)
          trainZeros++;
        else
          trainZeros = 0;
      }
    }

    attrsHfIt += (weight - 1) * numAttrs;
    out += weight * numAttrs;
  }

  // write-back reconstructed attributes
  assert(attrRec.size() == numAttrs * numPoints);
  auto attrOut = attributes;
  for (auto attr : attrRec) {
    attr += FixedPoint::kOneHalf;
    attr >>= FixedPoint::kFracBits;
    *attrOut++ = attr;
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
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int64_t* mortonCode,
  int* attributes,
  const int attribCount,
  const int voxelCount,
  int* coefficients)
{
  uraht_process<true>(
    rahtPredParams, qpset, pointQpOffsets, voxelCount, attribCount, mortonCode,
    attributes, coefficients);
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
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int64_t* mortonCode,
  int* attributes,
  const int attribCount,
  const int voxelCount,
  int* coefficients)
{
  uraht_process<false>(
    rahtPredParams, qpset, pointQpOffsets, voxelCount, attribCount, mortonCode,
    attributes, coefficients);
}

//============================================================================

}  // namespace pcc
