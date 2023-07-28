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
#include "TMC3.h"
#include <vector>

#include "PCCPointSet.h"
#include "PCCTMC3Encoder.h"
#include "entropy.h"
#include "hls.h"

namespace pcc {

//============================================================================
static const unsigned int motionParamPrec = 16;
static const unsigned int motionParamScale = 1 << motionParamPrec;
static const unsigned int motionParamOffset = 1 << (motionParamPrec - 1);
int plus1log2shifted4(int x);
struct PCCOctree3Node;

//============================================================================

int roundIntegerHalfInf(const double x);

//----------------------------------------- LOCAL MOTION -------------------
struct PUtree {
  std::vector<bool> popul_flags;
  std::vector<bool> split_flags;
  std::vector<Vec3<int>> MVs;
};

int deriveMotionMaxPrefixBits(const GeometryParameterSet::Motion& param);
int deriveMotionMaxSuffixBits(const GeometryParameterSet::Motion& param);

struct LPUwindow {
  Vec3<int> pos;
  Vec3<attr_t> color;
};

void buildActiveWindowAndBoundToBB(
  std::vector<std::vector<LPUwindow>>& lpuActiveWindow,
  int& LPUnumInAxis,
  const int maxBB,
  PCCPointSet3& predPointCloud,
  int th_dists,
  const int log2MotionBlockSize,
  Vec3<int> lvlNodeSizeLog2,
  point_t BBorig);

bool motionSearchForNode(
  const PCCPointSet3& pointCloud,
  const PCCOctree3Node* node0,
  const GeometryParameterSet::Motion& param,
  int nodeSizeLog2,
  EntropyEncoder* arithmeticEncoder,
  int8_t* bufferPoints,
  PUtree* local_PU_tree,
  const std::vector<std::vector<LPUwindow>>& lpuActiveWindow,
  int numLPUPerLine);


void encode_splitPU_MV_MC(
  PCCOctree3Node* node0,
  PUtree* local_PU_tree,
  const GeometryParameterSet::Motion& param,
  Vec3<int> nodeSizeLog2,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3* compensatedPointCloud,
  std::vector<std::vector<LPUwindow>>& lpuActiveWindow,
  int numLPUPerLine,
  int log2MotionBlkSize,
  std::vector<MotionVector>& motionVectors);

void extracPUsubtree(
  const GeometryParameterSet::Motion& param,
  PUtree* local_PU_tree,
  int block_size,
  int& pos_fs,
  int& pos_fp,
  int& pos_MV,
  PUtree* destination_tree);

// motion decoder

void decode_splitPU_MV_MC(
  PCCOctree3Node* node0,
  const GeometryParameterSet::Motion& param,
  Vec3<int> nodeSizeLog2,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3* compensatedPointCloud,
  std::vector<std::vector<LPUwindow>>& lpuActiveWindow,
  int numLPUPerLine,
  int log2MotionBlkSize,
  std::vector<MotionVector>& motionVectors);

  //============================================================================

}  // namespace pcc
