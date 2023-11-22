/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2019, ISO/IEC
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

#include "AttributeCommon.h"

#include "PCCTMC3Common.h"

#include "geometry_octree.h"

#include "PCCTMC3Encoder.h"

namespace pcc {

//============================================================================
// Attribute methods

void
AttributeInterPredParams::findMotion(
  const EncoderParams* params,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud
) {
  if (!gbh.interPredictionEnabledFlag) {
    mSOctreeCurr = MSOctree();
    motionPUTrees.clear();
    return;
  }

  MSOctree& mSOctree = mSOctreeRef;

  EntropyEncoder arithmeticEncoder;

  auto scaleFactor = params->codedGeomScale;
  int TriSoupSize = 1;
  int presetMode = params->motionPreset;
  if (gps.trisoup_enabled_flag)
    TriSoupSize = 1 << gbh.trisoupNodeSizeLog2(gps);
  auto motion = gps.motion;
  motion.motion_min_pu_size = motion.motion_min_pu_size_color;
  if (gps.interPredictionEnabledFlag) {
    switch (presetMode) {
    case 0:
      break;

    case 1:
      // search parameters
      motion.Amotion0 = motion.motion_window_size >> 2;
      motion.lambda = 0.5;
      motion.decimate = 6;
      break;

    case 2:
      // search parameters
      motion.Amotion0 = 1;
      motion.lambda = 0.5/std::sqrt(4.*scaleFactor)*4.0;
      motion.decimate = 7;
      motion.dgeom_color_factor = 2;
      motion.K = 10;
      break;

    case 3:
      if (params->gps.trisoup_enabled_flag)
        TriSoupSize = 1 << (params->trisoupNodeSizesLog2[0]);
      // search parameters
      motion.Amotion0 = 2;
      motion.lambda = 2.5 * TriSoupSize * 4 / 32;
      motion.decimate = 7;
      motion.dgeom_color_factor = 2;
      motion.K = 10;
      break;

    default:
      // search parameters
      motion.Amotion0 = motion.motion_window_size >> 2;
      motion.lambda = 0.5;
      motion.decimate = 6;
    }
  }

  mSOctreeCurr = MSOctree(
    &pointCloud, {}, std::min(
      ilog2(uint32_t(gps.motion.motion_min_pu_size_color)),
      gbh.trisoupNodeSizeLog2(gps)));
  motionPUTrees.clear();

  auto nodeLevelStart = mSOctreeCurr.nodes.begin();
  while (nodeLevelStart != mSOctreeCurr.nodes.end() && nodeLevelStart->sizeMinus1 != gps.motion.motion_block_size - 1)
    nodeLevelStart++;

  auto nodeLevelEnd = nodeLevelStart;
  while (nodeLevelEnd != mSOctreeCurr.nodes.end() &&  nodeLevelEnd->sizeMinus1 == gps.motion.motion_block_size - 1)
    nodeLevelEnd++;

  int nodeSizeLog2 = ilog2(uint32_t(gps.motion.motion_block_size));

  motionPUTrees.resize(std::distance(nodeLevelStart, nodeLevelEnd));

  // build root PU_trees
  for (auto it = nodeLevelStart; it != nodeLevelEnd; ++it) {
    PCCOctree3Node node0;
    node0.start = it->start;
    node0.end = it->end;
    node0.pos = it->pos0 >> nodeSizeLog2;
    node0.mSOctreeNodeIdx = mSOctreeCurr.nodeIdx(it->pos0, nodeSizeLog2);
    int i = std::distance(nodeLevelStart, it);
    node0.hasMotion = motionSearchForNode(mSOctreeCurr, mSOctree, &node0, motion, nodeSizeLog2,
        &arithmeticEncoder, &motionPUTrees[i].first, true);
    motionPUTrees[i].second = it;
  }
}

//----------------------------------------------------------------------------

void
AttributeInterPredParams::encodeMotionAndBuildCompensated(
  const GeometryParameterSet& gps,
  EntropyEncoder& arithmeticEncoder
) {
  MSOctree& mSOctree = mSOctreeRef;

  int log2MotionBlockSize = 0; // unused

  compensatedPointCloud.clear();

  // dirty hack
  auto motion = gps.motion;
  motion.motion_min_pu_size = motion.motion_min_pu_size_color;

  int nodeSizeLog2 = ilog2(uint32_t(gps.motion.motion_block_size));
  auto currentPUTrees = motionPUTrees;
  while (currentPUTrees.size()) {
    // coding (in morton order for simpler test)
    decltype(currentPUTrees) nextLevelPUTrees;
    nextLevelPUTrees.reserve(currentPUTrees.size());
    int childSizeLog2 = nodeSizeLog2 - 1;
    for (int i = 0; i < currentPUTrees.size(); ++i) {
      auto it = currentPUTrees[i].second;
      PCCOctree3Node node0;
      node0.start = it->start;
      node0.end = it->end;
      node0.pos = it->pos0 >> nodeSizeLog2;
      node0.isCompensated = false;
      assert((1 << nodeSizeLog2) - 1 == it->sizeMinus1);
      auto& local_PU_tree = currentPUTrees[i].first;

      encode_splitPU_MV_MC(mSOctree,
        &node0, &local_PU_tree, motion, nodeSizeLog2,
        &arithmeticEncoder, &compensatedPointCloud,
        log2MotionBlockSize);

      if (!node0.isCompensated && 1 << childSizeLog2 >= gps.motion.motion_min_pu_size_color) {
        node0.pos_fs = 1;
        node0.pos_fp = 0;
        node0.pos_MV = 0;
        for (int j = 0; j < 8; ++j) {
          //assert(local_PU_tree.popul_flags[node0.pos_fp] && it->child[j]
          //        || !local_PU_tree.popul_flags[node0.pos_fp] && !it->child[j]);
          if (local_PU_tree.popul_flags[node0.pos_fp++]) {
            assert(1 << childSizeLog2 >= gps.motion.motion_min_pu_size_color);
            // populated
            nextLevelPUTrees.emplace_back(
              std::make_pair(PUtree(), mSOctreeCurr.nodes.begin() + it->child[j]));
            extracPUsubtree(
              motion, &local_PU_tree, 1 << childSizeLog2, node0.pos_fs,
              node0.pos_fp, node0.pos_MV, &nextLevelPUTrees.back().first);
          }
        }
      }
    }
    currentPUTrees.clear();
    std::swap(currentPUTrees, nextLevelPUTrees);
    nodeSizeLog2--;
  }
}

//----------------------------------------------------------------------------

void
AttributeInterPredParams::prepareDecodeMotion(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud
) {
  if (!gbh.interPredictionEnabledFlag) {
    mSOctreeCurr = MSOctree();
    return;
  }
  mSOctreeCurr = MSOctree(
    &pointCloud, {}, std::min(
      ilog2(uint32_t(gps.motion.motion_min_pu_size_color)),
      gbh.trisoupNodeSizeLog2(gps)));
}

//----------------------------------------------------------------------------

void
AttributeInterPredParams::decodeMotionAndBuildCompensated(
  const GeometryParameterSet& gps,
  EntropyDecoder& arithmeticDecoder
) {
  MSOctree& mSOctree = mSOctreeRef;

  int log2MotionBlockSize = 0; // unused

  compensatedPointCloud.clear();

  // dirty hack
  auto motion = gps.motion;
  motion.motion_min_pu_size = motion.motion_min_pu_size_color;

  auto nodeLevelStart = mSOctreeCurr.nodes.begin();
  while (nodeLevelStart != mSOctreeCurr.nodes.end() && nodeLevelStart->sizeMinus1 != gps.motion.motion_block_size - 1)
    nodeLevelStart++;

  auto nodeLevelEnd = nodeLevelStart;
  while (nodeLevelEnd != mSOctreeCurr.nodes.end() &&  nodeLevelEnd->sizeMinus1 == gps.motion.motion_block_size - 1)
    nodeLevelEnd++;

  int nodeSizeLog2 = ilog2(uint32_t(gps.motion.motion_block_size));

  std::vector<decltype(nodeLevelStart)> currentPUTrees(std::distance(nodeLevelStart, nodeLevelEnd));
  // build root PU_trees
  for (auto it = nodeLevelStart; it != nodeLevelEnd; ++it) {
    int i = std::distance(nodeLevelStart, it);
    currentPUTrees[i] = it;
  }

  while (currentPUTrees.size()) {
    // coding (in morton order for simpler test)
    decltype(currentPUTrees) nextLevelPUTrees;
    nextLevelPUTrees.reserve(currentPUTrees.size());
    for (int i = 0; i < currentPUTrees.size(); ++i) {
      auto it = currentPUTrees[i];
      PCCOctree3Node node0;
      node0.start = it->start;
      node0.end = it->end;
      node0.pos = it->pos0 >> nodeSizeLog2;
      node0.isCompensated = false;
      assert((1 << nodeSizeLog2) - 1 == it->sizeMinus1);

      decode_splitPU_MV_MC(mSOctree,
        &node0, motion, nodeSizeLog2,
        &arithmeticDecoder, &compensatedPointCloud,
        log2MotionBlockSize);

      if (!node0.isCompensated && 1 << (nodeSizeLog2-1) >= gps.motion.motion_min_pu_size_color) {
        for (int j = 0; j < 8; ++j) {
          if (it->child[j]) {
            // populated
            nextLevelPUTrees.emplace_back(mSOctreeCurr.nodes.begin() + it->child[j]);
          }
        }
      }
    }
    currentPUTrees.clear();
    std::swap(currentPUTrees, nextLevelPUTrees);
    nodeSizeLog2--;
  }
}

//============================================================================

std::vector<int>
sortedPointCloud(
  const int attribCount,
  const PCCPointSet3& pointCloud,
  std::vector<int64_t>& mortonCode,
  std::vector<int>& attributes)
{
  const auto voxelCount = pointCloud.getPointCount();
  std::vector<MortonCodeWithIndex> packedVoxel;
  packedVoxel.reserve(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    packedVoxel.push_back({mortonAddr(pointCloud[n]), pointCloud[n], n});
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  std::vector<int> indexOrd;
  mortonCode.reserve(voxelCount);
  indexOrd.reserve(voxelCount);
  for (auto& voxel : packedVoxel) {
    mortonCode.push_back(voxel.mortonCode);
    indexOrd.push_back(voxel.index);
  }
  packedVoxel.clear();

  if (attribCount==3 && pointCloud.hasColors()) {
    attributes.reserve(voxelCount * 3);
    for (auto index : indexOrd) {
      const auto& color = pointCloud.getColor(index);
      attributes.push_back(color[0]);
      attributes.push_back(color[1]);
      attributes.push_back(color[2]);
    }
  } else if (attribCount == 1 && pointCloud.hasReflectances()) {
    attributes.reserve(voxelCount);
    for (auto index : indexOrd) {
      attributes.push_back(pointCloud.getReflectance(index));
    }
  }

  return indexOrd;
}

//============================================================================

}  // namespace pcc
