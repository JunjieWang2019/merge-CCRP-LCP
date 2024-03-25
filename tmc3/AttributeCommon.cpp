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
  const EncodeMotionSearchParams& msParams,
  const ParameterSetMotion& mvPS,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud
) {
  if (!gbh.interPredictionEnabledFlag) {
    mSOctreeCurr = MSOctree();
    motionPUTrees.clear();
    return;
  }

  const MSOctree& mSOctree = mSOctreeRef;

  EntropyEncoder arithmeticEncoder;

  mSOctreeCurr = MSOctree(
    &pointCloud, {}, ilog2(uint32_t(std::min(
      mvPS.motion_min_pu_size,
      gbh.trisoupNodeSize(gps)) - 1)) + 1);
  motionPUTrees.clear();

  auto& fifo = mSOctreeCurr.a;
  auto& fifo_next = mSOctreeCurr.b;
  fifo.clear();
  fifo_next.clear();

  // build node list, at the level of the prediction units,
  // in raster scan order, to follow node coding order.
  fifo.push(0);
  while(mSOctreeCurr.nodes[fifo.front()].sizeMinus1 > mvPS.motion_block_size - 1) {
    IterOneLevelSubnodesRSO(fifo.begin(), fifo.end(),
    [&](const decltype(fifo.begin())& it) -> const point_t& {
      return mSOctreeCurr.nodes[*it].pos0;
    },
    [&](const decltype(fifo.begin())& it, int childIdx) {
      auto& node = mSOctreeCurr.nodes[*it];
      if (node.child[childIdx])
        fifo_next.push_back(node.child[childIdx]);
    }
    );
    std::swap(fifo, fifo_next);
    fifo_next.clear();
  }

  motionPUTrees.resize(fifo.size());
  // build root PU_trees
  const int nodeSizeLog2 = ilog2(uint32_t(mvPS.motion_block_size - 1)) + 1;
  int i = 0;
  while (!fifo.empty()) {
    auto& node = mSOctreeCurr.nodes[fifo.front()];

    PCCOctree3Node node0;
    node0.start = node.start;
    node0.end = node.end;
    node0.pos = node.pos0 >> nodeSizeLog2;
    node0.mSOctreeNodeIdx = mSOctreeCurr.nodeIdx(node.pos0, nodeSizeLog2);//fifo.front();
    node0.hasMotion = motionSearchForNode(mSOctreeCurr, mSOctree, &node0,
      msParams, mvPS, 1 << nodeSizeLog2, &arithmeticEncoder,
      &motionPUTrees[i].first);
    motionPUTrees[i].second = fifo.front();

    fifo.pop_front();
    ++i;
  }
}

//----------------------------------------------------------------------------

void
AttributeInterPredParams::encodeMotionAndBuildCompensated(
  const ParameterSetMotion& mvPS,
  EntropyEncoder& arithmeticEncoder,
  bool mcap_to_rec_geom_flag
) {
  const MSOctree& mSOctree = mSOctreeRef;

  if (!mcap_to_rec_geom_flag) {
    compensatedPointCloud.clear();
    mortonCode_mc.clear();
    attributes_mc.clear();
  }

  int nodeSizeLog2 = ilog2(uint32_t(mvPS.motion_block_size - 1)) + 1;
  auto currentPUTrees = motionPUTrees;
  while (currentPUTrees.size()) {
    // coding (in morton order for simpler test)
    decltype(currentPUTrees) nextLevelPUTrees;
    nextLevelPUTrees.reserve(currentPUTrees.size());
    int childSizeLog2 = nodeSizeLog2 - 1;
    for (int i = 0; i < currentPUTrees.size(); ++i) {
      auto& node = mSOctreeCurr.nodes[currentPUTrees[i].second];
      PCCOctree3Node node0;
      node0.start = node.start;
      node0.end = node.end;
      node0.pos = node.pos0 >> nodeSizeLog2;
      node0.isCompensated = false;
      assert((1 << nodeSizeLog2) - 1 == node.sizeMinus1);
      auto& local_PU_tree = currentPUTrees[i].first;

      encode_splitPU_MV_MC(mSOctree,
        &node0, &local_PU_tree, mvPS, 1 << nodeSizeLog2,
        &arithmeticEncoder, &compensatedPointCloud,
        false, -1, -1, mcap_to_rec_geom_flag);

      if (!node0.isCompensated && 1 << nodeSizeLog2 > mvPS.motion_min_pu_size) {
        node0.pos_fs = 1;
        node0.pos_fp = 0;
        node0.pos_MV = 0;
        for (int j = 0; j < 8; ++j) {
          //assert(local_PU_tree.popul_flags[node0.pos_fp] && it->child[j]
          //        || !local_PU_tree.popul_flags[node0.pos_fp] && !it->child[j]);
          if (local_PU_tree.popul_flags[node0.pos_fp++]) {
            assert(1 << childSizeLog2 >= mvPS.motion_min_pu_size);
            // populated
            nextLevelPUTrees.emplace_back(
              std::make_pair(PUtree(), int(node.child[j])));
            extracPUsubtree(
              mvPS, &local_PU_tree, 1 << childSizeLog2, node0.pos_fs,
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
  const ParameterSetMotion& mvPS,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud
) {
  if (!gbh.interPredictionEnabledFlag) {
    mSOctreeCurr = MSOctree();
    return;
  }
  mSOctreeCurr = MSOctree(
    &pointCloud, {}, ilog2(uint32_t(std::min(
      mvPS.motion_min_pu_size,
      gbh.trisoupNodeSize(gps)) - 1)) + 1);
}

//----------------------------------------------------------------------------

void
AttributeInterPredParams::decodeMotionAndBuildCompensated(
  const ParameterSetMotion& mvPS,
  EntropyDecoder& arithmeticDecoder,
  bool mcap_to_rec_geom_flag
) {
  const MSOctree& mSOctree = mSOctreeRef;

  if (!mcap_to_rec_geom_flag) {
    compensatedPointCloud.clear();
    mortonCode_mc.clear();
    attributes_mc.clear();
  }

  auto& fifo = mSOctreeCurr.a;
  auto& fifo_next = mSOctreeCurr.b;
  fifo.clear();
  fifo_next.clear();

  fifo.push(0);
  while(mSOctreeCurr.nodes[fifo.front()].sizeMinus1 > mvPS.motion_block_size - 1) {
    IterOneLevelSubnodesRSO(fifo.begin(), fifo.end(),
    [&](const decltype(fifo.begin())& it) -> const point_t& {
      return mSOctreeCurr.nodes[*it].pos0;
    },
    [&](const decltype(fifo.begin())& it, int childIdx) {
      auto& node = mSOctreeCurr.nodes[*it];
      if (node.child[childIdx])
        fifo_next.push_back(node.child[childIdx]);
    }
    );
    std::swap(fifo, fifo_next);
    fifo_next.clear();
  }

  int nodeSizeLog2 = ilog2(uint32_t(mvPS.motion_block_size - 1)) + 1;

  while(1 << (nodeSizeLog2 + 1) > mvPS.motion_min_pu_size) {
    while (!fifo.empty()) {
      // coding (in morton order for simpler test)
      auto& node = mSOctreeCurr.nodes[fifo.front()];

      PCCOctree3Node node0;
      node0.start = node.start;
      node0.end = node.end;
      node0.pos = node.pos0 >> nodeSizeLog2;
      node0.isCompensated = false;
      assert((1 << nodeSizeLog2) - 1 == node.sizeMinus1);

      decode_splitPU_MV_MC(mSOctree,
        &node0, mvPS, 1 << nodeSizeLog2,
        &arithmeticDecoder, &compensatedPointCloud,
        false, -1, -1, mcap_to_rec_geom_flag);

      if (!node0.isCompensated && 1 << nodeSizeLog2 > mvPS.motion_min_pu_size) {
        for (int i = 0; i < 8; ++i) {
          if (node.child[i]) {
            // populated
            fifo_next.push_back(node.child[i]);
          }
        }
      }
      fifo.pop_front();
    }
    std::swap(fifo, fifo_next);
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
    packedVoxel.push_back({mortonAddr(pointCloud[n]), n});
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
