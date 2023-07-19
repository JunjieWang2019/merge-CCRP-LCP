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

namespace pcc {

//============================================================================
// Attribute methods

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
