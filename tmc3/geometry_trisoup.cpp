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

#include "geometry_trisoup.h"

#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"

#include "PCCTMC3Encoder.h"
#include "PCCTMC3Decoder.h"

namespace pcc {

//============================================================================

// instanciate for Encoder
template void RasterScanTrisoupEdges<true, TrisoupNodeEncoder>::processLocalAttributes(
  PCCPointSet3& recPointCloud, bool isLast);

// instanciate for Decoder
template void RasterScanTrisoupEdges<false, TrisoupNodeDecoder>::processLocalAttributes(
  PCCPointSet3& recPointCloud, bool isLast);

//---

template <bool isEncoder, typename TrisoupNode>
void
RasterScanTrisoupEdges<isEncoder, TrisoupNode>::processLocalAttributes(
  PCCPointSet3& recPointCloud, bool isLast)
{
  auto allocatedSizeLocal = localPointCloud.size();
  localPointCloud.resize(nRecPointsLocal);
  if (isEncoder)
    encoder->processNextSlabAttributes(localPointCloud, xStartLocalSlab, isLast);
  else
    decoder->processNextSlabAttributes(localPointCloud, xStartLocalSlab, isLast);
  localPointCloud.resize(allocatedSizeLocal);

  if (recPointCloud.getPointCount() < nRecPoints + nRecPointsLocal)
    recPointCloud.resize(nRecPoints + nRecPointsLocal + PC_PREALLOCATION_SIZE);

  recPointCloud.setFromPartition(localPointCloud, 0, nRecPointsLocal, nRecPoints);
  nRecPoints += nRecPointsLocal;
  nRecPointsLocal = 0; // point cloud buffer has been rendered
}

//============================================================================
}  // namespace pcc
