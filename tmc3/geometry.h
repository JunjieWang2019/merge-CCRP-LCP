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

#include <memory>
#include <vector>

#include "PCCPointSet.h"
#include "geometry_params.h"
#include "entropy.h"
#include "hls.h"
#include "partitioning.h"
#include "TMC3.h"
#include "attr_tools.h"
namespace pcc {

//============================================================================

struct GeometryOctreeContexts;
struct CloudFrame;
struct RasterScanTrisoupEdges;

//============================================================================

void encodeGeometryOctree(
  const OctreeEncOpts& opt,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoder,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  PCCPointSet3& compensatedPointCloud);

void decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const Vec3<int> minimum_position,
  PCCPointSet3& compensatedPointCloud);

void decodeGeometryOctreeScalable(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int minGeomNodeSizeLog2,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame);

//----------------------------------------------------------------------------

void encodeGeometryOctreeForTrisoup(
  const OctreeEncOpts& opt,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyEncoder* arithmeticEncoder,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  PCCPointSet3& compensatedPointCloud,
  RasterScanTrisoupEdges& rste);

void decodeGeometryOctreeForTrisoup(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const Vec3<int> minimum_position,
  PCCPointSet3& compensatedPointCloud,
  RasterScanTrisoupEdges& rste);

//----------------------------------------------------------------------------

void encodeGeometryTrisoup(
  const TrisoupEncOpts& opt,
  const OctreeEncOpts& optOctree,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoder,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  PCCPointSet3& compensatedPointCloud);

void decodeGeometryTrisoup(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const Vec3<int> minimum_position,
  PCCPointSet3& compensatedPointCloud);

//============================================================================

}  // namespace pcc
