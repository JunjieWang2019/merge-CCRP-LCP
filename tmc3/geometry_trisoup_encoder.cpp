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

namespace pcc {

//============================================================================
void
encodeGeometryTrisoup(
  const EncoderParams& encParams,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMemOctree,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Encoder3& encoder)
{
  bool isInter = gbh.interPredictionEnabledFlag;

  // prepare TriSoup parameters
  int blockWidth = gbh.trisoupNodeSize(gps);

  std::cout << "TriSoup QP = " << gbh.trisoup_QP << "\n";

  // get first encoder
  pcc::EntropyEncoder* arithmeticEncoder = arithmeticEncoders.begin()->get();

  // trisoup uses octree coding until reaching the triangulation level.
  std::vector<PCCOctree3Node> nodes;
  EntropyDecoder foo;
  RasterScanTrisoupEdges rste(nodes, blockWidth, pointCloud, true,
    1, isInter, interPredParams.compensatedPointCloud,
    gps, gbh, arithmeticEncoder, foo, ctxtMemOctree);
  rste.useLocalAttr = sps.localized_attributes_enabled_flag;
  if (rste.useLocalAttr) {
    rste.encoder = &encoder;
    rste.slabThickness = sps.localized_attributes_slab_thickness_minus1 + 1;
  }
  rste.init();
  rste.thVertexDetermination = encParams.trisoup.thVertexDetermination;

  // octree
  encodeGeometryOctree<true>(
    encParams, gps, gbh, pointCloud, ctxtMemOctree, arithmeticEncoders, &nodes,
    refFrame, sps, interPredParams, encoder, &rste);

  std::cout << "Size compensatedPointCloud for TriSoup = "
    << interPredParams.compensatedPointCloud.getPointCount() << "\n";
  std::cout << "Number of nodes for TriSoup = " << nodes.size() << "\n";
  std::cout << "TriSoup gives " << pointCloud.getPointCount() << " points \n";

}

//============================================================================
}  // namespace pcc
