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

namespace pcc {

//============================================================================
void
encodeGeometryTrisoup(
  const TrisoupEncOpts& opt,
  const OctreeEncOpts& optOctree,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMemOctree,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  PCCPointSet3& refPointCloud,
  MSOctree& mSOctree,
  PCCPointSet3& compensatedPointCloud)
{
  bool isInter = gbh.interPredictionEnabledFlag;

  // prepare TriSoup parameters
  int blockWidth = 1 << gbh.trisoupNodeSizeLog2(gps);
  const int maxVertexPrecisionLog2 = gbh.trisoup_vertex_quantization_bits ? gbh.trisoup_vertex_quantization_bits : gbh.trisoupNodeSizeLog2(gps);
  const int bitDropped = std::max(0, gbh.trisoupNodeSizeLog2(gps) - maxVertexPrecisionLog2);

  // TriSoup encoder parameters
  int distanceSearchEncoder = 1;
  /*if (opt.improvedVertexDetermination) { // !!!!!!!!!!!!!! impossible in one passe coding !!!!!!
    float estimatedSampling = float(nodes.size());
    estimatedSampling /= pointCloud.getPointCount();
    estimatedSampling = std::sqrt(estimatedSampling);
    estimatedSampling *= blockWidth;
    estimatedSampling = std::max(1.f, estimatedSampling);
    std::cout << "Estimation of sampling = " << estimatedSampling << "\n";

    distanceSearchEncoder = (1 << std::max(0, bitDropped - 2)) - 1;
    distanceSearchEncoder += int(std::round(estimatedSampling + 0.1f));
    distanceSearchEncoder = std::max(1, std::min(8, distanceSearchEncoder));
    std::cout << "distanceSearchEncoder = " << distanceSearchEncoder << "\n";
  } */

  // get first encoder
  pcc::EntropyEncoder* arithmeticEncoder = arithmeticEncoders.begin()->get();

  // trisoup uses octree coding until reaching the triangulation level.
  std::vector<PCCOctree3Node> nodes;
  EntropyDecoder foo;
  RasterScanTrisoupEdges rste(nodes, blockWidth, pointCloud, true, bitDropped,distanceSearchEncoder, isInter, compensatedPointCloud, gps, gbh, arithmeticEncoder, foo, ctxtMemOctree);
  rste.init();

  // octree
  encodeGeometryOctree(
    optOctree, gps, gbh, pointCloud, ctxtMemOctree, arithmeticEncoders, &nodes,
    refFrame, sps, refPointCloud, mSOctree, compensatedPointCloud, &rste);

  std::cout << "Size compensatedPointCloud for TriSoup = " << compensatedPointCloud.getPointCount() << "\n";
  std::cout << "Number of nodes for TriSoup = " << nodes.size() << "\n";
  std::cout << "TriSoup gives " << pointCloud.getPointCount() << " points \n";

  /*if (!(isInter
        && gps.gof_geom_entropy_continuation_enabled_flag)
      && !gbh.entropy_continuation_flag) {
    ctxtMemOctree.clearMap();
  }*/
}

//============================================================================
}  // namespace pcc
