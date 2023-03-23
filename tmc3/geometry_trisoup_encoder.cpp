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
  const SequenceParameterSet& sps)
{
  // trisoup uses octree coding until reaching the triangulation level.
  std::vector<PCCOctree3Node> nodes;
  PCCPointSet3 compensatedPointCloud;  // set of points after motion compensation
  encodeGeometryOctree(
    optOctree, gps, gbh, pointCloud, ctxtMemOctree, arithmeticEncoders, &nodes,
    refFrame, sps, compensatedPointCloud);

  std::cout << "\nSize compensatedPointCloud for TriSoup = " << compensatedPointCloud.getPointCount() << "\n";
  bool isInter = gbh.interPredictionEnabledFlag;

  // resume encoding with the last encoder
  pcc::EntropyEncoder* arithmeticEncoder = arithmeticEncoders.back().get();    

  int blockWidth = 1 << gbh.trisoupNodeSizeLog2(gps);
  const int maxVertexPrecisionLog2 = gbh.trisoup_vertex_quantization_bits
    ? gbh.trisoup_vertex_quantization_bits
    : gbh.trisoupNodeSizeLog2(gps);
  const int bitDropped =
    std::max(0, gbh.trisoupNodeSizeLog2(gps) - maxVertexPrecisionLog2);
  const bool isCentroidDriftActivated =
    gbh.trisoup_centroid_vertex_residual_flag;

  std::cout << "Number of points for TriSoup = " << pointCloud.getPointCount() << "\n";
  std::cout << "Number of nodes for TriSoup = " << nodes.size() << "\n";
  int distanceSearchEncoder = 1;
  if (opt.improvedVertexDetermination) {
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
  }

  // Determine neighbours
  std::vector<uint16_t> neighbNodes;
  std::vector<std::array<int, 18>> edgePattern;
  std::vector<int> segmentUniqueIndex;
  int Nunique;
  determineTrisoupNeighbours(nodes, neighbNodes, edgePattern, blockWidth, segmentUniqueIndex, Nunique);

  // Determine vertices
  std::vector<bool> segind;
  std::vector<uint8_t> vertices;
  determineTrisoupVertices(
    nodes, segind, vertices, pointCloud, pointCloud, gps, gbh, blockWidth,
    bitDropped, distanceSearchEncoder, false, segmentUniqueIndex, Nunique);

  // determine vertices from compensated point cloud
  std::vector<bool> segindPred;
  std::vector<uint8_t> verticesPred;
  if (isInter) {
    determineTrisoupVertices(
      nodes, segindPred, verticesPred, refFrame.cloud, compensatedPointCloud,
      gps, gbh, blockWidth, bitDropped, 1 /*distanceSearchEncoder*/, true, segmentUniqueIndex, Nunique);
  }

  gbh.num_unique_segments_minus1 = segind.size() - 1;
  gbh.num_unique_segments_bits_minus1 = numBits(gbh.num_unique_segments_minus1) - 1;

  // Encode vertex presence and position into bitstream
  assert(segind.size() > 0);
  encodeTrisoupVertices(segind, vertices, segindPred, verticesPred, neighbNodes, edgePattern, bitDropped,gps, gbh, arithmeticEncoder, ctxtMemOctree);

  // Decode vertices with certain sampling value
  bool haloFlag = gbh.trisoup_halo_flag;
  bool adaptiveHaloFlag = gbh.trisoup_adaptive_halo_flag;
  bool fineRayFlag = gbh.trisoup_fine_ray_tracing_flag;
  int thickness = gbh.trisoup_thickness;

  PCCPointSet3 recPointCloud;
  recPointCloud.addRemoveAttributes(pointCloud);

  std::vector<CentroidDrift> drifts;
  int subsample = 1;
  int32_t maxval = (1 << gbh.maxRootNodeDimLog2) - 1;
  std::cout << "Loop on sampling for max "
            << (gbh.footer.geom_num_points_minus1 + 1) << " points \n";
  if (gps.trisoup_sampling_value > 0) {
    subsample = gps.trisoup_sampling_value;
    decodeTrisoupCommon(
      nodes, segind, vertices, drifts, pointCloud, recPointCloud,
      compensatedPointCloud, gps, gbh, blockWidth,
      maxval, subsample, bitDropped, isCentroidDriftActivated, false,
      haloFlag, adaptiveHaloFlag, fineRayFlag, thickness, NULL, ctxtMemOctree, segmentUniqueIndex);
    std::cout << "Sub-sampling " << subsample << " gives "
              << recPointCloud.getPointCount() << " points \n";
  } else {
    int maxSubsample = 1 << gbh.trisoupNodeSizeLog2(gps);
    for (subsample = 1; subsample <= maxSubsample; subsample++) {
      decodeTrisoupCommon(
        nodes, segind, vertices, drifts, pointCloud, recPointCloud,
        compensatedPointCloud, gps, gbh, blockWidth,
        maxval, subsample, bitDropped, isCentroidDriftActivated, false,
        haloFlag, adaptiveHaloFlag, fineRayFlag, thickness, NULL, ctxtMemOctree, segmentUniqueIndex);

      std::cout << "Sub-sampling " << subsample << " gives "
                << recPointCloud.getPointCount() << " points \n";
      if (recPointCloud.getPointCount() <= gbh.footer.geom_num_points_minus1 + 1)
        break;
    }
  }

  pointCloud.resize(0);
  pointCloud = std::move(recPointCloud);

  gbh.trisoup_sampling_value_minus1 = subsample - 1;

  // encoder centroid residua into bitstream
  if (isCentroidDriftActivated)
    encodeTrisoupCentroidResidue(drifts, arithmeticEncoder, ctxtMemOctree);

  if (!(gps.localMotionEnabled && gps.gof_geom_entropy_continuation_enabled_flag) && !gbh.entropy_continuation_flag) {
    ctxtMemOctree.clearMap();
  }
}

//---------------------------------------------------------------------------
// Determine where the surface crosses each leaf
// (i.e., determine the segment indicators and vertices)
// from the set of leaves and the points in each leaf.
//
// @param leaves, list of blocks containing the surface
// @param segind, indicators for edges of blocks if they intersect the surface
// @param vertices, locations of intersections

void
determineTrisoupVertices(
  const std::vector<PCCOctree3Node>& leaves,
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  const PCCPointSet3& pointCloud,
  const PCCPointSet3& compensatedPointCloud,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const int defaultBlockWidth,
  const int bitDropped,
  int distanceSearchEncoder,
  bool isCompensated,
  std::vector<int>& segmentUniqueIndex,
  int Nunique)
{
  Box3<int32_t> sliceBB;
  sliceBB.min = gbh.slice_bb_pos << gbh.slice_bb_pos_log2_scale;
  sliceBB.max = sliceBB.min + ( gbh.slice_bb_width << gbh.slice_bb_width_log2_scale );

  // prepare accumulators
  std::vector<int> count(Nunique, 0);
  std::vector<int> count2(Nunique, 0);
  std::vector<int> distanceSum(Nunique, 0);
  std::vector<int> distanceSum2(Nunique, 0);

  // Put all leaves' edges into a list.
  const int bufferIdx0[27] = { 12, 12, 12, 12, 4, 7, 12, 5, 6, 12, 1, 3, 0, 4, 7, 2, 5, 6, 12, 9, 11, 8, 4, 7, 10, 5, 6 };
  const int voxelComp0[27] = { 0, 0, 0, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 2, 2 };
  const int bufferIdx1[27] = { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 0, 0, 12, 2, 2, 12, 12, 12, 12, 8, 8, 12, 10, 10 };
  const int voxelComp1[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  const int bufferIdx2[27] = { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 1, 3, 12, 1, 3, 12, 12, 12, 12, 9, 11, 12, 9, 11 };
  const int voxelComp2[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1 };

  // create segments
  int segIdx = 0; // for stable ordering
  for (int i = 0; i < leaves.size(); i++) {
    const auto& leaf = leaves[i];

    // Size of block
    Vec3<int32_t> newP, newW, corner[8];
    nonCubicNode( gps, gbh, leaf.pos, defaultBlockWidth, sliceBB, newP, newW, corner );

    // Each voxel votes for a position along each edge it is close to
    const int tmin = 1;
    const Vec3<int> tmax( newW.x() - tmin - 1,
                          newW.y() - tmin - 1,
                          newW.z() - tmin - 1 );
    const int tmin2 = distanceSearchEncoder;
    const Vec3<int> tmax2( newW.x() - tmin2 - 1,
                           newW.y() - tmin2 - 1,
                           newW.z() - tmin2 - 1 );

    // pick adequate point cloud
    int idxStart = leaf.start;
    int idxEnd = leaf.end;
    if (isCompensated) {
      idxStart = leaf.predStart;
      idxEnd = leaf.predEnd;
    }

    // loop on voxels
    for (int j = idxStart; j < idxEnd; j++) {
      Vec3<int> voxel = (isCompensated && leaf.isCompensated ? compensatedPointCloud[j] : pointCloud[j]) - newP;

      // parameter indicating threshold of how close voxels must be to edge
      // ----------- 1 -------------------
      int check0 = (voxel[0] < tmin) + 2 * (voxel[0] > tmax.x());
      int check1 = (voxel[1] < tmin) + 2 * (voxel[1] > tmax.y());
      int check2 = (voxel[2] < tmin) + 2 * (voxel[2] > tmax.z());
      if (check0 > 2) check0 = 1;
      if (check1 > 2) check1 = 1;
      if (check2 > 2) check2 = 1;

      int check = check0 + 3 * check1 + 9 * check2; // in [0,26]
      if (bufferIdx0[check] < 12) {
        int unique0 = segmentUniqueIndex[12 * i + bufferIdx0[check]];
        count[unique0]++;
        distanceSum[unique0] += voxel[voxelComp0[check]];
      }
      if (bufferIdx1[check] < 12) {
        int unique1 = segmentUniqueIndex[12 * i + bufferIdx1[check]];
        count[unique1]++;
        distanceSum[unique1] += voxel[voxelComp1[check]];
      }
      if (bufferIdx2[check] < 12) {
        int unique2 = segmentUniqueIndex[12 * i + bufferIdx2[check]];
        count[unique2]++;
        distanceSum[unique2] += voxel[voxelComp2[check]];
      }


      // parameter indicating threshold of how close voxels must be to edge
      // ----------- 2 -------------------
      if (distanceSearchEncoder > 1) {

        check0 = (voxel[0] < tmin2) + 2 * (voxel[0] > tmax2.x());
        check1 = (voxel[1] < tmin2) + 2 * (voxel[1] > tmax2.y());
        check2 = (voxel[2] < tmin2) + 2 * (voxel[2] > tmax2.z());
        if (check0 > 2) check0 = 1;
        if (check1 > 2) check1 = 1;
        if (check2 > 2) check2 = 1;

        check = check0 + 3 * check1 + 9 * check2; // in [0,26]

        if (bufferIdx0[check] < 12) {
          int unique0 = segmentUniqueIndex[12 * i + bufferIdx0[check]];
          count2[unique0]++;
          distanceSum2[unique0] += voxel[voxelComp0[check]];
        }
        if (bufferIdx1[check] < 12) {
          int unique1 = segmentUniqueIndex[12 * i + bufferIdx1[check]];
          count2[unique1]++;
          distanceSum2[unique1] += voxel[voxelComp1[check]];
        }
        if (bufferIdx2[check] < 12) {
          int unique2 = segmentUniqueIndex[12 * i + bufferIdx2[check]];
          count2[unique2]++;
          distanceSum2[unique2] += voxel[voxelComp2[check]];
        }
      }
    }
  }

  segind.resize(Nunique);
  for (int t = 0; t < Nunique; t++) {
    segind[t] = count[t] > 0 || count2[t] > 1;
    if (segind[t]) {
      int temp = ((2 * distanceSum[t] + distanceSum2[t]) << (10 - bitDropped))/ (2 * count[t] + count2[t]);
      int8_t vertex = (temp + (1 << 9 - bitDropped)) >> 10;
      vertices.push_back(vertex);
    }
  }
}

//------------------------------------------------------------------------------------
void
encodeTrisoupVertices(
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<bool>& segindPred,
  std::vector<uint8_t>& verticesPred,
  std::vector<uint16_t>& neighbNodes,
  std::vector<std::array<int, 18>>& edgePattern,
  int bitDropped,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  pcc::EntropyEncoder* arithmeticEncoder,
  GeometryOctreeContexts& ctxtMemOctree)
{
  const int nbitsVertices = gbh.trisoupNodeSizeLog2(gps) - bitDropped;
  const int max2bits = nbitsVertices > 1 ? 3 : 1;
  const int mid2bits = nbitsVertices > 1 ? 2 : 1;

  int iV = 0;
  int iVPred = 0;
  std::vector<int> correspondanceSegment2V(segind.size(), -1);

  for (int i = 0; i <= gbh.num_unique_segments_minus1; i++) {
    // reduced neighbour contexts
    int ctxE = (!!(neighbNodes[i] & 1)) + (!!(neighbNodes[i] & 2)) + (!!(neighbNodes[i] & 4)) + (!!(neighbNodes[i] & 8)) - 1; // at least one node is occupied
    int ctx0 = (!!(neighbNodes[i] & 16)) + (!!(neighbNodes[i] & 32)) + (!!(neighbNodes[i] & 64)) + (!!(neighbNodes[i] & 128));
    int ctx1 = (!!(neighbNodes[i] & 256)) + (!!(neighbNodes[i] & 512)) + (!!(neighbNodes[i] & 1024)) + (!!(neighbNodes[i] & 2048));
    int direction = neighbNodes[i] >> 13; // 0=x, 1=y, 2=z

    // construct pattern
    auto patternIdx = edgePattern[i];
    int pattern = 0;
    int patternClose  = 0;
    int patternClosest  = 0;
    int nclosestPattern = 0;

    int towardOrAway[18] = { // 0 = toward; 1 = away
      0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    int mapping18to9[3][9] = {
      { 0, 1, 2, 3,  4, 15, 14, 5,  7},
      { 0, 1, 2, 3,  9, 15, 14, 7, 12},
      { 0, 1, 2, 9, 10, 15, 14, 7, 12}
    };

    for (int v = 0; v < 9; v++) {
      int v18 = mapping18to9[direction][v];

      if (patternIdx[v18] != -1) {
        int idxEdge = patternIdx[v18];
        if (segind[idxEdge]) {
          pattern |= 1 << v;
          int vertexPos2bits = vertices[correspondanceSegment2V[idxEdge]] >> std::max(0, nbitsVertices - 2);
          if (towardOrAway[v18])
            vertexPos2bits = max2bits - vertexPos2bits; // reverses for away
          if (vertexPos2bits >= mid2bits)
            patternClose |= 1 << v;
          if (vertexPos2bits >= max2bits)
            patternClosest |= 1 << v;
          nclosestPattern += vertexPos2bits >= max2bits && v <= 4;
        }
      }
    }

    int missedCloseStart = /*!(pattern & 1)*/ + !(pattern & 2) + !(pattern & 4);
    int nclosestStart = !!(patternClosest & 1) + !!(patternClosest & 2) + !!(patternClosest & 4);
    if (direction == 0) {
      missedCloseStart +=  !(pattern & 8) + !(pattern & 16);
      nclosestStart +=  !!(patternClosest & 8) + !!(patternClosest & 16);
    }
    if (direction == 1) {
      missedCloseStart +=  !(pattern & 8);
      nclosestStart +=  !!(patternClosest & 8) - !!(patternClosest & 16) ;
    }
    if (direction == 2) {
      nclosestStart +=  - !!(patternClosest & 8) - !!(patternClosest & 16) ;
    }

    // reorganize neighbours of vertex /edge (endpoint) independently on xyz
    int neighbEdge = (neighbNodes[i] >> 0) & 15;
    int neighbEnd = (neighbNodes[i] >> 4) & 15;
    int neighbStart = (neighbNodes[i] >> 8) & 15;
    if (direction == 2) {
      neighbEdge = ((neighbNodes[i] >> 0 + 0) & 1);
      neighbEdge += ((neighbNodes[i] >> 0 + 3) & 1) << 1;
      neighbEdge += ((neighbNodes[i] >> 0 + 1) & 1) << 2;
      neighbEdge += ((neighbNodes[i] >> 0 + 2) & 1) << 3;

      neighbEnd = ((neighbNodes[i] >> 4 + 0) & 1);
      neighbEnd += ((neighbNodes[i] >> 4 + 3) & 1) << 1;
      neighbEnd += ((neighbNodes[i] >> 4 + 1) & 1) << 2;
      neighbEnd += ((neighbNodes[i] >> 4 + 2) & 1) << 3;

      neighbStart = ((neighbNodes[i] >> 8 + 0) & 1);
      neighbStart += ((neighbNodes[i] >> 8 + 3) & 1) << 1;
      neighbStart += ((neighbNodes[i] >> 8 + 1) & 1) << 2;
      neighbStart += ((neighbNodes[i] >> 8 + 2) & 1) << 3;
    }


    // encode flag vertex
    int ctxMap1 = std::min(nclosestPattern, 2) * 15 * 2 +  (neighbEdge-1) * 2 + ((ctx1 == 4));    // 2* 15 *3 = 90 -> 7 bits

    int ctxMap2 = neighbEnd << 11;
    ctxMap2 |= (patternClose & (0b00000110)) << 9 - 1 ; // perp that do not depend on direction = to start
    ctxMap2 |= direction << 7;
    ctxMap2 |= (patternClose & (0b00011000))<< 5-3; // perp that  depend on direction = to start or to end
    ctxMap2 |= (patternClose & (0b00000001))<< 4;  // before
    int orderedPclosePar = (((pattern >> 5) & 3) << 2) + (!!(pattern & 128) << 1) + !!(pattern & 256);
    ctxMap2 |= orderedPclosePar;

    bool isInter = gbh.interPredictionEnabledFlag  ;
    int ctxInter =  isInter ? 1 + segindPred[i] : 0;

    int ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[ctxInter][0].getEvolve(
      segind[i], ctxMap2, ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
      ctxtMemOctree._BufferOBUFleavesTrisoup);
    arithmeticEncoder->encode(
      (int)segind[i], ctxtMemOctree.ctxTriSoup[0][ctxInter][ctxTrisoup]);

    // encode position vertex
    if (segind[i]) {
      int v = 0;
      auto vertex = vertices[iV];
      correspondanceSegment2V[i] = iV;

      int ctxFullNbounds = (4 * (ctx0 <= 1 ? 0 : (ctx0 >= 3 ? 2 : 1)) + (std::max(1, ctx1) - 1)) * 2 + (ctxE == 3);
      int b = nbitsVertices - 1;

      // first bit
      ctxMap1 = ctxFullNbounds * 2 + (nclosestStart > 0);
      ctxMap2 = missedCloseStart << 8;
      ctxMap2 |= (patternClosest & 1) << 7;
      ctxMap2 |= direction << 5;
      ctxMap2 |= patternClose & (0b00011111);
      int orderedPclosePar = (((patternClose >> 5) & 3) << 2) + (!!(patternClose & 128) << 1) + !!(patternClose & 256);

      ctxInter = 0;
      if (isInter) {
        ctxInter = segindPred[i] ? 1 + ((verticesPred[iVPred] >> b-1) & 3) : 0;
      }

      int bit = (vertex >> b--) & 1;

      ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[ctxInter][1].getEvolve(
        bit, ctxMap2, ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
        ctxtMemOctree._BufferOBUFleavesTrisoup);
      arithmeticEncoder->encode(
        bit, ctxtMemOctree.ctxTriSoup[1][ctxInter][ctxTrisoup]);
      v = bit;

      // second bit
      if (b >= 0) {
        ctxMap1 = ctxFullNbounds * 2 + (nclosestStart > 0);
        ctxMap2 = missedCloseStart << 8;
        ctxMap2 |= (patternClose & 1) << 7;
        ctxMap2 |= (patternClosest & 1) << 6;
        ctxMap2 |= direction << 4;
        ctxMap2 |= (patternClose & (0b00011111)) >> 1;
        ctxMap2 = (ctxMap2 << 4) + orderedPclosePar;

        ctxInter = 0;
        if (isInter) {
          ctxInter = segindPred[i] ? 1 + ((verticesPred[iVPred] >> b) <= (v << 1)) : 0;
        }

        bit = (vertex >> b--) & 1;
        ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[ctxInter][2].getEvolve(
          bit, ctxMap2, (ctxMap1 << 1) + v,
          &ctxtMemOctree._OBUFleafNumberTrisoup,
          ctxtMemOctree._BufferOBUFleavesTrisoup);
        arithmeticEncoder->encode(
          bit, ctxtMemOctree.ctxTriSoup[2][ctxInter][ctxTrisoup]);
        v = (v << 1) | bit;
      }

      // third bit
      if (b >= 0) {
        int ctxFullNboundsReduced1 = (6 * (ctx0 >> 1) + missedCloseStart) * 2 + (ctxE == 3);
        bit = (vertex >> b--) & 1;
        arithmeticEncoder->encode(
          bit, ctxtMemOctree.ctxTempV2[4 * ctxFullNboundsReduced1 + v]);
        v = (v << 1) | bit;
      }

      // remaining bits are bypassed
      for (; b >= 0; b--)
        arithmeticEncoder->encode((vertex >> b) & 1);
      iV++;
    }

    if (isInter && segindPred[i])
      iVPred++;
  }

}

//-------------------------------------------------------------------------------------
void
encodeTrisoupCentroidResidue(
  std::vector<CentroidDrift>& drifts, pcc::EntropyEncoder* arithmeticEncoder, GeometryOctreeContexts& ctxtMemOctree)
{
  //AdaptiveBitModel ctxDrift0[9];
  //AdaptiveBitModel ctxDriftSign[3][8][8];
  //AdaptiveBitModel ctxDriftMag[4];
  for (int i = 0; i < drifts.size(); i++) {
    int driftQ = drifts[i].driftQ;
    int driftQPred = drifts[i].driftQPred;   
   
    if (driftQPred==-100) //intra 
      arithmeticEncoder->encode(driftQ == 0, ctxtMemOctree.ctxDrift0[drifts[i].ctxMinMax][0]);
    else //inter      
      arithmeticEncoder->encode(driftQ == 0, ctxtMemOctree.ctxDrift0[drifts[i].ctxMinMax][1+std::min(3,std::abs(driftQPred))]);    
    
    
    // if not 0
    // drift in [-lowBound; highBound]
    if (driftQ) {
      int lowBound = drifts[i].lowBound;
      int highBound = drifts[i].highBound;
      // code sign
      int lowS = std::min(7, drifts[i].lowBoundSurface);
      int highS = std::min(7, drifts[i].highBoundSurface);
      if (highBound && lowBound) {  // otherwise sign is known
        arithmeticEncoder->encode(driftQ > 0, ctxtMemOctree.ctxDriftSign[lowBound == highBound ? 0 : 1 + (lowBound < highBound)][lowS][highS][(driftQPred && driftQPred != -100) ? 1 + (driftQPred > 0) : 0]);
      }

      // code remaining bits 1 to 7 at most
      int magBound = (driftQ > 0 ? highBound : lowBound) - 1;
      bool sameSignPred = driftQPred != -100 && (driftQPred > 0 && driftQ > 0) || (driftQPred < 0 && driftQ < 0);      

      int magDrift = std::abs(driftQ) - 1;
      int ctx = 0;
      while (magBound > 0 && magDrift >= 0) {
        if (ctx < 4)
          arithmeticEncoder->encode(magDrift == 0, ctxtMemOctree.ctxDriftMag[ctx][driftQPred != -100 ? 1 + std::min(8, sameSignPred * std::abs(driftQPred)) : 0]);
        else
          arithmeticEncoder->encode(magDrift == 0);

        magDrift--;
        magBound--;
        ctx++;
      }
    }  // end if not 0

  }  // end loop on drifts
}

//============================================================================

}  // namespace pcc
