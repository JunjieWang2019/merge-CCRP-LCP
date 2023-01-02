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
  const InterGeomEncOpts& interParams)
{
  // trisoup uses octree coding until reaching the triangulation level.
  pcc::ringbuf<PCCOctree3Node> nodes;
  PCCPointSet3 compensatedPointCloud;  // set of points after motion compensation
  encodeGeometryOctree(
    optOctree, gps, gbh, pointCloud, ctxtMemOctree, arithmeticEncoders, &nodes,
    refFrame, sps, interParams, compensatedPointCloud);

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

  // Determine vertices
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

  std::vector<bool> segind;
  std::vector<uint8_t> vertices;
  determineTrisoupVertices(
    nodes, segind, vertices, pointCloud, blockWidth, bitDropped,
    distanceSearchEncoder, false);

  // determine vertices from compensated point cloud 
  std::vector<bool> segindPred;
  std::vector<uint8_t> verticesPred;
  if (isInter) {
    determineTrisoupVertices(
      nodes, segindPred, verticesPred, compensatedPointCloud, blockWidth, bitDropped,
      1 /*distanceSearchEncoder*/, true);    
  }

  
  // Determine neighbours
  std::vector<uint16_t> neighbNodes;
  std::vector<std::array<int, 18>> edgePattern;
  determineTrisoupNeighbours(nodes, neighbNodes, edgePattern, blockWidth);

  gbh.num_unique_segments_minus1 = segind.size() - 1;
  gbh.num_unique_segments_bits_minus1 = numBits(gbh.num_unique_segments_minus1) - 1;

  // Encode vertex presence and position into bitstream
  assert(segind.size() > 0);
  encodeTrisoupVertices(segind, vertices, segindPred, verticesPred, neighbNodes, edgePattern, bitDropped,gps, gbh, arithmeticEncoder, ctxtMemOctree);

  // Decode vertices with certain sampling value
  bool haloFlag = gbh.trisoup_halo_flag;
  bool adaptiveHaloFlag = gbh.trisoup_adaptive_halo_flag;
  bool fineRayFlag = gbh.trisoup_fine_ray_tracing_flag;

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
      nodes, segind, vertices, drifts, pointCloud, recPointCloud, compensatedPointCloud, blockWidth,
      maxval, subsample, bitDropped, isCentroidDriftActivated, false,
      haloFlag, adaptiveHaloFlag, fineRayFlag, NULL, ctxtMemOctree);
    std::cout << "Sub-sampling " << subsample << " gives "
              << recPointCloud.getPointCount() << " points \n";
  } else {
    int maxSubsample = 1 << gbh.trisoupNodeSizeLog2(gps);
    for (subsample = 1; subsample <= maxSubsample; subsample++) {
      decodeTrisoupCommon(
        nodes, segind, vertices, drifts, pointCloud, recPointCloud, compensatedPointCloud, blockWidth,
        maxval, subsample, bitDropped, isCentroidDriftActivated, false,
        haloFlag, adaptiveHaloFlag, fineRayFlag, NULL, ctxtMemOctree);

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
  const ringbuf<PCCOctree3Node>& leaves,
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  const PCCPointSet3& pointCloud,
  const int defaultBlockWidth,
  const int bitDropped,
  int distanceSearchEncoder,
  bool isCompensated)
{
  // Put all leaves' edges into a list.
  std::vector<TrisoupSegmentEnc> segments;
  segments.reserve(12 * leaves.size());
  for (int i = 0; i < leaves.size(); i++) {
    const auto& leaf = leaves[i];

    // Width of block.
    // in future, may override with leaf blockWidth
    const int32_t blockWidth = defaultBlockWidth;

    // Eight corners of block.
    const Vec3<int32_t> pos000({0, 0, 0});
    const Vec3<int32_t> posW00({blockWidth, 0, 0});
    const Vec3<int32_t> pos0W0({0, blockWidth, 0});
    const Vec3<int32_t> posWW0({blockWidth, blockWidth, 0});
    const Vec3<int32_t> pos00W({0, 0, blockWidth});
    const Vec3<int32_t> posW0W({blockWidth, 0, blockWidth});
    const Vec3<int32_t> pos0WW({0, blockWidth, blockWidth});
    const Vec3<int32_t> posWWW({blockWidth, blockWidth, blockWidth});

    // x: left to right; y: bottom to top; z: far to near
    TrisoupSegmentEnc seg000W00 =  // far bottom edge
      {leaf.pos + pos000,
       leaf.pos + posW00,
       12 * i + 0,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction x
    TrisoupSegmentEnc seg0000W0 =  // far left edge
      {leaf.pos + pos000,
       leaf.pos + pos0W0,
       12 * i + 1,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction y
    TrisoupSegmentEnc seg0W0WW0 =  // far top edge
      {leaf.pos + pos0W0,
       leaf.pos + posWW0,
       12 * i + 2,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction x
    TrisoupSegmentEnc segW00WW0 =  // far right edge
      {leaf.pos + posW00,
       leaf.pos + posWW0,
       12 * i + 3,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction y
    TrisoupSegmentEnc seg00000W =  // bottom left edge
      {leaf.pos + pos000,
       leaf.pos + pos00W,
       12 * i + 4,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction z
    TrisoupSegmentEnc seg0W00WW =  // top left edge
      {leaf.pos + pos0W0,
       leaf.pos + pos0WW,
       12 * i + 5,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction z
    TrisoupSegmentEnc segWW0WWW =  // top right edge
      {leaf.pos + posWW0,
       leaf.pos + posWWW,
       12 * i + 6,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction z
    TrisoupSegmentEnc segW00W0W =  // bottom right edge
      {leaf.pos + posW00,
       leaf.pos + posW0W,
       12 * i + 7,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction z
    TrisoupSegmentEnc seg00WW0W =  // near bottom edge
      {leaf.pos + pos00W,
       leaf.pos + posW0W,
       12 * i + 8,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction x
    TrisoupSegmentEnc seg00W0WW =  // near left edge
      {leaf.pos + pos00W,
       leaf.pos + pos0WW,
       12 * i + 9,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction y
    TrisoupSegmentEnc seg0WWWWW =  // near top edge
      {leaf.pos + pos0WW,
       leaf.pos + posWWW,
       12 * i + 10,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction x
    TrisoupSegmentEnc segW0WWWW =  // near right edge
      {leaf.pos + posW0W,
       leaf.pos + posWWW,
       12 * i + 11,
       -1,
       -1,
       0,
       0,
       0,
       0};  // // direction y

    // Each voxel votes for a position along each edge it is close to
    const int tmin = 1;
    const int tmax = blockWidth - tmin - 1;
    const int tmin2 = distanceSearchEncoder;
    const int tmax2 = blockWidth - tmin2 - 1;

    int idxStart = leaf.start;
    int idxEnd = leaf.end;
    if (isCompensated) {
      idxStart = leaf.predStart;
      idxEnd = leaf.predEnd;
    }


    for (int j = idxStart; j < idxEnd; j++) {
      Vec3<int> voxel = pointCloud[j] - leaf.pos;

      // parameter indicating threshold of how close voxels must be to edge ----------- 1 -------------------
      // to be relevant
      if (voxel[1] < tmin && voxel[2] < tmin) {
        seg000W00.count++;
        seg000W00.distanceSum += voxel[0];
      }  // far bottom edge
      if (voxel[0] < tmin && voxel[2] < tmin) {
        seg0000W0.count++;
        seg0000W0.distanceSum += voxel[1];
      }  // far left edge
      if (voxel[1] > tmax && voxel[2] < tmin) {
        seg0W0WW0.count++;
        seg0W0WW0.distanceSum += voxel[0];
      }  // far top edge
      if (voxel[0] > tmax && voxel[2] < tmin) {
        segW00WW0.count++;
        segW00WW0.distanceSum += voxel[1];
      }  // far right edge
      if (voxel[0] < tmin && voxel[1] < tmin) {
        seg00000W.count++;
        seg00000W.distanceSum += voxel[2];
      }  // bottom left edge
      if (voxel[0] < tmin && voxel[1] > tmax) {
        seg0W00WW.count++;
        seg0W00WW.distanceSum += voxel[2];
      }  // top left edge
      if (voxel[0] > tmax && voxel[1] > tmax) {
        segWW0WWW.count++;
        segWW0WWW.distanceSum += voxel[2];
      }  // top right edge
      if (voxel[0] > tmax && voxel[1] < tmin) {
        segW00W0W.count++;
        segW00W0W.distanceSum += voxel[2];
      }  // bottom right edge
      if (voxel[1] < tmin && voxel[2] > tmax) {
        seg00WW0W.count++;
        seg00WW0W.distanceSum += voxel[0];
      }  // near bottom edge
      if (voxel[0] < tmin && voxel[2] > tmax) {
        seg00W0WW.count++;
        seg00W0WW.distanceSum += voxel[1];
      }  // near left edge
      if (voxel[1] > tmax && voxel[2] > tmax) {
        seg0WWWWW.count++;
        seg0WWWWW.distanceSum += voxel[0];
      }  // near top edge
      if (voxel[0] > tmax && voxel[2] > tmax) {
        segW0WWWW.count++;
        segW0WWWW.distanceSum += voxel[1];
      }  // near right edge

      // parameter indicating threshold of how close voxels must be to edge ----------- 2 -------------------
      // to be relevant
      if (voxel[1] < tmin2 && voxel[2] < tmin2) {
        seg000W00.count2++;
        seg000W00.distanceSum2 += voxel[0];
      }  // far bottom edge
      if (voxel[0] < tmin2 && voxel[2] < tmin2) {
        seg0000W0.count2++;
        seg0000W0.distanceSum2 += voxel[1];
      }  // far left edge
      if (voxel[1] > tmax2 && voxel[2] < tmin2) {
        seg0W0WW0.count2++;
        seg0W0WW0.distanceSum2 += voxel[0];
      }  // far top edge
      if (voxel[0] > tmax2 && voxel[2] < tmin2) {
        segW00WW0.count2++;
        segW00WW0.distanceSum2 += voxel[1];
      }  // far right edge
      if (voxel[0] < tmin2 && voxel[1] < tmin2) {
        seg00000W.count2++;
        seg00000W.distanceSum2 += voxel[2];
      }  // bottom left edge
      if (voxel[0] < tmin2 && voxel[1] > tmax2) {
        seg0W00WW.count2++;
        seg0W00WW.distanceSum2 += voxel[2];
      }  // top left edge
      if (voxel[0] > tmax2 && voxel[1] > tmax2) {
        segWW0WWW.count2++;
        segWW0WWW.distanceSum2 += voxel[2];
      }  // top right edge
      if (voxel[0] > tmax2 && voxel[1] < tmin2) {
        segW00W0W.count2++;
        segW00W0W.distanceSum2 += voxel[2];
      }  // bottom right edge
      if (voxel[1] < tmin2 && voxel[2] > tmax2) {
        seg00WW0W.count2++;
        seg00WW0W.distanceSum2 += voxel[0];
      }  // near bottom edge
      if (voxel[0] < tmin2 && voxel[2] > tmax2) {
        seg00W0WW.count2++;
        seg00W0WW.distanceSum2 += voxel[1];
      }  // near left edge
      if (voxel[1] > tmax2 && voxel[2] > tmax2) {
        seg0WWWWW.count2++;
        seg0WWWWW.distanceSum2 += voxel[0];
      }  // near top edge
      if (voxel[0] > tmax2 && voxel[2] > tmax2) {
        segW0WWWW.count2++;
        segW0WWWW.distanceSum2 += voxel[1];
      }  // near right edge
    }

    // Push segments onto list.
    segments.push_back(seg000W00);  // far bottom edge
    segments.push_back(seg0000W0);  // far left edge
    segments.push_back(seg0W0WW0);  // far top edge
    segments.push_back(segW00WW0);  // far right edge
    segments.push_back(seg00000W);  // bottom left edge
    segments.push_back(seg0W00WW);  // top left edge
    segments.push_back(segWW0WWW);  // top right edge
    segments.push_back(segW00W0W);  // bottom right edge
    segments.push_back(seg00WW0W);  // near bottom edge
    segments.push_back(seg00W0WW);  // near left edge
    segments.push_back(seg0WWWWW);  // near top edge
    segments.push_back(segW0WWWW);  // near right edge
  }

  // Sort the list and find unique segments.
  std::sort(segments.begin(), segments.end());

  TrisoupSegmentEnc localSegment = segments[0];
  auto it = segments.begin() + 1;
  int i = 0;
  for (; it != segments.end(); it++) {
    if (
      localSegment.startpos != it->startpos
      || localSegment.endpos != it->endpos) {
      // Segment[i] is different from localSegment
      // Start a new uniqueSegment.
      segind.push_back(localSegment.count > 0 || localSegment.count2 > 1);
      if (segind.back()) {  // intersects the surface
        int temp = ((2 * localSegment.distanceSum + localSegment.distanceSum2)
                    << (10 - bitDropped))
          / (2 * localSegment.count + localSegment.count2);
        int8_t vertex = (temp + (1 << 9 - bitDropped)) >> 10;
        vertices.push_back(vertex);
      }
      localSegment = *it;  // unique segment
    } else {
      // Segment[i] is the same as localSegment
      // Accumulate
      localSegment.count += it->count;
      localSegment.distanceSum += it->distanceSum;
      localSegment.count2 += it->count2;
      localSegment.distanceSum2 += it->distanceSum2;
    }
  }
  segind.push_back(localSegment.count > 0 || localSegment.count2 > 1);
  if (segind.back()) {  // intersects the surface
    int temp = ((2 * localSegment.distanceSum + localSegment.distanceSum2)
                << (10 - bitDropped))
      / (2 * localSegment.count + localSegment.count2);
    int8_t vertex = (temp + (1 << 9 - bitDropped)) >> 10;
    vertices.push_back(vertex);
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

    int towardOrAway[18] = { 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; // 0 = toward; 1= away
    int mapping18to9[3][9] = { {0,1,2,3,4,15,14,5,7}, {0,1,2,3,9,15,14,7,12}, {0,1,2,9,10,15,14,7,12} };  

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

    arithmeticEncoder->encode((int)segind[i], ctxtMemOctree.ctxTriSoup[0][ctxInter][ctxtMemOctree.MapOBUFTriSoup[ctxInter][0].getEvolve(segind[i], ctxMap2, ctxMap1)]);

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
      arithmeticEncoder->encode(bit, ctxtMemOctree.ctxTriSoup[1][ctxInter][ctxtMemOctree.MapOBUFTriSoup[ctxInter][1].getEvolve(bit, ctxMap2, ctxMap1)]);
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
        arithmeticEncoder->encode(bit, ctxtMemOctree.ctxTriSoup[2][ctxInter][ctxtMemOctree.MapOBUFTriSoup[ctxInter][2].getEvolve(bit, ctxMap2, (ctxMap1 << 1) + v)]);
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
