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

#include <cstdio>

#include "geometry_trisoup.h"

#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"

namespace pcc {

//============================================================================

// The number of fractional bits used in trisoup triangle voxelisation
const int kTrisoupFpBits = 8;

// The value 1 in fixed-point representation
const int kTrisoupFpOne = 1 << (kTrisoupFpBits);
const int kTrisoupFpHalf = 1 << (kTrisoupFpBits - 1);
const int truncateValue = kTrisoupFpHalf;

//============================================================================

void
decodeGeometryTrisoup(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMemOctree,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const Vec3<int> minimum_position)
{
  // trisoup uses octree coding until reaching the triangulation level.
  // todo(df): pass trisoup node size rather than 0?
  std::vector<PCCOctree3Node> nodes;
  PCCPointSet3 compensatedPointCloud;  // set of points after motion compensation
  decodeGeometryOctree(
    gps, gbh, 0, pointCloud, ctxtMemOctree, arithmeticDecoder, &nodes,
    refFrame, minimum_position, compensatedPointCloud);

  std::cout << "\nSize compensatedPointCloud for TriSoup = " << compensatedPointCloud.getPointCount() << "\n";
  bool isInter = gbh.interPredictionEnabledFlag;

  int blockWidth = 1 << gbh.trisoupNodeSizeLog2(gps);
  const int maxVertexPrecisionLog2 = gbh.trisoup_vertex_quantization_bits ? gbh.trisoup_vertex_quantization_bits : gbh.trisoupNodeSizeLog2(gps);
  const int bitDropped =  std::max(0, gbh.trisoupNodeSizeLog2(gps) - maxVertexPrecisionLog2);
  const bool isCentroidDriftActivated = gbh.trisoup_centroid_vertex_residual_flag;

  // Determine neighbours
  std::vector<uint16_t> neighbNodes;
  std::vector<std::array<int, 18>> edgePattern;
  std::vector<int> segmentUniqueIndex;
  int Nunique;
  determineTrisoupNeighbours(nodes, neighbNodes, edgePattern, blockWidth, segmentUniqueIndex, Nunique);

  // determine vertices from compensated point cloud
  std::cout << "Number of points for TriSoup = " << pointCloud.getPointCount() << "\n";
  std::cout << "Number of nodes for TriSoup = " << nodes.size() << "\n";

  std::vector<bool> segindPred;
  std::vector<uint8_t> verticesPred;
  if (isInter) {
    determineTrisoupVertices(
      nodes, segindPred, verticesPred, refFrame->cloud, compensatedPointCloud,
      gps, gbh, blockWidth, bitDropped, 1 /*distanceSearchEncoder*/, true, segmentUniqueIndex, Nunique);
  }

  // Decode vertex presence and position into bitstream
  std::vector<bool> segind;
  std::vector<uint8_t> vertices;
  decodeTrisoupVertices(segind, vertices, segindPred, verticesPred, neighbNodes, edgePattern, bitDropped, gps, gbh, arithmeticDecoder, ctxtMemOctree);

  PCCPointSet3 recPointCloud;
  recPointCloud.addRemoveAttributes(pointCloud);

  // Compute refinedVertices.
  int32_t maxval = (1 << gbh.maxRootNodeDimLog2) - 1;
  bool haloFlag = gbh.trisoup_halo_flag;
  bool adaptiveHaloFlag = gbh.trisoup_adaptive_halo_flag;
  bool fineRayFlag = gbh.trisoup_fine_ray_tracing_flag;
  int thickness = gbh.trisoup_thickness;

  decodeTrisoupCommon(
    nodes, segind, vertices, pointCloud, recPointCloud,
    compensatedPointCloud, gps, gbh, blockWidth, maxval,
    bitDropped, isCentroidDriftActivated, true, haloFlag, adaptiveHaloFlag, fineRayFlag, thickness,
    &arithmeticDecoder,  NULL, ctxtMemOctree, segmentUniqueIndex);

  pointCloud.resize(0);
  pointCloud = std::move(recPointCloud);

  if (!(gps.localMotionEnabled && gps.gof_geom_entropy_continuation_enabled_flag) && !gbh.entropy_continuation_flag) {
    ctxtMemOctree.clearMap();
  }
}


//============================================================================

struct RasterScanTrisoupEdges {
  const int32_t blockWidth;
  // Eight corners of block.
  const Vec3<int32_t> pos000 { 0, 0, 0 };
  const Vec3<int32_t> posW00 { blockWidth, 0, 0 };
  const Vec3<int32_t> pos0W0 { 0, blockWidth, 0 };
  const Vec3<int32_t> posWW0 { blockWidth, blockWidth, 0 };
  const Vec3<int32_t> pos00W { 0, 0, blockWidth };
  const Vec3<int32_t> posW0W { blockWidth, 0, blockWidth };
  const Vec3<int32_t> pos0WW { 0, blockWidth, blockWidth };
  const Vec3<int32_t> posWWW { blockWidth, blockWidth, blockWidth };

  const Vec3<int32_t> offsets[8] = {
    { -blockWidth, -blockWidth, 0 },//-posWW0, // left-bottom
    { -blockWidth, 0, -blockWidth },//-posW0W, // left-front
    { -blockWidth, 0, 0 },//-posW00, // left
    { 0, -blockWidth, -blockWidth },//-pos0WW, // bottom-front
    { 0, -blockWidth, 0 },//-pos0W0, // bottom
    { 0, 0, -blockWidth },//-pos00W, // front
    { 0, 0, 0 },// pos000, // current
    { -blockWidth, -blockWidth, -blockWidth },//-posWWW, // left-bottom-front only useful for neighbors determination
  };

  std::array<int, 8> edgesNeighNodes; // neighboring nodes' index
  // The 7 firsts are used for unique segments generation/iteration
  // The 8-th is used for contextual information (mask)
  Vec3<int32_t> currWedgePos;

  const std::vector<PCCOctree3Node>& leaves;

  RasterScanTrisoupEdges(const std::vector<PCCOctree3Node>& leaves, int blockWidth)
  : leaves(leaves)
  , blockWidth(blockWidth)
  , currWedgePos(leaves.empty() ? Vec3<int32_t>{0,0,0} : leaves[0].pos)
  , edgesNeighNodes {0,0,0,0,0,0,0,0}
  {}

  void buildSegments(
      std::vector<uint16_t>& neighbNodes,
      std::vector<std::array<int, 18>>& edgePattern,
      std::vector<int>& segmentUniqueIndex,
      int& numUniqueIndexes)
  {
    neighbNodes.reserve(leaves.size() * 12); // at most 12 edges per node (to avoid reallocations)
    segmentUniqueIndex.clear();
    segmentUniqueIndex.resize(12 * leaves.size(), -1); // temporarily set to -1 to check everithing is working
    // TODO: set to -1 could be removed when everything will work properly

    int uniqueIndex = 0;

    while (nextIsAvailable()) {
      // process current wedge position

      std::array<bool, 8> processedNode;
      for (int i=0; i<8; ++i)
        processedNode[i] =
          edgesNeighNodes[i] < leaves.size()
          && currWedgePos + offsets[i] == leaves[edgesNeighNodes[i]].pos;

      // edge along z, then y, then x

      // index in edgesNeighNodes of neighboring nodes for each edge
      static const int edgesNeighNodesIdx[3][4] = {
        {6, 4, 0, 2}, // along z
        {6, 2, 5, 1}, // along y
        {6, 4, 5, 3}, // along x
      };

      // edge index in existing order
      static const int edgeIdx[3][4] = {
        {4, 5, 6, 7},
        {1, 3, 9, 11},
        {0, 2, 8, 10},
      };

      static const uint16_t neighMask[3][4] = {
        {0x4001, 0x4002, 0x4004, 0x4008},
        {0x2001, 0x2002, 0x2004, 0x2008},
        {0x0001, 0x0002, 0x0004, 0x0008},
      };

      // index in edgesNeighNodes of other neighboring nodes for the wedge
      static const int wedgeNeighNodesIdx[3][4] = {
        {1, 5, 3, 7}, // along z
        {0, 3, 4, 7}, // along y
        {0, 1, 2, 7}, // along x
      };

      // edge index in existing order
      static const int wedgeNeighNodesEdgeIdx[3][4] = {
        {7, 4, 5,  6},
        {3, 9, 1, 11},
        {2, 8, 0, 10},
      };

      static const uint16_t wedgeNeighMask[3][4] = {
        {0x0800, 0x0100, 0x0200, 0x0400},
        {0x0200, 0x0400, 0x0100, 0x0800},
        {0x0200, 0x0400, 0x0100, 0x0800},
      };

      static const uint16_t toPrevEdgeNeighMask[3][4] = {
        {0x0010, 0x0020, 0x0040, 0x0080},
        {0x0010, 0x0020, 0x0040, 0x0080},
        {0x0010, 0x0020, 0x0040, 0x0080},
      };

      // neighbourhood staic tables
      // ---------    8-bit pattern = 0 before, 1-4 perp, 5-12 others
      static const int localEdgeindex[12][11] = {
        { 4,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 0
        { 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 1
        { 1,  5,  4,  9,  0,  8, -1, -1, -1, -1, -1}, // vertex 2
        { 0,  7,  4,  8,  2, 10,  1,  9, -1, -1, -1}, // vertex 3
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 4
        { 1,  0,  9,  4, -1, -1, -1, -1, -1, -1, -1}, // vertex 5
        { 3,  2,  0, 10, 11,  9,  8,  7,  5,  4, -1}, // vertex 6
        { 0,  1,  2,  8, 10,  4,  5, -1, -1, -1, -1}, // vertex 7
        { 4,  9,  1,  0, -1, -1, -1, -1, -1, -1, -1}, // vertex 8
        { 4,  0,  1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 9
        { 5,  9,  1,  2,  8,  0, -1, -1, -1, -1, -1}, // vertex 10
        { 7,  8,  0, 10,  5,  2,  3,  9,  1, -1, -1}  // vertex 11
      };
      static const int patternIndex[12][11] = {
        { 3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 0
        { 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 1
        { 2,  3,  5,  8, 15, 17, -1, -1, -1, -1, -1}, // vertex 2
        { 2,  3,  5,  8,  9, 12, 15, 17, -1, -1, -1}, // vertex 3
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 4
        { 1,  7, 10, 14, -1, -1, -1, -1, -1, -1, -1}, // vertex 5
        { 1,  2,  6,  9, 10, 11, 13, 14, 15, 16, -1}, // vertex 6
        { 2,  5,  8,  9, 12, 15, 17, -1, -1, -1, -1}, // vertex 7
        { 1,  4,  7, 14, -1, -1, -1, -1, -1, -1, -1}, // vertex 8
        { 1,  7, 14, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 9
        { 1,  2,  6, 14, 15, 16, -1, -1, -1, -1, -1}, // vertex 10
        { 1,  2,  6,  9, 11, 13, 14, 15, 16, -1, -1}  // vertex 11
      };

      for (int dir=0; dir<3; ++dir) {
        bool processedEdge = false;
        uint16_t neighboursMask = 0;
        std::array<int, 18> pattern {
          -1, -1, -1, -1, -1, -1, -1, -1, -1,
          -1, -1, -1, -1, -1, -1, -1, -1, -1
        };
        for (int neighIdx=0; neighIdx<4; ++neighIdx) {
          int edgeNeighNodeIdx = edgesNeighNodesIdx[dir][neighIdx];
          if (processedNode[edgeNeighNodeIdx]) {
            processedEdge = true;
            int idx =
              edgesNeighNodes[edgeNeighNodeIdx] * 12
              + edgeIdx[dir][neighIdx];
            segmentUniqueIndex[idx] = uniqueIndex;

            // update mask from nodes touching edge
            neighboursMask |= neighMask[dir][neighIdx];

            int indexLow = edgeIdx[dir][neighIdx];
            for (int v = 0; v < 11; v++) {
              if (localEdgeindex[indexLow][v] == -1)
                break;

              int indexV = edgesNeighNodes[edgeNeighNodeIdx] * 12 + localEdgeindex[indexLow][v]; // index of segment
              int Vidx = segmentUniqueIndex[indexV];
              assert(Vidx != -1); // check if already coded
              pattern[patternIndex[indexLow][v]] = Vidx;
            }
          }
        }
        if (processedEdge) {
          int segmentUniqueIdxPrevEdge = -1;
          for (int prevNeighIdx=0; prevNeighIdx<4; ++prevNeighIdx) {
            int wedgeNeighNodeIdx = wedgeNeighNodesIdx[dir][prevNeighIdx];
            if (processedNode[wedgeNeighNodeIdx]) {
              // update current mask from nodes touching wedge
              neighboursMask |= wedgeNeighMask[dir][prevNeighIdx];
              if (segmentUniqueIdxPrevEdge == -1) {
                int idx = edgesNeighNodes[wedgeNeighNodeIdx] * 12
                  + wedgeNeighNodesEdgeIdx[dir][prevNeighIdx];
                segmentUniqueIdxPrevEdge = segmentUniqueIndex[idx];
                pattern[0] = segmentUniqueIdxPrevEdge;
                assert(segmentUniqueIdxPrevEdge != -1);
                for (int neighIdx=0; neighIdx<4; ++neighIdx) {
                  int edgeNeighNodeIdx = edgesNeighNodesIdx[dir][neighIdx];
                  if (processedNode[edgeNeighNodeIdx]) {
                    neighbNodes[segmentUniqueIdxPrevEdge] |= toPrevEdgeNeighMask[dir][neighIdx];
                  }
                }
              }
            }
          }
          ++uniqueIndex;
          neighbNodes.push_back(neighboursMask);
          edgePattern.push_back(pattern);
        }
      }

      goNextWedge(processedNode);
    }
    numUniqueIndexes = uniqueIndex;
  }
private:
  bool nextIsAvailable() const { return edgesNeighNodes[0] < leaves.size(); }

  void goNextWedge(const std::array<bool, 8>& processedNode) {
    bool nextWedgeIsNext = false;
    if (processedNode[0])
      edgesNeighNodes[7] = edgesNeighNodes[0];
    // handle neighboring nodes with same z as wedge
    for (int i=0; i<7; i+=2)
      if (processedNode[i]) {
        nextWedgeIsNext = true;
        edgesNeighNodes[i]++;
      }
    for (int i=1; i<7; i+=2)
      if (processedNode[i]) {
        edgesNeighNodes[i]++;
      }
    if (nextWedgeIsNext) {
      currWedgePos += pos00W;
    }
    else {
      currWedgePos = leaves[edgesNeighNodes[0]].pos - offsets[0];
      for (int i=1; i<7; ++i) {
        if (edgesNeighNodes[i] >= leaves.size())
          break;
        auto wedgePos = leaves[edgesNeighNodes[i]].pos - offsets[i];
        if (currWedgePos > wedgePos) {
          currWedgePos = wedgePos;
        }
      }
    }
  }

};

void determineTrisoupNeighbours(
  const std::vector<PCCOctree3Node>& leaves,
  std::vector<uint16_t>& neighbNodes,
  std::vector<std::array<int, 18>>& edgePattern,
  const int defaultBlockWidth,
  std::vector<int>& segmentUniqueIndex,
  int& Nunique) {

  // Width of block.
  // in future, may override with leaf blockWidth
  const int32_t blockWidth = defaultBlockWidth;

  RasterScanTrisoupEdges rste(leaves, blockWidth);
  rste.buildSegments(neighbNodes, edgePattern, segmentUniqueIndex, Nunique);
}
//============================================================================

template<typename T>
Vec3<T>
crossProduct(const Vec3<T> a, const Vec3<T> b)
{
  Vec3<T> ret;
  ret[0] = a[1] * b[2] - a[2] * b[1];
  ret[1] = a[2] * b[0] - a[0] * b[2];
  ret[2] = a[0] * b[1] - a[1] * b[0];
  return ret;
}

//---------------------------------------------------------------------------

Vec3<int32_t>
truncate(const Vec3<int32_t> in, const int32_t offset)
{
  Vec3<int32_t> out = in + offset;  
  if (out[0] < 0)
    out[0] = 0;
  if (out[1] < 0)
    out[1] = 0;
  if (out[2] < 0)
    out[2] = 0;

  return out;
}

//---------------------------------------------------------------------------
bool
boundaryinsidecheck(const Vec3<int32_t> a, const int bbsize)
{
  return a[0] >= 0 && a[0] <= bbsize && a[1] >= 0 && a[1] <= bbsize
    && a[2] >= 0 && a[2] <= bbsize;
}

//---------------------------------------------------------------------------

void nonCubicNode
(
 const GeometryParameterSet& gps,
 const GeometryBrickHeader& gbh,
 const Vec3<int32_t>& leafpos,
 const int32_t blockWidth,
 const Box3<int32_t>& bbox,
 Vec3<int32_t>& newp,
 Vec3<int32_t>& neww,
 Vec3<int32_t>* corner )
{
  bool flag_n = gps.non_cubic_node_start_edge && ( gbh.slice_bb_pos_bits   > 0 );
  bool flag_f = gps.non_cubic_node_end_edge   && ( gbh.slice_bb_width_bits > 0 );
  for( int k=0; k<3; k++ ) {
    newp[k] = ( ( flag_n ) && ( leafpos[k] < bbox.min[k] ) ) ? bbox.min[k] : leafpos[k];
    neww[k] = ( ( flag_n ) && ( leafpos[k] < bbox.min[k] ) ) ?
      (blockWidth-(bbox.min[k]-leafpos[k])) :
      ( flag_f ) ? std::min(bbox.max[k]-leafpos[k]+1, blockWidth) : blockWidth;
  }
  corner[POS_000] = {       0,       0,       0 };
  corner[POS_W00] = { neww[0],       0,       0 };
  corner[POS_0W0] = {       0, neww[1],       0 };
  corner[POS_WW0] = { neww[0], neww[1],       0 };
  corner[POS_00W] = {       0,       0, neww[2] };
  corner[POS_W0W] = { neww[0],       0, neww[2] };
  corner[POS_0WW] = {       0, neww[1], neww[2] };
  corner[POS_WWW] = { neww[0], neww[1], neww[2] };
  return;
}


// --------------------------------------------------------------------------
void
determineCentroidAndDominantAxis(
  Vec3<int32_t> &blockCentroid,
  int &dominantAxis,
  std::vector<Vertex> &leafVertices,
  Vec3<int32_t> nodew)

{
  // compute centroid
  int triCount = (int)leafVertices.size();
  blockCentroid = 0;
  for (int j = 0; j < triCount; j++) {
    blockCentroid += leafVertices[j].pos;
  }
  blockCentroid /= triCount;

  // order vertices along a dominant axis only if more than three (otherwise only one triangle, whatever...)
  dominantAxis = findDominantAxis(leafVertices, nodew, blockCentroid);

  // compute centroid
  std::vector<int> Weigths(leafVertices.size(), 0);
  int Wtotal = 0;
  for (int k = 0; k < triCount; k++) {
    int k2 = k + 1;
    if (k2 >= triCount)
       k2 -= triCount;
    Vec3<int32_t> segment = (leafVertices[k].pos - leafVertices[k2].pos).abs();
    int weight = segment[0] + segment[1] + segment[2];

    Weigths[k] += weight;
    Weigths[k2] += weight;
    Wtotal += 2 * weight;
  }

  Vec3<int64_t> blockCentroid2 = 0;
  for (int j = 0; j < triCount; j++) {
    blockCentroid2 += int64_t(Weigths[j]) * leafVertices[j].pos;
  }
  blockCentroid2 /= int64_t(Wtotal);
  blockCentroid = { int(blockCentroid2[0]),int(blockCentroid2[1]), int(blockCentroid2[2]) };
}

// --------------------------------------------------------------------------
//  encoding of centroid residual in TriSOup node
Vec3<int32_t>
determineCentroidNormalAndBounds(
  int& lowBound,
  int& highBound,
  int& lowBoundSurface,
  int& highBoundSurface,
  int & ctxMinMax,
  int bitDropped,
  int bitDropped2,
  int triCount,
  Vec3<int32_t> blockCentroid,
  int dominantAxis,
  std::vector<Vertex>& leafVertices,
  int nodewDominant)
{
  int halfDropped2 = bitDropped2 == 0 ? 0 : 1 << bitDropped2 - 1;

  // contextual information  for drift coding
  int minPos = leafVertices[0].pos[dominantAxis];
  int maxPos = leafVertices[0].pos[dominantAxis];
  for (int k = 1; k < triCount; k++) {
    if (leafVertices[k].pos[dominantAxis] < minPos)
      minPos = leafVertices[k].pos[dominantAxis];
    if (leafVertices[k].pos[dominantAxis] > maxPos)
      maxPos = leafVertices[k].pos[dominantAxis];
  }

  // find normal vector
  Vec3<int64_t> accuNormal = 0;
  for (int k = 0; k < triCount; k++) {
    int k2 = k + 1;
    if (k2 >= triCount)
      k2 -= triCount;
    accuNormal += crossProduct(leafVertices[k].pos - blockCentroid, leafVertices[k2].pos - blockCentroid);
  }
  int64_t normN = isqrt(accuNormal[0] * accuNormal[0] + accuNormal[1] * accuNormal[1] + accuNormal[2] * accuNormal[2]);
  Vec3<int32_t> normalV = ((accuNormal << kTrisoupFpBits) / normN);

  // drift bounds
  ctxMinMax = std::min(8, (maxPos - minPos) >> (kTrisoupFpBits + bitDropped));
  int bound = (nodewDominant - 1) << kTrisoupFpBits;
  int m = 1;
  for (; m < nodewDominant; m++) {
    Vec3<int32_t> temp = blockCentroid + m * normalV;
    if (temp[0]<0 || temp[1]<0 || temp[2]<0 || temp[0]>bound || temp[1]>bound || temp[2]> bound)
      break;
  }
  highBound = (m - 1) + halfDropped2 >> bitDropped2;

  m = 1;
  for (; m < nodewDominant; m++) {
    Vec3<int32_t> temp = blockCentroid - m * normalV;
    if (temp[0]<0 || temp[1]<0 || temp[2]<0 || temp[0]>bound || temp[1]>bound || temp[2]> bound)
      break;
  }
  lowBound = (m - 1) + halfDropped2 >> bitDropped2;
  lowBoundSurface = std::max(0, ((blockCentroid[dominantAxis] - minPos) + kTrisoupFpHalf >> kTrisoupFpBits) + halfDropped2 >> bitDropped2);
  highBoundSurface = std::max(0, ((maxPos - blockCentroid[dominantAxis]) + kTrisoupFpHalf >> kTrisoupFpBits) + halfDropped2 >> bitDropped2);

  return normalV;
}

// --------------------------------------------------------------------------
int
determineCentroidPredictor(
  int bitDropped2,
  Vec3<int32_t> normalV,
  Vec3<int32_t> blockCentroid,
  Vec3<int32_t> nodepos,
  PCCPointSet3& compensatedPointCloud,
  int start,
  int end,
  int lowBound,
  int  highBound)
{

  int driftQPred = -100;
  // determine quantized drift for predictor
  if (end > start) {
    int driftPred = 0;
    driftQPred = 0;
    int counter = 0;
    int maxD = bitDropped2;

    for (int p = start; p < end; p++) {
      auto point = (compensatedPointCloud[p] - nodepos) << kTrisoupFpBits;

      Vec3<int32_t> CP = crossProduct(normalV, point - blockCentroid) >> kTrisoupFpBits;
      int dist = std::max(std::max(std::abs(CP[0]), std::abs(CP[1])), std::abs(CP[2]));
      dist >>= kTrisoupFpBits;

      if (dist <= maxD) {
        int w = 1 + 4 * (maxD - dist);
        counter += w;
        driftPred += w * ((normalV * (point - blockCentroid)) >> kTrisoupFpBits);
      }
    }

    if (counter) { // drift is shift by kTrisoupFpBits
      driftPred = (driftPred >> kTrisoupFpBits - 6) / counter; // drift is shift by 6 bits
    }

    int half = 1 << 5 + bitDropped2;
    int DZ = 2 * half / 3;

    if (abs(driftPred) >= DZ) {
      driftQPred = (abs(driftPred) - DZ + 2 * half) >> 6 + bitDropped2 - 1;
      if (driftPred < 0)
        driftQPred = -driftQPred;
    }
    driftQPred = std::min(std::max(driftQPred, -2 * lowBound), 2 * highBound);  // drift in [-lowBound; highBound] but quantization is twice better

  }

  return driftQPred;
}

// --------------------------------------------------------------------------
int
determineCentroidResidual(
  int bitDropped2,
  Vec3<int32_t> normalV,
  Vec3<int32_t> blockCentroid,
  Vec3<int32_t> nodepos,
  PCCPointSet3& pointCloud,
  int start,
  int end,
  int lowBound,
  int  highBound)
{
  // determine quantized drift
  int counter = 0;
  int drift = 0;
  int maxD = bitDropped2;

  // determine quantized drift
  for (int p = start; p < end; p++) {
    auto point = (pointCloud[p] - nodepos) << kTrisoupFpBits;

    Vec3<int64_t> CP = crossProduct(normalV, point - blockCentroid) >> kTrisoupFpBits;
    int64_t dist = isqrt(CP[0] * CP[0] + CP[1] * CP[1] + CP[2] * CP[2]);
    dist >>= kTrisoupFpBits;

    if ((dist << 10) <= 1774 * maxD) {
      int32_t w = (1 << 10) + 4 * (1774 * maxD - ((1 << 10) * dist));
      counter += w >> 10;
      drift += (w >> 10) * ((normalV * (point - blockCentroid)) >> kTrisoupFpBits);
    }
  }

  if (counter) { // drift is shift by kTrisoupFpBits
    drift = (drift >> kTrisoupFpBits - 6) / counter; // drift is shift by 6 bits
  }

  int half = 1 << 5 + bitDropped2;
  int DZ = 2 * half / 3;

  int driftQ = 0;
  if (abs(drift) >= DZ) {
    //driftQ = (abs(drift) - DZ + 2 * half) >> 6 + bitDropped2;
    driftQ = (abs(drift) + 2 * half) >> 6 + bitDropped2;
    if (drift < 0)
      driftQ = -driftQ;
  }
  driftQ = std::min(std::max(driftQ, -lowBound), highBound);  // drift in [-lowBound; highBound]

  return driftQ;
}


// --------------------------------------------------------------------------
//  encoding of centroid residual in TriSOup node
void
encodeCentroidResidual(
  int driftQ,
  pcc::EntropyEncoder* arithmeticEncoder,
  GeometryOctreeContexts & ctxtMemOctree,
  int driftQPred,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound)
{
  if (driftQPred == -100) //intra
    arithmeticEncoder->encode(driftQ == 0, ctxtMemOctree.ctxDrift0[ctxMinMax][0]);
  else //inter
    arithmeticEncoder->encode(driftQ == 0, ctxtMemOctree.ctxDrift0[ctxMinMax][1 + std::min(3, std::abs(driftQPred))]);

  // if not 0, drift in [-lowBound; highBound]
  if (driftQ) {
    // code sign
    int lowS = std::min(7, lowBoundSurface);
    int highS = std::min(7, highBoundSurface);
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
}

//---------------------------------------------------------------------------
//  decoding of centroid residual in TriSOup node
int
decodeCentroidResidual(
  pcc::EntropyDecoder* arithmeticDecoder,
  GeometryOctreeContexts& ctxtMemOctree,
  int driftQPred,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound)
{
  // decode drift
  int driftQ = 0;
  if (driftQPred == -100) //intra
    driftQ = arithmeticDecoder->decode(ctxtMemOctree.ctxDrift0[ctxMinMax][0]) ? 0 : 1;
  else //inter
    driftQ = arithmeticDecoder->decode(ctxtMemOctree.ctxDrift0[ctxMinMax][1 + std::min(3, std::abs(driftQPred))]) ? 0 : 1;

  // if not 0, drift in [-lowBound; highBound]
  if (driftQ) {
    // code sign
    int lowS = std::min(7, lowBoundSurface);
    int highS = std::min(7, highBoundSurface);

    int sign = 1;
    if (highBound && lowBound) // otherwise sign is knwow
      sign = arithmeticDecoder->decode(ctxtMemOctree.ctxDriftSign[lowBound == highBound ? 0 : 1 + (lowBound < highBound)][lowS][highS][(driftQPred && driftQPred != -100) ? 1 + (driftQPred > 0) : 0]);
    else if (!highBound) // highbound is 0 , so sign is negative; otherwise sign is already set to positive
      sign = 0;

    // code remaining bits 1 to 7 at most
    int magBound = (sign ? highBound : lowBound) - 1;
    bool sameSignPred = driftQPred != -100 && (driftQPred > 0 && sign) || (driftQPred < 0 && !sign);


    int ctx = 0;
    while (magBound > 0) {
      int bit;
      if (ctx < 4)
        bit = arithmeticDecoder->decode(ctxtMemOctree.ctxDriftMag[ctx][driftQPred != -100 ? 1 + std::min(8, sameSignPred * std::abs(driftQPred)) : 0]);
      else
        bit = arithmeticDecoder->decode();

      if (bit) // magDrift==0 and magnitude coding is finished
        break;

      driftQ++;
      magBound--;
      ctx++;
    }

    if (!sign)
      driftQ = -driftQ;
  }
  return driftQ;
}

//---------------------------------------------------------------------------
// Trisoup geometry decoding, at both encoder and decoder.
// Compute from leaves, segment indicators, and vertices
// a set of triangles, refine the triangles, and output their vertices.
//
// @param leaves  list of blocks containing the surface
// @param segind, indicators for edges of blocks if they intersect the surface
// @param vertices, locations of intersections

void
decodeTrisoupCommon(
  const std::vector<PCCOctree3Node>& leaves,
  const std::vector<bool>& segind,
  const std::vector<uint8_t>& vertices,
  PCCPointSet3& pointCloud,
  PCCPointSet3& recPointCloud,
  PCCPointSet3& compensatedPointCloud,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int defaultBlockWidth,
  int poistionClipValue,
  const int bitDropped,
  const bool isCentroidDriftActivated,
  bool isDecoder,
  bool haloFlag,
  bool adaptiveHaloFlag,
  bool fineRayflag,
  int thickness,
  pcc::EntropyDecoder* arithmeticDecoder,
  pcc::EntropyEncoder* arithmeticEncoder,
  GeometryOctreeContexts& ctxtMemOctree,
  std::vector<int>& segmentUniqueIndex)
{
  recPointCloud.resize(0);

  Box3<int32_t> sliceBB;
  sliceBB.min = gbh.slice_bb_pos << gbh.slice_bb_pos_log2_scale;
  sliceBB.max = sliceBB.min + ( gbh.slice_bb_width << gbh.slice_bb_width_log2_scale );

  // Width of block. In future, may override with leaf blockWidth
  const int32_t blockWidth = defaultBlockWidth;

  std::vector<int> vertexList(segind.size());
  int vertexCount = 0;
  for (int i = 0; i < segind.size(); i++) {

    if (segind[i]) {  // intersects the surface
      vertexList[i]= vertices[vertexCount++];
    }
    else {  // does not intersect the surface
      vertexList[i] = -1;
    }
  }

  const int startCorner[12] = { POS_000, POS_000, POS_0W0, POS_W00, POS_000, POS_0W0, POS_WW0, POS_W00, POS_00W, POS_00W, POS_0WW, POS_W0W };
  const int endCorner[12] =   { POS_W00, POS_0W0, POS_WW0, POS_WW0, POS_00W, POS_0WW, POS_WWW, POS_W0W, POS_W0W, POS_0WW, POS_WWW, POS_WWW };

  // ----------- loop on leaf nodes ----------------------
  int idxSegment = 0;
  for (int i = 0; i < leaves.size(); i++) {
    Vec3<int32_t> nodepos, nodew, corner[8];
    nonCubicNode( gps, gbh, leaves[i].pos, defaultBlockWidth, sliceBB, nodepos, nodew, corner );

    // Find up to 12 vertices for this leaf.
    std::vector<Vertex> leafVertices;
    std::vector<Vec3<int32_t>> refinedVerticesBlock;
    refinedVerticesBlock.reserve(blockWidth * blockWidth * 4);

    for (int j = 0; j < 12; j++) {
      int uniqueIndex = segmentUniqueIndex[idxSegment++];
      int vertex = vertexList[uniqueIndex];

      if (vertex < 0)
        continue;  // skip segments that do not intersect the surface
      auto startSegment = corner[startCorner[j]];
      auto endSegment = corner[endCorner[j]];

      // Get distance along edge of vertex.
      // Vertex code is the index of the voxel along the edge of the block
      // of surface intersection./ Put decoded vertex at center of voxel,
      // unless voxel is first or last along the edge, in which case put the
      // decoded vertex at the start or endpoint of the segment.
      Vec3<int32_t> direction = endSegment - startSegment;
      uint32_t segment_len = direction.max();

      // Get 3D position of point of intersection.
      Vec3<int32_t> point = startSegment << kTrisoupFpBits;
      point -= kTrisoupFpHalf; // the volume is [-0.5; B-0.5]^3

      // points on edges are located at integer values
      int32_t distance = (vertex << (kTrisoupFpBits + bitDropped)) + (kTrisoupFpHalf << bitDropped);
      if (direction[0])
        point[0] += distance; // in {0,1,...,B-1}
      else if (direction[1])
        point[1] += distance;
      else  // direction[2] 
        point[2] += distance;

      // Add vertex to list of vertices.
      leafVertices.push_back({ point, 0, 0 });

      // vertex to list of points
      if (bitDropped) {
        Vec3<int32_t> foundvoxel = (point + truncateValue) >> kTrisoupFpBits;
        if (boundaryinsidecheck(foundvoxel, blockWidth - 1))
           refinedVerticesBlock.push_back(nodepos + foundvoxel);
      }
    }

    // Skip leaves that have fewer than 3 vertices.
    int triCount = (int)leafVertices.size();
    if (triCount < 3) {
      std::sort(refinedVerticesBlock.begin(), refinedVerticesBlock.end());
      refinedVerticesBlock.erase(std::unique(refinedVerticesBlock.begin(), refinedVerticesBlock.end()), refinedVerticesBlock.end());

      // Move list of points to pointCloud
      int nPointInCloud = recPointCloud.getPointCount();
      recPointCloud.resize(nPointInCloud + refinedVerticesBlock.size());
      for (int i = 0; i < refinedVerticesBlock.size(); i++)
        recPointCloud[nPointInCloud + i] = refinedVerticesBlock[i];

      continue;
    }

    // compute centroid
    Vec3<int32_t> blockCentroid;
    int dominantAxis;
    determineCentroidAndDominantAxis(blockCentroid, dominantAxis, leafVertices, nodew);

    // Refinement of the centroid along the domiannt axis
    if (triCount > 3 && isCentroidDriftActivated) {
      int bitDropped2 = bitDropped;

      int lowBound, highBound, lowBoundSurface, highBoundSurface, ctxMinMax;
      Vec3<int32_t> normalV = determineCentroidNormalAndBounds(lowBound, highBound, lowBoundSurface, highBoundSurface, ctxMinMax, bitDropped, bitDropped2, triCount, blockCentroid, dominantAxis, leafVertices, nodew[dominantAxis]);

      int driftQPred = determineCentroidPredictor(bitDropped2, normalV, blockCentroid, nodepos, compensatedPointCloud, leaves[i].predStart, leaves[i].predEnd, lowBound, highBound);

      int driftQ = 0;
      if (!isDecoder) { // encode centroid residual
        driftQ = determineCentroidResidual(bitDropped2, normalV, blockCentroid, nodepos, pointCloud, leaves[i].start, leaves[i].end, lowBound, highBound);
        encodeCentroidResidual(driftQ, arithmeticEncoder, ctxtMemOctree, driftQPred, ctxMinMax, lowBoundSurface, highBoundSurface, lowBound, highBound);
      }
      else { // decode centroid residual
        driftQ = decodeCentroidResidual(arithmeticDecoder, ctxtMemOctree, driftQPred, ctxMinMax, lowBoundSurface, highBoundSurface, lowBound, highBound);
      }

      // dequantize and apply drift 
      int driftDQ = 0;
      if (driftQ) {
        driftDQ = std::abs(driftQ) << bitDropped2 + 6;
        int half = 1 << 5 + bitDropped2;
        int DZ = 2*half/3; 
        driftDQ += DZ - half; 
        if (driftQ < 0)
          driftDQ = -driftDQ;
      }

      blockCentroid += (driftDQ * normalV) >> 6;
      blockCentroid[0] = std::max(-kTrisoupFpHalf, blockCentroid[0]);
      blockCentroid[1] = std::max(-kTrisoupFpHalf, blockCentroid[1]);
      blockCentroid[2] = std::max(-kTrisoupFpHalf, blockCentroid[2]);
      blockCentroid[0] = std::min(((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1, blockCentroid[0]);
      blockCentroid[1] = std::min(((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1, blockCentroid[1]);
      blockCentroid[2] = std::min(((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1, blockCentroid[2]);
    } // end refinement of the centroid


    // Divide vertices into triangles around centroid
    // and upsample each triangle by an upsamplingFactor.
    if (triCount > 3) {
      Vec3<int32_t> foundvoxel = (blockCentroid + truncateValue) >> kTrisoupFpBits;
      if (boundaryinsidecheck(foundvoxel, blockWidth - 1))
        refinedVerticesBlock.push_back(nodepos + foundvoxel);
    }

    int haloTriangle = (((1 << bitDropped) - 1) << kTrisoupFpBits) / blockWidth;
    haloTriangle = (haloTriangle * 28) >> 5; // / 32;
    haloTriangle = haloTriangle > 36 ? 36 : haloTriangle;

    Vec3<int32_t> v2 = triCount == 3 ? leafVertices[2].pos : blockCentroid;
    Vec3<int32_t> v1 = leafVertices[0].pos;
    for (int triIndex = 0; triIndex < (triCount == 3 ? 1 : triCount); triIndex++) {
      int j1 = triIndex + 1;
      if (j1 >= triCount)
        j1 -= triCount;

      Vec3<int32_t> v0 = v1;
      v1 = leafVertices[j1].pos;

      // range      
      int minRange[3];
      int maxRange[3];
      for (int k = 0; k < 3; k++) {
        minRange[k] = std::max(0, std::min(std::min(v0[k], v1[k]), v2[k]) + truncateValue >> kTrisoupFpBits);
        maxRange[k] = std::min(blockWidth - 1, std::max(std::max(v0[k], v1[k]), v2[k]) + truncateValue >> kTrisoupFpBits);
      }

      // choose ray direction
      Vec3<int32_t> edge1 = v1 - v0;
      Vec3<int32_t> edge2 = v2 - v0;
      Vec3<int32_t> h = crossProduct(edge1, edge2) >> kTrisoupFpBits;
      int minDir = std::abs(h[0]);
      int directionOk = 0;
      if (std::abs(h[1]) >= minDir) {
        minDir = std::abs(h[1]);
        directionOk = 1;
      }
      if (std::abs(h[2]) >= minDir) {
        directionOk = 2;
      }

      // applying ray tracing along direction
      rayTracingAlongdirection_samp1_optim(
        refinedVerticesBlock, directionOk, blockWidth, nodepos, minRange,
        maxRange, edge1, edge2, v0, haloTriangle, thickness);      

    }  // end loop on triangles

    // remove points present twice or more for node
    std::sort(refinedVerticesBlock.begin(), refinedVerticesBlock.end());
    refinedVerticesBlock.erase(std::unique(refinedVerticesBlock.begin(), refinedVerticesBlock.end()), refinedVerticesBlock.end());

    // Move list of points to pointCloud
    int nPointInCloud = recPointCloud.getPointCount();
    recPointCloud.resize(nPointInCloud + refinedVerticesBlock.size());
    for (int i = 0; i < refinedVerticesBlock.size(); i++)
      recPointCloud[nPointInCloud+i] = refinedVerticesBlock[i];

  }// end loop on leaves

}


// ---------------------------------------------------------------------------
void decodeTrisoupVertices(
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<bool>& segindPred,
  std::vector<uint8_t>& verticesPred,
  std::vector<uint16_t>& neighbNodes,
  std::vector<std::array<int, 18>>& edgePattern,
  int bitDropped,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  pcc::EntropyDecoder& arithmeticDecoder,
  GeometryOctreeContexts& ctxtMemOctree)
{
  const int nbitsVertices = gbh.trisoupNodeSizeLog2(gps) - bitDropped;
  const int max2bits = nbitsVertices > 1 ? 3 : 1;
  const int mid2bits = nbitsVertices > 1 ? 2 : 1;

  int iV = 0;
  int iVPred = 0;
  std::vector<int> correspondanceSegment2V;

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

    bool c = ctxtMemOctree.MapOBUFTriSoup[ctxInter][0].decodeEvolve(
      &arithmeticDecoder, ctxtMemOctree.ctxTriSoup[0][ctxInter], ctxMap2,
      ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
      ctxtMemOctree._BufferOBUFleavesTrisoup);

    segind.push_back(c);
    correspondanceSegment2V.push_back(-1);

    // encode position vertex 
    if (c) {
      correspondanceSegment2V.back() = iV;

      uint8_t v = 0;
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

      int bit = ctxtMemOctree.MapOBUFTriSoup[ctxInter][1].decodeEvolve(
        &arithmeticDecoder, ctxtMemOctree.ctxTriSoup[1][ctxInter], ctxMap2,
        ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
        ctxtMemOctree._BufferOBUFleavesTrisoup);
      v = (v << 1) | bit;
      b--;

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

        bit = ctxtMemOctree.MapOBUFTriSoup[ctxInter][2].decodeEvolve(
          &arithmeticDecoder, ctxtMemOctree.ctxTriSoup[2][ctxInter], ctxMap2,
          (ctxMap1 << 1) + v, &ctxtMemOctree._OBUFleafNumberTrisoup,
          ctxtMemOctree._BufferOBUFleavesTrisoup);
        v = (v << 1) | bit;
        b--;
      }


      // third bit
      if (b >= 0) {
        int ctxFullNboundsReduced1 = (6 * (ctx0 >> 1) + missedCloseStart) * 2 + (ctxE == 3);
        v = (v << 1) | arithmeticDecoder.decode(ctxtMemOctree.ctxTempV2[4 * ctxFullNboundsReduced1 + v]);
        b--;
      }

      // remaining bits are bypassed
      for (; b >= 0; b--)
        v = (v << 1) | arithmeticDecoder.decode();
      vertices.push_back(v);
      iV++;
    }

    if (isInter && segindPred[i])
      iVPred++;

  }

}


//-----------------------
// Project vertices along dominant axis (i.e., into YZ, XZ, or XY plane).
// Sort projected vertices by decreasing angle in [-pi,+pi] around center
// of block (i.e., clockwise) breaking ties in angle by
// increasing distance along the dominant axis.   

int findDominantAxis(
  std::vector<Vertex>& leafVertices,
  Vec3<uint32_t> blockWidth,
  Vec3<int32_t> blockCentroid ) {

  int dominantAxis = 0;
  int triCount = leafVertices.size();

  auto leafVerticesTemp = leafVertices;

  if (triCount > 3) {
    Vertex vertex;
    Vec3<int32_t> Width = blockWidth << kTrisoupFpBits;

    const int sIdx1[3] = { 2,2,1 };
    const int sIdx2[3] = { 1,0,0 };

    int maxNormTri = 0;
    for (int axis = 0; axis <= 2; axis++) {
      int axis1 = sIdx1[axis];
      int axis2 = sIdx2[axis];
      // order along axis
      for (int j = 0; j < triCount; j++) {
        // compute score closckwise
        int x = leafVerticesTemp[j].pos[axis1] + kTrisoupFpHalf; // back to [0,B]^3 for ordering
        int y = leafVerticesTemp[j].pos[axis2] + kTrisoupFpHalf; // back to [0,B]^3 for ordering
        int Width_x = Width[axis1];
        int Width_y = Width[axis2];

        int flag3 = x <= 0;
        int score3 = Width_y - flag3 * y + (!flag3) * x;
        int flag2 = y >= Width_y;
        int score2 = Width_y + Width_x - flag2 * x + (!flag2) * score3;
        int flag1 = x >= Width_x;
        int score = flag1 * y + (!flag1) * score2;
        leafVerticesTemp[j].theta = score;

        // stable sort if same score
        leafVerticesTemp[j].tiebreaker = leafVerticesTemp[j].pos[axis] + kTrisoupFpHalf;
      }
      std::sort(leafVerticesTemp.begin(), leafVerticesTemp.end(), vertex);

      // compute sum normal
      int32_t accuN = 0;
      for (int k = 0; k < triCount; k++) {
        int k2 = k == triCount-1 ? 0 : k + 1;
        int32_t h = (leafVerticesTemp[k].pos[axis1] - blockCentroid[axis1])* (leafVerticesTemp[k2].pos[axis2] - blockCentroid[axis2]);
        h -= (leafVerticesTemp[k].pos[axis2] - blockCentroid[axis2]) * (leafVerticesTemp[k2].pos[axis1] - blockCentroid[axis1]);
        accuN += std::abs(h);
      }

      // if sumnormal is bigger , this is dominantAxis
      if (accuN > maxNormTri) {
        maxNormTri = accuN;
        dominantAxis = axis;
        leafVertices = leafVerticesTemp;
      }
    }
  } // end find dominant axis

  return dominantAxis;
}





// -------------------------------------------
void rayTracingAlongdirection_samp1_optim(
  std::vector<Vec3<int32_t>>& refinedVerticesBlock,
  int direction,
  int blockWidth,
  Vec3<int32_t> nodepos,
  int minRange[3],
  int maxRange[3],
  Vec3<int32_t> edge1,
  Vec3<int32_t> edge2,
  Vec3<int32_t> Ver0,
  int haloTriangle,
  int thickness) {

  // check if ray tracing is valid; if not skip the direction
  Vec3<int32_t> rayVector = 0;
  rayVector[direction] = 1;// << kTrisoupFpBits;
  Vec3<int32_t> h = crossProduct(rayVector, edge2);// >> kTrisoupFpBits;
  int32_t a = (edge1 * h) >> kTrisoupFpBits; // max is node size square, shifted left by kTrisoupFpBits; max bits = 2*log22Nodesize + kTrisoupFpBits <=2*6 +8 = 20 bits 
  if (std::abs(a) <= kTrisoupFpOne)
    return;

  const int precDivA = 30;
  int64_t inva = (int64_t(1) << precDivA) / a;

  //bounds
  const int g1pos[3] = { 1, 0, 0 };
  const int g2pos[3] = { 2, 2, 1 };
  const int i1 = g1pos[direction];
  const int i2 = g2pos[direction];

  const int32_t startposG1 = minRange[i1];
  const int32_t startposG2 = minRange[i2];
  const int32_t endposG1 = maxRange[i1];
  const int32_t endposG2 = maxRange[i2];

  Vec3<int32_t>  rayOrigin0 = minRange[direction] << kTrisoupFpBits;;
  rayOrigin0[i1] = startposG1 << kTrisoupFpBits;
  rayOrigin0[i2] = startposG2 << kTrisoupFpBits;

  Vec3<int32_t> s0 = rayOrigin0 - Ver0;
  int32_t u0 =  ((s0 * h) * inva) >> precDivA;
  Vec3<int32_t> q0 = crossProduct(s0, edge1);
  int32_t v0 = (q0[direction] * inva) >> precDivA;
  int32_t t0 = ((edge2 * (q0 >> kTrisoupFpBits)) * inva) >> precDivA;

  Vec3<int32_t>  ray1 = { 0,0,0 };
  ray1[i1] = kTrisoupFpOne;
  int32_t u1 = (h[i1] * inva) >> (precDivA - kTrisoupFpBits); //(ray1 * h) / a;
  Vec3<int32_t> q1 = crossProduct(ray1, edge1);
  int32_t v1 = (q1[direction] * inva) >> precDivA;
  int32_t t1 = ((edge2 * (q1 >> kTrisoupFpBits)) * inva) >> precDivA;

  Vec3<int32_t>  ray2 = { 0,0,0 };
  ray2[i2] = kTrisoupFpOne;
  int32_t u2 = (h[i2] * inva) >> (precDivA - kTrisoupFpBits); //(ray2 * h) / a;
  Vec3<int32_t> q2 = crossProduct(ray2, edge1);
  int32_t v2 = (q2[direction] * inva) >> precDivA;
  int32_t t2 = ((edge2 * (q2 >> kTrisoupFpBits)) * inva) >> precDivA;

  for (int32_t g1 = startposG1;
      g1 <= endposG1;
      g1++, u0 += u1, v0 += v1, t0 += t1, rayOrigin0[i1] += kTrisoupFpOne) {

    Vec3<int32_t> rayOrigin = rayOrigin0;
    int32_t u = u0;
    int32_t v = v0;
    int32_t t = t0;

    for (int32_t g2 = startposG2;
        g2 <= endposG2;
        g2++, u += u2, v += v2, t += t2, rayOrigin[i2] += kTrisoupFpOne) {

      int w = kTrisoupFpOne - u - v;
      if (u >= -haloTriangle && v >= -haloTriangle && w >= -haloTriangle) {
        Vec3<int32_t>  intersection = rayOrigin;
        intersection[direction] += t;
        Vec3<int32_t> foundvoxel =
          (intersection + truncateValue) >> kTrisoupFpBits;
        if (foundvoxel[direction]>=0 && foundvoxel[direction] < blockWidth) {
          refinedVerticesBlock.push_back(nodepos + foundvoxel);
        }

        intersection[direction] += thickness;
        Vec3<int32_t> foundvoxelUp =
          (intersection + truncateValue) >> kTrisoupFpBits;
        if (foundvoxelUp != foundvoxel
            && foundvoxelUp[direction] >= 0
            && foundvoxelUp[direction] < blockWidth) {
          refinedVerticesBlock.push_back(nodepos + foundvoxelUp);
        }

        intersection[direction] -= 2 * thickness;
        Vec3<int32_t> foundvoxelDown =
          (intersection + truncateValue) >> kTrisoupFpBits;
        if (foundvoxelDown != foundvoxel
            && foundvoxelDown[direction] >= 0
            && foundvoxelDown[direction] < blockWidth) {
          refinedVerticesBlock.push_back(nodepos + foundvoxelDown);
        }
      }
    }// loop g2
  }//loop g1

}

//============================================================================

}  // namespace pcc
