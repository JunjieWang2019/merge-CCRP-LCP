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

#include <cstdint>
#include <vector>
#include <cstring>

#include "PCCPointSet.h"
#include "geometry_octree.h"

namespace pcc {
//============================================================================
  struct codeVertexCtxInfo {
    int ctxE;
    int ctx0;
    int ctx1;
    int direction;
    int pattern = 0;
    int patternClose = 0;
    int patternClosest = 0;
    int nclosestPattern = 0;
    int missedCloseStart;
    int nclosestStart;
    int neighbEdge;
    int neighbEnd;
    int neighbStart;
    int orderedPclosePar;
    int orderedPcloseParPos;
  };

  static const int towardOrAway[18] = { // 0 = toward; 1 = away
   0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

  static const int mapping18to9[3][9] = {
    { 0, 1, 2, 3,  4, 15, 14, 5,  7},
    { 0, 1, 2, 3,  9, 15, 14, 7, 12},
    { 0, 1, 2, 9, 10, 15, 14, 7, 12}
  };

  void constructCtxInfo(
    codeVertexCtxInfo& ctxInfo,
    int neigh, std::array<int, 18>& patternIdx,
    std::vector<int8_t>& TriSoupVertices,
    int nbitsVertices,
    int max2bits,
    int mid2bits);

  void constructCtxPresence(
    int& ctxMap1,
    int& ctxMap2,
    int& ctxInter,
    codeVertexCtxInfo& ctxInfo,
    bool isInter,
    int8_t TriSoupVerticesPred);

  void constructCtxPos1(
    int& ctxMap1,
    int& ctxMap2,
    int& ctxInter,
    codeVertexCtxInfo& ctxInfo,
    bool isInter,
    int8_t TriSoupVerticesPred,
    int b);

  void constructCtxPos2(
    int& ctxMap1,
    int& ctxMap2,
    int& ctxInter,
    codeVertexCtxInfo& ctxInfo,
    bool isInter,
    int8_t TriSoupVerticesPred,
    int b,
    int v);

  //============================================================================
  // Representation for a vertex in preparation for sorting.
  struct Vertex {
    Vec3<int32_t> pos;  // position of vertex
    int32_t theta;      // angle of vertex when projected along dominant axis
    int32_t tiebreaker;  // coordinate of vertex along dominant axis
    bool operator()(Vertex v1, Vertex v2)
    {
      if (v1.theta > v2.theta)
        return true;  // sort in decreasing order of theta
      if (v1.theta == v2.theta && v1.tiebreaker < v2.tiebreaker)
        return true;
      return false;
    }

    bool operator==(Vertex v1)
    {
      return this->pos == v1.pos;
    }

    bool operator<(Vertex v1)
    {
      return this->pos < v1.pos;
    }

  } ;


 //============================================================================
  // index in edgesNeighNodes of neighboring nodes for each edge
  static const int edgesNeighNodesIdx[3][4] = {
    {6, 4, 0, 2}, // along z
    {6, 2, 5, 1}, // along y
    {6, 4, 5, 3}, // along x
  };

  // axis direction
  static const int axisdirection[3][3] = { {2,0,1}, {1, 0,2}, {0, 1, 2 } }; // z, y ,x

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


void determineTrisoupNeighbours(
  const std::vector<PCCOctree3Node>& leaves,
  const int defaultBlockWidth,
  PCCPointSet3& pointCloud,
  bool isEncoder,
  int bitDropped,
  int distanceSearchEncoder,
  bool isInter,
  const PCCPointSet3& refPointCloud,
  const PCCPointSet3& compensatedPointCloud,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  pcc::EntropyEncoder* arithmeticEncoder,
  pcc::EntropyDecoder& arithmeticDecoder,
  GeometryOctreeContexts& ctxtMemOctree,
  const bool isCentroidDriftActivated,
  bool haloFlag,
  bool adaptiveHaloFlag,
  int thickness,
  int &nSegments);


//============================================================================
void encodeCentroidResidual(
  int driftQ,
  pcc::EntropyEncoder* arithmeticEncoder,
  GeometryOctreeContexts& ctxtMemOctree,
  int driftQPred,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound);

int decodeCentroidResidual(
  pcc::EntropyDecoder* arithmeticDecoder,
  GeometryOctreeContexts& ctxtMemOctree,
  int driftQPred,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound);

int findDominantAxis(
  std::vector<Vertex>& leafVertices,
  Vec3<uint32_t> blockWidth,
  Vec3<int32_t> blockCentroid);


void rayTracingAlongdirection_samp1_optim(
  std::vector<Vec3<int32_t>>& refinedVerticesBlock,
  int direction,
  int blockWidth,
  Vec3<int32_t> posNnodeposode,
  int minRange[3],
  int maxRange[3],
  Vec3<int32_t> edge1,
  Vec3<int32_t> edge2,
  Vec3<int32_t> v0,
  int haloTriangle,
  int thickness);

//============================================================================

enum{
  POS_000 = 0,
  POS_W00 = 1,
  POS_0W0 = 2,
  POS_WW0 = 3,
  POS_00W = 4,
  POS_W0W = 5,
  POS_0WW = 6,
  POS_WWW = 7
};

void nonCubicNode
(
 const GeometryParameterSet& gps,
 const GeometryBrickHeader& gbh,
 const Vec3<int32_t>& leafpos,
 const int32_t blockWidth,
 const Box3<int32_t>& bbox,
 Vec3<int32_t>& newp,
 Vec3<int32_t>& neww,
 Vec3<int32_t>* corner
 );

//============================================================================
bool
boundaryinsidecheck(const Vec3<int32_t> a, const int bbsize);

template<typename T>
Vec3<T>
crossProduct(const Vec3<T> a, const Vec3<T> b);

void
determineCentroidAndDominantAxis(
  Vec3<int32_t>& blockCentroid,
  int& dominantAxis,
  std::vector<Vertex>& leafVertices,
  Vec3<int32_t> nodew);

Vec3<int32_t>
determineCentroidNormalAndBounds(
  int& lowBound,
  int& highBound,
  int& lowBoundSurface,
  int& highBoundSurface,
  int& ctxMinMax,
  int bitDropped,
  int bitDropped2,
  int triCount,
  Vec3<int32_t> blockCentroid,
  int dominantAxis,
  std::vector<Vertex>& leafVertices,
  int nodewDominant);

int
determineCentroidPredictor(
  int bitDropped2,
  Vec3<int32_t> normalV,
  Vec3<int32_t> blockCentroid,
  Vec3<int32_t> nodepos,
  const PCCPointSet3& compensatedPointCloud,
  int start,
  int end,
  int lowBound,
  int  highBound);

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
  int  highBound);

//============================================================================

}  // namespace pcc
