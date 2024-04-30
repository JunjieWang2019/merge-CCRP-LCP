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

#define PC_PREALLOCATION_SIZE 200000

static const int precDivA = 30;
static const int LUTdivtriCount[13] = { 32767, 1024, 512, 341, 256, 205, 171, 146, 128, 114, 102, 93, 85 }; // 10 bits precision

namespace pcc {

//============================================================================
// The number of fractional bits used in trisoup triangle voxelisation
  const int kTrisoupFpBits = 8;

  // The value 1 in fixed-point representation
  const int kTrisoupFpOne = 1 << (kTrisoupFpBits);
  const int kTrisoupFpHalf = 1 << (kTrisoupFpBits - 1);
  const int truncateValue = kTrisoupFpHalf;

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

    int nBadPredRef = 0;
    int nBadPredRef1 = 0;
    int nBadPredRef2 = 0;

    int nBadPredComp = 0;
    int nBadPredComp1 = 0;
    int nBadPredComp2 = 0;
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
    int mid2bits,
    std::vector<int8_t>& qualityRef,
    std::vector<int8_t>& qualityComp);

  void constructCtxPresence(
    int& ctxMap1,
    int& ctxMap2,
    int& ctxInter,
    codeVertexCtxInfo& ctxInfo,
    bool isInter,
    int8_t TriSoupVerticesPred,
    int8_t colocatedVertex);

  void constructCtxPos1(
    int& ctxMap1,
    int& ctxMap2,
    int& ctxInter,
    codeVertexCtxInfo& ctxInfo,
    bool isInter,
    int8_t TriSoupVerticesPred,
    int b,
    int8_t colocatedVertex);

  void constructCtxPos2(
    int& ctxMap1,
    int& ctxMap2,
    int& ctxInter,
    codeVertexCtxInfo& ctxInfo,
    bool isInter,
    int8_t TriSoupVerticesPred,
    int b,
    int v,
    int8_t colocatedVertex);

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

//============================================================================
struct CentroidInfo{
  int driftQPred = 0;
  bool possibleSKIP;
  int driftSKIP;
  int qualitySKIP;
};


//============================================================================
void encodeCentroidResidual(
  int driftQ,
  pcc::EntropyEncoder* arithmeticEncoder,
  GeometryOctreeContexts& ctxtMemOctree,
  CentroidInfo centroidInfo,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound);

int decodeCentroidResidual(
  pcc::EntropyDecoder* arithmeticDecoder,
  GeometryOctreeContexts& ctxtMemOctree,
  CentroidInfo centroidInfo,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound);

int findDominantAxis(
  std::vector<Vertex>& leafVertices,
  Vec3<uint32_t> blockWidth,
  Vec3<int32_t> blockCentroid);

  void rayTracingAlongdirection_samp1_optimX(
    std::vector<int64_t>& renderedBlock,
    int& nPointsInBlock,
    int blockWidth,
    Vec3<int32_t>& nodepos,
    int minRange[3],
    int maxRange[3],
    Vec3<int32_t>& edge1,
    Vec3<int32_t>& edge2,
    Vec3<int32_t>& s0,
    int64_t inva,
    int haloTriangle,
    int thickness);

  void rayTracingAlongdirection_samp1_optimY(
    std::vector<int64_t>& renderedBlock,
    int& nPointsInBlock,
    int blockWidth,
    Vec3<int32_t>& nodepos,
    int minRange[3],
    int maxRange[3],
    Vec3<int32_t>& edge1,
    Vec3<int32_t>& edge2,
    Vec3<int32_t>& s0,
    int64_t inva,
    int haloTriangle,
    int thickness);


  void rayTracingAlongdirection_samp1_optimZ(
    std::vector<int64_t>& renderedBlock,
    int& nPointsInBlock,
    int blockWidth,
    Vec3<int32_t>& nodepos,
    int minRange[3],
    int maxRange[3],
    Vec3<int32_t>& edge1,
    Vec3<int32_t>& edge2,
    Vec3<int32_t>& s0,
    int64_t inva,
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
crossProduct(const Vec3<T> a, const Vec3<T> b)
{
  Vec3<T> ret;
  ret[0] = a[1] * b[2] - a[2] * b[1];
  ret[1] = a[2] * b[0] - a[0] * b[2];
  ret[2] = a[0] * b[1] - a[1] * b[0];
  return ret;
}

void
determineCentroidAndDominantAxis(
  Vec3<int32_t>& blockCentroid,
  int& dominantAxis,
  std::vector<Vertex>& leafVertices,
  Vec3<int32_t> nodew,
  bool &flagCentroOK,
  int bitDropped,
  int &scaleQ);

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
  int nodewDominant,
  int blockWidth,
  int scaleQ);

void
determineCentroidPredictor(
  CentroidInfo& centroidInfo,
  int bitDropped2,
  Vec3<int32_t> normalV,
  Vec3<int32_t> blockCentroid,
  Vec3<int32_t> nodepos,
  const PCCPointSet3& compensatedPointCloud,
  int start,
  int end,
  int lowBound,
  int  highBound,
  int badQualityComp,
  int badQualityRef,
  int driftRef,
  bool possibleSKIPRef,
  int scaleQ);

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
  int  highBound,
  int scaleQ,
  int &drift);

bool nodeBoundaryInsideCheck(Vec3<int32_t> bw, Vec3<int32_t> pt);


struct CentroidDrift {
  int driftQ;
  int lowBound;
  int highBound;
  int ctxMinMax;
  int lowBoundSurface;
  int highBoundSurface;
};

struct TrisoupNodeEdgeVertex {
  int dominantAxis;
  // ( x - 0.5 ) << kTrisoupFpBits ( [-0.5, W-0.5] x256 )
  std::vector<Vertex> vertices;
};


struct TrisoupCentroidVertex {
  bool valid;  // this represents centroid existence
  Vec3<int32_t> pos;
  int32_t drift;
  bool boundaryInside;  // true if pos is inside of node boundary
};


struct TrisoupNodeFaceVertex {
  std::vector<Vertex> vertices;
  std::vector<int> formerEdgeVertexIdx;
};


struct TrisoupFace {
  bool connect;

  TrisoupFace(const bool cn)
  : connect(cn)
  {}

  TrisoupFace()
  { this->clear(); }

  void clear(void)
  { connect = false; }
};


//============================================================================
struct RasterScanTrisoupEdges {
  const int32_t blockWidth;
  // Eight corners of block.
  const Vec3<int32_t> pos000{ 0, 0, 0 };
  const Vec3<int32_t> posW00{ blockWidth, 0, 0 };
  const Vec3<int32_t> pos0W0{ 0, blockWidth, 0 };
  const Vec3<int32_t> posWW0{ blockWidth, blockWidth, 0 };
  const Vec3<int32_t> pos00W{ 0, 0, blockWidth };
  const Vec3<int32_t> posW0W{ blockWidth, 0, blockWidth };
  const Vec3<int32_t> pos0WW{ 0, blockWidth, blockWidth };
  const Vec3<int32_t> posWWW{ blockWidth, blockWidth, blockWidth };

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

  const int startCorner[12] = { POS_000, POS_000, POS_0W0, POS_W00, POS_000, POS_0W0, POS_WW0, POS_W00, POS_00W, POS_00W, POS_0WW, POS_W0W };
  const int endCorner[12] = { POS_W00, POS_0W0, POS_WW0, POS_WW0, POS_00W, POS_0WW, POS_WWW, POS_W0W, POS_W0W, POS_0WW, POS_WWW, POS_WWW };
  const int LUTsegmentDirection[12] = { 0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1 };

  std::array<int, 8> edgesNeighNodes; // neighboring nodes' index
  // The 7 firsts are used for unique segments generation/iteration
  // All the 8 are used for contextual information to be used by edge entropy coder
  Vec3<int32_t> currWedgePos = { INT32_MIN,INT32_MIN,INT32_MIN };
  int lastWedgex = 0;

  const std::vector<PCCOctree3Node>& leaves;
  PCCPointSet3& pointCloud;
  const bool isEncoder;
  const int bitDropped;
  const int distanceSearchEncoder;
  const bool isInter;
  const bool interSkipEnabled;
  const PCCPointSet3& compensatedPointCloud;

  // for coding
  const GeometryParameterSet& gps;
  const GeometryBrickHeader& gbh;
  pcc::EntropyEncoder* arithmeticEncoder;
  pcc::EntropyDecoder& arithmeticDecoder;
  GeometryOctreeContexts& ctxtMemOctree;

  AdaptiveBitModel ctxFaces;
  std::vector<TrisoupNodeEdgeVertex> eVerts;
  std::vector<TrisoupCentroidVertex> cVerts;
  std::vector<TrisoupNodeFaceVertex> fVerts;
  std::vector<Vec3<int32_t>> gravityCenter;
  std::vector<std::array<int, 3>> neiNodeIdxVec;

  // Box
  point_t BBorig;
  point_t keyshift;
  Box3<int32_t> sliceBB;

  // local variables for loop on wedges
  std::vector<int8_t> TriSoupVertices;
  std::vector<uint16_t> neighbNodes;
  std::queue<std::array<int, 18>> edgePattern;
  std::queue<int8_t> TriSoupVerticesPred;
  std::vector<int> segmentUniqueIndex;

  // for slice tracking
  int uniqueIndex = 0;
  int firstVertexToCode = 0;
  int nodeIdxC = 0;
  int firstNodeToRender = 0;

  // for edge coding
  std::queue<int> xForedgeOfVertex;
  int nbitsVertices ;
  int max2bits;
  int mid2bits;

  // for colocated edge tracking
  std::vector<int64_t> currentFrameEdgeKeys;
  int colocatedEdgeIdx = 0;
  std::vector<int8_t> qualityRef;
  std::vector<int8_t> qualityComp;

  // for colocated centroid tracking
  std::vector<int64_t> currentFrameNodeKeys;
  std::vector<int8_t> CentroValue;
  int colocatedNodeIdx = 0;

  // for rendering
  PCCPointSet3 recPointCloud;
  int nRecPoints = 0;
  int idxSegment = 0;
  std::vector<int64_t> renderedBlock;

  bool haloFlag;
  int thickness ;
  bool isCentroidDriftActivated;
  bool isFaceVertexActivated;

  int haloTriangle = 0;

  // flag to indicate last TriSoup pass
  bool isFinalPass = true;


  // constructor
  RasterScanTrisoupEdges(const std::vector<PCCOctree3Node>& leaves, int blockWidth, PCCPointSet3& pointCloud, bool isEncoder,
    int bitDropped, int distanceSearchEncoder, bool isInter, const PCCPointSet3& compensatedPointCloud,
    const GeometryParameterSet& gps, const GeometryBrickHeader& gbh, pcc::EntropyEncoder* arithmeticEncoder, pcc::EntropyDecoder& arithmeticDecoder, GeometryOctreeContexts& ctxtMemOctree)
    : leaves(leaves)
    , blockWidth(blockWidth)
    , pointCloud(pointCloud)
    , isEncoder(isEncoder)
    , bitDropped(bitDropped)
    , distanceSearchEncoder(distanceSearchEncoder)
    , isInter(isInter)
    , interSkipEnabled(isInter&& gps.trisoup_skip_mode_enabled_flag)
    , compensatedPointCloud(compensatedPointCloud)
    , gps(gps)
    , gbh(gbh)
    , arithmeticEncoder(arithmeticEncoder)
    , arithmeticDecoder(arithmeticDecoder)
    , ctxtMemOctree(ctxtMemOctree)
    , currWedgePos( Vec3<int32_t>{INT32_MIN, INT32_MIN, INT32_MIN})
    , edgesNeighNodes{ 0,0,0,0,0,0,0,0 }
  {
  }

  //---------------------------------------------------------------------------
  void  encodeOneTriSoupVertexRasterScan(
    int8_t vertex,
    pcc::EntropyEncoder* arithmeticEncoder,
    GeometryOctreeContexts& ctxtMemOctree,
    std::vector<int8_t>& TriSoupVertices,
    int neigh,
    std::array<int, 18>& patternIdx,
    int8_t interPredictor,
    int8_t colocatedVertex,
    std::vector<int8_t>& qualityRef,
    std::vector<int8_t>& qualityComp,
    int nbitsVertices,
    int max2bits,
    int mid2bits,
    int i) {

    codeVertexCtxInfo ctxInfo;
    constructCtxInfo(ctxInfo, neigh, patternIdx, TriSoupVertices, nbitsVertices, max2bits, mid2bits, qualityRef, qualityComp);

    // encode vertex presence
    int ctxMap1, ctxMap2, ctxInter;
    constructCtxPresence(ctxMap1, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, colocatedVertex);

    int ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[ctxInter][0].getEvolve(
      vertex >= 0, ctxMap2, ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
      ctxtMemOctree._BufferOBUFleavesTrisoup);

    arithmeticEncoder->encode(
      (int)(vertex >= 0), ctxTrisoup >> 3,
      ctxtMemOctree.ctxTriSoup[0][ctxInter][ctxTrisoup],
      ctxtMemOctree.ctxTriSoup[0][ctxInter].obufSingleBound);


    // quality ref edge
    if (interSkipEnabled && (vertex >= 0) == (colocatedVertex >= 0)) {
      qualityRef[i] = 1 + (vertex >= 0);
      if (qualityRef[i] >= 2) {
        qualityRef[i] += (vertex >> nbitsVertices - 1) == (colocatedVertex >> nbitsVertices - 1);
        qualityRef[i] += (vertex >> nbitsVertices - 2) == (colocatedVertex >> nbitsVertices - 2);
      }
    }
    // quality comp edge
    if (isInter && (vertex >= 0) == (interPredictor >= 0 ? 1 : 0)) {
      qualityComp[i] = 1 + (vertex >= 0);
      if (qualityComp[i] >= 2) {
        qualityComp[i] += (vertex >> nbitsVertices - 1) == (interPredictor >> nbitsVertices - 1);
        qualityComp[i] += (vertex >> nbitsVertices - 2) == (interPredictor >> nbitsVertices - 2);
      }
    }

    // encode  vertex position
    if (vertex >= 0) {
      int v = 0;
      int b = nbitsVertices - 1;

      // first position bit
      constructCtxPos1(ctxMap1, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, b, colocatedVertex);
      int bit = (vertex >> b--) & 1;

      ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[ctxInter][1].getEvolve(
        bit, ctxMap2, ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
        ctxtMemOctree._BufferOBUFleavesTrisoup);
      arithmeticEncoder->encode(
        bit, ctxTrisoup >> 3,
        ctxtMemOctree.ctxTriSoup[1][ctxInter][ctxTrisoup],
        ctxtMemOctree.ctxTriSoup[1][ctxInter].obufSingleBound);
      v = bit;

      // second position bit
      if (b >= 0) {
        constructCtxPos2(ctxMap1, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, b, v, colocatedVertex);

        bit = (vertex >> b--) & 1;
        ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[ctxInter][2].getEvolve(
          bit, ctxMap2, (ctxMap1 << 1) + v,
          &ctxtMemOctree._OBUFleafNumberTrisoup,
          ctxtMemOctree._BufferOBUFleavesTrisoup);
        arithmeticEncoder->encode(
          bit, ctxTrisoup >> 3,
          ctxtMemOctree.ctxTriSoup[2][ctxInter][ctxTrisoup],
          ctxtMemOctree.ctxTriSoup[2][ctxInter].obufSingleBound);
        v = (v << 1) | bit;
      }

      // third bit
      if (b >= 0) {
        int ctxFullNboundsReduced1 = (6 * (ctxInfo.ctx0 >> 1) + ctxInfo.missedCloseStart) * 2 + (ctxInfo.ctxE == 3);
        bit = (vertex >> b--) & 1;
        arithmeticEncoder->encode(
          bit, ctxtMemOctree.ctxTempV2[4 * ctxFullNboundsReduced1 + v]);
        v = (v << 1) | bit;
      }

      // remaining bits are bypassed
      for (; b >= 0; b--)
        arithmeticEncoder->encode((vertex >> b) & 1);
    }
  }

  //---------------------------------------------------------------------------
  void  decodeOneTriSoupVertexRasterScan(
    pcc::EntropyDecoder& arithmeticDecoder,
    GeometryOctreeContexts& ctxtMemOctree,
    std::vector<int8_t>& TriSoupVertices,
    int neigh,
    std::array<int, 18>& patternIdx,
    int8_t interPredictor,
    int8_t colocatedVertex,
    std::vector<int8_t>& qualityRef,
    std::vector<int8_t>& qualityComp,
    int nbitsVertices,
    int max2bits,
    int mid2bits,
    int i) {

    codeVertexCtxInfo ctxInfo;
    constructCtxInfo(ctxInfo, neigh, patternIdx, TriSoupVertices, nbitsVertices, max2bits, mid2bits, qualityRef, qualityComp);

    // decode vertex presence
    int ctxMap1, ctxMap2, ctxInter;
    constructCtxPresence(ctxMap1, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, colocatedVertex);

    bool c = ctxtMemOctree.MapOBUFTriSoup[ctxInter][0].decodeEvolve(
      &arithmeticDecoder, ctxtMemOctree.ctxTriSoup[0][ctxInter], ctxMap2,
      ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
      ctxtMemOctree._BufferOBUFleavesTrisoup);

    if (!c)
      TriSoupVertices.push_back(-1);

    // quality ref edge
    if (interSkipEnabled && c == (colocatedVertex >= 0 ? 1 : 0)) {
      qualityRef[i] = 1 + (c != 0);
    }
    // quality comp edge
    if (isInter && c == (interPredictor >= 0 ? 1 : 0)) {
      qualityComp[i] = 1 + (c != 0);
    }


    // decode vertex position
    if (c) {
      uint8_t v = 0;
      int b = nbitsVertices - 1;

      // first position bit
      constructCtxPos1(ctxMap1, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, b, colocatedVertex);
      int bit = ctxtMemOctree.MapOBUFTriSoup[ctxInter][1].decodeEvolve(
        &arithmeticDecoder, ctxtMemOctree.ctxTriSoup[1][ctxInter], ctxMap2,
        ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
        ctxtMemOctree._BufferOBUFleavesTrisoup);
      v = (v << 1) | bit;
      b--;

      // second position bit
      if (b >= 0) {
        constructCtxPos2(ctxMap1, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, b, v, colocatedVertex);
        bit = ctxtMemOctree.MapOBUFTriSoup[ctxInter][2].decodeEvolve(
          &arithmeticDecoder, ctxtMemOctree.ctxTriSoup[2][ctxInter], ctxMap2,
          (ctxMap1 << 1) + v, &ctxtMemOctree._OBUFleafNumberTrisoup,
          ctxtMemOctree._BufferOBUFleavesTrisoup);
        v = (v << 1) | bit;
        b--;
      }

      // third bit
      if (b >= 0) {
        int ctxFullNboundsReduced1 = (6 * (ctxInfo.ctx0 >> 1) + ctxInfo.missedCloseStart) * 2 + (ctxInfo.ctxE == 3);
        v = (v << 1) | arithmeticDecoder.decode(ctxtMemOctree.ctxTempV2[4 * ctxFullNboundsReduced1 + v]);
        b--;
      }

      // remaining bits are bypassed
      for (; b >= 0; b--)
        v = (v << 1) | arithmeticDecoder.decode();

      TriSoupVertices.push_back(v);

      // quality ref edge
      if (interSkipEnabled && qualityRef[i] >= 2) {
        qualityRef[i] += (v >> nbitsVertices - 1) == (colocatedVertex >> nbitsVertices - 1);
        qualityRef[i] += (v >> nbitsVertices - 2) == (colocatedVertex >> nbitsVertices - 2);
      }

      // quality comp edge
      if (isInter && qualityComp[i] >= 2) {
        qualityComp[i] += (v >> nbitsVertices - 1) == (interPredictor >> nbitsVertices - 1);
        qualityComp[i] += (v >> nbitsVertices - 2) == (interPredictor >> nbitsVertices - 2);
      }

    }
  }

  //---------------------------------------------------------------------------
  void flush2PointCloud(
    int& nRecPoints,
    std::vector<int64_t>::iterator itBegingBlock,
    int nPointsInBlock,
    PCCPointSet3& recPointCloud
  )
  {
    std::sort(itBegingBlock, itBegingBlock + nPointsInBlock);
    auto last = std::unique(itBegingBlock, itBegingBlock + nPointsInBlock);

    // Move list of points to pointCloud
    int nPointInCloud = recPointCloud.getPointCount();

    int nPointInNode = last - itBegingBlock;
    if (nPointInCloud < nRecPoints + nPointInNode)
      recPointCloud.resize(nRecPoints + nPointInNode + PC_PREALLOCATION_SIZE);

    for (auto it = itBegingBlock; it != last; it++)
      recPointCloud[nRecPoints++] = { int(*it >> 40), int(*it >> 20) & 0xFFFFF, int(*it) & 0xFFFFF };
  }

  void generateCentroidsInNodeRasterScan(
    const PCCOctree3Node& leaf,
    const std::vector<int8_t>& TriSoupVertices,
    int& idxSegment,
    Box3<int32_t>& sliceBB,
    const bool isCentroidDriftActivated,
    std::vector<int>& segmentUniqueIndex,
    std::vector<int8_t>& qualityRef,
    std::vector<int8_t>& qualityComp,
    bool nodeRefExist,
    int8_t colocatedCentroid,
    std::vector<int8_t>& CentroValue)
  {
    Vec3<int32_t> nodepos, nodew, corner[8];
    nonCubicNode(gps, gbh, leaf.pos, blockWidth, sliceBB, nodepos, nodew, corner);

    // Find up to 12 vertices for this leaf.
    std::vector<Vertex> leafVertices;

    TrisoupNodeEdgeVertex neVertex;

    // inter quality for centroid
    int badQualityRef = 0;
    int badQualityComp = 0;

    for (int j = 0; j < 12; j++) {
      int uniqueIndex = segmentUniqueIndex[idxSegment++];
      int vertex = TriSoupVertices[uniqueIndex];

      int qualityR = qualityRef[uniqueIndex];
      int qualityC = qualityComp[uniqueIndex];
      badQualityRef += qualityR == 0 || qualityR == 2 || qualityR == 3;
      badQualityComp += qualityC == 0 || qualityC == 2 || qualityC == 3;

      if (vertex < 0)
        continue;  // skip segments that do not intersect the surface

      // Get 3D position of point of intersection.
      Vec3<int32_t> point = corner[startCorner[j]] << kTrisoupFpBits;
      point -= kTrisoupFpHalf; // the volume is [-0.5; B-0.5]^3

      // points on edges are located at integer values
      int32_t dequantizedPosition = (vertex << (kTrisoupFpBits + bitDropped)) + (kTrisoupFpHalf << bitDropped);
      point[LUTsegmentDirection[j]] += dequantizedPosition;

      // Add vertex to list of vertices.
      leafVertices.push_back({ point, 0, 0 });

      neVertex.vertices.push_back({ point, 0, 0 });
    }

    int vtxCount = (int)neVertex.vertices.size();
    Vec3<int32_t> gCenter = 0;
    for (int j = 0; j < vtxCount; j++) {
      gCenter += neVertex.vertices[j].pos;
    }
    if (vtxCount) {
      gCenter /= vtxCount;
    }
    neVertex.dominantAxis =
      findDominantAxis(neVertex.vertices, nodew, gCenter);
    eVerts.push_back(neVertex);

    // Skip leaves that have fewer than 3 vertices.
    if (vtxCount < 3) {
      cVerts.push_back({ false, { 0, 0, 0}, 0, true });
      gravityCenter.push_back({ 0, 0, 0 });
      return;
    }

    // compute centroid
    Vec3<int32_t> blockCentroid;
    int dominantAxis;
    int driftDQ = 0;
    bool flagCentroOK = isCentroidDriftActivated;
    int scaleQ = 256; // scaled on 8 bits; 256 is one
    determineCentroidAndDominantAxis(
      blockCentroid, dominantAxis, leafVertices, nodew, flagCentroOK,
      bitDropped, scaleQ);

    gCenter = blockCentroid;

    if (!isCentroidDriftActivated) {
      cVerts.push_back({ false, gCenter, 0, true });
      gravityCenter.push_back(gCenter);
      return;
    }

    // Refinement of the centroid along the domiannt axis
    if (vtxCount >= 3 && isCentroidDriftActivated) {
      int bitDropped2 = bitDropped;

      // colocated centroid
      bool possibleSKIPRef = nodeRefExist && badQualityRef <= 1;
      int driftRef = possibleSKIPRef ? colocatedCentroid : 0;

      int driftQ = 0;

      int half = 1 << 5 + bitDropped2;
      int DZ = 682 * half >> 10; // 2 * half / 3;

      int lowBound, highBound, lowBoundSurface, highBoundSurface, ctxMinMax;
      Vec3<int32_t> normalV =
        determineCentroidNormalAndBounds(
          lowBound, highBound, lowBoundSurface, highBoundSurface, ctxMinMax,
          bitDropped, bitDropped2, vtxCount, blockCentroid, dominantAxis,
          leafVertices, nodew[dominantAxis], blockWidth, scaleQ);

      flagCentroOK = lowBound != 0 || highBound != 0;
      if (flagCentroOK || possibleSKIPRef) {
        CentroidInfo centroidInfo;
        determineCentroidPredictor(
          centroidInfo, bitDropped2, normalV, blockCentroid, nodepos,
          compensatedPointCloud, leaf.predStart, leaf.predEnd, lowBound,
          highBound, badQualityComp, badQualityRef, driftRef, possibleSKIPRef, scaleQ);

        if (isEncoder) { // encode centroid residual

          int drift = 0;
          driftQ = determineCentroidResidual(
            bitDropped2, normalV, blockCentroid, nodepos, pointCloud,
            leaf.start, leaf.end, lowBound, highBound, scaleQ, drift);

          // naive RDO
          if (centroidInfo.possibleSKIP) {
            int driftDQSKip = 0;
            if (centroidInfo.driftSKIP) {
              driftDQSKip =
                std::abs(centroidInfo.driftSKIP) << bitDropped2 + 6;
              driftDQSKip += DZ - half;
              if (centroidInfo.driftSKIP < 0)
                driftDQSKip = -driftDQSKip;
            }

            if (std::abs(drift - driftDQSKip) <= DZ + half)
              driftQ = centroidInfo.driftSKIP;
          }

          encodeCentroidResidual(
            driftQ, arithmeticEncoder, ctxtMemOctree, centroidInfo, ctxMinMax,
            lowBoundSurface, highBoundSurface, lowBound, highBound);
        }
        else { // decode centroid residual
          driftQ =
            decodeCentroidResidual(
              &arithmeticDecoder, ctxtMemOctree, centroidInfo, ctxMinMax,
              lowBoundSurface, highBoundSurface, lowBound, highBound);
        }
      }

      // store centroid residual for next frame
      CentroValue.back() = driftQ;

      // dequantize and apply drift
      if (driftQ) {
        driftDQ = std::abs(driftQ) << bitDropped2 + 6;
        driftDQ += DZ - half;
        if (driftQ < 0)
          driftDQ = -driftDQ;
      }

      driftDQ =
        driftDQ > 0 ? (driftDQ * scaleQ >> 8) : -((-driftDQ) * scaleQ >> 8);

      blockCentroid += (driftDQ * normalV) >> 6;
      blockCentroid[0] =
        std::min(
          ((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1,
          std::max(-kTrisoupFpHalf, blockCentroid[0]));
      blockCentroid[1] =
        std::min(
          ((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1,
          std::max(-kTrisoupFpHalf, blockCentroid[1]));
      blockCentroid[2] =
        std::min(
          ((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1,
          std::max(-kTrisoupFpHalf, blockCentroid[2]));

      bool boundaryInside = true;
      if (!nodeBoundaryInsideCheck(nodew << kTrisoupFpBits, blockCentroid)) {
        boundaryInside = false;
      }

      cVerts.push_back({ true, blockCentroid, driftDQ, boundaryInside });
      gravityCenter.push_back(gCenter);
    } // end refinement of the centroid

    return;
  }



  // find face vertex position from connection of centroid vertices of
  // current and neighbour nodes one by one.
  // (totally this function is called up to three times per node.)
  void findTrisoupFaceVertex(
    const TrisoupCentroidVertex& cVerts0,
    const TrisoupCentroidVertex& cVerts1,
    const int axis, // order : -x,-y,-z
    const Vec3<int32_t>& neiNodeW,
    Vertex* fVert)
  {
    int32_t c0face = neiNodeW[axis] - kTrisoupFpHalf;
    Vec3<int32_t> c0 = cVerts0.pos;
    Vec3<int32_t> c1 = cVerts1.pos;
    c0[axis] += neiNodeW[axis];
    int32_t denom = c0[axis] - c1[axis];
    // denom=0 means that c0 and c1 are on the same grid(t=0).
    int32_t t = denom ? ((c0face - c1[axis] << kTrisoupFpBits) / denom) : 0;
    Vertex faceVertex(
      { c1 + (t * (c0 - c1) + kTrisoupFpHalf >> kTrisoupFpBits), 0, 0 });
    fVert[0] = faceVertex;
    fVert[0].pos[axis] = -kTrisoupFpHalf;
    fVert[1] = faceVertex;
    fVert[1].pos[axis] = neiNodeW[axis] - kTrisoupFpHalf;
  }


  int countTrisoupEdgeVerticesOnFace(
    const TrisoupNodeEdgeVertex& eVerts, Vec3<int32_t>& nodeW, int axis)
  {
    int neVtxBoundaryFace = 0;
    // count edge vertices included in the face
    for (int k = 0; k < eVerts.vertices.size(); k++) {
      Vec3<int32_t> vtxC = eVerts.vertices[k].pos + kTrisoupFpHalf;
      if (nodeW[axis] == vtxC[axis]) {
        neVtxBoundaryFace++;
      }
    }
    return neVtxBoundaryFace;
  }



  void determineTrisoupEdgeBoundaryLine(
    int i,
    const TrisoupNodeEdgeVertex& eVerts,
    const TrisoupCentroidVertex& cVerts,
    const Vec3<int32_t>& gravityCenter,
    Vec3<int32_t>& nodeW, int axis, Vertex& fvert, int* eIdx)
  {
    // if there were two or three edge vertices on the face,
    // then to select the nearest bridge segment
    // between two edge vertices from tentative face vertex of current node,
    // edge vertices within current node are already sorted,
    // and two vertices are selected which make the nearest segment from temtative face vertex.
    // if the surface has a hole within the current node,
    // the sequential number of edge vertex couldn't be found,
    // false is returned without creating a face vertex.
    int evCnt = eVerts.vertices.size();
    int dist = 0, distMin = 0x7fffffff; // initial value must be larger than any case.
    int evIdxMin[2] = { -1, -1 };
    nodeW -= kTrisoupFpHalf;
    for (int evIdx = 0; evIdx < (evCnt == 3 ? 1 : evCnt); evIdx++) {
      int ev0 = evIdx;
      int ev1 = evIdx + 1;
      if (ev1 >= evCnt) { ev1 -= evCnt; }

      Vec3<int32_t> evCoord0 = eVerts.vertices[ev0].pos;
      Vec3<int32_t> evCoord1 = eVerts.vertices[ev1].pos;
      if (!(nodeW[axis] == evCoord0[axis] && nodeW[axis] == evCoord1[axis])) {
        continue;
      }
      Vec3<int32_t> middlePoint = (evCoord0 + evCoord1) / 2;
      Vec3<int32_t> distVec = (middlePoint - fvert.pos) >> kTrisoupFpBits;
      dist =
        distVec[0] * distVec[0]
        + distVec[1] * distVec[1]
        + distVec[2] * distVec[2];
      if (distMin > dist) {
        evIdxMin[0] = ev0;
        evIdxMin[1] = ev1;
        distMin = dist;
      }
    }
    eIdx[0] = evIdxMin[0];
    eIdx[1] = evIdxMin[1];
  }



  // finding vector:
  //   1. gravityCenter to Centroid of current node
  //   2. gravityCenter to Centroid of current neighbour node
  //   3. face vertex vector on face of current node
  // and confirm directions of these three vectors are
  // not invert with inner product.
  bool determineTrisoupDirectionOfCentroidsAndFvert(
    const TrisoupNodeEdgeVertex& eVerts,
    const TrisoupCentroidVertex& cVerts0,
    const TrisoupCentroidVertex& cVerts1,
    const Vec3<int32_t>& gravityCenter0,
    const Vec3<int32_t>& gravityCenter1,
    int w, int e0, int e1, Vertex* fVert)
  {
    // unit vector between two edge vertices on boundary face
    Vec3<int64_t> euv = eVerts.vertices[e1].pos - eVerts.vertices[e0].pos;
    int64_t euvNorm =
      isqrt(euv[0] * euv[0] + euv[1] * euv[1] + euv[2] * euv[2]);
    euv = euvNorm ? ((euv << kTrisoupFpBits) / euvNorm) : 0;
    Vec3<int32_t> c0 = cVerts0.pos;
    Vec3<int32_t> c1 = cVerts1.pos;
    Vec3<int32_t> g0 = gravityCenter0;
    Vec3<int32_t> g1 = gravityCenter1;
    Vec3<int32_t> ef = fVert[0].pos - eVerts.vertices[e0].pos;
    int64_t en = ef * euv >> kTrisoupFpBits;
    int32_t dp0 = (c0 - g0) * (ef - (en * euv >> kTrisoupFpBits));
    int32_t dp1 = (c1 - g1) * (ef - (en * euv >> kTrisoupFpBits));
    bool judge = dp0 > 0 && dp1 > 0;
    return judge;
  }




  void generateFaceVerticesInNodeRasterScan(
    const std::vector<PCCOctree3Node>& leaves, const int nodeIdx, bool isEncoder)
  {
    const int32_t tmin1 = 2 * 4;
    const int32_t tmin2 = distanceSearchEncoder * 4;
    Box3<int32_t> sliceBB;
    sliceBB.min = gbh.slice_bb_pos << gbh.slice_bb_pos_log2_scale;
    sliceBB.max = sliceBB.min + (gbh.slice_bb_width << gbh.slice_bb_width_log2_scale);
    int32_t w = blockWidth;
    int i = nodeIdx;
    // For 6-neighbour nodes of three which have smaller coordinates than current node,
    // if current node and its neighbour node have refined centroid vertex each other,
    // and when the centroids are connected, if there is no bite on the current surface,
    // then the intersection of centroid connection segment and
    // node boundary face is defined as the temporary face vertex.
    // And if original points are distributed around the temporary face vertex,
    // it is defined as a true face vertex and determine to connect these centroids,
    // and then face-flag becomes true.
    // to generate 3 faces per node, neighbour direction loop must be placed at the most outer loop.
    // x,y,z-axis order 0,1,2
    for (int axis = 0; axis < 3; axis++) {
      if (cVerts[i].valid && cVerts[i].boundaryInside) {
        int ii = neiNodeIdxVec[i][axis];
        // neighbour-node exists on this direction
        if (-1 != ii) {
          Vec3<int32_t> nodepos, neiNodew, corner[8];
          nonCubicNode(
            gps, gbh, leaves[ii].pos, blockWidth, sliceBB, nodepos, neiNodew,
            corner);
          // centroid of the neighbour-node exists and inside of the boundary
          if (cVerts[ii].valid && cVerts[ii].boundaryInside) {
            int eIdx[2][2] = { -1 };
            Vec3<int32_t> neiNodeW = neiNodew << kTrisoupFpBits;
            Vec3<int32_t> zeroW = 0 << kTrisoupFpBits;
            int neVtxBoundaryFace =
              countTrisoupEdgeVerticesOnFace(eVerts[i], zeroW, axis);
            if (2 == neVtxBoundaryFace || 3 == neVtxBoundaryFace) {
              Vertex fVert[2];
              findTrisoupFaceVertex(
                cVerts[i], cVerts[ii], axis, neiNodeW, fVert);
              determineTrisoupEdgeBoundaryLine(
                i, eVerts[i], cVerts[i], gravityCenter[i], zeroW, axis,
                fVert[0], eIdx[0]);
              determineTrisoupEdgeBoundaryLine(
                ii, eVerts[ii], cVerts[ii], gravityCenter[ii], neiNodeW,
                axis, fVert[1], eIdx[1]);
              if (-1 != eIdx[0][0] && -1 != eIdx[0][1]) {
                bool judge =
                  determineTrisoupDirectionOfCentroidsAndFvert(
                    eVerts[i], cVerts[i], cVerts[ii], gravityCenter[i],
                    gravityCenter[ii], w, eIdx[0][0], eIdx[0][1], fVert);
                if (judge) {
                  TrisoupFace face(false);
                  if (isEncoder) {
                    // c0, c1, and face vertex is on the same side of the surface
                    int32_t weight1 = 0, weight2 = 0;
                    uint32_t st[2] = { leaves[i].start, leaves[ii].start };
                    uint32_t ed[2] = { leaves[i].end, leaves[ii].end };
                    // 0:current-node  1:nei-node
                    for (int n = 0; n < 2; n++) {
                      int _nodeidx = !n ? i : ii;
                      for (int k = st[n]; k < ed[n]; k++) {
                        Vec3<int32_t> dist =
                          fVert[n].pos - (
                            pointCloud[k] - leaves[_nodeidx].pos
                            << kTrisoupFpBits);
                        int32_t d =
                          dist.abs().max() + kTrisoupFpHalf >> kTrisoupFpBits;
                        if (d < tmin1) { weight1++; }
                        if (d < tmin2) { weight2++; }
                      }
                    }
                    if (weight1 > 0 || weight2 > 1) {
                      face.connect = true;
                    }
                    arithmeticEncoder->encode((int)face.connect, ctxFaces);
                  }
                  else {
                    face.connect = !!(arithmeticDecoder.decode(ctxFaces));
                  }
                  if (face.connect) {
                    fVerts[i].formerEdgeVertexIdx.push_back(eIdx[0][0]);
                    fVerts[i].vertices.push_back(fVert[0]);
                    fVerts[ii].formerEdgeVertexIdx.push_back(eIdx[1][0]);
                    fVerts[ii].vertices.push_back(fVert[1]);
                  }
                }
              }
              // if (2 or 3 == neVtxBoundaryFace)
            }
          }
        }
      }
    }
    return;
  }



  //---------------------------------------------------------------------------
  int generateTrianglesInNodeRasterScan(
    const PCCOctree3Node& leaf, int i,
    std::vector<int64_t>& renderedBlock,
    Box3<int32_t>& sliceBB,
    int haloTriangle,
    int thickness,
    bool isFaceVertexActivated)
  {
    Vec3<int32_t> nodepos, nodew, corner[8];
    nonCubicNode(
      gps, gbh, leaf.pos, blockWidth, sliceBB, nodepos, nodew, corner);

    int nPointsInBlock = 0;

    for (int j = 0; j < eVerts[i].vertices.size(); j++) {
      Vec3<int32_t> point =
        eVerts[i].vertices[j].pos + kTrisoupFpHalf >> kTrisoupFpBits;
      // vertex to list of points
      if (bitDropped) {
        if (boundaryinsidecheck(point, blockWidth - 1)) {
          Vec3<int64_t> renderedPoint = nodepos + point;
          renderedBlock[nPointsInBlock++] =
            (renderedPoint[0] << 40)
            + (renderedPoint[1] << 20)
            + renderedPoint[2];
        }
      }
    }
    // Skip leaves that have fewer than 3 vertices.
    if (eVerts[i].vertices.size() < 3) {
      return nPointsInBlock;
    }

    if (eVerts[i].vertices.size() >= 3) {
      Vec3<int32_t> foundvoxel =
        cVerts[i].pos + truncateValue >> kTrisoupFpBits;
      if (boundaryinsidecheck(foundvoxel, blockWidth - 1)) {
        Vec3<int64_t> renderedPoint = nodepos + foundvoxel;
        renderedBlock[nPointsInBlock++] =
          (renderedPoint[0] << 40)
          + (renderedPoint[1] << 20)
          + renderedPoint[2];
      }
    }

    std::vector<Vertex> nodeVertices;
    for (int j = 0; j < eVerts[i].vertices.size(); j++) {
      nodeVertices.push_back(eVerts[i].vertices[j]);
      if (isFaceVertexActivated) {
        for (int k = 0; k < fVerts[i].vertices.size(); k++) {
          if (j == fVerts[i].formerEdgeVertexIdx[k]) {
            nodeVertices.push_back(fVerts[i].vertices[k]);
          }
        }
      }
    }

    // Divide vertices into triangles around centroid
    // and upsample each triangle by an upsamplingFactor.
    int vtxCount = nodeVertices.size();
    Vec3<int32_t> blockCentroid = cVerts[i].pos;
    Vec3<int32_t> v2 = blockCentroid;
    Vec3<int32_t> v1 = nodeVertices[0].pos;
    Vec3<int32_t> posNode = nodepos << kTrisoupFpBits;

    for (int vtxIndex = 0; vtxIndex < vtxCount; vtxIndex++) {
      int j1 = vtxIndex + 1;
      if (j1 >= vtxCount)
        j1 -= vtxCount;

      Vec3<int32_t> v0 = v1;
      v1 = nodeVertices[j1].pos;


      // choose ray direction
      Vec3<int32_t> edge1 = v1 - v0;
      Vec3<int32_t> edge2 = v2 - v0;
      Vec3<int32_t> a = crossProduct(edge2, edge1) >> kTrisoupFpBits;
      Vec3<int32_t> h = a.abs();
      int directionOk = (h[0] > h[1] && h[0] > h[2]) ? 0 : h[1] > h[2] ? 1 : 2;

      // check if ray tracing is valid; if not skip triangle which is too small
      if (h[directionOk] <= kTrisoupFpOne) // < 2*kTrisoupFpOne should be ok
        continue;

      int64_t inva =
        divApprox(int64_t(1) << precDivA, std::abs(a[directionOk]), 0);
      inva = a[directionOk] > 0 ? inva : -inva;

      // range
      int minRange[3];
      int maxRange[3];
      for (int k = 0; k < 3; k++) {
        minRange[k] =
          std::max(
            0,
            std::min(std::min(v0[k], v1[k]), v2[k])
            + truncateValue >> kTrisoupFpBits);
        maxRange[k] =
          std::min(
            blockWidth - 1,
            std::max(std::max(v0[k], v1[k]), v2[k])
            + truncateValue >> kTrisoupFpBits);
      }
      Vec3<int32_t> s0 = {
        (minRange[0] << kTrisoupFpBits) - v0[0],
        (minRange[1] << kTrisoupFpBits) - v0[1],
        (minRange[2] << kTrisoupFpBits) - v0[2]
      };

      // ensure there is enough space in the block buffer
      if (renderedBlock.size() <= nPointsInBlock + blockWidth * blockWidth)
        renderedBlock.resize(renderedBlock.size() + blockWidth * blockWidth);

      // applying ray tracing along direction
      if (directionOk == 0)
        rayTracingAlongdirection_samp1_optimX(
          renderedBlock, nPointsInBlock, blockWidth, nodepos, minRange,
          maxRange, edge1, edge2, s0, inva, haloTriangle, thickness);

      if (directionOk == 1)
        rayTracingAlongdirection_samp1_optimY(
          renderedBlock, nPointsInBlock, blockWidth, nodepos, minRange,
          maxRange, edge1, edge2, s0, inva, haloTriangle, thickness);

      if (directionOk == 2)
        rayTracingAlongdirection_samp1_optimZ(
          renderedBlock, nPointsInBlock, blockWidth, nodepos, minRange,
          maxRange, edge1, edge2, s0, inva, haloTriangle, thickness);

    }  // end loop on triangles

    return nPointsInBlock;
  }




  void clearTrisoupElements(void)
  {
    eVerts.clear();
    cVerts.clear();
    fVerts.clear();
    gravityCenter.clear();
    neiNodeIdxVec.clear();
    return;
  }


  //---------------------------------------------------------------------------
  void callTriSoupSlice(bool flagFinal)
  {
    isFinalPass = flagFinal;
    buildSegments();
  }


  //---------------------------------------------------------------------------
  void init() {
    BBorig = gbh.geomBoxOrigin;
    keyshift = BBorig + (1 << 18);

    lastWedgex = currWedgePos[0];

    nbitsVertices = gbh.trisoupNodeSizeLog2(gps) - bitDropped;
    max2bits = nbitsVertices > 1 ? 3 : 1;
    mid2bits = nbitsVertices > 1 ? 2 : 1;

    if (isEncoder)
      recPointCloud.resize(pointCloud.getPointCount());

    sliceBB.min = gbh.slice_bb_pos << gbh.slice_bb_pos_log2_scale;
    sliceBB.max = sliceBB.min + (gbh.slice_bb_width << gbh.slice_bb_width_log2_scale);

    renderedBlock.resize(blockWidth * blockWidth * 16, 0);

    haloFlag = gbh.trisoup_halo_flag;
    thickness = gbh.trisoup_thickness;
    isCentroidDriftActivated = gbh.trisoup_centroid_vertex_residual_flag;
    isFaceVertexActivated = gbh.trisoup_face_vertex_flag;

    if (haloFlag) {
      haloTriangle = (((1 << bitDropped) - 1) << kTrisoupFpBits) / blockWidth; // this division is a shift if width is a power of 2
      haloTriangle = (haloTriangle * 28) >> 5; // / 32;
      haloTriangle = haloTriangle > 36 ? 36 : haloTriangle;
    }
  }


  //---------------------------------------------------------------------------
  void buildSegments()
  {
    neighbNodes.reserve(leaves.size() * 12); // at most 12 edges per node (to avoid reallocations)
    segmentUniqueIndex.resize(12 * leaves.size(), -1); // temporarily set to -1 to check everything is working
    // TODO: set to -1 could be removed when everything will work properly

    if (isFaceVertexActivated) {
      fVerts.resize(leaves.size());
    }

    if (currWedgePos == Vec3<int32_t>{INT32_MIN, INT32_MIN, INT32_MIN}) ;
      currWedgePos = leaves[0].pos;

    int lastAcceptable = INT_MAX;
    if (!isFinalPass)
      lastAcceptable = leaves.back().pos[0];

    while (nextIsAvailable() && currWedgePos[0] < lastAcceptable) { // this a loop on start position of edges; 3 edges along x,y and z per start position
      // process current wedge position

      std::array<bool, 8> isNeigbourSane;
      for (int i = 0; i < 8; ++i) // sanity of neighbouring nodes
        isNeigbourSane[i] = edgesNeighNodes[i] < leaves.size() && currWedgePos + offsets[i] == leaves[edgesNeighNodes[i]].pos;


      if (isFaceVertexActivated) {
        std::array<int, 3> nei3idx = { -1, -1, -1 };
        if (isNeigbourSane[6]) {
          nei3idx[2] = isNeigbourSane[5] ? edgesNeighNodes[5] : -1; // -z
          nei3idx[1] = isNeigbourSane[4] ? edgesNeighNodes[4] : -1; // -y
          nei3idx[0] = isNeigbourSane[2] ? edgesNeighNodes[2] : -1; // -x
          neiNodeIdxVec.push_back(nei3idx);
        }
      }

      for (int dir = 0; dir < 3; ++dir) { // this the loop on the 3 edges along z, then y, then x
        bool processedEdge = false;
        uint16_t neighboursMask = 0;
        std::array<int, 18> pattern{ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

        // for TriSoup Vertex by the encoder
        int countNearPoints = 0;
        int distanceSum = 0;
        int countNearPoints2 = 0;
        int distanceSum2 = 0;
        int dir0 = axisdirection[dir][0];
        int dir1 = axisdirection[dir][1];
        int dir2 = axisdirection[dir][2];

        // for TriSoup Vertex inter prediction
        int countNearPointsPred = 0;
        int distanceSumPred = 0;

        for (int neighIdx = 0; neighIdx < 4; ++neighIdx) { // this the loop on the 4 nodes to interesect the edge
          int edgeNeighNodeIdx = edgesNeighNodesIdx[dir][neighIdx]; // gives the neighbour inde in the list of 8 neighbours

          if (isNeigbourSane[edgeNeighNodeIdx]) { // test for sanity of the neighbour, process only if sane
            processedEdge = true; // at least one neighbour is sane, so the edge become a TriSoup edge and will be added to the list
            int neighbNodeIndex = edgesNeighNodes[edgeNeighNodeIdx];

            int idx = neighbNodeIndex * 12 + edgeIdx[dir][neighIdx];
            segmentUniqueIndex[idx] = uniqueIndex;

            // update mask from nodes touching edge
            neighboursMask |= neighMask[dir][neighIdx];

            int indexLow = edgeIdx[dir][neighIdx];
            for (int v = 0; v < 11; v++) {
              if (localEdgeindex[indexLow][v] == -1)
                break;

              int indexV = neighbNodeIndex * 12 + localEdgeindex[indexLow][v]; // index of segment
              int Vidx = segmentUniqueIndex[indexV];
              //assert(Vidx != -1); // check if already coded
              pattern[patternIndex[indexLow][v]] = Vidx;
            }

            // determine TriSoup Vertex by the encoder
            const int offset1 = leaves[neighbNodeIndex].pos[dir1] < currWedgePos[dir1];
            const int offset2 = leaves[neighbNodeIndex].pos[dir2] < currWedgePos[dir2];
            const int pos1 = currWedgePos[dir1] - offset1;
            const int pos2 = currWedgePos[dir2] - offset2;
            if (isEncoder) {
              for (int j = leaves[neighbNodeIndex].start; j < leaves[neighbNodeIndex].end; j++) {
                Vec3<int> voxel = pointCloud[j];
                if (voxel[dir1] == pos1 && voxel[dir2] == pos2) {
                  countNearPoints++;
                  distanceSum += voxel[dir0] - currWedgePos[dir0];
                }
                if (distanceSearchEncoder > 1 && std::abs(voxel[dir1] - pos1) < distanceSearchEncoder && std::abs(voxel[dir2] - pos2) < distanceSearchEncoder) {
                  countNearPoints2++;
                  distanceSum2 += voxel[dir0] - currWedgePos[dir0];
                }
              }
            }

            // determine TriSoup Vertex inter prediction
            if (isInter) {
              const PCCPointSet3& PC = compensatedPointCloud;
              for (int j = leaves[neighbNodeIndex].predStart; j < leaves[neighbNodeIndex].predEnd; j++) {
                Vec3<int> voxel = PC[j];
                if (voxel[dir1] == pos1 && voxel[dir2] == pos2) {
                  countNearPointsPred++;
                  distanceSumPred += voxel[dir0] - currWedgePos[dir0];
                }
              }
            }
          }
        } // end loop on 4 neighbouring nodes


        if (processedEdge) { // the TriSoup edge is added to the list
          int segmentUniqueIdxPrevEdge = -1;
          for (int prevNeighIdx = 0; prevNeighIdx < 4; ++prevNeighIdx) {
            int wedgeNeighNodeIdx = wedgeNeighNodesIdx[dir][prevNeighIdx];
            if (isNeigbourSane[wedgeNeighNodeIdx]) {
              // update current mask from nodes touching wedge
              neighboursMask |= wedgeNeighMask[dir][prevNeighIdx];
              if (segmentUniqueIdxPrevEdge == -1) {
                int idx = edgesNeighNodes[wedgeNeighNodeIdx] * 12 + wedgeNeighNodesEdgeIdx[dir][prevNeighIdx];
                segmentUniqueIdxPrevEdge = segmentUniqueIndex[idx];
                pattern[0] = segmentUniqueIdxPrevEdge;
                //assert(segmentUniqueIdxPrevEdge != -1);
                for (int neighIdx = 0; neighIdx < 4; ++neighIdx) {
                  int edgeNeighNodeIdx = edgesNeighNodesIdx[dir][neighIdx];
                  if (isNeigbourSane[edgeNeighNodeIdx]) {
                    neighbNodes[segmentUniqueIdxPrevEdge] |= toPrevEdgeNeighMask[dir][neighIdx];
                  }
                }
              }
            }
          }

          ++uniqueIndex;
          neighbNodes.push_back(neighboursMask);
          edgePattern.push(pattern);
          xForedgeOfVertex.push(currWedgePos[0]);

          // for colocated edges
          int64_t key = (int64_t(currWedgePos[0] + keyshift[0]) << 42) + (int64_t(currWedgePos[1] + keyshift[1]) << 22) + (int64_t(currWedgePos[2] + keyshift[2]) << 2) + dir;
          currentFrameEdgeKeys.push_back(key);
          qualityRef.push_back(0);
          qualityComp.push_back(0);

          // determine TriSoup Vertex by the encoder
          if (isEncoder)
          {
            int8_t vertexPos = -1;
            if (countNearPoints > 0 || countNearPoints2 > 1) {
              int temp = ((2 * distanceSum + distanceSum2) << (10 - bitDropped)) / (2 * countNearPoints + countNearPoints2); // division on the encoder side
              vertexPos = (temp + (1 << 9 - bitDropped)) >> 10;
            }
            TriSoupVertices.push_back(vertexPos);
          }

          // determine TriSoup Vertex inter prediction
          if (isInter) {
            int8_t vertexPos = -1;
            if (countNearPointsPred > 0) {
              int temp = divApprox(distanceSumPred << (10 - bitDropped), countNearPointsPred, 0);
              vertexPos = (temp + (1 << 9 - bitDropped)) >> 10;
            }
            TriSoupVerticesPred.push(vertexPos);
          }
        } // end if on is a TriSoup wedge
      } // end loop on three deirections

      // move to next wedge
      goNextWedge(isNeigbourSane);


      // code vertices and rendering of preceding slices in case the loop has moved up one slice or if finished
      if (changeSlice()) {
        // coding vertices
        int upperxForCoding = !nextIsAvailable() ? INT32_MAX : currWedgePos[0] - blockWidth;
        while (!xForedgeOfVertex.empty() && xForedgeOfVertex.front() < upperxForCoding) {
          // spatial neighbour and inter comp predictors
          int8_t  interPredictor = isInter ? TriSoupVerticesPred.front() : 0;
          auto pattern = edgePattern.front();

          // colocated edge predictor
          int8_t colocatedVertex = -1;
          if (interSkipEnabled) {
            auto keyCurrent = currentFrameEdgeKeys[firstVertexToCode];
            while (colocatedEdgeIdx < ctxtMemOctree.refFrameEdgeKeys.size() - 1 && ctxtMemOctree.refFrameEdgeKeys[colocatedEdgeIdx] < keyCurrent)
              colocatedEdgeIdx++;

            if (ctxtMemOctree.refFrameEdgeKeys[colocatedEdgeIdx] == keyCurrent)
              colocatedVertex = ctxtMemOctree.refFrameEdgeValue[colocatedEdgeIdx];
          }

          // code edge
          if (isEncoder) { // encode vertex
            auto vertex = TriSoupVertices[firstVertexToCode];
            encodeOneTriSoupVertexRasterScan(vertex, arithmeticEncoder, ctxtMemOctree, TriSoupVertices, neighbNodes[firstVertexToCode], pattern, interPredictor, colocatedVertex, qualityRef, qualityComp, nbitsVertices, max2bits, mid2bits, firstVertexToCode);
          }
          else
            decodeOneTriSoupVertexRasterScan(arithmeticDecoder, ctxtMemOctree, TriSoupVertices, neighbNodes[firstVertexToCode], pattern, interPredictor, colocatedVertex, qualityRef, qualityComp, nbitsVertices, max2bits, mid2bits, firstVertexToCode);

          xForedgeOfVertex.pop();
          edgePattern.pop();
          if (isInter)
            TriSoupVerticesPred.pop();
          firstVertexToCode++;
        }

        // centroid processing
        int upperxForVertex = !nextIsAvailable() ? INT32_MAX : currWedgePos[0] - 2 * blockWidth;
        while (nodeIdxC < leaves.size() && leaves[nodeIdxC].pos[0] < upperxForVertex) {
          auto leaf = leaves[nodeIdxC];

          // centroid predictor
          int8_t colocatedCentroid = 0;
          bool nodeRefExist = false;
          auto nodePos = leaf.pos;
          auto keyCurrent =
            (int64_t(nodePos[0] + keyshift[0]) << 40)
            + (int64_t(nodePos[1] + keyshift[1]) << 20)
            + int64_t(nodePos[2] + keyshift[2]);
          currentFrameNodeKeys.push_back(keyCurrent);
          CentroValue.push_back(0);

          if (interSkipEnabled) {
            while (
              colocatedNodeIdx
              < ctxtMemOctree.refFrameNodeKeys.size() - 1
              && ctxtMemOctree.refFrameNodeKeys[colocatedNodeIdx]
              < keyCurrent)
              colocatedNodeIdx++;

            nodeRefExist =
              ctxtMemOctree.refFrameNodeKeys[colocatedNodeIdx] == keyCurrent;
            if (nodeRefExist)
              colocatedCentroid =
              ctxtMemOctree.refFrameCentroValue[colocatedNodeIdx];
          }

          generateCentroidsInNodeRasterScan(
            leaf, TriSoupVertices, idxSegment, sliceBB,
            isCentroidDriftActivated, segmentUniqueIndex, qualityRef,
            qualityComp, nodeRefExist, colocatedCentroid, CentroValue);
          if (isFaceVertexActivated) {
            generateFaceVerticesInNodeRasterScan(leaves, nodeIdxC, isEncoder);
          }
          nodeIdxC++;
        }

        // rendering by TriSoup triangles
        int upperxForRendering = !nextIsAvailable() ?  INT32_MAX  : currWedgePos[0] - 3 * blockWidth;
        while (firstNodeToRender < leaves.size()
          && leaves[firstNodeToRender].pos[0] < upperxForRendering) {
          auto leaf = leaves[firstNodeToRender];

          int nPointsInBlock =
            generateTrianglesInNodeRasterScan(
              leaf, firstNodeToRender, renderedBlock, sliceBB, haloTriangle,
              thickness, isFaceVertexActivated);
          flush2PointCloud(
            nRecPoints, renderedBlock.begin(), nPointsInBlock,
            isEncoder ? recPointCloud : pointCloud);
          firstNodeToRender++;
        }
      } // end if on slice chnage

      lastWedgex = currWedgePos[0];
    } // end while loop on wedges
  }

  // ------------------------------------------------------------------------
  void finishSlice() {
    // store edges for colocated next frame
    ctxtMemOctree.refFrameEdgeKeys = currentFrameEdgeKeys;
    ctxtMemOctree.refFrameEdgeValue = TriSoupVertices;

    // store centroids for colocated next frame
    ctxtMemOctree.refFrameNodeKeys = currentFrameNodeKeys;
    ctxtMemOctree.refFrameCentroValue = CentroValue;

    if (isEncoder) {
      // copy reconstructed point cloud to point cloud
      recPointCloud.resize(nRecPoints);
      pointCloud.resize(0);
      pointCloud.swap(recPointCloud);
    }
    else {
      assert(nRecPoints == pointCloud.getPointCount());
    }
    clearTrisoupElements();
  }

private:
  //---------------------------------------------------------------------------
  bool nextIsAvailable() const { return edgesNeighNodes[0] < leaves.size(); }

  //---------------------------------------------------------------------------
  bool changeSlice() const {
    return  (!nextIsAvailable()) || currWedgePos[0] > lastWedgex;
  }

  //---------------------------------------------------------------------------
  void goNextWedge(const std::array<bool, 8>& isNeigbourSane) {
    if (isNeigbourSane[0])
      edgesNeighNodes[7] = edgesNeighNodes[0];

    // move ++ sane neigbours
    for (int i = 0; i < 7; i++)
      edgesNeighNodes[i] += isNeigbourSane[i];

    if (edgesNeighNodes[0] >= leaves.size())
      return;

    currWedgePos = leaves[edgesNeighNodes[0]].pos - offsets[0];
    for (int i = 1; i < 7; ++i) {
      if (edgesNeighNodes[i] >= leaves.size())
        break;
      auto wedgePos = leaves[edgesNeighNodes[i]].pos - offsets[i];
      if (currWedgePos > wedgePos) {
        currWedgePos = wedgePos;
      }
    }
  }
};


//============================================================================

}  // namespace pcc
