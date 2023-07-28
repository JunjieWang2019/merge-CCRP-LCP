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
#include <queue>

#include "geometry_trisoup.h"

#include "pointset_processing.h"
#include "geometry.h"
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
void
decodeGeometryTrisoup(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMemOctree,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const Vec3<int> minimum_position,
  PCCPointSet3& compensatedPointCloud,
  std::vector<MotionVector>& motionVectors)
{
  // trisoup uses octree coding until reaching the triangulation level.
  // todo(df): pass trisoup node size rather than 0?
  std::vector<PCCOctree3Node> nodes;
  decodeGeometryOctree(
    gps, gbh, 0, pointCloud, ctxtMemOctree, arithmeticDecoder, &nodes,
    refFrame, minimum_position, compensatedPointCloud, motionVectors);

  std::cout << "\nSize compensatedPointCloud for TriSoup = " << compensatedPointCloud.getPointCount() << "\n";
  bool isInter = gbh.interPredictionEnabledFlag;

  int blockWidth = 1 << gbh.trisoupNodeSizeLog2(gps);
  const int maxVertexPrecisionLog2 = gbh.trisoup_vertex_quantization_bits ? gbh.trisoup_vertex_quantization_bits : gbh.trisoupNodeSizeLog2(gps);
  const int bitDropped =  std::max(0, gbh.trisoupNodeSizeLog2(gps) - maxVertexPrecisionLog2);

  std::cout << "Number of nodes for TriSoup = " << nodes.size() << "\n";

  // Determine neighbours
  int nSegments = 0;
  codeAndRenderTriSoupRasterScan(nodes, blockWidth, pointCloud, false, bitDropped, 1 /*distanceSearchEncoder*/,
    isInter, refFrame->cloud, compensatedPointCloud,  gps, gbh, NULL, arithmeticDecoder, ctxtMemOctree, nSegments);

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

  const int startCorner[12] = { POS_000, POS_000, POS_0W0, POS_W00, POS_000, POS_0W0, POS_WW0, POS_W00, POS_00W, POS_00W, POS_0WW, POS_W0W };
  const int endCorner[12] = { POS_W00, POS_0W0, POS_WW0, POS_WW0, POS_00W, POS_0WW, POS_WWW, POS_W0W, POS_W0W, POS_0WW, POS_WWW, POS_WWW };
  const int LUTsegmentDirection[12] = { 0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1 };

  std::array<int, 8> edgesNeighNodes; // neighboring nodes' index
  // The 7 firsts are used for unique segments generation/iteration
  // All the 8 are used for contextual information to be used by edge entropy coder
  Vec3<int32_t> currWedgePos;
  int lastWedgex;

  const std::vector<PCCOctree3Node>& leaves;
  PCCPointSet3& pointCloud;
  const bool isEncoder;
  const int bitDropped;
  const int distanceSearchEncoder;
  const bool isInter;
  const bool interSkipEnabled;
  const PCCPointSet3& refPointCloud;
  const PCCPointSet3& compensatedPointCloud;

  // for coding
  const GeometryParameterSet& gps;
  const GeometryBrickHeader& gbh;
  pcc::EntropyEncoder* arithmeticEncoder;
  pcc::EntropyDecoder& arithmeticDecoder;
  GeometryOctreeContexts& ctxtMemOctree;

  RasterScanTrisoupEdges(const std::vector<PCCOctree3Node>& leaves, int blockWidth, PCCPointSet3& pointCloud, bool isEncoder,
    int bitDropped, int distanceSearchEncoder, bool isInter, const PCCPointSet3& refPointCloud, const PCCPointSet3& compensatedPointCloud,
    const GeometryParameterSet& gps, const GeometryBrickHeader& gbh, pcc::EntropyEncoder* arithmeticEncoder, pcc::EntropyDecoder& arithmeticDecoder, GeometryOctreeContexts& ctxtMemOctree)
  : leaves(leaves)
  , blockWidth(blockWidth)
  , pointCloud(pointCloud)
  , isEncoder(isEncoder)
  , bitDropped(bitDropped)
  , distanceSearchEncoder(distanceSearchEncoder)
  , isInter(isInter)
  , interSkipEnabled(isInter && gps.trisoup_skip_mode_enabled_flag)
  , refPointCloud(refPointCloud)
  , compensatedPointCloud(compensatedPointCloud)
  , gps(gps)
  , gbh(gbh)
  , arithmeticEncoder(arithmeticEncoder)
  , arithmeticDecoder(arithmeticDecoder)
  , ctxtMemOctree(ctxtMemOctree)
  , currWedgePos(leaves.empty() ? Vec3<int32_t>{0,0,0} : leaves[0].pos)
  , edgesNeighNodes {0,0,0,0,0,0,0,0}
  {}

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
    if(interSkipEnabled && (vertex >= 0) == (colocatedVertex >= 0 ? 1 : 0)) {
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
    if (nPointInCloud <= nRecPoints + nPointInNode)
      recPointCloud.resize(nRecPoints + nPointInNode + PC_PREALLOCATION_SIZE);

    for (auto it = itBegingBlock; it != last; it++)
      recPointCloud[nRecPoints++] = { int(*it >> 40), int(*it >> 20) & 0xFFFFF, int(*it) & 0xFFFFF };
  }

  //---------------------------------------------------------------------------
  int generateTrianglesInNodeRasterScan(
    const PCCOctree3Node& leaf,
    std::vector<int64_t>& renderedBlock,
    const std::vector<int8_t>& TriSoupVertices,
    int& idxSegment,
    Box3<int32_t>& sliceBB,
    const bool isCentroidDriftActivated,
    int haloTriangle,
    int thickness,
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
    int nPointsInBlock = 0;

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

      // vertex to list of points
      if (bitDropped) {
        Vec3<int32_t> foundvoxel = (point + truncateValue) >> kTrisoupFpBits;
        if (boundaryinsidecheck(foundvoxel, blockWidth - 1)) {
          Vec3<int64_t> renderedPoint = nodepos + foundvoxel;
          renderedBlock[nPointsInBlock++] = (renderedPoint[0] << 40) + (renderedPoint[1] << 20) + renderedPoint[2];
        }
      }
    }

    // Skip leaves that have fewer than 3 vertices.
    int triCount = (int)leafVertices.size();
    if (triCount < 3)
      return nPointsInBlock;

    // compute centroid
    Vec3<int32_t> blockCentroid;
    int dominantAxis;
    determineCentroidAndDominantAxis(blockCentroid, dominantAxis, leafVertices, nodew);

    // Refinement of the centroid along the domiannt axis
    if (triCount > 3 && isCentroidDriftActivated) {
      int bitDropped2 = bitDropped;

      // colocated centroid
      bool possibleSKIPRef = nodeRefExist && badQualityRef <= 1;
      int driftRef = possibleSKIPRef ? colocatedCentroid : 0;

      int lowBound, highBound, lowBoundSurface, highBoundSurface, ctxMinMax;
      Vec3<int32_t> normalV = determineCentroidNormalAndBounds(lowBound, highBound, lowBoundSurface, highBoundSurface, ctxMinMax, bitDropped, bitDropped2, triCount, blockCentroid, dominantAxis, leafVertices, nodew[dominantAxis], blockWidth);

      CentroidInfo centroidInfo;
      determineCentroidPredictor(centroidInfo, bitDropped2, normalV, blockCentroid, nodepos, leaf.isCompensated ? compensatedPointCloud : refPointCloud, leaf.predStart, leaf.predEnd, lowBound, highBound, badQualityComp, badQualityRef, driftRef, possibleSKIPRef);

      int driftQ = 0;
      if (isEncoder) { // encode centroid residual
        driftQ = determineCentroidResidual(bitDropped2, normalV, blockCentroid, nodepos, pointCloud, leaf.start, leaf.end, lowBound, highBound);
        // naive RDO
        if (centroidInfo.possibleSKIP && std::abs(driftQ - centroidInfo.driftSKIP) <= 1)
          driftQ = centroidInfo.driftSKIP;

        encodeCentroidResidual(driftQ, arithmeticEncoder, ctxtMemOctree, centroidInfo, ctxMinMax, lowBoundSurface, highBoundSurface, lowBound, highBound);
      }
      else { // decode centroid residual
        driftQ = decodeCentroidResidual(&arithmeticDecoder, ctxtMemOctree, centroidInfo, ctxMinMax, lowBoundSurface, highBoundSurface, lowBound, highBound);
      }

      // store centroid residual for next frame
      CentroValue.back() = driftQ;

      // dequantize and apply drift
      int driftDQ = 0;
      if (driftQ) {
        driftDQ = std::abs(driftQ) << bitDropped2 + 6;
        int half = 1 << 5 + bitDropped2;
        int DZ = 2 * half * 341 >> 10; // 2 * half / 3;
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
      if (boundaryinsidecheck(foundvoxel, blockWidth - 1)) {
        Vec3<int64_t> renderedPoint = nodepos + foundvoxel;
        renderedBlock[nPointsInBlock++]= (renderedPoint[0] << 40) + (renderedPoint[1] << 20) + renderedPoint[2];
      }
    }

    Vec3<int32_t> v2 = triCount == 3 ? leafVertices[2].pos : blockCentroid;
    Vec3<int32_t> v1 = leafVertices[0].pos;
    for (int triIndex = 0; triIndex < (triCount == 3 ? 1 : triCount); triIndex++) {
      int j1 = triIndex == triCount - 1 ? 0 : triIndex + 1;
      Vec3<int32_t> v0 = v1;
      v1 = leafVertices[j1].pos;

      // choose ray direction
      Vec3<int32_t> edge1 = v1 - v0;
      Vec3<int32_t> edge2 = v2 - v0;
      Vec3<int32_t> a = crossProduct(edge2, edge1) >> kTrisoupFpBits;
      Vec3<int32_t> h = a.abs();
      int directionOk = (h[0] > h[1] && h[0] > h[2]) ? 0 : h[1] > h[2] ? 1 : 2;

      // check if ray tracing is valid; if not skip triangle which is too small
      if (h[directionOk] <= kTrisoupFpOne) // < 2*kTrisoupFpOne should be ok
        continue;

      int64_t inva = divApprox(int64_t(1) << precDivA,  std::abs(a[directionOk]), 0);
      inva = a[directionOk] > 0 ? inva : -inva;

      // range
      int minRange[3];
      int maxRange[3];
      for (int k = 0; k < 3; k++) {
        minRange[k] = std::max(0, std::min(std::min(v0[k], v1[k]), v2[k]) + truncateValue >> kTrisoupFpBits);
        maxRange[k] = std::min(blockWidth - 1, std::max(std::max(v0[k], v1[k]), v2[k]) + truncateValue >> kTrisoupFpBits);
      }
      Vec3<int32_t> s0 = { (minRange[0] << kTrisoupFpBits) - v0[0], (minRange[1] << kTrisoupFpBits) - v0[1], (minRange[2] << kTrisoupFpBits) - v0[2] };

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


  //---------------------------------------------------------------------------
  void buildSegments(int& nSegments)
  {
    point_t BBorig = gbh.geomBoxOrigin;
    point_t keyshift = BBorig + (1 << 18);

    std::vector<int8_t> TriSoupVertices;
    std::vector<uint16_t> neighbNodes;
    std::queue<std::array<int, 18>> edgePattern;
    std::queue<int8_t> TriSoupVerticesPred;
    std::vector<int> segmentUniqueIndex;

    neighbNodes.reserve(leaves.size() * 12); // at most 12 edges per node (to avoid reallocations)
    segmentUniqueIndex.resize(12 * leaves.size(), -1); // temporarily set to -1 to check everything is working
    // TODO: set to -1 could be removed when everything will work properly

    // for slice tracking
    int uniqueIndex = 0;
    lastWedgex = currWedgePos[0];
    int firstVertexToCode = 0;
    int firstNodeToRender = 0;

    // for edge coding
    std::queue<int> xForedgeOfVertex;
    const int nbitsVertices = gbh.trisoupNodeSizeLog2(gps) - bitDropped;
    const int max2bits = nbitsVertices > 1 ? 3 : 1;
    const int mid2bits = nbitsVertices > 1 ? 2 : 1;

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
    recPointCloud.addRemoveAttributes(pointCloud);
    if (isEncoder)
      recPointCloud.resize(pointCloud.getPointCount());
    else
      recPointCloud.resize(PC_PREALLOCATION_SIZE);
    int nRecPoints = 0;
    Box3<int32_t> sliceBB;
    sliceBB.min = gbh.slice_bb_pos << gbh.slice_bb_pos_log2_scale;
    sliceBB.max = sliceBB.min + (gbh.slice_bb_width << gbh.slice_bb_width_log2_scale);

    int idxSegment = 0;
    std::vector<int64_t> renderedBlock(blockWidth * blockWidth * 16, 0) ;

    bool haloFlag = gbh.trisoup_halo_flag;
    bool adaptiveHaloFlag = gbh.trisoup_adaptive_halo_flag;
    int thickness = gbh.trisoup_thickness;
    bool isCentroidDriftActivated = gbh.trisoup_centroid_vertex_residual_flag;

    // halo
    int haloTriangle = 0;
    if (haloFlag) {
      haloTriangle = (((1 << bitDropped) - 1) << kTrisoupFpBits) / blockWidth; // this division is a shift if width is a power of 2
      haloTriangle = (haloTriangle * 28) >> 5; // / 32;
      haloTriangle = haloTriangle > 36 ? 36 : haloTriangle;
    }

    while (nextIsAvailable()) { // this a loop on start position of edges; 3 edges along x,y and z per start position
      // process current wedge position

      std::array<bool, 8> isNeigbourSane;
      for (int i=0; i<8; ++i) // sanity of neighbouring nodes
        isNeigbourSane[i] = edgesNeighNodes[i] < leaves.size() && currWedgePos + offsets[i] == leaves[edgesNeighNodes[i]].pos;

      for (int dir=0; dir<3; ++dir) { // this the loop on the 3 edges along z, then y, then x
        bool processedEdge = false;
        uint16_t neighboursMask = 0;
        std::array<int, 18> pattern {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

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

        for (int neighIdx=0; neighIdx<4; ++neighIdx) { // this the loop on the 4 nodes to interesect the edge
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
                if (voxel[dir1] ==  pos1 && voxel[dir2] == pos2) {
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
              const PCCPointSet3& PC = leaves[neighbNodeIndex].isCompensated ? compensatedPointCloud : refPointCloud;
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
          for (int prevNeighIdx=0; prevNeighIdx<4; ++prevNeighIdx) {
            int wedgeNeighNodeIdx = wedgeNeighNodesIdx[dir][prevNeighIdx];
            if (isNeigbourSane[wedgeNeighNodeIdx]) {
              // update current mask from nodes touching wedge
              neighboursMask |= wedgeNeighMask[dir][prevNeighIdx];
              if (segmentUniqueIdxPrevEdge == -1) {
                int idx = edgesNeighNodes[wedgeNeighNodeIdx] * 12 + wedgeNeighNodesEdgeIdx[dir][prevNeighIdx];
                segmentUniqueIdxPrevEdge = segmentUniqueIndex[idx];
                pattern[0] = segmentUniqueIdxPrevEdge;
                //assert(segmentUniqueIdxPrevEdge != -1);
                for (int neighIdx=0; neighIdx<4; ++neighIdx) {
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

      // code vertices adn rfendering of preceding slices in case the loop has moved up one slice or if finished
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

        // rendering by TriSoup triangles
        int upperxForRendering = !nextIsAvailable() ? INT32_MAX : currWedgePos[0] - 2*blockWidth;
        while (firstNodeToRender < leaves.size() && leaves[firstNodeToRender].pos[0] < upperxForRendering) {
          auto leaf = leaves[firstNodeToRender];

          // centroid predictor
          int8_t colocatedCentroid = 0;
          bool nodeRefExist = false;
          auto nodePos = leaf.pos;
          auto keyCurrent = (int64_t(nodePos[0] + keyshift[0]) << 40) + (int64_t(nodePos[1] + keyshift[1]) << 20) + int64_t(nodePos[2] + keyshift[2]);
          currentFrameNodeKeys.push_back(keyCurrent);
          CentroValue.push_back(0);

          if (interSkipEnabled) {
            while (colocatedNodeIdx < ctxtMemOctree.refFrameNodeKeys.size() - 1 && ctxtMemOctree.refFrameNodeKeys[colocatedNodeIdx] < keyCurrent)
              colocatedNodeIdx++;

            nodeRefExist = ctxtMemOctree.refFrameNodeKeys[colocatedNodeIdx] == keyCurrent;
            if (nodeRefExist)
              colocatedCentroid = ctxtMemOctree.refFrameCentroValue[colocatedNodeIdx];
          }

          int nPointsInBlock = generateTrianglesInNodeRasterScan(leaf, renderedBlock, TriSoupVertices, idxSegment, sliceBB, isCentroidDriftActivated, haloTriangle, thickness, segmentUniqueIndex, qualityRef, qualityComp, nodeRefExist, colocatedCentroid, CentroValue);

          flush2PointCloud(nRecPoints, renderedBlock.begin(), nPointsInBlock, recPointCloud);
          firstNodeToRender++;
        }
      } // end if on slice chnage

      lastWedgex = currWedgePos[0];
    } // end while loop on wedges

    // store edges for colocated next frame
    ctxtMemOctree.refFrameEdgeKeys = currentFrameEdgeKeys;
    ctxtMemOctree.refFrameEdgeValue = TriSoupVertices;

    // store centroids for colocated next frame
    ctxtMemOctree.refFrameNodeKeys = currentFrameNodeKeys;
    ctxtMemOctree.refFrameCentroValue = CentroValue;

    // copy reconstructed point cloud to point cloud
    recPointCloud.resize(nRecPoints);
    pointCloud.resize(0);
    pointCloud = std::move(recPointCloud);
    nSegments = TriSoupVertices.size();
  }

private:
  //---------------------------------------------------------------------------
  bool nextIsAvailable() const { return edgesNeighNodes[0] < leaves.size(); }

  //---------------------------------------------------------------------------
  bool changeSlice() const {
    return !nextIsAvailable() || currWedgePos[0] > lastWedgex;}

  //---------------------------------------------------------------------------
  void goNextWedge(const std::array<bool, 8>& isNeigbourSane) {
    if (isNeigbourSane[0])
      edgesNeighNodes[7] = edgesNeighNodes[0];

    // move ++ sane neigbours
    for (int i=0; i<7; i++)
        edgesNeighNodes[i] += isNeigbourSane[i];

    if (edgesNeighNodes[0] >= leaves.size())
      return;

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
};


//---------------------------------------------------------------------------
void codeAndRenderTriSoupRasterScan(
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
  int &nSegments) {

  const int32_t blockWidth = defaultBlockWidth; // Width of block. In future, may override with leaf blockWidth
  RasterScanTrisoupEdges rste(leaves, blockWidth, pointCloud, isEncoder, bitDropped, distanceSearchEncoder, isInter, refPointCloud, compensatedPointCloud, gps, gbh, arithmeticEncoder, arithmeticDecoder, ctxtMemOctree);
  rste.buildSegments(nSegments);
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
  blockCentroid = blockCentroid * LUTdivtriCount[triCount] >> 10;

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
  int nodewDominant,
  int blockWidth)
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
  int64_t invNormN = irsqrt(accuNormal[0] * accuNormal[0] + accuNormal[1] * accuNormal[1] + accuNormal[2] * accuNormal[2]);
  Vec3<int32_t> normalV = accuNormal  * invNormN >> 40 - kTrisoupFpBits;

  // drift bounds
  ctxMinMax = std::min(8, (maxPos - minPos) >> (kTrisoupFpBits + bitDropped));
  int boundL = -kTrisoupFpHalf;
  int boundH = ((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1;
  int m = 1;
  int half = 1 << 5 + bitDropped2;
  int DZ = 2 * half * 341 >> 10; ; // 2 * half / 3
  for (; m < nodewDominant; m++) {
    int driftDm = m << bitDropped2 + 6;
    driftDm += DZ - half;

    Vec3<int32_t> temp = blockCentroid + (driftDm * normalV >> 6);
    if (temp[0]<boundL || temp[1]<boundL || temp[2]<boundL || temp[0]>boundH || temp[1]>boundH || temp[2]> boundH)
      break;
  }
  highBound = m - 1;

  m = 1;
  for (; m < nodewDominant; m++) {
    int driftDm = m << bitDropped2 + 6;
    driftDm += DZ - half;

    Vec3<int32_t> temp = blockCentroid + (-driftDm * normalV >> 6);
    if (temp[0]<boundL || temp[1]<boundL || temp[2]<boundL || temp[0]>boundH || temp[1]>boundH || temp[2]> boundH)
      break;
  }
  lowBound = m - 1;
  lowBoundSurface = std::max(0, ((blockCentroid[dominantAxis] - minPos) + kTrisoupFpHalf >> kTrisoupFpBits) + halfDropped2 >> bitDropped2);
  highBoundSurface = std::max(0, ((maxPos - blockCentroid[dominantAxis]) + kTrisoupFpHalf >> kTrisoupFpBits) + halfDropped2 >> bitDropped2);

  return normalV;
}

// --------------------------------------------------------------------------
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
  bool possibleSKIPRef)
{
  int driftQPred = -100;
  int driftQComp = -100;
  // determine quantized drift for predictor
  if (end > start) {
    int driftPred = 0;
    driftQPred = 0;
    driftQComp = 0;
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
      driftPred = divApprox(driftPred >> kTrisoupFpBits - 6, counter, 0); // drift is shift by 6 bits
    }

    int half = 1 << 5 + bitDropped2;
    int DZ = 2 * half * 341 >> 10; //2 * half / 3;

    if (abs(driftPred) >= DZ) {
      driftQPred = (abs(driftPred) - DZ + 2 * half) >> 6 + bitDropped2 - 1;
      driftQComp = (abs(driftPred) - DZ + 2 * half) >> 6 + bitDropped2;
      if (driftPred < 0) {
        driftQPred = -driftQPred;
        driftQComp = -driftQComp;
      }
    }
    driftQPred = std::min(std::max(driftQPred, -2 * lowBound), 2 * highBound);  // drift in [-lowBound; highBound] but quantization is twice better
    driftQComp = std::min(std::max(driftQComp, -lowBound), highBound);  // drift in [-lowBound; highBound]
  }

  bool possibleSKIPComp = driftQComp != -100 && badQualityComp <= 0;
  int driftComp = possibleSKIPComp ? driftQComp : 0;

  bool possibleSKIP = possibleSKIPRef || possibleSKIPComp;
  int driftSKIP = 0;
  int qualitySKIP = badQualityRef;
  if (possibleSKIP) {
    if (possibleSKIPRef)
      driftSKIP = driftRef;

    if (possibleSKIPComp && (!possibleSKIPRef || badQualityComp < badQualityRef)) {
      driftSKIP = driftComp;
      qualitySKIP = badQualityComp;
    }
  }

  centroidInfo.driftQPred = driftQPred;
  centroidInfo.possibleSKIP = possibleSKIP;
  centroidInfo.driftSKIP = driftSKIP;
  centroidInfo.qualitySKIP = qualitySKIP;
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
    drift = divApprox(drift >> kTrisoupFpBits - 6, counter, 0); // drift is shift by 6 bits
  }

  int half = 1 << 5 + bitDropped2;
  int DZ = 2 * half * 341 >> 10; //2 * half / 3;

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
  CentroidInfo centroidInfo,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound)
{
  int driftQPred = centroidInfo.driftQPred;
  bool possibleSKIP = centroidInfo.possibleSKIP;
  int driftSKIP = centroidInfo.driftSKIP;

  if (possibleSKIP) {
    arithmeticEncoder->encode(driftQ == driftSKIP, ctxtMemOctree.ctxDriftSKIP[centroidInfo.qualitySKIP][driftSKIP == 0]);
    if (driftQ == driftSKIP)
      return;
  }

  bool isNotZero = false;
  isNotZero = possibleSKIP && driftSKIP == 0;
  if (!isNotZero) {
    if (driftQPred == -100) //intra
       arithmeticEncoder->encode(driftQ == 0, ctxtMemOctree.ctxDrift0[ctxMinMax][0]);
    else { //inter
      if (possibleSKIP)
        arithmeticEncoder->encode(driftQ == 0, ctxtMemOctree.ctxDrift0Skip[std::min(3, std::abs(driftSKIP))]);
      else
        arithmeticEncoder->encode(driftQ == 0, ctxtMemOctree.ctxDrift0[ctxMinMax][1 + std::min(3, std::abs(driftQPred))]);
    }
  }


  // if not 0, drift in [-lowBound; highBound]
  if (driftQ) {
    // code sign
    if (highBound && lowBound) {  // otherwise sign is known
      int lowS = std::min(7, lowBoundSurface);
      int highS = std::min(7, highBoundSurface);
      if (possibleSKIP)
        arithmeticEncoder->encode(driftQ > 0, ctxtMemOctree.ctxDriftSignSkip[driftSKIP == 0 ? 0 : 1 + (driftSKIP > 0) + 2 * (std::abs(driftSKIP) > 1)]);
      else
        arithmeticEncoder->encode(driftQ > 0, ctxtMemOctree.ctxDriftSign[lowBound == highBound ? 0 : 1 + (lowBound < highBound)][lowS][highS][(driftQPred && driftQPred != -100) ? 1 + (driftQPred > 0) : 0]);
    }

    // code remaining bits 1 to 7 at most
    int magBound = (driftQ > 0 ? highBound : lowBound) - 1;
    bool sameSignPred = driftQPred != -100 && (driftQPred > 0 && driftQ > 0) || (driftQPred < 0 && driftQ < 0);

    int magDrift = std::abs(driftQ) - 1;
    int ctx = 0;
    int ctx2 = driftQPred != -100 ? 1 + std::min(8, sameSignPred * std::abs(driftQPred)) : 0;
    while (magBound > 0 && magDrift >= 0) {
      if (ctx < 4)
        arithmeticEncoder->encode(magDrift == 0, ctxtMemOctree.ctxDriftMag[ctx][ctx2]);
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
  CentroidInfo centroidInfo,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound)
{
  int driftQPred = centroidInfo.driftQPred;
  bool possibleSKIP = centroidInfo.possibleSKIP;
  int driftSKIP = centroidInfo.driftSKIP;

  bool isSKIP = false;
  if (possibleSKIP) {
    isSKIP = arithmeticDecoder->decode(ctxtMemOctree.ctxDriftSKIP[centroidInfo.qualitySKIP][driftSKIP == 0]);
    if (isSKIP)
      return driftSKIP;
  }

  bool isNotZero = false;
  isNotZero = possibleSKIP && driftSKIP == 0;
  int driftQ = 1;
  if (!isNotZero) {
    if (driftQPred == -100) //intra
      driftQ = arithmeticDecoder->decode(ctxtMemOctree.ctxDrift0[ctxMinMax][0]) ? 0 : 1;
    else {//inter
      if (possibleSKIP)
        driftQ = arithmeticDecoder->decode(ctxtMemOctree.ctxDrift0Skip[std::min(3, std::abs(driftSKIP))]) ? 0 : 1;
      else
        driftQ = arithmeticDecoder->decode(ctxtMemOctree.ctxDrift0[ctxMinMax][1 + std::min(3, std::abs(driftQPred))]) ? 0 : 1;
    }
    if (driftQ == 0)
      return 0;
  }

  // if not 0, drift in [-lowBound; highBound]
  // code sign
  int sign = 1;
  if (highBound && lowBound) {// otherwise sign is knwow
    int lowS = std::min(7, lowBoundSurface);
    int highS = std::min(7, highBoundSurface);
    if (possibleSKIP)
      sign = arithmeticDecoder->decode(ctxtMemOctree.ctxDriftSignSkip[driftSKIP == 0 ? 0 : 1 + (driftSKIP > 0) + 2 * (std::abs(driftSKIP) > 1)]);
    else
      sign = arithmeticDecoder->decode(ctxtMemOctree.ctxDriftSign[lowBound == highBound ? 0 : 1 + (lowBound < highBound)][lowS][highS][(driftQPred && driftQPred != -100) ? 1 + (driftQPred > 0) : 0]);
  }
  else if (!highBound) // highbound is 0 , so sign is negative; otherwise sign is already set to positive
    sign = 0;

  // code remaining bits 1 to 7 at most
  int magBound = (sign ? highBound : lowBound) - 1;
  bool sameSignPred = driftQPred != -100 && (driftQPred > 0 && sign) || (driftQPred < 0 && !sign);

  int ctx = 0;
  int ctx2 = driftQPred != -100 ? 1 + std::min(8, sameSignPred * std::abs(driftQPred)) : 0;
  while (magBound > 0) {
    int bit;
    if (ctx < 4)
      bit = arithmeticDecoder->decode(ctxtMemOctree.ctxDriftMag[ctx][ctx2]);
    else
      bit = arithmeticDecoder->decode();

    if (bit) // magDrift==0 and magnitude coding is finished
      break;

    driftQ++;
    magBound--;
    ctx++;
  }

  return sign ? driftQ : -driftQ;
}

// ---------------------------------------------------------------------------
void
constructCtxInfo(
  codeVertexCtxInfo& ctxInfo,
  int neigh,
  std::array<int, 18>& patternIdx,
  std::vector<int8_t>& TriSoupVertices,
  int nbitsVertices,
  int max2bits,
  int mid2bits,
  std::vector<int8_t>& qualityRef,
  std::vector<int8_t>& qualityComp) {

  // node info
  ctxInfo.ctxE = (!!(neigh & 1)) + (!!(neigh & 2)) + (!!(neigh & 4)) + (!!(neigh & 8)) - 1; // at least one node is occupied
  ctxInfo.ctx0 = (!!(neigh & 16)) + (!!(neigh & 32)) + (!!(neigh & 64)) + (!!(neigh & 128));
  ctxInfo.ctx1 = (!!(neigh & 256)) + (!!(neigh & 512)) + (!!(neigh & 1024)) + (!!(neigh & 2048));
  int direction = neigh >> 13; // 0=x, 1=y, 2=z
  ctxInfo.direction = direction;

  // colocated info
  for (int v = 0; v < 18; v++) {
    if (patternIdx[v] != -1) {
      int idxEdge = patternIdx[v];
      ctxInfo.nBadPredRef += qualityRef[idxEdge] == 0; // presence badly predicted
      ctxInfo.nBadPredRef1 += qualityRef[idxEdge] == 0 || qualityRef[idxEdge] == 2; // first bit pos badly predicted
      ctxInfo.nBadPredRef2 += qualityRef[idxEdge] == 0 || qualityRef[idxEdge] == 2 || qualityRef[idxEdge] == 3; // second bit pos badly predicted

      ctxInfo.nBadPredComp += qualityComp[idxEdge] == 0; // presence badly predicted
      ctxInfo.nBadPredComp1 += qualityComp[idxEdge] == 0 || qualityComp[idxEdge] == 2; // first bit pos badly predicted
      ctxInfo.nBadPredComp2 += qualityComp[idxEdge] == 0 || qualityComp[idxEdge] == 2 || qualityComp[idxEdge] == 3; // second bit pos badly predicted
    }
  }


  //neighbours info
  for (int v = 0; v < 9; v++) {
    int v18 = mapping18to9[direction][v];

    if (patternIdx[v18] != -1) {
      int vertex = TriSoupVertices[patternIdx[v18]];
      if (vertex >= 0) {
        ctxInfo.pattern |= 1 << v;
        int vertexPos2bits = vertex >> std::max(0, nbitsVertices - 2);
        if (towardOrAway[v18])
          vertexPos2bits = max2bits - vertexPos2bits; // reverses for away
        if (vertexPos2bits >= mid2bits)
          ctxInfo.patternClose |= 1 << v;
        if (vertexPos2bits >= max2bits)
          ctxInfo.patternClosest |= 1 << v;
        ctxInfo.nclosestPattern += vertexPos2bits >= max2bits && v <= 4;
      }
    }
  }

  ctxInfo.missedCloseStart = /*!(ctxInfo.pattern & 1) + */  !(ctxInfo.pattern & 2) + !(ctxInfo.pattern & 4);
  ctxInfo.nclosestStart = !!(ctxInfo.patternClosest & 1) + !!(ctxInfo.patternClosest & 2) + !!(ctxInfo.patternClosest & 4);
  if (direction == 0) {
    ctxInfo.missedCloseStart += !(ctxInfo.pattern & 8) + !(ctxInfo.pattern & 16);
    ctxInfo.nclosestStart += !!(ctxInfo.patternClosest & 8) + !!(ctxInfo.patternClosest & 16);
  }
  if (direction == 1) {
    ctxInfo.missedCloseStart += !(ctxInfo.pattern & 8);
    ctxInfo.nclosestStart += !!(ctxInfo.patternClosest & 8) - !!(ctxInfo.patternClosest & 16);
  }
  if (direction == 2) {
    ctxInfo.nclosestStart += -!!(ctxInfo.patternClosest & 8) - !!(ctxInfo.patternClosest & 16);
  }

  // reorganize neighbours of vertex /edge (endpoint) independently on xyz
  ctxInfo.neighbEdge = (neigh >> 0) & 15;
  ctxInfo.neighbEnd = (neigh >> 4) & 15;
  ctxInfo.neighbStart = (neigh >> 8) & 15;
  if (direction == 2) {
    ctxInfo.neighbEdge = ((neigh >> 0 + 0) & 1);
    ctxInfo.neighbEdge += ((neigh >> 0 + 3) & 1) << 1;
    ctxInfo.neighbEdge += ((neigh >> 0 + 1) & 1) << 2;
    ctxInfo.neighbEdge += ((neigh >> 0 + 2) & 1) << 3;

    ctxInfo.neighbEnd = ((neigh >> 4 + 0) & 1);
    ctxInfo.neighbEnd += ((neigh >> 4 + 3) & 1) << 1;
    ctxInfo.neighbEnd += ((neigh >> 4 + 1) & 1) << 2;
    ctxInfo.neighbEnd += ((neigh >> 4 + 2) & 1) << 3;

    ctxInfo.neighbStart = ((neigh >> 8 + 0) & 1);
    ctxInfo.neighbStart += ((neigh >> 8 + 3) & 1) << 1;
    ctxInfo.neighbStart += ((neigh >> 8 + 1) & 1) << 2;
    ctxInfo.neighbStart += ((neigh >> 8 + 2) & 1) << 3;
  }

  ctxInfo.orderedPclosePar = (((ctxInfo.pattern >> 5) & 3) << 2) + (!!(ctxInfo.pattern & 128) << 1) + !!(ctxInfo.pattern & 256);
  ctxInfo.orderedPcloseParPos = (((ctxInfo.patternClose >> 5) & 3) << 2) + (!!(ctxInfo.patternClose & 128) << 1) + !!(ctxInfo.patternClose & 256);
}

// -------------------------------------------------------------------------- -
void
constructCtxPresence(
  int& ctxMap1,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int8_t colocatedVertex) {

  ctxMap1 = std::min(ctxInfo.nclosestPattern, 2) * 15 * 2 + (ctxInfo.neighbEdge - 1) * 2 + ((ctxInfo.ctx1 == 4));    // 2* 15 *3 = 90 -> 7 bits
  ctxMap2 = ctxInfo.neighbEnd << 11;
  ctxMap2 |= (ctxInfo.patternClose & (0b00000110)) << 9 - 1; // perp that do not depend on direction = to start
  ctxMap2 |= ctxInfo.direction << 7;
  ctxMap2 |= (ctxInfo.patternClose & (0b00011000)) << 5 - 3; // perp that  depend on direction = to start or to end
  ctxMap2 |= (ctxInfo.patternClose & (0b00000001)) << 4;  // before
  ctxMap2 |= ctxInfo.orderedPclosePar;

  //ctxInter = isInter ? 1 + (TriSoupVerticesPred >= 0) : 0;
  bool isInterGood = isInter && (ctxInfo.nBadPredRef <= 0 || ctxInfo.nBadPredComp <= 3);
  ctxInter = 0;
  if (isInterGood) {
    ctxInter = 1 + (TriSoupVerticesPred >= 0);

    bool goodRef = ctxInfo.nBadPredRef <= 0;
    if (goodRef) {
      ctxMap2 |= (colocatedVertex >= 0 ? 1 : 0) << 15;
    }
    else
       ctxMap2 <<= 1;
    ctxMap2 |= goodRef << 16;
  }
  else {
    ctxMap2 <<= 2;
  }

}

// -------------------------------------------------------------------------- -
void
constructCtxPos1(
  int& ctxMap1,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int b,
  int8_t colocatedVertex) {

  int ctxFullNbounds = (4 * (ctxInfo.ctx0 <= 1 ? 0 : (ctxInfo.ctx0 >= 3 ? 2 : 1)) + (std::max(1, ctxInfo.ctx1) - 1)) * 2 + (ctxInfo.ctxE == 3);
  ctxMap1 = ctxFullNbounds * 2 + (ctxInfo.nclosestStart > 0);
  ctxMap2 = ctxInfo.missedCloseStart << 8;
  ctxMap2 |= (ctxInfo.patternClosest & 1) << 7;
  ctxMap2 |= ctxInfo.direction << 5;
  ctxMap2 |= ctxInfo.patternClose & (0b00011111);
  ctxMap2 = (ctxMap2 << 4) + ctxInfo.orderedPcloseParPos;

  bool isGoodRef = isInter && colocatedVertex >= 0;
  bool isInterGood = isInter && ((isGoodRef && ctxInfo.nBadPredRef1 <= 4) || ctxInfo.nBadPredComp1 <= 4);

  ctxInter = 0;
  if (isInterGood) {
    ctxInter = TriSoupVerticesPred >= 0 ? 1 + ((TriSoupVerticesPred >> b - 1) & 3) : 0;

    int goodPresence = colocatedVertex >= 0 && ctxInfo.nBadPredRef1 <= 0;
    ctxMap2 |= goodPresence << 16;
    if (goodPresence)
      ctxMap2 |= ((colocatedVertex >> b) & 1) << 15;
  }
  else {
    ctxMap2 <<= 2;
  }
}

// ------------------------------------------------------------------------
void
constructCtxPos2(
  int& ctxMap1,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int b,
  int v,
  int8_t colocatedVertex){

  int ctxFullNbounds = (4 * (ctxInfo.ctx0 <= 1 ? 0 : (ctxInfo.ctx0 >= 3 ? 2 : 1)) + (std::max(1, ctxInfo.ctx1) - 1)) * 2 + (ctxInfo.ctxE == 3);
  ctxMap1 = ctxFullNbounds * 2 + (ctxInfo.nclosestStart > 0);
  ctxMap2 = ctxInfo.missedCloseStart << 8;
  ctxMap2 |= (ctxInfo.patternClose & 1) << 7;
  ctxMap2 |= (ctxInfo.patternClosest & 1) << 6;
  ctxMap2 |= ctxInfo.direction << 4;
  ctxMap2 |= (ctxInfo.patternClose & (0b00011111)) >> 1;
  ctxMap2 = (ctxMap2 << 4) + ctxInfo.orderedPcloseParPos;

  bool isGoodRef = isInter && colocatedVertex >= 0 && ((colocatedVertex >> b + 1) == v);
  bool isInterGood = isInter && ((isGoodRef && ctxInfo.nBadPredRef1 <= 4) || ctxInfo.nBadPredComp2 <= 4);
  ctxInter = 0;
  if (isInter) {
    ctxInter = TriSoupVerticesPred >= 0 ? 1 + ((TriSoupVerticesPred >> b) <= (v << 1)) : 0;

    int goodPresence = colocatedVertex >= 0 ? 1 : 0;
    goodPresence = goodPresence && ((colocatedVertex >> b + 1) == v) && ctxInfo.nBadPredRef2 <= 1;
    ctxMap2 |= goodPresence << 16;
    if (goodPresence)
      ctxMap2 |= ((colocatedVertex >> b) & 1) << 15;
  }
  else {
    ctxMap2 <<= 2;
  }
}


// ---------------------------------------------------------------------------
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

  if (triCount > 3) {
    Vertex vertex;
    const int sIdx1[3] = { 2,2,1 };
    const int sIdx2[3] = { 1,0,0 };
    auto leafVerticesTemp = leafVertices;

    int maxNormTri = 0;
    for (int axis = 0; axis <= 2; axis++) {
      int axis1 = sIdx1[axis];
      int axis2 = sIdx2[axis];
      const int Width_x = blockWidth[axis1] << kTrisoupFpBits;
      const int Width_y = blockWidth[axis2] << kTrisoupFpBits;

      // order along axis
      for (int j = 0; j < triCount; j++) {
        // compute score closckwise
        int x = leafVerticesTemp[j].pos[axis1] + kTrisoupFpHalf; // back to [0,B]^3 for ordering
        int y = leafVerticesTemp[j].pos[axis2] + kTrisoupFpHalf; // back to [0,B]^3 for ordering

        int flag3 = x <= 0;
        int score3 = Width_y - flag3 * y + (!flag3) * x;
        int flag2 = y >= Width_y;
        int score2 = Width_y + Width_x - flag2 * x + (!flag2) * score3;
        int flag1 = x >= Width_x;
        int score = flag1 * y + (!flag1) * score2;
        leafVerticesTemp[j].theta = score;
        leafVerticesTemp[j].tiebreaker = leafVerticesTemp[j].pos[axis]; // stable sort if same score
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
  }

  return dominantAxis;
}

// --------------------------------------------------------
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
  int thickness)
{
  int32_t u0 = ((-s0[1] * edge2[2] + s0[2] * edge2[1]) * inva) >> precDivA;
  Vec3<int32_t> q0 = crossProduct(s0, edge1);
  int32_t v0 = (q0[0] * inva) >> precDivA;
  int32_t t0 = ((edge2 * (q0 >> kTrisoupFpBits)) * inva) >> precDivA;

  int32_t u1 = ((-edge2[2] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t v1 = ((edge1[2] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t t1 = ((edge2[0]* edge1[2] - edge2[2] * edge1[0]) * inva) >> precDivA;

  int32_t u2 = ((edge2[1] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t v2 = ((-edge1[1] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t t2 = ((-edge2[0]* edge1[1] + edge2[1] * edge1[0]) * inva) >> precDivA;

  int64_t renderedPoint1D0 = (int64_t(nodepos[0]) << 40) + (int64_t(nodepos[1]) << 20) + int64_t(nodepos[2]);
  for (int32_t g1 = minRange[1]; g1 <= maxRange[1]; g1++, u0 += u1, v0 += v1, t0 += t1) {
    int32_t u = u0, v = v0, t = t0;
    for (int32_t g2 = minRange[2];  g2 <= maxRange[2]; g2++, u += u2, v += v2, t += t2) {
      int w = kTrisoupFpOne - u - v;
      if (u >= -haloTriangle && v >= -haloTriangle && w >= -haloTriangle) {
        int32_t foundvoxel = minRange[0] + (t + truncateValue >> kTrisoupFpBits);
        if (foundvoxel >= 0 && foundvoxel < blockWidth) {
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(foundvoxel) << 40) + (int64_t(g1) << 20) + int64_t(g2);
        }
        int32_t foundvoxelUp = minRange[0] + (t + thickness + truncateValue >> kTrisoupFpBits);
        if (foundvoxelUp != foundvoxel && foundvoxelUp >= 0 && foundvoxelUp < blockWidth) {
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(foundvoxelUp) << 40) + (int64_t(g1) << 20) + int64_t(g2);
        }
        int32_t foundvoxelDown = minRange[0] + (t - thickness + truncateValue >> kTrisoupFpBits);
        if (foundvoxelDown != foundvoxel && foundvoxelDown >= 0 && foundvoxelDown < blockWidth) {
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(foundvoxelDown) << 40) + (int64_t(g1) << 20) + int64_t(g2);
        }
      }
    }// loop g2
  }//loop g1
}

// --------------------------------------------------------
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
  int thickness)
{
  int32_t u0 = ((s0[0] * edge2[2] - s0[2] * edge2[0]) * inva) >> precDivA;
  Vec3<int32_t> q0 = crossProduct(s0, edge1);
  int32_t v0 = (q0[1] * inva) >> precDivA;
  int32_t t0 = ((edge2 * (q0 >> kTrisoupFpBits)) * inva) >> precDivA;

  int32_t u1 = ((edge2[2] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t v1 = ((-edge1[2] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t t1 = ((-edge2[1] * edge1[2] + edge2[2] * edge1[1]) * inva) >> precDivA;

  int32_t u2 = ((-edge2[0] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t v2 = ((edge1[0] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t t2 = ((-edge2[0] * edge1[1] + edge2[1] * edge1[0]) * inva) >> precDivA;

  int64_t renderedPoint1D0 = (int64_t(nodepos[0]) << 40) + (int64_t(nodepos[1]) << 20) + int64_t(nodepos[2]);
  for (int32_t g1 = minRange[0]; g1 <= maxRange[0]; g1++, u0 += u1, v0 += v1, t0 += t1) {
    int32_t u = u0, v = v0, t = t0;
    for (int32_t g2 = minRange[2]; g2 <= maxRange[2]; g2++, u += u2, v += v2, t += t2) {
      int w = kTrisoupFpOne - u - v;
      if (u >= -haloTriangle && v >= -haloTriangle && w >= -haloTriangle) {
        int32_t foundvoxel = minRange[1] + (t + truncateValue >> kTrisoupFpBits);
        if (foundvoxel >= 0 && foundvoxel < blockWidth) {
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(g1) << 40) + (int64_t(foundvoxel) << 20) + int64_t(g2);
        }
        int32_t foundvoxelUp = minRange[1] + (t + thickness + truncateValue >> kTrisoupFpBits);
        if (foundvoxelUp >= 0 && foundvoxelUp < blockWidth) {
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(g1) << 40) + (int64_t(foundvoxelUp) << 20) + int64_t(g2);
        }
        int32_t foundvoxelDown = minRange[1] + (t - thickness + truncateValue >> kTrisoupFpBits);
        if (foundvoxelDown >= 0 && foundvoxelDown < blockWidth) {
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(g1) << 40) + (int64_t(foundvoxelDown) << 20) + int64_t(g2);
        }
      }
    }// loop g2
  }//loop g1
}

// --------------------------------------------------------
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
  int thickness)
{
  int32_t u0 = ((-s0[0] * edge2[1] + s0[1] * edge2[0]) * inva) >> precDivA;
  Vec3<int32_t> q0 = crossProduct(s0, edge1);
  int32_t v0 = (q0[2] * inva) >> precDivA;
  int32_t t0 = ((edge2 * (q0 >> kTrisoupFpBits)) * inva) >> precDivA;

  int32_t u1 = ((-edge2[1] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t v1 = ((edge1[1] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t t1 = ((edge2[2] * edge1[1] - edge2[1] * edge1[2]) * inva) >> precDivA;

  int32_t u2 = ((edge2[0] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t v2 = ((-edge1[0] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t t2 = ((edge2[0] * edge1[2] - edge2[2] * edge1[0]) * inva) >> precDivA;

  int64_t renderedPoint1D0 = (int64_t(nodepos[0]) << 40) + (int64_t(nodepos[1]) << 20) + int64_t(nodepos[2]);
  for (int32_t g1 = minRange[0]; g1 <= maxRange[0]; g1++, u0 += u1, v0 += v1, t0 += t1) {
    int32_t u = u0, v = v0, t = t0;
    for (int32_t g2 = minRange[1]; g2 <= maxRange[1]; g2++, u += u2, v += v2, t += t2) {
      int w = kTrisoupFpOne - u - v;
      if (u >= -haloTriangle && v >= -haloTriangle && w >= -haloTriangle) {
        int32_t foundvoxel = minRange[2] + (t + truncateValue >> kTrisoupFpBits);
        if (foundvoxel >= 0 && foundvoxel < blockWidth) {
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(g1) << 40) + (int64_t(g2) << 20) + int64_t(foundvoxel);
        }
        int32_t foundvoxelUp = minRange[2] + (t + thickness + truncateValue >> kTrisoupFpBits);
        if (foundvoxelUp != foundvoxel && foundvoxelUp >= 0 && foundvoxelUp < blockWidth) {
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(g1) << 40) + (int64_t(g2) << 20) + int64_t(foundvoxelUp);
        }
        int32_t foundvoxelDown = minRange[2] + (t - thickness + truncateValue >> kTrisoupFpBits);
        if (foundvoxelDown != foundvoxel && foundvoxelDown >= 0 && foundvoxelDown < blockWidth) {
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(g1) << 40) + (int64_t(g2) << 20) + int64_t(foundvoxelDown);
        }
      }
    }// loop g2
  }//loop g1
}

//============================================================================
}  // namespace pcc
