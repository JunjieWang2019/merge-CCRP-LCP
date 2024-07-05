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

#include "PCCTMC3Decoder.h"

namespace pcc {

//============================================================================
void
decodeGeometryTrisoup(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMemOctree,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Decoder3& decoder)
{
  const Vec3<int> minimum_position = sps.seqBoundingBoxOrigin;
  bool isInter = gbh.interPredictionEnabledFlag;

  // prepare TriSoup parameters
  int blockWidth = gbh.trisoupNodeSize(gps);
  std::cout << "TriSoup QP = " << gbh.trisoup_QP << "\n";

  // trisoup uses octree coding until reaching the triangulation level.
  RasterScanTrisoupEdgesDecoder rste(blockWidth, pointCloud,
    1 /*distanceSearchEncoder*/, isInter, interPredParams.compensatedPointCloud,
    gps, gbh, NULL, arithmeticDecoder, ctxtMemOctree);
  rste.useLocalAttr = sps.localized_attributes_enabled_flag;
  if (rste.useLocalAttr) {
    rste.decoder = &decoder;
    rste.slabThickness = sps.localized_attributes_slab_thickness_minus1 + 1;
  }
  rste.init();

  // octree
  decodeGeometryOctree<true>(
    gps, gbh, 0, pointCloud, ctxtMemOctree, arithmeticDecoder, nullptr, refFrame,
    sps, minimum_position, interPredParams, decoder, &rste);

  std::cout << "\nSize compensatedPointCloud for TriSoup = "
    << interPredParams.compensatedPointCloud.getPointCount() << "\n";
  std::cout << "Number of nodes for TriSoup = " << rste.leaves.size() << "\n";
  std::cout << "TriSoup gives " << pointCloud.getPointCount() << " points \n";
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


bool nodeBoundaryInsideCheck(Vec3<int32_t> bw, Vec3<int32_t> pt)
{
  if (0 <= pt[0] && pt[0] <= bw[0]
      && 0 <= pt[1] && pt[1] <= bw[1]
      && 0 <= pt[2] && pt[2] <= bw[2]) {
    return true;
  }
  return false;
}


// --------------------------------------------------------------------------
static const int LUTsqrt[13] = { 0, 256, 362, 443, 512, 572, 627, 677, 724, 768, 810, 849, 887 };

void
determineCentroidAndDominantAxis(
  Vec3<int32_t> &blockCentroid,
  int &dominantAxis,
  std::vector<Vertex> &leafVertices,
  Vec3<int32_t> nodew,
  bool &flagCentroOK,
  int &scaleQ,
  int stepQcentro)

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

  // >> bitDropped  is like / (step/256)
  int64_t a = int64_t(7680) * nodew[0] * 256;  // cst * node size
  int64_t b = int64_t(LUTsqrt[triCount]) * Wtotal * stepQcentro;  // sqrt(#vertices) * triangle area * Qstep
  int ratio = divApprox(a << 16, b, 0); // cst * node size / sqrt(#vertices) / triangle area / Qstep
  scaleQ = std::max(190, ratio); // precision on 8 bits

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
  int triCount,
  Vec3<int32_t> blockCentroid,
  int dominantAxis,
  std::vector<Vertex>& leafVertices,
  int nodewDominant,
  int blockWidth,
  int stepQcentro,
  int scaleQ)
{
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
  ctxMinMax = std::min(8, (maxPos - minPos) / stepQcentro);
  int boundL = -kTrisoupFpHalf;
  int boundH = ((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1;
  int m = 1;

  int half = stepQcentro >> 1;
  int DZ = 682 * half >> 10; // 2 * half / 3;

  for (; m < nodewDominant; m++) {
    int driftDm = m * stepQcentro;
    driftDm += DZ - half;
    driftDm = driftDm * scaleQ  >> 8;

    Vec3<int32_t> temp = blockCentroid + (driftDm * normalV >> 8); // >>6
    if (temp[0]<boundL || temp[1]<boundL || temp[2]<boundL || temp[0]>boundH || temp[1]>boundH || temp[2]> boundH)
      break;
  }
  highBound = m - 1;

  m = 1;
  for (; m < nodewDominant; m++) {
    int driftDm = m * stepQcentro;
    driftDm += DZ - half;
    driftDm = driftDm * scaleQ >> 8;

    Vec3<int32_t> temp = blockCentroid + (-driftDm * normalV >> 8); // >>6
    if (temp[0]<boundL || temp[1]<boundL || temp[2]<boundL || temp[0]>boundH || temp[1]>boundH || temp[2]> boundH)
      break;
  }
  lowBound = m - 1;
  lowBoundSurface = std::max(0, ((blockCentroid[dominantAxis] - minPos) + kTrisoupFpHalf >> kTrisoupFpBits));
  highBoundSurface = std::max(0, ((maxPos - blockCentroid[dominantAxis]) + kTrisoupFpHalf >> kTrisoupFpBits));

  return normalV;
}

// --------------------------------------------------------------------------
void
determineCentroidPredictor(
  CentroidInfo& centroidInfo,
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
  int stepQcentro,
  int scaleQ,
  int8_t colocatedCentroidQP,
  int qpNode)
{
  int driftQPred = -100;
  int driftQComp = -100;
  // determine quantized drift for predictor
  if (end > start) {
    int driftPred = 0;
    driftQPred = 0;
    driftQComp = 0;
    int counter = 0;
    int maxD = std::max(1, stepQcentro >> 9);

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

    if (counter) {
      driftPred = divApprox(int64_t(driftPred) << 8, counter * scaleQ, 0); // drift is shift by 8 bit // drift is shift by +8 due to scaleQ,
    }

    int half = stepQcentro >> 1;
    int DZ = 682 * half >> 10; //2 * half / 3;

    if (abs(driftPred) >= DZ) {
      driftQPred = (1 + abs(driftPred - DZ) << 1) / stepQcentro;
      driftQComp = 1 + abs(driftPred - DZ) / stepQcentro;

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
    if (possibleSKIPRef) {
      driftSKIP = driftRef;
      centroidInfo.QP = colocatedCentroidQP;
    }

    if (possibleSKIPComp && (!possibleSKIPRef || badQualityComp < badQualityRef)) {
      driftSKIP = driftComp;
      qualitySKIP = badQualityComp;
      centroidInfo.QP = qpNode;
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
  Vec3<int32_t> normalV,
  Vec3<int32_t> blockCentroid,
  Vec3<int32_t> nodepos,
  PCCPointSet3& pointCloud,
  int start,
  int end,
  int lowBound,
  int  highBound,
  int stepQcentro,
  int  scaleQ,
  int& drift)
{
  // determine quantized drift
  int counter = 0;
  int maxD = std::max(1, stepQcentro >> 9);

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

  if (counter) {
    drift = divApprox(int64_t(drift) << 8, counter * scaleQ, 0);// drift is shift by 8 bit // drift is shift by +8 due to scaleQ
  }

  int half = stepQcentro >> 1;
  int DZ = 682 * half >> 10; //2 * half / 3;

  int driftQ = 0;
  if (abs(drift) >= DZ) {
    driftQ = 1 + (abs(drift) ) / stepQcentro;
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
  int highBound,
  bool isSKIP)
{
  int driftQPred = centroidInfo.driftQPred;
  bool possibleSKIP = centroidInfo.possibleSKIP;
  int driftSKIP = centroidInfo.driftSKIP;

  if (possibleSKIP) {
    arithmeticEncoder->encode(
      isSKIP,
      ctxtMemOctree.ctxDriftSKIP[centroidInfo.qualitySKIP][driftSKIP == 0]);
    if (isSKIP)
      return;
  }

  if (driftQPred == -100) //intra
     arithmeticEncoder->encode(
      driftQ == 0, ctxtMemOctree.ctxDrift0[ctxMinMax][0]);
  else { //inter
    if (possibleSKIP)
      arithmeticEncoder->encode(
        driftQ == 0,
        ctxtMemOctree.ctxDrift0Skip[std::min(3, std::abs(driftSKIP))]);
    else
      arithmeticEncoder->encode(
        driftQ == 0,
        ctxtMemOctree.ctxDrift0[ctxMinMax][1 + std::min(3, std::abs(driftQPred))]);
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
  int highBound,
  bool& isSKIP)
{
  int driftQPred = centroidInfo.driftQPred;
  bool possibleSKIP = centroidInfo.possibleSKIP;
  int driftSKIP = centroidInfo.driftSKIP;

  isSKIP = false;

  if (possibleSKIP) {
    isSKIP = arithmeticDecoder->decode(
      ctxtMemOctree.ctxDriftSKIP[centroidInfo.qualitySKIP][driftSKIP == 0]);
    if (isSKIP)
      return driftSKIP;
  }

  int driftQ = 1;
  if (driftQPred == -100) //intra
    driftQ = arithmeticDecoder->decode(
      ctxtMemOctree.ctxDrift0[ctxMinMax][0]) ? 0 : 1;
  else {//inter
    if (possibleSKIP)
      driftQ = arithmeticDecoder->decode(
        ctxtMemOctree.ctxDrift0Skip[std::min(3, std::abs(driftSKIP))]) ? 0 : 1;
    else
      driftQ = arithmeticDecoder->decode(
        ctxtMemOctree.ctxDrift0[ctxMinMax][1 + std::min(3, std::abs(driftQPred))]) ? 0 : 1;
  }
  if (driftQ == 0)
    return 0;

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
  std::vector<int8_t>& TriSoupVertices2bits,
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
      int vertexPos2bits = TriSoupVertices2bits[patternIdx[v18]];
      if (vertexPos2bits >= 0) {
        ctxInfo.pattern |= 1 << v;
        if (towardOrAway[v18])
          vertexPos2bits = /*max2bits*/ 3 - vertexPos2bits; // reverses for away
        if (vertexPos2bits >= /*mid2bits*/ 2)
          ctxInfo.patternClose |= 1 << v;
        if (vertexPos2bits >= /*max2bits*/ 3)
          ctxInfo.patternClosest |= 1 << v;
        ctxInfo.nclosestPattern += vertexPos2bits >=/* max2bits*/ 3 && v <= 4;
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
  int8_t colocatedVertex,
  int pos2Pred,
  int& mapIdx) {


  ctxMap1 = std::min(ctxInfo.nclosestPattern, 2) * 15 * 2 + (ctxInfo.neighbEdge - 1) * 2 + ((ctxInfo.ctx1 == 4));    // 2* 15 *3 = 90 -> 7 bits
  ctxMap2 = ctxInfo.neighbEnd << 11;
  ctxMap2 |= (ctxInfo.patternClose & (0b00000110)) << 9 - 1; // perp that do not depend on direction = to start
  ctxMap2 |= ctxInfo.direction << 7;
  ctxMap2 |= (ctxInfo.patternClose & (0b00011000)) << 5 - 3; // perp that  depend on direction = to start or to end
  ctxMap2 |= (ctxInfo.patternClose & (0b00000001)) << 4;  // before
  ctxMap2 |= ctxInfo.orderedPclosePar;

  bool isInterGood = isInter && (ctxInfo.nBadPredRef2 <= 0 || ctxInfo.nBadPredComp <= 3);
  ctxInter = 0;
  mapIdx = 0;
  if (isInterGood) {
    ctxInter = 1 + (TriSoupVerticesPred != 0) * 3;//1 4;
    mapIdx = 1 + (TriSoupVerticesPred != 0);

    bool goodRef = ctxInfo.nBadPredRef2 <= 0;
    if (goodRef) {
      ctxInter += (colocatedVertex != 0 ? 2 : 1);
    }
  }
}

// ---------------------------------------------------------------------------
void
constructCtxPos1(
  int& ctxMap1,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int8_t colocatedVertex,
  int pos2Pred,
  int& mapIdx) {

  int ctxFullNbounds = (4 * (ctxInfo.ctx0 <= 1 ? 0 : (ctxInfo.ctx0 >= 3 ? 2 : 1)) + (std::max(1, ctxInfo.ctx1) - 1)) * 2 + (ctxInfo.ctxE == 3);
  ctxMap1 = ctxFullNbounds * 2 + (ctxInfo.nclosestStart > 0);
  ctxMap2 = ctxInfo.missedCloseStart << 8;
  ctxMap2 |= (ctxInfo.patternClosest & 1) << 7;
  ctxMap2 |= ctxInfo.direction << 5;
  ctxMap2 |= ctxInfo.patternClose & (0b00011111);
  ctxMap2 = (ctxMap2 << 4) + ctxInfo.orderedPcloseParPos;

  bool isGoodRef = isInter && colocatedVertex != 0;
  bool isInterGood = isInter && ((isGoodRef && ctxInfo.nBadPredRef2 <= 0) || ctxInfo.nBadPredComp1 <= 4);

  ctxInter = 0;
  mapIdx = 0;
  if (isInterGood) {
    ctxInter =
      TriSoupVerticesPred != 0 ? 1 + (TriSoupVerticesPred > 0 ? 1 : 0) : 0;
    if (ctxInter > 0) {
      ctxInter += 2 * (pos2Pred == 1 || pos2Pred == 2);
      mapIdx = ctxInter;
      ctxInter = ctxInter * 3 - 2;

      int goodPresence = colocatedVertex != 0 && ctxInfo.nBadPredRef2 <= 0;
      if (goodPresence)
        ctxInter += (colocatedVertex > 0 ? 2 : 1);
    }
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
  int Nshift4Mag,
  int v,
  int8_t colocatedVertex,
  int blockWidthLog2,
  int& mapIdx){

  int ctxFullNbounds = (4 * (ctxInfo.ctx0 <= 1 ? 0 : (ctxInfo.ctx0 >= 3 ? 2 : 1)) + (std::max(1, ctxInfo.ctx1) - 1)) * 2 + (ctxInfo.ctxE == 3);
  ctxMap1 = ctxFullNbounds * 2 + (ctxInfo.nclosestStart > 0);
  ctxMap2 = ctxInfo.missedCloseStart << 8;
  ctxMap2 |= (ctxInfo.patternClose & 1) << 7;
  ctxMap2 |= (ctxInfo.patternClosest & 1) << 6;
  ctxMap2 |= ctxInfo.direction << 4;
  ctxMap2 |= (ctxInfo.patternClose & (0b00011111)) >> 1;
  ctxMap2 = (ctxMap2 << 4) + ctxInfo.orderedPcloseParPos;

  ctxInter = 0;
  mapIdx = 0;
  if (isInter) {
    ctxInter = TriSoupVerticesPred != 0
      ? 1 + (!!((std::abs(TriSoupVerticesPred) - 1) >> Nshift4Mag)) : 0;
    mapIdx = ctxInter;
    if (ctxInter > 0) {
      ctxInter = ctxInter * 3 - 2;//1,4
      int goodPresence = colocatedVertex != 0 ? 1 : 0;
      goodPresence = goodPresence
        && ((colocatedVertex > 0 ? 1 : 0) == v)
        && ctxInfo.nBadPredRef2 <= blockWidthLog2 - 2;
      if (goodPresence)
        ctxInter += (!!((std::abs(colocatedVertex) - 1) >> Nshift4Mag)) ? 2 : 1;
    }
  }
}

void
constructCtxPos3(
  int& ctxMap1,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int Nshift4Mag,
  int v,
  int8_t colocatedVertex,
  int blockWidthLog2) {

  int ctxFullNbounds = (4 * (ctxInfo.ctx0 <= 1 ? 0 : (ctxInfo.ctx0 >= 3 ? 2 : 1)) + (std::max(1, ctxInfo.ctx1) - 1)) * 2 + (ctxInfo.ctxE == 3);
  ctxMap1 = ctxFullNbounds * 2 + (ctxInfo.nclosestStart > 0);
  ctxMap2 = ctxInfo.missedCloseStart << 8;
  ctxMap2 |= (ctxInfo.patternClose & 1) << 7;
  ctxMap2 |= (ctxInfo.patternClosest & 1) << 6;
  ctxMap2 |= ctxInfo.direction << 4;
  ctxMap2 |= (ctxInfo.patternClose & (0b00011111)) >> 1;
  ctxMap2 = (ctxMap2 << 2) + (ctxInfo.orderedPcloseParPos >> 2);

  ctxInter = 0;
  if (isInter) {
    int temp = ((std::abs(TriSoupVerticesPred) - 1) >> std::max(0,Nshift4Mag - 1)) & 1;
    ctxInter = TriSoupVerticesPred != 0 ? 1 + temp : 0;

    int goodPresence = colocatedVertex != 0 ? 1 : 0;
    goodPresence = goodPresence && ((colocatedVertex > 0 ? 1 : 0) == (v >> 1));
    goodPresence = goodPresence && (!!((std::abs(colocatedVertex) - 1) >> Nshift4Mag)) == (v & 1) && ctxInfo.nBadPredRef2 <= blockWidthLog2 - 2;
    ctxMap2 |= goodPresence << 14;
    if (goodPresence) {
      int temp2 = ((std::abs(colocatedVertex) - 1) >> std::max(0, Nshift4Mag - 1)) & 1;
      ctxMap2 |= temp2 << 13;
    }
  }
  else {
    ctxMap2 <<= 2;
  }
}



void
constructCtxPos4(
  int& ctxMap1,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int Nshift4Mag,
  int v,
  int8_t colocatedVertex,
  int blockWidthLog2) {

  int ctxFullNbounds = (4 * (ctxInfo.ctx0 <= 1 ? 0 : (ctxInfo.ctx0 >= 3 ? 2 : 1)) + (std::max(1, ctxInfo.ctx1) - 1)) * 2 + (ctxInfo.ctxE == 3);
  ctxMap1 = ctxFullNbounds * 2 + (ctxInfo.nclosestStart > 0);
  ctxMap2 = ctxInfo.missedCloseStart << 8;
  ctxMap2 |= (ctxInfo.patternClose & 1) << 7;
  ctxMap2 |= (ctxInfo.patternClosest & 1) << 6;
  ctxMap2 |= ctxInfo.direction << 4;
  ctxMap2 |= (ctxInfo.patternClose & (0b00011111)) >> 1;

  ctxInter = 0;
  if (isInter) {
    int temp = ((std::abs(TriSoupVerticesPred) - 1) >> std::max(0, Nshift4Mag - 2)) & 1;
    ctxInter = TriSoupVerticesPred != 0 ? 1 + temp : 0;

    int goodPresence = colocatedVertex != 0 ? 1 : 0;
    goodPresence = goodPresence && ((colocatedVertex > 0 ? 1 : 0) == (v >> 2));
    goodPresence = goodPresence && (!!((std::abs(colocatedVertex) - 1) >> Nshift4Mag)) == ((v >> 1) & 1);
    goodPresence = goodPresence && (((std::abs(colocatedVertex) - 1) >> std::max(0, Nshift4Mag - 1)) & 1) == (v & 1) && ctxInfo.nBadPredRef2 <= blockWidthLog2 - 2;
    ctxMap2 |= goodPresence << 12;
    if (goodPresence) {
      int temp2 = ((std::abs(colocatedVertex) - 1) >> std::max(0, Nshift4Mag - 2)) & 1;
      ctxMap2 |= temp2 << 11;
    }
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

  if (triCount >= 3) {
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
