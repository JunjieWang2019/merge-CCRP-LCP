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

#include "geometry.h"

#include "OctreeNeighMap.h"
#include "geometry_octree.h"
//#include "geometry_intra_pred.h"
#include "io_hls.h"
#include "tables.h"
#include "quantization.h"
#include "motionWip.h"
#include <unordered_map>

namespace pcc {

//============================================================================

class GeometryOctreeDecoder : protected GeometryOctreeContexts {
public:
  GeometryOctreeDecoder(
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    const GeometryOctreeContexts& ctxMem,
    EntropyDecoder* arithmeticDecoder);

  GeometryOctreeDecoder(const GeometryOctreeDecoder&) = default;
  GeometryOctreeDecoder(GeometryOctreeDecoder&&) = default;
  GeometryOctreeDecoder& operator=(const GeometryOctreeDecoder&) = default;
  GeometryOctreeDecoder& operator=(GeometryOctreeDecoder&&) = default;

  void beginOctreeLevel(const Vec3<int>& planarDepth);

  // dynamic OBUF
  void resetMap() { GeometryOctreeContexts::resetMap(); }
  void clearMap() { GeometryOctreeContexts::clearMap(); };

  int decodePositionLeafNumPoints();

  /*int decodePlanarMode(
    OctreeNodePlanar& planar,
    int planeZ,
    int dist,
    int adjPlanes,
    int planeId,
    bool* multiPlanarFlag,
    bool* multiPlanarEligible, 
    OctreeNodePlanar& planarRef

  );*/

  /*void derivePlanarPCMContextBuffer(
    OctreeNodePlanar& planar,
    OctreeNodePlanar& planarRef,
    OctreePlanarBuffer& planeBuffer,
    int xx,
    int yy,
    int zz);*/

  /*void determinePlanarMode(
    bool adjacent_child_contextualization_enabled_flag,
    int planeId,
    OctreeNodePlanar& child,
    OctreePlanarBuffer::Row* planeBuffer,
    int coord1,
    int coord2,
    int coord3,
    int posInParent,
    const GeometryNeighPattern& gnp,
    uint8_t siblingOccupancy,
    int planarRate[3],
    bool* multiPlanarFlag,
    bool* multiPlanarEligible,
    OctreeNodePlanar& planarRef

  );*/

  /*void determinePlanarMode(
    bool adjacent_child_contextualization_enabled_flag,
    const bool planarEligible[3],
    int posInParent,
    const GeometryNeighPattern& gnp,
    PCCOctree3Node& child,
    OctreeNodePlanar& planar,
    OctreeNodePlanar& planarRef
  );*/

  uint32_t decodeOccupancyFullNeihbourgsNZ(
    const GeometryNeighPattern& gnp,
    int planarMaskX,
    int planarMaskY,
    int planarMaskZ,
    bool planarPossibleX,
    bool planarPossibleY,
    bool planarPossibleZ,
    const MortonMap3D& occupancyAtlas,
    Vec3<int32_t> &pos,
    const int atlasShift,
    int predOcc,
    bool isInter,
    bool flagNoSingle);

  uint32_t decodeOccupancyFullNeihbourgs(
    const GeometryNeighPattern& gnp,
    int planarMaskX,
    int planarMaskY,
    int planarMaskZ,
    bool planarPossibleX,
    bool planarPossibleY,
    bool planarPossibleZ,
    const MortonMap3D& occupancyAtlas,
    Vec3<int32_t> pos,
    const int atlasShift,
    bool flagWord4,
    bool adjacent_child_contextualization_enabled_flag
    , int predOcc, bool isInter
  );

  Vec3<int32_t> decodePointPosition(
    const Vec3<int>& nodeSizeLog2, Vec3<int32_t>& deltaPlanar);

  void decodeOrdered2ptPrefix(
    Vec3<bool> directIdcm,
    Vec3<int>& nodeSizeLog2AfterUnordered,
    Vec3<int32_t> deltaUnordered[2]);

  bool decodeNodeQpOffsetsPresent();
  int decodeQpOffset();

  bool decodeIsIdcm();

  template<class OutputIt>
  int decodeDirectPosition(
    bool geom_unique_points_flag,
    bool joint_2pt_idcm_enabled_flag,
    const Vec3<int>& nodeSizeLog2,
    const Vec3<int>& posQuantBitMasks,
    const PCCOctree3Node& node,
    const OctreeNodePlanar& planar,
    OutputIt outputPoints);

  const GeometryOctreeContexts& getCtx() const { return *this; }

public:
  const uint8_t* _neighPattern64toR1;

  EntropyDecoder* _arithmeticDecoder;

  // Planar state
  //OctreePlanarState _planar;
};

//============================================================================

GeometryOctreeDecoder::GeometryOctreeDecoder(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const GeometryOctreeContexts& ctxtMem,
  EntropyDecoder* arithmeticDecoder)
  : GeometryOctreeContexts(ctxtMem)
  , _neighPattern64toR1(neighPattern64toR1(gps))
  , _arithmeticDecoder(arithmeticDecoder)
  /*, _planar(gps)*/
{
}
//============================================================================

void
GeometryOctreeDecoder::beginOctreeLevel(const Vec3<int>& planarDepth)
{
  //_planar.initPlanes(planarDepth);
}


//============================================================================
// Decode the number of points in a leaf node of the octree.

int
GeometryOctreeDecoder::decodePositionLeafNumPoints()
{
  int val = _arithmeticDecoder->decode(_ctxDupPointCntGt0);
  if (val)
    val += _arithmeticDecoder->decodeExpGolomb(0, _ctxDupPointCntEgl);

  return val + 1;
}

//============================================================================

//int
//GeometryOctreeDecoder::decodePlanarMode(
//  OctreeNodePlanar& planar,
//  int planeZ,
//  int dist,
//  int adjPlanes,
//  int planeId,
//  bool* multiPlanarFlag,
//  bool* multiPlanarEligible,
//  OctreeNodePlanar& planarRef)
//{
//  const int mask0 = (1 << planeId);
//  const int mask1[3] = {6, 5, 3};
//
//  // decode planar mode
//  bool isPlanarRef = planarRef.planarMode & mask0;
//  int planeBitRef = (planarRef.planePosBits & mask0) == 0 ? 0 : 1;
//
//  int ctxIdx_Planar_flag = planeId;
//  if (isPlanarRef)
//    ctxIdx_Planar_flag += 3 * (planeBitRef + 1);
//
//  bool isPlanar = isPlanarRef;
//  
//
//  if (!planar.isPCM) {
//    if (_planar._geom_multiple_planar_mode_enable_flag) {
//      bool multiPlanarFlagFalse = true;
//      static const int planeId2Index[3][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}};
//      for (int i = 0; i < 3; i++) {
//        multiPlanarFlagFalse &= !(multiPlanarFlag[planeId2Index[planeId][i]]);
//      }
//      bool inferredPlanarFalse =
//        multiPlanarFlagFalse;  // todo: consider renaming inferredPlaneFalse
//
//      if (multiPlanarFlagFalse) {
//        if (planeId == 2) {
//          if (multiPlanarEligible[0])  //xyz
//            inferredPlanarFalse =
//              !((planar.planarMode & 2) && (planar.planarMode & 1));
//          else if (multiPlanarEligible[2])  //xz
//            inferredPlanarFalse = !(planar.planarMode & 1);
//          else if (multiPlanarEligible[3])  //yz
//            inferredPlanarFalse = !(planar.planarMode & 2);
//
//        } else if (planeId == 1) {
//          if (multiPlanarEligible[1])  //xy
//            inferredPlanarFalse = !(planar.planarMode & 1);
//        }
//      }
//
//      if (inferredPlanarFalse)
//        isPlanar =
//          _arithmeticDecoder->decode(_ctxPlanarMode[ctxIdx_Planar_flag]);
//      else if (!multiPlanarFlagFalse)
//        isPlanar = true;
//      else
//        isPlanar = false;
//    } else {
//      isPlanar =
//        _arithmeticDecoder->decode(_ctxPlanarMode[ctxIdx_Planar_flag]);
//    }
//  }
//
//  planar.planarMode |= isPlanar ? mask0 : 0;
//
//  if (!isPlanar) {
//    planar.planarPossible &= mask1[planeId];
//    return -1;
//  }
//
//  // decode the plane index // encode the plane index
//  int planeBit;
//
//  if (planar.isPCM) {
//    planeBit = planeBitRef;
//    planar.planePosBits |= (planeBit << planeId);
//    return planeBit;
//  }
//  // Not PCM and signal the plane position bit
//  if (
//    planeId == planar.lastDirIdx && planar.isPreDirMatch && planar.allowPCM
//    && isPlanarRef) {
//    planeBit = (planeBitRef == 1) ? 0 : 1;
//    planar.planePosBits |= (planeBit << planeId);
//    return planeBit;
//  }
//
//  {
//    static const int kAdjPlaneCtx[4] = {0, 1, 2, 0};
//    int planePosCtx = kAdjPlaneCtx[adjPlanes];
//    if (planeZ < 0) {
//      int planePostCtxTmp = planePosCtx;
//      if (isPlanarRef) {
//        planePostCtxTmp += 3 * (planeBitRef + 1);
//      }
//      planeBit =
//        _arithmeticDecoder->decode(_ctxPlanarPlaneLastIndexZ[planePostCtxTmp]);
//
//    } else {
//      int discreteDist = dist > (8 >> OctreePlanarBuffer::shiftAb);
//      int lastIndexPlane2d = planeZ + (discreteDist << 1);
//      int refPlane = 0;
//      if (isPlanarRef) {
//        refPlane = 1 + planeBitRef;
//      }
//      planeBit = _arithmeticDecoder->decode(
//        _ctxPlanarPlaneLastIndex[refPlane][planeId][planePosCtx]
//                                [lastIndexPlane2d]);
//    }
//  }
//
//  planar.planePosBits |= (planeBit << planeId);
//  return planeBit;
//}

//============================================================================
//void
//GeometryOctreeDecoder::derivePlanarPCMContextBuffer(
//  OctreeNodePlanar& planar,
//  OctreeNodePlanar& planarRef,
//  OctreePlanarBuffer& planeBuffer,
//  int xx,
//  int yy,
//  int zz)
//{
//  int matchedDir = 0;
//
//  planarRef.ctxBufPCM = 4
//    * (int(planar.eligible[0]) + int(planar.eligible[1])
//       + int(planar.eligible[2]) - 1);
//  assert(planarRef.ctxBufPCM >= 0);
//
//  for (int planeId = 0; planeId < 3; planeId++) {
//    if (planar.eligible[planeId]) {
//      const int mask0 = (1 << planeId);
//
//      bool isPlanarRef = planarRef.planarMode & mask0;
//      int planeBitRef = (planarRef.planePosBits & mask0) == 0 ? 0 : 1;
//
//      // Get the buffer information
//      OctreePlanarBuffer::Row* planeBufferDir = planeBuffer.getBuffer(planeId);
//
//      int coord3 = (planeId == 2) ? zz : (planeId == 1 ? yy : xx);
//
//      const int rowLen = OctreePlanarBuffer::rowSize;
//
//      if (planeBufferDir) {
//        coord3 &= OctreePlanarBuffer::maskC;
//
//        const auto row = planeBufferDir[coord3];
//
//        const int idxMinDist = rowLen - 1;
//        const int closestPlanarFlag = row[idxMinDist].planeIdx;
//        const bool closestPL = (closestPlanarFlag > -1) ? true : false;
//        const int closestPlane = closestPL ? closestPlanarFlag : 0;
//
//        matchedDir +=
//          int(closestPL == isPlanarRef && closestPlane == planeBitRef);
//
//      }
//    }
//  }
//  planarRef.ctxBufPCM += matchedDir;
//}


//============================================================================

//void
//GeometryOctreeDecoder::determinePlanarMode(
//  bool adjacent_child_contextualization_enabled_flag,
//  int planeId,
//  OctreeNodePlanar& planar,
//  OctreePlanarBuffer::Row* planeBuffer,
//  int coord1,
//  int coord2,
//  int coord3,
//  int posInParent,
//  const GeometryNeighPattern& gnp,
//  uint8_t siblingOccupancy,
//  int planarRate[3],
//  bool* multiPlanarFlag,
//  bool* multiPlanarEligible,
//  OctreeNodePlanar& planarRef)
//{
//  const int kPlanarChildThreshold = 63;
//  const int kAdjNeighIdxFromPlanePos[3][2] = {1, 0, 2, 3, 4, 5};
//  const int planeSelector = 1 << planeId;
//  static const uint8_t KAdjNeighIdxMask[3][2] = {0x0f, 0xf0, 0x33,
//                                                 0xcc, 0x55, 0xaa};
//  OctreePlanarBuffer::Elmt* row;
//  int rowLen = OctreePlanarBuffer::rowSize;
//  int closestPlanarFlag;
//  int closestDist;
//  int maxCoord;
//
//  if (!planeBuffer) {
//    // angular: buffer disabled
//    closestPlanarFlag = -1;
//    closestDist = 0;
//  } else {
//    coord1 =
//      (coord1 & OctreePlanarBuffer::maskAb) >> OctreePlanarBuffer::shiftAb;
//    coord2 =
//      (coord2 & OctreePlanarBuffer::maskAb) >> OctreePlanarBuffer::shiftAb;
//    coord3 = coord3 & OctreePlanarBuffer::maskC;
//
//    row = planeBuffer[coord3];
//
//    maxCoord = std::max(coord1, coord2);
//    closestDist = std::abs(maxCoord - int(row[rowLen - 1].pos));
//    int idxMinDist = rowLen - 1;
//
//    // push closest point front
//    row[rowLen - 1] = row[idxMinDist];
//
//    closestPlanarFlag = row[idxMinDist].planeIdx;
//  }
//
//  // The relative plane position (0|1) along the planeId axis.
//  int pos = !(KAdjNeighIdxMask[planeId][0] & (1 << posInParent));
//
//  // Determine which adjacent planes are occupied
//  // The low plane is at position axis - 1
//  bool lowAdjPlaneOccupied = adjacent_child_contextualization_enabled_flag
//    ? KAdjNeighIdxMask[planeId][1] & gnp.adjNeighOcc[planeId]
//    : (gnp.neighPattern >> kAdjNeighIdxFromPlanePos[planeId][0]) & 1;
//
//  // The high adjacent plane is at position axis + 1
//  bool highAdjPlaneOccupied = !pos
//    ? KAdjNeighIdxMask[planeId][1] & siblingOccupancy
//    : (gnp.neighPattern >> kAdjNeighIdxFromPlanePos[planeId][1]) & 1;
//
//  int adjPlanes = (highAdjPlaneOccupied << 1) | lowAdjPlaneOccupied;
//
//  int planeBit = decodePlanarMode(
//    planar, closestPlanarFlag, closestDist, adjPlanes, planeId,
//    multiPlanarFlag, multiPlanarEligible, planarRef);
//  bool isPlanar = (planar.planarMode & planeSelector);
//
//  planarRate[planeId] =
//    (255 * planarRate[planeId] + (isPlanar ? 256 * 8 : 0) + 128) >> 8;
//
//  if (planeBuffer)
//    row[rowLen - 1] = {unsigned(maxCoord), planeBit};
//
//  bool isPlanarRef = (planarRef.planarMode & planeSelector);
//  int planeBitRef = (planarRef.planePosBits & planeSelector) == 0 ? 0 : 1;
//
//  if (!((isPlanar == isPlanarRef) && (planeBit == planeBitRef))) {
//    planar.isPreDirMatch = false;
//  }
//}

//============================================================================

//void
//GeometryOctreeDecoder::determinePlanarMode(
//  bool adjacent_child_contextualization_enabled_flag,
//  const bool planarEligible[3],
//  int posInParent,
//  const GeometryNeighPattern& gnp,
//  PCCOctree3Node& child,
//  OctreeNodePlanar& planar,
//  OctreeNodePlanar& planarRef
//)
//{
//  int xx = child.pos[0];
//  int yy = child.pos[1];
//  int zz = child.pos[2];
//
//  auto& planeBuffer = _planar._planarBuffer;
//
//  uint8_t planarEligibleMask = 0;
//  planarEligibleMask |= planarEligible[2] << 2;
//  planarEligibleMask |= planarEligible[1] << 1;
//  planarEligibleMask |= planarEligible[0] << 0;
//  planarRef.planarMode &= planarEligibleMask;
//  planarRef.planePosBits &= planarEligibleMask;
//
//  if (planar.allowPCM) {
//    derivePlanarPCMContextBuffer(planar, planarRef, planeBuffer, xx, yy, zz);
//  }
//
//  if (!planar.isRead && planar.allowPCM) {
//    planar.isPCM = _arithmeticDecoder->decode(
//      _ctxPlanarCopyMode[planarRef.ctxBufPCM][planarRef.planarMode]);
//    planar.isRead = true;
//  }
//  bool multiPlanarFlag[4] = {false, false, false, false};
//  bool multiPlanarEligible[4] = {false, false, false, false};
//  if (_planar._geom_multiple_planar_mode_enable_flag) {
//    if (!planar.isPCM) {
//      if (planarEligible[2] && planarEligible[1] && planarEligible[0]) {
//        multiPlanarEligible[0] = true;
//        multiPlanarFlag[0] = _arithmeticDecoder->decode(_ctxMultiPlanarMode);
//      } else if (
//        (!planarEligible[2]) && planarEligible[1] && planarEligible[0]) {  //xy
//        multiPlanarEligible[1] = true;
//        multiPlanarFlag[1] = _arithmeticDecoder->decode(_ctxMultiPlanarMode);
//      } else if (
//        planarEligible[2] && (!planarEligible[1]) && planarEligible[0]) {  //xz
//        multiPlanarEligible[2] = true;
//        multiPlanarFlag[2] = _arithmeticDecoder->decode(_ctxMultiPlanarMode);
//      } else if (
//        planarEligible[2] && planarEligible[1] && (!planarEligible[0])) {  //yz
//        multiPlanarEligible[3] = true;
//        multiPlanarFlag[3] = _arithmeticDecoder->decode(_ctxMultiPlanarMode);
//      }
//    }
//  }
//
//
//
//  // planar x
//  if (planarEligible[0]) {
//    determinePlanarMode(
//      adjacent_child_contextualization_enabled_flag, 0, planar,
//      planeBuffer.getBuffer(0), yy, zz, xx, posInParent, gnp,
//      child.siblingOccupancy, _planar._rate.data(),
//      multiPlanarFlag, multiPlanarEligible, planarRef);
//  }
//  // planar y
//  if (planarEligible[1]) {
//    determinePlanarMode(
//      adjacent_child_contextualization_enabled_flag, 1, planar,
//      planeBuffer.getBuffer(1), xx, zz, yy, posInParent, gnp,
//      child.siblingOccupancy, _planar._rate.data(),
//      multiPlanarFlag, multiPlanarEligible, planarRef
//
//    );
//  }
//  // planar z
//  if (planarEligible[2]) {
//    determinePlanarMode(
//      adjacent_child_contextualization_enabled_flag, 2, planar,
//      planeBuffer.getBuffer(2), xx, yy, zz, posInParent, gnp,
//      child.siblingOccupancy, _planar._rate.data(),
//      multiPlanarFlag, multiPlanarEligible, planarRef);
//  }
//}

//-------------------------------------------------------------------------
// decode node occupancy bits
//

static void (*pointer2FunctionContext[8])(
  OctreeNeighours&, int, int&, int&, bool&) = {
  makeGeometryAdvancedNeighPattern0,
  makeGeometryAdvancedNeighPattern1,
  makeGeometryAdvancedNeighPattern2,
  makeGeometryAdvancedNeighPattern3,
  makeGeometryAdvancedNeighPattern4,
  makeGeometryAdvancedNeighPattern5,
  makeGeometryAdvancedNeighPattern6,
  makeGeometryAdvancedNeighPattern7};

static const int LUTinitCoded0[27][6] = {
  {0, 0, 0, 0, 0, 0}, {4, 0, 2, 2, 2, 2}, {0, 4, 2, 2, 2, 2},
  {2, 2, 4, 0, 2, 2}, {4, 2, 4, 2, 3, 3}, {2, 4, 4, 2, 3, 3},
  {2, 2, 0, 4, 2, 2}, {4, 2, 2, 4, 3, 3}, {2, 4, 2, 4, 3, 3},
  {2, 2, 2, 2, 4, 0}, {4, 2, 3, 3, 4, 2}, {2, 4, 3, 3, 4, 2},
  {3, 3, 4, 2, 4, 2}, {4, 3, 4, 3, 4, 3}, {3, 4, 4, 3, 4, 3},
  {3, 3, 2, 4, 4, 2}, {4, 3, 3, 4, 4, 3}, {3, 4, 3, 4, 4, 3},
  {2, 2, 2, 2, 0, 4}, {4, 2, 3, 3, 2, 4}, {2, 4, 3, 3, 2, 4},
  {3, 3, 4, 2, 2, 4}, {4, 3, 4, 3, 3, 4}, {3, 4, 4, 3, 3, 4},
  {3, 3, 2, 4, 2, 4}, {4, 3, 3, 4, 3, 4}, {3, 4, 3, 4, 3, 4}};

uint32_t
GeometryOctreeDecoder::decodeOccupancyFullNeihbourgsNZ(
  const GeometryNeighPattern& gnp,
  int planarMaskX,
  int planarMaskY,
  int planarMaskZ,
  bool planarPossibleX,
  bool planarPossibleY,
  bool planarPossibleZ,
  const MortonMap3D& occupancyAtlas,
  Vec3<int32_t> &pos,
  const int atlasShift,
  int predOcc,
  bool isInter,
  bool flagNoSingle)
{
  bool sure_planarityX = planarMaskX || !planarPossibleX;
  bool sure_planarityY = planarMaskY || !planarPossibleY;
  bool sure_planarityZ = planarMaskZ || !planarPossibleZ;
  const int maxPerPlaneX = planarMaskX && flagNoSingle ? 2 : 3;
  const int maxPerPlaneY = planarMaskY && flagNoSingle ? 2 : 3;
  const int maxPerPlaneZ = planarMaskZ && flagNoSingle ? 2 : 3;
  const int maxAll = flagNoSingle ? 6 : 7;

  //int  MaskConfig = !planarMaskX ? 0 : planarMaskX == 15 ? 1 : 2;
  //MaskConfig += !planarMaskY ? 0 : planarMaskY == 51 ? 3 : 6;
  //MaskConfig += !planarMaskZ ? 0 : planarMaskZ == 85 ? 9 : 18;
  int MaskConfig = (!!planarMaskX) * (1 + (planarMaskX != 0x0F));
  MaskConfig += (!!planarMaskY) * 3 * (1 + (planarMaskY != 0x33));
  MaskConfig += (!!planarMaskZ) * 9 * (1 + (planarMaskZ != 0x55));

  int coded0[6] = {0, 0, 0, 0, 0, 0};  // mask x0 x1 y0 y1 z0 z1
  if (MaskConfig) {
    memcpy(coded0, LUTinitCoded0[MaskConfig], 6 * sizeof(int));
  }

  OctreeNeighours octreeNeighours;
  prepareGeometryAdvancedNeighPattern(
    octreeNeighours, gnp, pos, atlasShift, occupancyAtlas);

  // loop on occupancy bits from occupancy map
  uint32_t partialOccupancy = 0;
  uint32_t occupancy = 0;
  int maskedOccupancy = planarMaskX | planarMaskY | planarMaskZ;
  for (int i = 0; i < 8; i++) {
    if ((maskedOccupancy >> i) & 1) {
      //  bit is 0 because masked by QTBT or planar
      partialOccupancy <<= 1;
      continue;
    }

    int mask0X = (0xf0 >> i) & 1;
    int mask0Y = 2 + ((0xcc >> i) & 1);
    int mask0Z = 4 + ((0xaa >> i) & 1);

    bool bitIsOne =
      (sure_planarityX && coded0[mask0X] >= maxPerPlaneX)
      || (coded0[0] + coded0[1] >= maxAll)
      || (sure_planarityY && coded0[mask0Y] >= maxPerPlaneY)
      || (coded0[2] + coded0[3] >= maxAll)
      || (sure_planarityZ && coded0[mask0Z] >= maxPerPlaneZ)
      || (coded0[4] + coded0[5] >= maxAll);

    if (bitIsOne) {
      occupancy += 1 << i;
      partialOccupancy <<= 1;
      partialOccupancy |= 1;
      continue;
    }

    int bitPred = (predOcc >> i) & 1;
    int interCtx = bitPred;

    // OBUF contexts
    int ctx1, ctx2;
    bool Sparse;
    (*pointer2FunctionContext[i])(
      octreeNeighours, occupancy,  ctx1, ctx2, Sparse);

    bool isInter2 = isInter && predOcc;
    if (isInter2) {
      int bitPred = (predOcc >> i) & 1;
      int bitPred2 = (predOcc >> i + 8) & 1;
      ctx1 = (ctx1 << 2) | bitPred | (bitPred2 << 1);
    }

    // decode

    int bit;
    if (Sparse) {
      bit = _MapOccupancySparse[isInter2][i].decodeEvolve(
        _arithmeticDecoder, _CtxMapDynamicOBUF[isInter2], ctx2, ctx1,
        &_OBUFleafNumber, _BufferOBUFleaves);
    }
    else {
      bit = _MapOccupancy[isInter2][i].decodeEvolve(
        _arithmeticDecoder, _CtxMapDynamicOBUF[2+ isInter2], ctx2, ctx1,
        &_OBUFleafNumber, _BufferOBUFleaves);
    }

    // update partial occupancy of current node
    occupancy += bit << i;
    coded0[mask0X] += !bit;
    coded0[mask0Y] += !bit;
    coded0[mask0Z] += !bit;
    partialOccupancy <<= 1;
    partialOccupancy |= bit;
  }

  return occupancy;
}

//-------------------------------------------------------------------------
// decode node occupancy bits
//
uint32_t
GeometryOctreeDecoder::decodeOccupancyFullNeihbourgs(
  const GeometryNeighPattern& gnp,
  int planarMaskX,
  int planarMaskY,
  int planarMaskZ,
  bool planarPossibleX,
  bool planarPossibleY,
  bool planarPossibleZ,
  const MortonMap3D& occupancyAtlas,
  Vec3<int32_t> pos,
  const int atlasShift,
  bool flagWord4,
  bool adjacent_child_contextualization_enabled_flag,
  int predOcc,
  bool isInter)
{
  // decode occupancy pattern
  uint32_t occupancy;

  // single child and we know its position
  if (planarMaskX && planarMaskY && planarMaskZ) {
    uint32_t cnt = (planarMaskZ & 1);
    cnt |= (planarMaskY & 1) << 1;
    cnt |= (planarMaskX & 1) << 2;
    occupancy = 1 << cnt;
    return occupancy;
  }

  // neighbour empty and only one point => decode index, not pattern
  //------ Z occupancy decoding from here ----------------

  bool flagNoSingle = false;

  if (
    gnp.neighPattern == 0
    && (!predOcc || (planarMaskX | planarMaskY | planarMaskZ))) {
    bool singleChild = false;
    if (planarPossibleX && planarPossibleY && planarPossibleZ) {
      singleChild = _arithmeticDecoder->decode(_ctxSingleChild) == 1;
    }

    if (singleChild) {
      uint32_t cnt;
      if (!planarMaskZ)
        cnt = _arithmeticDecoder->decode();
      else
        cnt = (planarMaskZ & 1);

      if (!planarMaskY)
        cnt |= _arithmeticDecoder->decode() << 1;
      else
        cnt |= (planarMaskY & 1) << 1;

      if (!planarMaskX)
        cnt |= _arithmeticDecoder->decode() << 2;
      else
        cnt |= (planarMaskX & 1) << 2;

      occupancy = 1 << cnt;
      return occupancy;
    }

    flagNoSingle = true;
    // at least two child nodes occupied and two planars => we know the occupancy
    if (planarMaskX && planarMaskY) {
      uint32_t cnt = ((planarMaskX & 1) << 2) | ((planarMaskY & 1) << 1);
      occupancy = (1 << cnt) | (1 << (cnt + 1));
      return occupancy;
    }

    if (planarMaskY && planarMaskZ) {
      uint32_t cnt = ((planarMaskY & 1) << 1) | (planarMaskZ & 1);
      occupancy = (1 << cnt) | (1 << (cnt + 4));
      return occupancy;
    }

    if (planarMaskX && planarMaskZ) {
      uint32_t cnt = ((planarMaskX & 1) << 2) | (planarMaskZ & 1);
      occupancy = (1 << cnt) | (1 << (cnt + 2));
      return occupancy;
    }
  }

  return decodeOccupancyFullNeihbourgsNZ(
    gnp, planarMaskX, planarMaskY, planarMaskZ,
    planarPossibleX, planarPossibleY, planarPossibleZ,
    occupancyAtlas, pos, atlasShift, predOcc, isInter, flagNoSingle);
}
//-------------------------------------------------------------------------

bool
GeometryOctreeDecoder::decodeNodeQpOffsetsPresent()
{
  return _arithmeticDecoder->decode();
}

//-------------------------------------------------------------------------

int
GeometryOctreeDecoder::decodeQpOffset()
{
  if (!_arithmeticDecoder->decode(_ctxQpOffsetAbsGt0))
    return 0;

  int dqp = _arithmeticDecoder->decodeExpGolomb(0, _ctxQpOffsetAbsEgl) + 1;
  int dqp_sign = _arithmeticDecoder->decode(_ctxQpOffsetSign);
  return dqp_sign ? -dqp : dqp;
}

//-------------------------------------------------------------------------
// Decode a position of a point in a given volume.
Vec3<int32_t>
GeometryOctreeDecoder::decodePointPosition(
  const Vec3<int>& nodeSizeLog2, Vec3<int32_t>& deltaPlanar)
{
  Vec3<int32_t> delta = deltaPlanar;
  for (int k = 0; k < 3; k++) {
    if (nodeSizeLog2[k] <= 0)
      continue;

    for (int i = nodeSizeLog2[k]; i > 0; i--) {
      delta[k] <<= 1;
      delta[k] |= _arithmeticDecoder->decode();
    }
  }

  return delta;
}

//-------------------------------------------------------------------------
// Decode part of the position of two unordred points  point in a given volume.
void
GeometryOctreeDecoder::decodeOrdered2ptPrefix(
  Vec3<bool> directIdcm, Vec3<int>& nodeSizeLog2, Vec3<int32_t> pointPrefix[2])
{
  if (nodeSizeLog2[0] >= 1 && directIdcm[0]) {
    int ctxIdx = 0;
    bool sameBit = true;
    while (nodeSizeLog2[0] && sameBit) {
      pointPrefix[0][0] <<= 1;
      pointPrefix[1][0] <<= 1;
      nodeSizeLog2[0]--;

      sameBit = _arithmeticDecoder->decode(_ctxSameBitHighx[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      if (sameBit) {
        int bit = _arithmeticDecoder->decode();
        pointPrefix[0][0] |= bit;
        pointPrefix[1][0] |= bit;
      } else {
        pointPrefix[1][0] |= 1;
      }
    }
  }

  if (nodeSizeLog2[1] >= 1 && directIdcm[1]) {
    int ctxIdx = 0;
    bool sameBit = true;
    bool sameX = !directIdcm[0] || pointPrefix[0][0] == pointPrefix[1][0];

    while (nodeSizeLog2[1] && sameBit) {
      pointPrefix[0][1] <<= 1;
      pointPrefix[1][1] <<= 1;
      nodeSizeLog2[1]--;

      sameBit = _arithmeticDecoder->decode(_ctxSameBitHighy[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      int bit = 0;
      if (!(sameX && !sameBit))
        bit = _arithmeticDecoder->decode();
      pointPrefix[0][1] |= bit;
      pointPrefix[1][1] |= sameBit ? bit : !bit;
    }
  }

  if (nodeSizeLog2[2] >= 1 && directIdcm[2]) {
    int ctxIdx = 0;
    bool sameBit = true;
    bool sameXy = (!directIdcm[0] || pointPrefix[0][0] == pointPrefix[1][0])
      && (!directIdcm[1] || pointPrefix[0][1] == pointPrefix[1][1]);

    while (nodeSizeLog2[2] && sameBit) {
      pointPrefix[0][2] <<= 1;
      pointPrefix[1][2] <<= 1;
      nodeSizeLog2[2]--;

      sameBit = _arithmeticDecoder->decode(_ctxSameBitHighz[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      int bit = 0;
      if (!(sameXy && !sameBit))
        bit = _arithmeticDecoder->decode();
      pointPrefix[0][2] |= bit;
      pointPrefix[1][2] |= sameBit ? bit : !bit;
    }
  }
}

//-------------------------------------------------------------------------

bool
GeometryOctreeDecoder::decodeIsIdcm()
{
  return _arithmeticDecoder->decode(_ctxBlockSkipTh);
}

//-------------------------------------------------------------------------
// Direct coding of position of points in node (early tree termination).
// Decoded points are written to @outputPoints
// Returns the number of points emitted.

template<class OutputIt>
int
GeometryOctreeDecoder::decodeDirectPosition(
  bool geom_unique_points_flag,
  bool joint_2pt_idcm_enabled_flag,
  const Vec3<int>& nodeSizeLog2,
  const Vec3<int>& posQuantBitMask,
  const PCCOctree3Node& node,
  const OctreeNodePlanar& planar,
  OutputIt outputPoints)
{
  int numPoints = 1;
  bool numPointsGt1 = _arithmeticDecoder->decode(_ctxNumIdcmPointsGt1);
  numPoints += numPointsGt1;

  int numDuplicatePoints = 0;
  if (!geom_unique_points_flag && !numPointsGt1) {
    numDuplicatePoints = _arithmeticDecoder->decode(_ctxDupPointCntGt0);
    if (numDuplicatePoints) {
      numDuplicatePoints += _arithmeticDecoder->decode(_ctxDupPointCntGt1);
      if (numDuplicatePoints == 2)
        numDuplicatePoints +=
          _arithmeticDecoder->decodeExpGolomb(0, _ctxDupPointCntEgl);
    }
  }

  // nodeSizeLog2Rem indicates the number of bits left to decode
  // the first bit may be inferred from the planar information
  Vec3<int32_t> deltaPlanar{0, 0, 0};
  Vec3<int> nodeSizeLog2Rem = nodeSizeLog2;
  for (int k = 0; k < 3; k++)
    if (nodeSizeLog2Rem[k] > 0 && (planar.planarMode & (1 << k))) {
      deltaPlanar[k] |= (planar.planePosBits & (1 << k) ? 1 : 0);
      nodeSizeLog2Rem[k]--;
    }

  // quantised partial positions must be scaled for angular coding
  // nb, the decoded position remains quantised.
  OctreeAngPosScaler quant(node.qp, posQuantBitMask);

  // Indicates which components are directly coded
  Vec3<bool> directIdcm = true;

  // decode (ordred) two points
  Vec3<int32_t> deltaPos[2] = {deltaPlanar, deltaPlanar};
  if (numPoints == 2 && joint_2pt_idcm_enabled_flag)
    decodeOrdered2ptPrefix(directIdcm, nodeSizeLog2Rem, deltaPos);

  Vec3<int32_t> pos;
  for (int i = 0; i < numPoints; i++) {
    *(outputPoints++) = pos =
      decodePointPosition(nodeSizeLog2Rem, deltaPos[i]);
  }

  for (int i = 0; i < numDuplicatePoints; i++)
    *(outputPoints++) = pos;

  return numPoints + numDuplicatePoints;
}

//-------------------------------------------------------------------------
// Helper to inverse quantise positions

Vec3<int32_t>
invQuantPosition(int qp, Vec3<uint32_t> quantMasks, const Vec3<int32_t>& pos)
{
  // pos represents the position within the coded tree as follows:
  //     |pppppqqqqqq|00
  //  - p = unquantised bit
  //  - q = quantised bit
  //  - 0 = bits that were not coded (MSBs of q)
  // The reconstruction is:
  //   |ppppp00qqqqqq| <- just prior to scaling
  //   |pppppssssssss| <  after scaling (s = scale(q))

  QuantizerGeom quantizer(qp);
  int shiftBits = QuantizerGeom::qpShift(qp);
  Vec3<int32_t> recon;
  for (int k = 0; k < 3; k++) {
    int lowPart = pos[k] & (quantMasks[k] >> shiftBits);
    int highPart = pos[k] ^ lowPart;
    int lowPartScaled = PCCClip(quantizer.scale(lowPart), 0, quantMasks[k]);
    recon[k] = (highPart << shiftBits) | lowPartScaled;
  }

  return recon;
}

//-------------------------------------------------------------------------

void
decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int skipLastLayers,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining,
  const CloudFrame* refFrame,
  const Vec3<int> minimum_position,
  PCCPointSet3& compensatedPointCloud
)
{
  const bool isInter = gbh.interPredictionEnabledFlag;

  PCCPointSet3 predPointCloud;

  if (isInter)
    predPointCloud = refFrame->cloud;

  // init main fifo
  //  -- worst case size is the last level containing every input poit
  //     and each point being isolated in the previous level.
  // NB: some trisoup configurations can generate fewer points than
  //     octree nodes.  Blindly trusting the number of points to guide
  //     the ringbuffer size is problematic.
  // todo(df): derive buffer size from level limit
  size_t ringBufferSize = gbh.footer.geom_num_points_minus1 + 1;
  if (gbh.trisoupNodeSizeLog2(gps))
    ringBufferSize = 5000000;
  pcc::ringbuf<PCCOctree3Node> fifo(ringBufferSize + 1);

  size_t processedPointCount = 0;
  std::vector<uint32_t> values;

  // rotating mask used to enable idcm
  uint32_t idcmEnableMaskInit = /*mkIdcmEnableMask(gps)*/ 0; //NOTE[FT] : set to 0 by construction

  MortonMap3D occupancyAtlas;
  if (gps.neighbour_avail_boundary_log2_minus1) {
    occupancyAtlas.resize(
      gps.adjacent_child_contextualization_enabled_flag,
      gps.neighbour_avail_boundary_log2_minus1 + 1);
    occupancyAtlas.clear();
  }

  Vec3<uint32_t> posQuantBitMasks = 0xffffffff;
  int idcmQp = 0;
  int sliceQp = gbh.sliceQp(gps);
  int nodeQpOffsetsSignalled = !gps.geom_scaling_enabled_flag;

  // generate the list of the node size for each level in the tree
  //  - starts with the smallest node and works up
  std::vector<Vec3<int>> lvlNodeSizeLog2{gbh.trisoupNodeSizeLog2(gps)};
  for (auto split : inReverse(gbh.tree_lvl_coded_axis_list)) {
    Vec3<int> splitStv = {!!(split & 4), !!(split & 2), !!(split & 1)};
    lvlNodeSizeLog2.push_back(lvlNodeSizeLog2.back() + splitStv);
  }
  std::reverse(lvlNodeSizeLog2.begin(), lvlNodeSizeLog2.end());

  // Derived parameter used by trisoup.
  gbh.maxRootNodeDimLog2 = lvlNodeSizeLog2[0].max();

  // the termination depth of the octree phase
  // NB: minNodeSizeLog2 is only non-zero for partial decoding (not trisoup)
  int maxDepth = lvlNodeSizeLog2.size() - skipLastLayers - 1;

  // append a dummy entry to the list so that depth+2 access is always valid
  lvlNodeSizeLog2.emplace_back(lvlNodeSizeLog2.back());

  // NB: this needs to be after the root node size is determined to
  //     allocate the planar buffer
  GeometryOctreeDecoder decoder(gps, gbh, ctxtMem, &arithmeticDecoder);

  // saved state for use with parallel bistream coding.
  // the saved state is restored at the start of each parallel octree level
  std::unique_ptr<GeometryOctreeDecoder> savedState;

  int LPUnumInAxis = 0;
  int log2MotionBlockSize = 0;

  // local motion prediction structure -> LPUs from predPointCloud
  std::vector<std::vector<Vec3<int>>> firstLpuActiveWindow;
  if (isInter && gps.localMotionEnabled) {
    for (int i = 0; i < predPointCloud.getPointCount(); i++) {
      predPointCloud[i] -= gbh.geomBoxOrigin;
    }
    const int extended_window = gps.motion.motion_window_size;

    log2MotionBlockSize = int(log2(gps.motion.motion_block_size));
    const int maxBB = (1 << gbh.maxRootNodeDimLog2) - 1;

    if (gbh.maxRootNodeDimLog2 < log2MotionBlockSize) { // LPU is bigger than root note, must adjust if possible
      int log2MotionBlockSizeMin = int(log2(gps.motion.motion_min_pu_size));
      if (log2MotionBlockSizeMin <= gbh.maxRootNodeDimLog2)
        log2MotionBlockSize = gbh.maxRootNodeDimLog2;
    }

    LPUnumInAxis = (maxBB) >> log2MotionBlockSize;
    if ((LPUnumInAxis << log2MotionBlockSize) != maxBB)
      LPUnumInAxis++;

    // N.B. after this, predPointCloud need to be in same slice boundaries as current slice
    firstLpuActiveWindow = buildActiveWindow(LPUnumInAxis, predPointCloud, extended_window, log2MotionBlockSize);

    // Remove all points outside boundaries
    Box3<int> sliceBox (0, 1 << lvlNodeSizeLog2[0]);
    int numPoints = 0;
    for (int i = 0; i < predPointCloud.getPointCount(); i++) {
      const auto & point = predPointCloud[i];
      if (sliceBox.contains(point)) {
        predPointCloud[numPoints++] = point;
      }
    }
    predPointCloud.resize(numPoints);

    std::cout << "Predictor size = " << predPointCloud.getPointCount() << std::endl;
  }


  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(0);
  node00.pos = int32_t(0);
  node00.predStart = uint32_t(0);
  node00.predEnd = isInter && gps.localMotionEnabled ? predPointCloud.getPointCount() : uint32_t(0);
  node00.numSiblingsMispredicted = 0;
  node00.numSiblingsPlus1 = 8;
  node00.siblingOccupancy = 0;
  node00.qp = 0;
  //node00.idcmEligible = 0; NOTE[FT]: idcmEligible is already set to false at construction

  // local motion
  node00.hasMotion = 0;
  node00.isCompensated = 0;

  // The number of nodes to wait before updating the planar rate.
  // This is to match the prior behaviour where planar is updated once
  // per coded occupancy.
  int nodesBeforePlanarUpdate = 1;

  if (!(isInter && gps.localMotionEnabled && gps.gof_geom_entropy_continuation_enabled_flag) && !gbh.entropy_continuation_flag) {
    decoder.clearMap();
    decoder.resetMap();
  }

  bool planarEligibleKOctreeDepth = 0;
  int numPointsCodedByIdcm = 0;
  const bool checkPlanarEligibilityBasedOnOctreeDepth = false; //NOTE[FT] : FORCING geom_planar_mode_enabled_flag=false
    //gps.geom_planar_mode_enabled_flag
    //&& gps.geom_octree_depth_planar_eligibiity_enabled_flag;

  for (int depth = 0; depth < maxDepth; depth++) {
    int numSubnodes = 0;
    // setup at the start of each level
    auto fifoCurrLvlEnd = fifo.end();
    int numNodesNextLvl = 0;
    Vec3<int32_t> occupancyAtlasOrigin = 0xffffffff;

    // derive per-level node size related parameters
    auto nodeSizeLog2 = lvlNodeSizeLog2[depth];
    auto childSizeLog2 = lvlNodeSizeLog2[depth + 1];
    // represents the largest dimension of the current node
    int nodeMaxDimLog2 = nodeSizeLog2.max();

    auto pointSortMask = qtBtChildSize(nodeSizeLog2, childSizeLog2);

    // if one dimension is not split, atlasShift[k] = 0
    int codedAxesPrevLvl = depth ? gbh.tree_lvl_coded_axis_list[depth - 1] : 7;
    int codedAxesCurLvl = gbh.tree_lvl_coded_axis_list[depth];

    // Determine if this is the level where node QPs are sent
    bool nodeQpOffsetsPresent =
      !nodeQpOffsetsSignalled && decoder.decodeNodeQpOffsetsPresent();

    // record the node size when quantisation is signalled -- all subsequnt
    // coded occupancy bits are quantised
    // after the qp offset, idcm nodes do not receive special treatment
    if (nodeQpOffsetsPresent) {
      nodeQpOffsetsSignalled = true;
      idcmQp = 0;
      posQuantBitMasks = Vec3<uint32_t>((1 << nodeSizeLog2) - 1);
    }

    // Idcm quantisation applies to child nodes before per node qps
    if (!nodeQpOffsetsSignalled) {
      // If planar is enabled, the planar bits are not quantised (since
      // the planar mode is determined before quantisation)
      auto quantNodeSizeLog2 = nodeSizeLog2;
      /*if (gps.geom_planar_mode_enabled_flag) //NOTE[FT] : FORCING geom_planar_mode_enabled_flag=false
        quantNodeSizeLog2 -= 1;*/

      for (int k = 0; k < 3; k++)
        quantNodeSizeLog2[k] = std::max(0, quantNodeSizeLog2[k]);

      // limit the idcmQp such that it cannot overquantise the node
      auto minNs = quantNodeSizeLog2.min();
      idcmQp = gps.geom_base_qp + gps.geom_idcm_qp_offset;
      idcmQp <<= gps.geom_qp_multiplier_log2;
      idcmQp = std::min(idcmQp, minNs * 8);

      posQuantBitMasks = Vec3<uint32_t>((1 << quantNodeSizeLog2) - 1);
    }

    // save context state for parallel coding
    if (depth == maxDepth - 1 - gbh.geom_stream_cnt_minus1)
      if (gbh.geom_stream_cnt_minus1)
        savedState.reset(new GeometryOctreeDecoder(decoder));

    // a new entropy stream starts one level after the context state is saved.
    // restore the saved state and flush the arithmetic decoder
    if (depth > maxDepth - 1 - gbh.geom_stream_cnt_minus1) {
      decoder = *savedState;
      arithmeticDecoder.flushAndRestart();
    }

    // reset the idcm eligibility mask at the start of each level to
    // support multiple streams
    auto idcmEnableMask = /*rotateRight(idcmEnableMaskInit, depth)*/ 0; //NOTE[FT]: still 0

    auto planarDepth = lvlNodeSizeLog2[0] - nodeSizeLog2;
    decoder.beginOctreeLevel(planarDepth);


    // process all nodes within a single level
    for (; fifo.begin() != fifoCurrLvlEnd; fifo.pop_front()) {
      PCCOctree3Node& node0 = fifo.front();

      // decode local motion PU tree
      if (isInter && gps.localMotionEnabled) {

        if (nodeSizeLog2[0] == log2MotionBlockSize) {
          const Vec3<int32_t> pos = node0.pos << nodeSizeLog2;
          const int lpuX = pos[0] >> log2MotionBlockSize;
          const int lpuY = pos[1] >> log2MotionBlockSize;
          const int lpuZ = pos[2] >> log2MotionBlockSize;
          const int lpuIdx = (lpuX * LPUnumInAxis + lpuY) * LPUnumInAxis + lpuZ;

          // no local motion if not enough points
          const bool isLocalEnabled = firstLpuActiveWindow[lpuIdx].size() > 50;
          node0.hasMotion = isLocalEnabled;
        }

        // decode LPU/PU/MV
        if (node0.hasMotion && !node0.isCompensated) {
          decode_splitPU_MV_MC(
            &node0, gps.motion, nodeSizeLog2,
            &arithmeticDecoder, &compensatedPointCloud,
            firstLpuActiveWindow, LPUnumInAxis, log2MotionBlockSize);
        }
      }


      // sort the predictor into eight child partitions
      //  - perform an 8-way counting sort of the current node's points
      //  - (later) map to child nodes

      auto sortPredicate = [=](const PCCPointSet3::Proxy& proxy) {
        const auto& point = *proxy;
        return !!(int(point[2]) & pointSortMask[2])
          | (!!(int(point[1]) & pointSortMask[1]) << 1)
          | (!!(int(point[0]) & pointSortMask[0]) << 2);
      };

      // sort and partition the predictor...
      std::array<int, 8> predCounts = {};

      // ...for local motion
      if (isInter && gps.localMotionEnabled) {
        if (node0.isCompensated) {
          countingSort(
            PCCPointSet3::iterator(&compensatedPointCloud, node0.predStart),  // Need to update the predStar
            PCCPointSet3::iterator(&compensatedPointCloud, node0.predEnd),
            predCounts, [=](const PCCPointSet3::Proxy& proxy) {
            const auto & point = *proxy;
            return !!(int(point[2]) & pointSortMask[2])
              | (!!(int(point[1]) & pointSortMask[1]) << 1)
              | (!!(int(point[0]) & pointSortMask[0]) << 2);
          });
        }
        else {
          countingSort(
            PCCPointSet3::iterator(&predPointCloud, node0.predStart),
            PCCPointSet3::iterator(&predPointCloud, node0.predEnd),
            predCounts, [=](const PCCPointSet3::Proxy& proxy) {
            const auto & point = *proxy;
            return !!(int(point[2]) & pointSortMask[2])
              | (!!(int(point[1]) & pointSortMask[1]) << 1)
              | (!!(int(point[0]) & pointSortMask[0]) << 2);
          });
        }
      }


      // generate the bitmap of child occupancy and count
      // the number of occupied children in node0.
      int predOccupancy = 0;
      int predOccupancyStrong = 0;

      for (int i = 0; i < 8; i++) {
        if (predCounts[i]) {
          predOccupancy |= 1 << i;
        }
        if (predCounts[i] > 2) {
          predOccupancyStrong |= 1 << i;
        }
      }

      bool occupancyIsPredictable =
        predOccupancy && node0.numSiblingsMispredicted <= 5;
      // The predictor may be cleared for the purpose of context
      // selection if the prediction is unlikely to be good.
      // NB: any other tests should use the original prediction.
      int predOccupancyReal = predOccupancy;

      int occupancyIsPredicted = 0;
      int occupancyPrediction = 0;

      if (nodeQpOffsetsPresent) {
        node0.qp = sliceQp;
        node0.qp += decoder.decodeQpOffset() << gps.geom_qp_multiplier_log2;
      }

      /*OctreeNodePlanar planarRef;
      if (isInter)
        setPlanesFromOccupancy(predOccupancy, planarRef);*/

      int shiftBits = QuantizerGeom::qpShift(node0.qp);
      auto effectiveNodeSizeLog2 = nodeSizeLog2 - shiftBits;
      auto effectiveChildSizeLog2 = childSizeLog2 - shiftBits;

      // make quantisation work with qtbt and planar.
      auto codedAxesCurNode = codedAxesCurLvl;
      if (shiftBits != 0) {
        for (int k = 0; k < 3; k++) {
          if (effectiveChildSizeLog2[k] < 0)
            codedAxesCurNode &= ~(4 >> k);
        }
      }

      GeometryNeighPattern gnp{};
      // The position of the node in the parent's occupancy map
      int posInParent = 0;
      posInParent |= (node0.pos[0] & 1) << 2;
      posInParent |= (node0.pos[1] & 1) << 1;
      posInParent |= (node0.pos[2] & 1) << 0;
      posInParent &= codedAxesPrevLvl;

      if (gps.neighbour_avail_boundary_log2_minus1) {
        updateGeometryOccupancyAtlas(
          node0.pos, codedAxesPrevLvl, fifo, fifoCurrLvlEnd, &occupancyAtlas,
          &occupancyAtlasOrigin);

        gnp = makeGeometryNeighPattern(
          gps.adjacent_child_contextualization_enabled_flag,
          node0.pos, codedAxesPrevLvl, occupancyAtlas);

      } else {
        gnp.neighPattern =
          neighPatternFromOccupancy(posInParent, node0.siblingOccupancy);
      }

      bool isDirectMode = false;
      // At the scaling depth, it is possible for a node that has previously
      // been marked as being eligible for idcm to be fully quantised due
      // to the choice of QP.  There is therefore nothing to code with idcm.
      /*if (isLeafNode(effectiveNodeSizeLog2))
        node0.idcmEligible = false;*/ //NOTE[FT]: always false

      OctreeNodePlanar planar;
      //if (!isLeafNode(effectiveNodeSizeLog2)) {
      //  // planar eligibility
      //  bool planarEligible[3] = {false, false, false};

      //  planar.allowPCM = /*(isInter) && (occupancyIsPredictable)
      //    && (planarEligible[0] || planarEligible[1] || planarEligible[2])*/ false;
      //  planar.isPreDirMatch = true;
      //  planar.eligible[0] = /*planarEligible[0]*/ false;
      //  planar.eligible[1] = /*planarEligible[1]*/ false;
      //  planar.eligible[2] = /*planarEligible[2]*/ false;
      //  planar.lastDirIdx =
      //    /*planarEligible[2] ? 2 : (planarEligible[1] ? 1 : 0)*/ 0;

      //  /*decoder.determinePlanarMode(
      //    gps.adjacent_child_contextualization_enabled_flag, planarEligible,
      //    posInParent, gnp, node0, planar, planarRef);*/
      //}

      //if (node0.idcmEligible) {
      //  if (true) //NOTE[FT]: FORCING geom_planar_disabled_idcm_angular_flag to false
      //    isDirectMode = decoder.decodeIsIdcm();
      //  if (isDirectMode) {
      //    auto idcmSize = effectiveNodeSizeLog2;
      //    if (idcmQp) {
      //      node0.qp = idcmQp;
      //      idcmSize = nodeSizeLog2 - QuantizerGeom::qpShift(idcmQp);
      //    }

      //    int numPoints = decoder.decodeDirectPosition(
      //      gps.geom_unique_points_flag, gps.joint_2pt_idcm_enabled_flag,
      //      idcmSize, posQuantBitMasks,
      //      node0, planar,
      //      &pointCloud[processedPointCount]);

      //    // calculating the number of points coded by IDCM for determining eligibility of planar mode
      //    if (checkPlanarEligibilityBasedOnOctreeDepth)
      //      numPointsCodedByIdcm += numPoints;

      //    for (int j = 0; j < numPoints; j++) {
      //      auto& point = pointCloud[processedPointCount++];
      //      for (int k = 0; k < 3; k++)
      //        point[k] += rotateLeft(node0.pos[k], idcmSize[k]);

      //      point = invQuantPosition(node0.qp, posQuantBitMasks, point);
      //    }

      //    // NB: no further siblings to decode by definition of IDCM
      //    if (gps.inferred_direct_coding_mode <= 1)
      //      assert(node0.numSiblingsPlus1 == 1);

      //    // This node has no children, ensure that future nodes avoid
      //    // accessing stale child occupancy data.
      //    if (gps.adjacent_child_contextualization_enabled_flag)
      //      updateGeometryOccupancyAtlasOccChild(
      //        node0.pos, 0, &occupancyAtlas);

      //    continue;
      //  }
      //}

      uint8_t occupancy = 1;
      if (!isLeafNode(effectiveNodeSizeLog2)) {
        // planar mode for current node
        // mask to be used for the occupancy coding
        // (bit =1 => occupancy bit not coded due to not belonging to the plane)
        int planarMask[3] = {0, 0, 0};
        maskPlanar(planar, planarMask, codedAxesCurNode);

        bool flagWord4 =
          gps.neighbour_avail_boundary_log2_minus1 > 0;  //&& intraPredUsed;
        occupancy = decoder.decodeOccupancyFullNeihbourgs(
          gnp, planarMask[0], planarMask[1], planarMask[2],
          /*planar.planarPossible & 1*/ true, /*planar.planarPossible & 2*/ true,
          /*planar.planarPossible & 4*/ true, occupancyAtlas, node0.pos,
          codedAxesPrevLvl, flagWord4,
          gps.adjacent_child_contextualization_enabled_flag, predOccupancy | (predOccupancyStrong <<8) , isInter);
      }

      assert(occupancy > 0);

      // update atlas for child neighbours
      // NB: the child occupancy atlas must be updated even if the current
      //     node has no occupancy coded in order to clear any stale state in
      //     the atlas.
      if (gps.adjacent_child_contextualization_enabled_flag)
        updateGeometryOccupancyAtlasOccChild(
          node0.pos, occupancy, &occupancyAtlas);

      // population count of occupancy for IDCM
      int numOccupied = popcnt(occupancy);

      // calculating the number of subnodes for determining eligibility of planar mode
      if (checkPlanarEligibilityBasedOnOctreeDepth)
        numSubnodes += numOccupied;

      int predFailureCount = popcnt(uint8_t(occupancy ^ predOccupancyReal));
      int predPointsStartIdx = node0.predStart;
      // nodeSizeLog2 > 1: for each child:
      //  - determine elegibility for IDCM
      //  - directly decode point positions if IDCM allowed and selected
      //  - otherwise, insert split children into fifo while updating neighbour state
      for (int i = 0; i < 8; i++) {
        uint32_t mask = 1 << i;
        if (!(occupancy & mask)) {
          predPointsStartIdx += predCounts[i];
          // child is empty: skip
          continue;
        }

        int x = !!(i & 4);
        int y = !!(i & 2);
        int z = !!(i & 1);

        // point counts for leaf nodes are coded immediately upon
        // encountering the leaf node.
        if (isLeafNode(effectiveChildSizeLog2)) {
          int numPoints = 1;

          if (!gps.geom_unique_points_flag) {
            numPoints = decoder.decodePositionLeafNumPoints();
          }

          // the final bits from the leaf:
          Vec3<int32_t> point{
            (node0.pos[0] << !!(codedAxesCurLvl & 4)) + x,
            (node0.pos[1] << !!(codedAxesCurLvl & 2)) + y,
            (node0.pos[2] << !!(codedAxesCurLvl & 1)) + z};

          // remove any padding bits that were not coded
          for (int k = 0; k < 3; k++)
            point[k] = rotateLeft(point[k], effectiveChildSizeLog2[k]);

          point = invQuantPosition(node0.qp, posQuantBitMasks, point);

          for (int i = 0; i < numPoints; ++i)
            pointCloud[processedPointCount++] = point;

          // do not recurse into leaf nodes
          continue;
        }

        // create & enqueue new child.
        fifo.emplace_back();
        auto& child = fifo.back();

        child.qp = node0.qp;
        // only shift position if an occupancy bit was coded for the axis
        child.pos[0] = (node0.pos[0] << !!(codedAxesCurLvl & 4)) + x;
        child.pos[1] = (node0.pos[1] << !!(codedAxesCurLvl & 2)) + y;
        child.pos[2] = (node0.pos[2] << !!(codedAxesCurLvl & 1)) + z;
        child.numSiblingsPlus1 = numOccupied;
        child.siblingOccupancy = occupancy;
        child.numSiblingsMispredicted = predFailureCount;
        child.predStart = predPointsStartIdx;
        child.predEnd = predPointsStartIdx + predCounts[i];
        predPointsStartIdx = child.predEnd;

        //local motion PU inheritance
        child.hasMotion = node0.hasMotion;
        child.isCompensated = node0.isCompensated;

        //if (isInter)
        //  child.idcmEligible = isDirectModeEligible_Inter(
        //    /*gps.inferred_direct_coding_mode*/ 0, nodeMaxDimLog2, gnp.neighPattern,
        //    node0, child, occupancyIsPredictable);
        //else
        //  child.idcmEligible = isDirectModeEligible(
        //    /*gps.inferred_direct_coding_mode*/ 0, nodeMaxDimLog2, gnp.neighPattern,
        //    node0, child, occupancyIsPredictable);

        //if (child.idcmEligible) {
        //  child.idcmEligible &= idcmEnableMask & 1; //NOTE[FT] : stays at 0, whatever the childidcmEligible status
        //  idcmEnableMask = rotateRight(idcmEnableMask, 1);
        //}

        numNodesNextLvl++;
      }
    }
    if (checkPlanarEligibilityBasedOnOctreeDepth)
      planarEligibleKOctreeDepth =
        (ringBufferSize - numPointsCodedByIdcm) * 10 < numSubnodes * 13;

    // Check that one level hasn't produced too many nodes
    // todo(df): this check is too weak to spot overflowing the fifo
    assert(numNodesNextLvl <= ringBufferSize);
  }
  if (!(gps.localMotionEnabled && gps.gof_geom_entropy_continuation_enabled_flag) && !(gps.trisoup_enabled_flag || gbh.entropy_continuation_flag))
    decoder.clearMap();

  // save the context state for re-use by a future slice if required
  ctxtMem = decoder.getCtx();

  // NB: the point cloud needs to be resized if partially decoded
  // OR: if geometry quantisation has changed the number of points
  pointCloud.resize(processedPointCount);

  // return partial coding result
  //  - add missing levels to node positions and inverse quantise
  if (nodesRemaining) {
    auto nodeSizeLog2 = lvlNodeSizeLog2[maxDepth];
    for (auto& node : fifo) {
      node.pos <<= nodeSizeLog2 - QuantizerGeom::qpShift(node.qp);
      node.pos = invQuantPosition(node.qp, posQuantBitMasks, node.pos);
    }
    *nodesRemaining = std::move(fifo);
    return;
  }
}

//-------------------------------------------------------------------------

void
decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const Vec3<int> minimum_position,
  PCCPointSet3& compensatedPointCloud
)
{
  decodeGeometryOctree(
    gps, gbh, 0, pointCloud, ctxtMem, arithmeticDecoder, nullptr,
    refFrame, minimum_position, compensatedPointCloud);
}

//-------------------------------------------------------------------------

void
decodeGeometryOctreeScalable(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int minGeomNodeSizeLog2,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame
)
{
  pcc::ringbuf<PCCOctree3Node> nodes;
  PCCPointSet3 compensatedPointCloud;
  decodeGeometryOctree(
    gps, gbh, minGeomNodeSizeLog2, pointCloud, ctxtMem, arithmeticDecoder,
    &nodes, refFrame, { 0, 0, 0 }, compensatedPointCloud);

  if (minGeomNodeSizeLog2 > 0) {
    size_t size =
      pointCloud.removeDuplicatePointInQuantizedPoint(minGeomNodeSizeLog2);

    pointCloud.resize(size + nodes.size());
    size_t processedPointCount = size;

    if (minGeomNodeSizeLog2 > 1) {
      uint32_t mask = uint32_t(-1) << minGeomNodeSizeLog2;
      for (auto node0 : nodes) {
        for (int k = 0; k < 3; k++)
          node0.pos[k] &= mask;
        node0.pos += 1 << (minGeomNodeSizeLog2 - 1);
        pointCloud[processedPointCount++] = node0.pos;
      }
    } else {
      for (const auto& node0 : nodes)
        pointCloud[processedPointCount++] = node0.pos;
    }
  }
}

//============================================================================

}  // namespace pcc
