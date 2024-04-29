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

#include "geometry_octree.h"

#include <algorithm>
#include <iterator>
#include <climits>

#include "PCCMisc.h"
#include "geometry_params.h"
#include "quantization.h"
#include "tables.h"

namespace pcc {

//============================================================================

Vec3<int>
oneQtBtDecision(
  const QtBtParameters& qtbt,
  Vec3<int> nodeSizeLog2,
  int maxNumQtbtBeforeOt,
  int minDepthQtbt)
{
  int nodeMinDimLog2 = nodeSizeLog2.min();

  if (maxNumQtbtBeforeOt || nodeMinDimLog2 == minDepthQtbt) {
    int nodeMaxDimLog2 = nodeSizeLog2.max();
    for (int k = 0; k < 3; k++) {
      if (nodeSizeLog2[k] == nodeMaxDimLog2)
        nodeSizeLog2[k]--;
    }
  } else // octree partition
    nodeSizeLog2 = nodeSizeLog2 - 1;

  return nodeSizeLog2;
}

//---------------------------------------------------------------------------

void
updateQtBtParameters(
  const Vec3<int>& nodeSizeLog2,
  bool trisoup_enabled_flag,
  int* maxNumQtbtBeforeOt,
  int* minSizeQtbt)
{
  int nodeMinDimLog2 = nodeSizeLog2.min();
  int nodeMaxDimLog2 = nodeSizeLog2.max();

  // max number of qtbt partitions before ot is bounded by difference between
  // max and min node size
  if (*maxNumQtbtBeforeOt > (nodeMaxDimLog2 - nodeMinDimLog2))
    *maxNumQtbtBeforeOt = nodeMaxDimLog2 - nodeMinDimLog2;
  // min depth of qtbt partition is bounded by min node size
  if (*minSizeQtbt > nodeMinDimLog2)
    *minSizeQtbt = nodeMinDimLog2;
  // if all dimensions have same size, min depth of qtbt should be 0
  if (nodeMaxDimLog2 == nodeMinDimLog2) {
    *minSizeQtbt = 0;
  }

  // if trisoup is enabled, perform qtbt first before ot
  if (trisoup_enabled_flag) {
    *maxNumQtbtBeforeOt = nodeMaxDimLog2 - nodeMinDimLog2;
    *minSizeQtbt = 0;
  }
}

//---------------------------------------------------------------------------

std::vector<Vec3<int>>
mkQtBtNodeSizeList(
  const GeometryParameterSet& gps,
  const QtBtParameters& qtbt,
  const GeometryBrickHeader& gbh)
{
  std::vector<Vec3<int>> nodeSizeLog2List;

  // size of the current node (each dimension can vary due to qtbt)
  Vec3<int> nodeSizeLog2 = gbh.rootNodeSizeLog2;
  nodeSizeLog2List.push_back(nodeSizeLog2);

  // update qtbt parameters
  int maxNumQtbtBeforeOt = qtbt.maxNumQtBtBeforeOt;
  int minSizeQtbt = qtbt.minQtbtSizeLog2;
  updateQtBtParameters(
    nodeSizeLog2, qtbt.trisoupEnabled, &maxNumQtbtBeforeOt, &minSizeQtbt);

  while (!isLeafNode(nodeSizeLog2)) {
    if (!gps.qtbt_enabled_flag)
      nodeSizeLog2 -= 1;
    else
      nodeSizeLog2 =
        oneQtBtDecision(qtbt, nodeSizeLog2, maxNumQtbtBeforeOt, minSizeQtbt);

    nodeSizeLog2List.push_back(nodeSizeLog2);

    if (maxNumQtbtBeforeOt)
      maxNumQtbtBeforeOt--;

    // if all dimensions have same size, then use octree for remaining nodes
    if (
      nodeSizeLog2[0] == minSizeQtbt && nodeSizeLog2[0] == nodeSizeLog2[1]
      && nodeSizeLog2[1] == nodeSizeLog2[2])
      minSizeQtbt = -1;
  }

  return nodeSizeLog2List;
}

//============================================================================

const int CtxModelDynamicOBUF::kContextsInitProbability[] = {
  65461, 65160, 64551, 63637, 62426, 60929, 59163, 57141, 54884, 52413, 49753,
  46929, 43969, 40899, 37750, 34553, 31338, 28135, 24977, 21893, 18914, 16067,
  13382, 10883, 8596,  6542,  4740,  3210,  1967,  1023,  388,   75
};

//============================================================================

uint32_t
mkIdcmEnableMask(const GeometryParameterSet& gps)
{
  if (/*!gps.inferred_direct_coding_mode*/ true)  //NOTE[FT]: FORCING idcm to 0, hence !0 is true
    return 0;

  //// intense IDCM requires idcm to be enabled all the time
  //if (gps.inferred_direct_coding_mode != 1)
  //  return 0xffffffff;

  //// if planar is disabled, there is no control over the rate
  //if (/*!gps.geom_planar_mode_enabled_flag*/ true)
  //  return 0xffffffff;
}

//============================================================================
// directional mask depending on the planarity

int
maskPlanarX(const OctreeNodePlanar& planar)
{
  if ((planar.planarMode & 1) == 0)
    return 0;

  return (planar.planePosBits & 1) ? 0x0f : 0xf0;
}

//----------------------------------------------------------------------------

int
maskPlanarY(const OctreeNodePlanar& planar)
{
  if ((planar.planarMode & 2) == 0)
    return 0;

  return (planar.planePosBits & 2) ? 0x33 : 0xcc;
}

//----------------------------------------------------------------------------

int
maskPlanarZ(const OctreeNodePlanar& planar)
{
  if ((planar.planarMode & 4) == 0)
    return 0;

  return (planar.planePosBits & 4) ? 0x55 : 0xaa;
}

//----------------------------------------------------------------------------

// three direction mask
void
maskPlanar(OctreeNodePlanar& planar, int mask[3], int codedAxes)
{
  for (int k = 0; k <= 2; k++) {
    // QTBT does not split in this direction
    //   => infer the mask low for occupancy bit coding
    if (!(codedAxes & (4 >> k))) {
      planar.planePosBits &= ~(1 << k);
      planar.planarMode |= 1 << k;
    }
  }

  mask[0] = maskPlanarX(planar);
  mask[1] = maskPlanarY(planar);
  mask[2] = maskPlanarZ(planar);
}

//============================================================================

void
GeometryOctreeContexts::resetMap()
{
  for (int i = 0; i < 2; i++) {
    int isInter = 2 * (i > 0);
    const int n2 = 6;
    _MapOccupancy[i][0].reset(6 + n2 + 1, 18 - 6 - n2 + isInter);
    _MapOccupancy[i][1].reset(6 + n2 + 1, 18 - 6 - n2 + isInter);
    _MapOccupancy[i][2].reset(6 + n2 + 1, 18 - 6 - n2 + isInter);
    _MapOccupancy[i][3].reset(4 + n2 + 1, 18 - 6 - n2 + isInter);
    _MapOccupancy[i][4].reset(6 + n2 + 1, 18 - 6 - n2 + isInter);
    _MapOccupancy[i][5].reset(6 + n2 + 1, 18 - 6 - n2 + isInter);
    _MapOccupancy[i][6].reset(6 + n2 + 1, 18 - 6 - n2 + isInter);
    _MapOccupancy[i][7].reset(4 + n2 + 1, 18 - 6 - n2 + isInter);

    const int n3 = 5;
    _MapOccupancySparse[i][0].reset(6 + n3 + 1, 9 - n3 + isInter);
    _MapOccupancySparse[i][1].reset(6 + n3 + 1, 12 - n3 + isInter);
    _MapOccupancySparse[i][2].reset(6 + n3 + 1, 12 - n3 + isInter);
    _MapOccupancySparse[i][3].reset(6 + n3 + 1, 11 - n3 + isInter);
    _MapOccupancySparse[i][4].reset(6 + n3 + 1, 9 - n3 + isInter);
    _MapOccupancySparse[i][5].reset(6 + n3 + 1, 12 - n3 + isInter);
    _MapOccupancySparse[i][6].reset(6 + n3 + 1, 12 - n3 + isInter);
    _MapOccupancySparse[i][7].reset(6 + n3 + 1, 11 - n3 + isInter);
  }

  for (int i = 0; i < 5; i++) {
    MapOBUFTriSoup[i][0].reset(14 + 1 + 2, 7);      // flag
    MapOBUFTriSoup[i][1].reset(10 + 1 + 3 + 1 + 2, 6);      // first bit position
    MapOBUFTriSoup[i][2].reset(10 + 1 + 3 + 1 + 2, 6 + 1);  // second bit position
    MapOBUFTriSoup[i][3].reset(10 + 1 + 3 + 1 + 2, 6 + 2);  // third bit position
    MapOBUFTriSoup[i][4].reset(10 + 1 + 3 + 1 + 2, 6 + 3);  // fourth bit position
  }

  const int bufferSize =
    CtxMapDynamicOBUF::kLeafBufferSize * (1 << CtxMapDynamicOBUF::kLeafDepth);
  std::fill_n(_BufferOBUFleaves, bufferSize, 0);
  _OBUFleafNumber = 0;
  std::fill_n(_BufferOBUFleavesTrisoup, bufferSize, 0);
  _OBUFleafNumberTrisoup = 0;

  // octree intra
  const uint8_t initValueOcc0[64] = {  127, 17, 82, 38, 127, 105, 141, 81, 127, 15, 45, 43, 116, 105, 152, 115, 127, 53, 21, 20, 127, 127, 127, 37, 127, 127, 127, 127, 127, 127, 127, 127, 171, 186, 170, 240, 182, 209, 223, 240, 44, 101, 101, 74, 65, 66, 134, 199, 47, 27, 141, 113, 126, 61, 240, 151, 45, 68, 113, 101, 47, 84, 153, 234, };
  const uint8_t initValueOcc1[64] = {  240, 240, 222, 240, 175, 181, 127, 127, 120, 152, 132, 116, 57, 127, 127, 127, 105, 185, 127, 87, 105, 116, 65, 69, 66, 105, 58, 43, 44, 49, 18, 15, 228, 240, 138, 240, 178, 198, 114, 152, 173, 240, 204, 127, 70, 141, 127, 127, 184, 192, 105, 116, 121, 181, 35, 46, 58, 87, 114, 73, 51, 15, 101, 40, };
  const uint8_t initValueOcc2[64] = {  194, 240, 173, 190, 115, 129, 87, 87, 168, 161, 116, 92, 127, 127, 26, 96, 160, 106, 96, 127, 86, 109, 105, 127, 116, 68, 80, 27, 116, 116, 46, 19, 240, 240, 205, 114, 215, 194, 134, 78, 225, 182, 191, 141, 122, 127, 58, 127, 200, 214, 124, 89, 188, 161, 91, 59, 126, 126, 74, 152, 80, 96, 59, 127, };
  const uint8_t initValueOcc3[64] = {  59, 121, 160, 210, 171, 211, 240, 231, 127, 56, 149, 125, 127, 115, 230, 204, 55, 127, 78, 192, 127, 182, 197, 218, 35, 39, 15, 72, 96, 87, 151, 139, 46, 141, 152, 240, 114, 162, 240, 240, 87, 69, 127, 96, 44, 67, 129, 155, 53, 105, 141, 73, 96, 105, 198, 128, 15, 35, 96, 57, 127, 96, 127, 96, };
  const uint8_t initValueOcc4[64] = {  23, 30, 130, 66, 139, 127, 30, 105, 113, 127, 87, 127, 127, 127, 127, 127, 166, 146, 70, 15, 209, 116, 141, 90, 114, 138, 71, 15, 127, 127, 127, 127, 204, 240, 198, 219, 232, 240, 142, 240, 151, 139, 87, 127, 209, 190, 43, 141, 141, 181, 116, 127, 240, 210, 88, 127, 73, 170, 65, 61, 140, 194, 48, 65, };
  const uint8_t initValueOcc5[64] = {  240, 99, 240, 69, 189, 96, 105, 80, 154, 233, 152, 141, 127, 152, 127, 127, 166, 48, 57, 15, 97, 41, 43, 15, 127, 116, 127, 127, 127, 85, 127, 127, 235, 214, 177, 154, 240, 240, 161, 61, 219, 185, 152, 208, 157, 90, 127, 127, 117, 138, 69, 30, 154, 80, 62, 15, 141, 121, 127, 127, 127, 41, 127, 105, };
  const uint8_t initValueOcc6[64] = {  227, 199, 188, 103, 212, 141, 205, 55, 240, 240, 210, 141, 178, 70, 127, 127, 240, 84, 139, 73, 139, 60, 127, 59, 161, 127, 127, 127, 80, 65, 127, 127, 201, 195, 127, 69, 175, 80, 87, 39, 115, 240, 127, 175, 116, 168, 127, 127, 115, 96, 42, 23, 65, 65, 49, 15, 96, 141, 127, 127, 105, 127, 127, 127, };
  const uint8_t initValueOcc7[64] = {  141, 141, 139, 146, 127, 144, 177, 218, 127, 63, 127, 115, 127, 164, 240, 194, 127, 127, 73, 97, 127, 190, 186, 128, 73, 16, 15, 88, 116, 127, 80, 161, 127, 116, 116, 240, 42, 166, 161, 230, 96, 47, 127, 127, 58, 88, 116, 109, 105, 116, 15, 61, 15, 80, 73, 155, 15, 15, 15, 45, 36, 73, 57, 121,};
  _MapOccupancy[0][0].init(initValueOcc0);
  _MapOccupancy[0][1].init(initValueOcc1);
  _MapOccupancy[0][2].init(initValueOcc2);
  _MapOccupancy[0][3].init(initValueOcc3);
  _MapOccupancy[0][4].init(initValueOcc4);
  _MapOccupancy[0][5].init(initValueOcc5);
  _MapOccupancy[0][6].init(initValueOcc6);
  _MapOccupancy[0][7].init(initValueOcc7);


  // octree inter
  const uint8_t initValueOcc0_inter[256] = { 127,127,127,127,15,127,127,168,15,127,127,175,26,127,127,175,127,127,127,127,61,127,127,127,116,127,127,127,15,99,127,208,127,127,127,127,15,127,127,96,15,116,127,116,15,105,127,141,127,127,127,127,116,116,127,127,116,127,127,168,36,183,127,233,127,127,127,127,17,127,127,168,15,127,127,127,15,127,127,127,127,127,127,127,80,127,127,127,80,127,127,127,41,162,127,134,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,45,167,127,189,55,113,127,237,94,136,127,240,103,141,127,240,46,154,127,218,74,93,127,215,87,186,127,240,129,212,127,240,28,116,127,130,75,136,127,139,44,105,127,212,62,127,127,158,41,116,127,116,44,140,127,125,39,168,127,199,74,112,127,240,15,105,127,127,15,66,127,108,142,127,127,141,52,93,127,196,110,141,127,211,31,113,127,212,198,152,127,240,61,169,127,228,15,73,127,73,17,119,127,115,22,108,127,176,52,116,127,172,15,96,127,92,40,113,127,194,70,151,127,234,93,186,127,240, };
  const uint8_t initValueOcc1_inter[256] = { 100,186,127,240,114,185,127,240,80,127,127,240,161,127,127,240,58,127,127,159,74,152,127,240,96,127,127,127,127,127,127,127,16,136,127,240,127,141,127,204,55,137,127,217,127,127,127,175,33,127,127,127,116,127,127,127,105,127,127,127,127,127,127,127,69,141,127,127,45,106,127,197,105,127,127,127,69,127,127,116,65,127,127,116,77,136,127,218,37,127,127,127,31,127,127,127,49,127,127,141,69,127,127,127,15,85,127,143,15,127,127,116,15,127,127,127,15,127,127,127,15,127,127,110,15,127,127,116,103,168,127,240,141,130,127,240,74,153,127,220,175,152,127,240,73,172,127,240,116,127,127,232,52,169,127,189,105,127,127,200,77,190,127,214,127,141,127,240,141,127,127,228,127,127,127,127,42,101,127,230,105,127,127,152,127,127,127,141,127,127,127,127,57,157,127,240,80,127,127,240,42,159,127,187,106,127,127,141,59,114,127,199,57,136,127,218,15,81,127,141,15,127,127,141,57,110,127,171,41,127,127,141,60,128,127,194,53,127,127,152,15,127,127,168,30,141,127,152,60,127,127,152,15,127,127,141, };
  const uint8_t initValueOcc2_inter[256] = { 68,73,127,212,57,172,127,240,87,141,127,216,87,127,127,225,53,116,127,164,53,105,127,205,96,127,127,127,69,127,127,116,82,161,127,182,70,141,127,191,27,129,127,192,50,170,127,199,116,127,127,127,116,127,127,127,15,127,127,122,61,127,127,127,45,116,127,186,27,127,127,161,96,127,127,127,87,116,127,116,81,141,127,180,61,147,127,126,69,127,127,127,80,127,127,127,116,127,127,127,56,145,127,185,20,116,127,116,15,70,127,146,80,127,127,127,51,127,127,127,15,127,127,127,15,127,127,127,149,157,127,240,105,196,127,240,113,193,127,225,59,129,127,203,62,130,127,240,76,127,127,240,55,116,127,188,15,96,127,181,142,165,127,240,141,127,127,230,61,139,127,240,127,127,127,141,69,127,127,161,127,127,127,127,42,136,127,186,127,127,127,127,114,161,127,221,135,152,127,240,42,96,127,171,28,105,127,158,62,156,127,233,73,127,127,186,15,141,127,141,15,127,127,116,56,116,127,182,65,142,127,227,53,116,127,130,130,127,127,130,39,127,127,127,51,127,127,127,34,127,127,127,127,127,127,127, };
  const uint8_t initValueOcc3_inter[256] = { 44,131,127,200,83,127,127,127,125,140,127,224,61,140,127,240,91,165,127,208,66,136,127,202,87,165,127,240,87,176,127,240,57,127,127,127,15,112,127,110,61,127,127,214,42,128,127,221,127,127,127,161,44,152,127,210,119,187,127,207,53,181,127,240,22,116,127,96,105,127,127,127,32,86,127,165,90,170,127,221,116,127,127,116,105,127,127,237,40,137,127,225,40,129,127,240,15,127,127,127,15,127,127,116,15,87,127,142,15,104,127,150,65,127,127,127,55,127,127,127,37,135,127,199,29,148,127,212,38,127,127,85,96,127,127,127,116,127,127,198,135,127,127,240,63,188,127,199,67,126,127,215,111,149,127,240,61,190,127,236,44,127,127,116,31,127,127,127,127,127,127,127,80,127,127,141,15,127,127,106,20,125,127,182,74,127,127,212,43,161,127,205,15,116,127,105,116,127,127,127,73,116,127,142,65,127,127,127,53,127,127,168,53,105,127,96,69,167,127,240,30,110,127,202,15,127,127,127,15,127,127,127,55,127,127,127,15,127,127,127,61,127,127,127,40,127,127,127,96,127,127,127,38,127,127,152, };
  const uint8_t initValueOcc4_inter[256] = { 41,130,127,184,15,127,127,127,92,141,127,213,21,127,127,196,54,127,127,210,127,127,127,116,15,116,127,105,116,127,127,127,43,140,127,192,127,127,127,127,40,127,127,141,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,46,115,127,235,34,134,127,150,25,122,127,190,15,105,127,175,92,204,127,240,75,127,127,195,78,127,127,198,56,127,127,152,96,159,127,240,38,153,127,209,15,152,127,182,15,106,127,164,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,88,168,127,240,78,127,127,240,65,127,127,240,116,127,127,240,101,185,127,234,194,141,127,240,77,143,127,211,202,168,127,240,59,127,127,211,51,140,127,209,51,127,127,127,127,127,127,127,61,154,127,237,87,127,127,232,15,92,127,124,141,127,127,200,87,127,127,161,106,116,127,240,105,127,127,127,116,127,127,127,125,170,127,240,111,141,127,235,53,115,127,188,78,127,127,121,36,127,127,130,69,141,127,206,42,127,127,127,30,127,127,127,28,141,127,188,82,127,127,240,15,81,127,115,58,127,127,127, };
  const uint8_t initValueOcc5_inter[256] = { 93,180,127,240,59,168,127,240,127,141,127,240,32,127,127,112,89,116,127,183,59,116,127,116,127,116,127,116,53,127,127,116,126,175,127,240,155,178,127,240,127,141,127,141,127,127,127,152,127,127,127,127,116,116,127,141,127,127,127,127,127,127,127,127,44,78,127,200,25,91,127,181,40,127,127,141,15,127,127,80,41,127,127,147,15,127,127,127,15,127,127,127,15,127,127,127,80,127,127,127,96,127,127,141,127,127,127,127,127,127,127,127,127,127,127,127,61,127,127,116,127,127,127,127,127,127,127,127,138,151,127,240,99,107,127,240,63,169,127,240,47,115,127,197,85,136,127,240,105,127,127,240,43,126,127,182,46,127,127,116,105,127,127,240,112,175,127,240,127,127,127,168,127,127,127,240,52,141,127,190,34,152,127,212,127,127,127,127,127,127,127,127,62,129,127,240,20,101,127,163,72,100,127,188,25,116,127,139,44,130,127,240,57,127,127,127,15,92,127,210,15,116,127,87,105,127,127,141,56,149,127,175,127,127,127,127,127,127,127,127,105,127,127,127,15,96,127,143,127,127,127,127,55,127,127,127, };
  const uint8_t initValueOcc6_inter[256] = { 103,168,127,240,60,116,127,229,116,158,127,192,53,161,127,159,60,128,127,240,46,122,127,195,105,127,127,240,32,116,127,127,90,168,127,240,118,165,127,240,127,127,127,240,127,127,127,141,65,141,127,212,38,155,127,198,127,127,127,127,61,127,127,136,88,127,127,240,97,105,127,202,68,116,127,196,15,127,127,127,34,127,127,195,30,141,127,194,105,127,127,127,15,127,127,127,88,116,127,207,96,127,127,141,127,127,127,127,127,127,127,127,59,141,127,118,28,96,127,116,127,127,127,127,116,127,127,127,37,147,127,240,74,127,127,226,15,108,127,159,15,127,127,141,39,116,127,236,45,127,127,127,42,127,127,127,15,116,127,127,77,127,127,200,131,204,127,240,127,127,127,127,127,127,127,194,87,127,127,127,127,152,127,168,116,127,127,127,127,127,127,127,63,161,127,219,116,127,127,127,15,106,127,169,15,127,127,116,15,116,127,141,43,127,127,127,15,116,127,127,15,127,127,127,87,127,127,127,127,127,127,141,127,127,127,127,127,127,127,127,65,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127, };
  const uint8_t initValueOcc7_inter[256] = { 127,127,127,127,127,127,127,161,61,127,127,217,106,181,127,227,116,127,127,152,84,161,127,240,115,140,127,240,81,168,127,240,127,127,127,127,15,130,127,177,116,127,127,127,60,117,127,224,127,127,127,127,37,161,127,226,127,127,127,240,53,182,127,240,127,127,127,127,116,127,127,127,15,155,127,159,53,154,127,215,127,127,127,127,127,127,127,224,26,114,127,216,25,96,127,229,40,127,127,127,15,116,127,105,15,116,127,80,15,94,127,216,127,127,127,127,65,116,127,127,57,127,127,127,33,169,127,209,116,127,127,116,127,127,127,127,73,127,127,141,127,127,127,240,41,142,127,170,71,160,127,202,47,114,127,222,98,199,127,240,87,127,127,127,15,105,127,127,127,127,127,127,80,127,127,152,15,116,127,116,15,142,127,207,65,127,127,127,51,121,127,226,43,127,127,127,127,127,127,127,15,105,127,65,49,127,127,127,15,127,127,130,45,127,127,127,19,146,127,200,26,116,127,194,15,127,127,127,15,127,127,127,15,127,127,127,15,127,127,127,15,127,127,127,34,127,127,127,29,127,127,127,55,127,127,105, };
  _MapOccupancy[1][0].init(initValueOcc0_inter);
  _MapOccupancy[1][1].init(initValueOcc1_inter);
  _MapOccupancy[1][2].init(initValueOcc2_inter);
  _MapOccupancy[1][3].init(initValueOcc3_inter);
  _MapOccupancy[1][4].init(initValueOcc4_inter);
  _MapOccupancy[1][5].init(initValueOcc5_inter);
  _MapOccupancy[1][6].init(initValueOcc6_inter);
  _MapOccupancy[1][7].init(initValueOcc7_inter);


  // Intra
  const uint8_t initValue0[128] = { 15,15,15,15,15,15,15,15,15,15,42,96,71,37,15,15,22,51,15,15,30,27,15,15,64,15,48,15,224,171,127,24,127,34,80,46,141,44,66,49,127,116,140,116,105,39,127,116,114,46,172,109,60,73,181,161,112,65,240,159,127,127,127,87,183,127,116,116,195,88,152,141,228,141,127,80,127,127,160,92,224,167,129,135,240,183,240,184,240,240,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127 };
  const uint8_t initValue1[64] = { 116,127,118,15,104,56,97,15,96,15,29,15,95,15,46,15,196,116,182,53,210,104,163,69,169,15,114,15,121,15,167,63,240,127,184,92,240,163,197,77,239,73,179,59,213,48,185,108,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127, };
  const uint8_t initValue2[128] = { 141,127,127,127,189,81,36,127,143,105,103,116,201,60,38,116,116,127,15,127,153,59,15,116,69,105,15,127,158,93,36,79,141,161,116,127,197,102,53,127,177,125,88,79,209,75,102,28,95,74,72,56,189,62,78,18,88,116,28,45,237,100,152,35,141,240,127,127,208,133,101,141,186,210,168,98,201,124,138,15,195,194,103,94,229,82,167,23,92,197,112,59,185,87,156,79,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127 };
  MapOBUFTriSoup[0][0].init(initValue0);
  MapOBUFTriSoup[0][1].init(initValue1);
  MapOBUFTriSoup[0][2].init(initValue2);

  // inter flag
  const uint8_t initValue0_inter0[128] = { 15,15,15,15,15,15,15,15,15,15,15,45,26,15,15,15,15,38,15,15,15,15,15,15,59,15,24,15,90,62,116,15,80,15,47,75,80,15,71,51,127,87,155,51,69,15,116,80,51,15,136,66,84,34,141,57,60,42,184,118,127,87,105,65,174,61,127,87,196,48,141,127,204,121,105,61,127,127,91,71,196,132,153, 87, 198, 156, 155, 115, 228, 145, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 };
  const uint8_t initValue0_inter1[128] = { 17,53,23,43,81,56,15,80,134,62,80,96,177,62,15,28,80,105,63,15,115,72,101,31,185,97,154,66,240,212,127,73,127,51,102,81,127,73,152,121,127,127,175,176,121,53,127,116,123,117,149,70,119,124,200,142,178,138,240,197,127,116,127,116,182,154,127,105,235,141,186,127,240,185,116, 87, 161, 116, 167, 128, 240, 219, 167, 183, 229, 171, 240, 214, 240, 234, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 };
  MapOBUFTriSoup[1][0].init(initValue0_inter0);
  MapOBUFTriSoup[2][0].init(initValue0_inter1);

  // bit1 pos
  const uint8_t initValue1_inter0[64] = { 141,127,39,15,96,44,15,15,111,15,15,15,47,15,21,15,127,127,64,42,123,57,36,17,76,15,20,15,57,15,132,15,175,127,105,49,229,110,132,54,207,51,73,32,119,38,100,61,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,};
  const uint8_t initValue1_inter1[64] = { 127,127,39,46,110,80,28,15,80,15,48,15,96,15,62,28,152,127,106,73,139,61,96,72,100,15,61,43,56,15,104,49,203,127,69,96,207,133,96,82,130,85,141,66,173,78,192,121,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,};
  const uint8_t initValue1_inter2[64] = { 127,127,153,105,127,105,149,46,121,55,144,25,127,47,182,81,175,116,197,105,181,153,214,133,116,132,200,135,152,75,134,120,240,116,187,116,240,210,202,86,240,140,203,137,194,68,203,161,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,};
  const uint8_t initValue1_inter3[64] = { 152,127,180,127,127,127,193,116,141,87,178,58,141,92,217,80,223,127,240,116,238,180,240,99,151,97,200,96,140,66,224,102,240,127,240,130,240,170,240,157,240,142,240,129,240,113,234,116,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,};
  MapOBUFTriSoup[1][1].init(initValue1_inter0);
  MapOBUFTriSoup[2][1].init(initValue1_inter1);
  MapOBUFTriSoup[3][1].init(initValue1_inter2);
  MapOBUFTriSoup[4][1].init(initValue1_inter3);

}

//============================================================================

void
GeometryOctreeContexts::clearMap()
{
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 8; i++) {
      _MapOccupancy[j][i].clear();
      _MapOccupancySparse[j][i].clear();
    }

  std::cout << "Size used buffer OBUF LEAF = " << _OBUFleafNumber << "\n";

  for (int i = 0; i < 5; i++) {
    MapOBUFTriSoup[i][0].clear();
    MapOBUFTriSoup[i][1].clear();
    MapOBUFTriSoup[i][2].clear();
    MapOBUFTriSoup[i][3].clear();
    MapOBUFTriSoup[i][4].clear();
  }

  // colocated edges
  refFrameEdgeKeys.clear();
  refFrameEdgeValue.clear();
  refFrameNodeKeys.clear();
  refFrameCentroValue.clear();
}

//============================================================================

}  // namespace pcc
