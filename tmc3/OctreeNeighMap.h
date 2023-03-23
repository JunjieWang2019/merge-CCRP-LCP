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

#include "PCCMath.h"
#include "geometry_octree.h"
#include "tables.h"

#include <cassert>
#include <memory>
#include <vector>

namespace pcc {

//============================================================================

struct GeometryNeighPattern {
  // Mask indicating presence of neigbours of the corresponding tree node
  //    32 8 (y)
  //     |/
  //  2--n--1 (x)
  //    /|
  //   4 16 (z)
  uint8_t neighPattern;

  // mask indicating the number of external child neighbours
  uint8_t adjacencyGt0;
  uint8_t adjacencyGt1;

  // mask indicating unoccupied external child neighbours
  uint8_t adjacencyUnocc;

  // occupancy map of {x-1, y-1, z-1} neighbours.
  uint8_t adjNeighOcc[3];
};


//============================================================================
struct OctreeNeighours {

  int occLeft = 0;
  int occFront = 0;
  int occBottom = 0;

  int occL = 0;
  int occF = 0;
  int occB = 0;
  int occOrLFBfb = 0;

  int edgeBits = 0;

  int N3 = 0;
  int N2 = 0;

  int neighPatternLFB = 0;

  int neighb20 = 0;
};




//============================================================================

void prepareGeometryAdvancedNeighPattern(
  const RasterScanContext::occupancy& occ,
  OctreeNeighours& octreeNeighours
  );

void makeGeometryAdvancedNeighPattern0(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int &ctx1,
  int &ctx2,
  bool& Sparse);

void makeGeometryAdvancedNeighPattern1(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int &ctx1,
  int &ctx2,
  bool& Sparse);

void makeGeometryAdvancedNeighPattern2(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int &ctx1,
  int &ctx2,
  bool& Sparse);

void makeGeometryAdvancedNeighPattern3(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int &ctx1,
  int &ctx2,
  bool& Sparse);

void makeGeometryAdvancedNeighPattern4(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int &ctx1,
  int &ctx2,
  bool& Sparse);

void makeGeometryAdvancedNeighPattern5(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int &ctx1,
  int &ctx2,
  bool& Sparse);

void makeGeometryAdvancedNeighPattern6(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int &ctx1,
  int &ctx2,
  bool& Sparse);

void makeGeometryAdvancedNeighPattern7(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int &ctx1,
  int &ctx2,
  bool& Sparse);

}  // namespace pcc
