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
#include "TMC3.h"

#include "motionWip.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <set>
#include <vector>
#include <unordered_map>

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "geometry_octree.h"

namespace pcc {

//============================================================================
struct MotionEntropy {
  // global
  AdaptiveBitModelFast _ctxIsWorld;

  //local
  AdaptiveBitModel splitPu;
  StaticBitModel mvSign;
  AdaptiveBitModel mvIsZero;
  AdaptiveBitModel mvIsOne;
  AdaptiveBitModel mvIsTwo;
  AdaptiveBitModel mvIsThree;
  StaticBitModel expGolombV0;
  AdaptiveBitModel expGolombV[6];
  AdaptiveBitModel _ctxLocalMV;

  MotionEntropy();
};

//----------------------------------------------------------------------------
class MotionEntropyEncoder : public MotionEntropy {
public:
  MotionEntropyEncoder(EntropyEncoder* arithmeticEncoder)
    : _arithmeticEncoder(arithmeticEncoder)
  {}

  // local
  void encodeSplitPu(int symbol);
  void encodeVector(const Vec3<int>& mv);

private:
  EntropyEncoder* _arithmeticEncoder;
};

//----------------------------------------------------------------------------
class MotionEntropyDecoder : public MotionEntropy {
public:
  MotionEntropyDecoder(EntropyDecoder* arithmeticDecoder)
    : _arithmeticDecoder(arithmeticDecoder)
  {}

  //local
  bool decodeSplitPu();
  void decodeVector(Vec3<int>* mv);

private:
  EntropyDecoder* _arithmeticDecoder;
};

//----------------------------------------------------------------------------

struct MotionEntropyEstimate {
  double hMvIsZero[2];
  double hMvIsOne[2];
  double hMvIsTwo[2];
  double hMvIsThree[2];
  double hExpGolombV[6][2];
  double hSplitPu[2];

  MotionEntropyEstimate(const MotionEntropy& codec);


  void prepareEstimateVector(
    double* LUT_MVestimate,
    int wSize,
    int boundPrefix,
    int boundSuffix) const;

  double estimateVector(
    const Vec3<int>& mv,
    double* LUT_MVestimate) const;
};

//============================================================================
MotionEntropy::MotionEntropy()
{
  //local
  splitPu.reset(true);
  mvIsZero.reset(true);
  mvIsOne.reset(true);
  mvIsTwo.reset(true);
  mvIsThree.reset(true);
  for (int i = 0; i < 6; i++)
    expGolombV[i].reset(true);
}

//----------------------------------------------------------------------------
MotionEntropyEstimate::MotionEntropyEstimate(const MotionEntropy& codec)
{
  codec.splitPu.getEntropy(hSplitPu);
  codec.mvIsZero.getEntropy(hMvIsZero);
  codec.mvIsOne.getEntropy(hMvIsOne);
  codec.mvIsTwo.getEntropy(hMvIsTwo);
  codec.mvIsThree.getEntropy(hMvIsThree);
  for (int i = 0; i < 6; i++)
    codec.expGolombV[i].getEntropy(hExpGolombV[i]);
}

//----------------------------------------------------------------------------
inline void
MotionEntropyEncoder::encodeSplitPu(int symbol)
{
  _arithmeticEncoder->encode(symbol, splitPu);
}

//----------------------------------------------------------------------------

inline bool
MotionEntropyDecoder::decodeSplitPu()
{
  return _arithmeticDecoder->decode(splitPu);
}

//----------------------------------------------------------------------------

inline void
MotionEntropyEncoder::encodeVector(
  const Vec3<int>& mv)
{
  for (int comp = 0; comp < 3; comp++) {
    int v = mv[comp];
    if (v == 0) {
      _arithmeticEncoder->encode(1, mvIsZero);
    }
    else {
      _arithmeticEncoder->encode(0, mvIsZero);

      _arithmeticEncoder->encode(v < 0, mvSign);
      if (v < 0)
        v = -v;
      v--;
      _arithmeticEncoder->encode(v == 0, mvIsOne);

      if (!v) {
        continue;
      }
      v--;

      _arithmeticEncoder->encode(v == 0, mvIsTwo);
      if (!v) {
        continue;
      }
      v--;

      _arithmeticEncoder->encode(v == 0, mvIsThree);
      if (!v) {
        continue;
      }
      v--;

      // expGolomb on |v|-1 with truncation
      _arithmeticEncoder->encodeExpGolomb(uint32_t(v), 1, _ctxLocalMV);
    }
  }
}

//----------------------------------------------------------------------------

inline void
MotionEntropyDecoder::decodeVector(Vec3<int>* mv)
{
  for (int comp = 0; comp < 3; comp++) {
    if (_arithmeticDecoder->decode(mvIsZero)) {
      (*mv)[comp] = 0;
      continue;
    }
    bool sign = _arithmeticDecoder->decode(mvSign);

    if (_arithmeticDecoder->decode(mvIsOne)) {
      (*mv)[comp] = sign ? -1 : 1;
      continue;
    }
    if (_arithmeticDecoder->decode(mvIsTwo)) {
      (*mv)[comp] = sign ? -2 : 2;
      continue;
    }
    if (_arithmeticDecoder->decode(mvIsThree)) {
      (*mv)[comp] = sign ? -3 : 3;
      continue;
    }

    int v = 4 + _arithmeticDecoder->decodeExpGolomb(1, _ctxLocalMV);
    if (sign)
      v = -v;
    (*mv)[comp] = v;
  }
}



//----------------------------------------------------------------------------

void
MotionEntropyEstimate::prepareEstimateVector(
  double* LUT_MVestimate,
  int wSize,
  int boundPrefix,
  int boundSuffix) const
{
  for (int v0 = 0; v0 <= wSize; v0++) {
    int v = v0;

    if (!v) {
      LUT_MVestimate[v0] = hMvIsZero[1];
    }
    else {
      double r = hMvIsZero[0] + 1.;  // unpredictable sign
      v--;

      r += hMvIsOne[v == 0];
      if (!v) {
        LUT_MVestimate[v0] = r;
        continue;
      }
      v--;

      r += hMvIsTwo[v == 0];
      if (!v) {
        LUT_MVestimate[v0] = r;
        continue;
      }
      v--;

      r += hMvIsThree[v == 0];
      if (!v) {
        LUT_MVestimate[v0] = r;
        continue;
      }
      v--;

      int k = 1;  // initially, expgolomb order
      int ctx_number = 0;
      int num_bit_prefix = 0;
      while (1) {
        if (boundPrefix && v >= static_cast<unsigned int>(1 << k)) {
          // unary part
          r += hExpGolombV[ctx_number][1];
          v = v - (1 << k);
          k++;

          ctx_number++;
          if (ctx_number > 5)
            ctx_number = 5;
          num_bit_prefix++;
          continue;
        }

        // terminated zero of unary part + suffix
        if (num_bit_prefix < boundPrefix)
          r += hExpGolombV[ctx_number][0];
        if (num_bit_prefix == boundPrefix)
          k = boundSuffix;
        r += k;
        break;
      }
      LUT_MVestimate[v0] = r;
    }
  }
}


//----------------------------------------------------------------------------

double
MotionEntropyEstimate::estimateVector(
  const Vec3<int>& mv,
  double* LUT_MVestimate) const
{
  return LUT_MVestimate[std::abs(mv[0])] + LUT_MVestimate[std::abs(mv[1])] + LUT_MVestimate[std::abs(mv[2])];
}

//============================================================================
static const int LUT_LOG2[64]{
  INT_MIN, 0,  16, 25, 32, 37, 41, 45, 48, 51, 53, 55, 57, 59, 61, 63,
  64,      65, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 79,
  80,      81, 81, 82, 83, 83, 84, 85, 85, 86, 86, 87, 87, 88, 88, 89,
  89,      90, 90, 91, 91, 92, 92, 93, 93, 93, 94, 94, 95, 95, 95, 96};

//----------------------------------------------------------------------------
inline int
plus1log2shifted4(int x)
{
  if (x < 62)
    return LUT_LOG2[x + 1];

  x++;
  int result = 0;
  while (x >= 64) {
    x >>= 1;
    result += 16;
  }

  return result + LUT_LOG2[x];
}


//----------------------------------------------------------------------------
int
roundIntegerHalfInf(const double x)
{
  return (x >= 0) ? int(x + 0.5) : -int(-x + 0.5);
}

//----------------------------------------- LOCAL MOTION -------------------


//============================================================================

int
deriveMotionMaxPrefixBits(const GeometryParameterSet::Motion& param)
{
  int max_MV =
    std::max(0, param.motion_window_size - 1);
  if (max_MV >= 256)
    return 31;

  static int LUT_bound_prefix[256] = {
    0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7};

  return LUT_bound_prefix[max_MV];
}

int
deriveMotionMaxSuffixBits(const GeometryParameterSet::Motion& param)
{

  int max_MV =
    std::max(0, param.motion_window_size  - 1);
  if (max_MV >= 256)
    return 31;

  static int LUT_bound_suffix[256] = {
    0, 1, 0, 1, 2, 2, 0, 1, 2, 2, 3, 3, 3, 3, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 1};

  return LUT_bound_suffix[max_MV];
}

//----------------------------------------------------------------------------

static double
find_motion(
  const GeometryParameterSet::Motion& param,
  const MotionEntropyEstimate& motionEntropy,
  const std::vector<Vec3<int>>& Window,
  const std::vector<Vec3<int>>& Block0,
  int x0,
  int y0,
  int z0,
  int local_size,
  PUtree* local_PU_tree,
  int8_t* bufferPoints,
  double LUT_MVestimate[64])
{
  if (!Window.size())
    return DBL_MAX;

  // ---------------------------- test no split --------------------
  // test V=0 and  dynamic planing at Block0-Window +/- window_size
  double cost_NoSplit = DBL_MAX;
  int wSize = param.motion_window_size;

  int max_distance = 3 * wSize;
  Vec3<int> bestV_NoSplit = Vec3<int>(0, 0, 0);
  Vec3<int> bestV_NoSplit_2nd = Vec3<int>(0, 0, 0);

  std::vector<int> blockEnds;
  blockEnds.push_back(0);
  int blockPos = 0;
  int8_t* pBuffer = bufferPoints;

  const auto* itB = Block0.data();
  int jumpBlock = 1 + (Block0.size() >> param.decimate);  // (kind of) random sampling of the original block to code
  int Nsamples = 0;
  int t0, t1, t2;
  int a0, a1, a2;

  int wSize1 = std::max(-62, std::min(62, wSize));
  int wSize2 = wSize1 >> 1;
  int wSize3 = 3 * wSize1 >> 2;
  int cntFar2 = 0;
  int cntFar3 = 0;

  int Dist = 0;
  for (int Nb = 0; Nb < Block0.size(); Nb += jumpBlock, itB += jumpBlock, Nsamples++) {
    Vec3<int> b = *itB;
    int min_d = max_distance;
    //for (const auto w : Window) {
    for (int idxW = 0; idxW < Window.size(); idxW++) {
      auto w = Window[idxW];
      t0 = b[0] - w[0];
      t1 = b[1] - w[1];
      t2 = b[2] - w[2];
      a0 = std::abs(t0);
      a1 = std::abs(t1);
      a2 = std::abs(t2);

      if (a0 <= wSize1 && a1 <= wSize1 && a2 <= wSize1) {

        bool flag = true;
        if (a0 > wSize2 || a1 > wSize2 || a2 > wSize2) {
          flag = (cntFar2 & 3) == 0;
          cntFar2++;

        }
        if (a0 > wSize3 || a1 > wSize3 || a2 > wSize3) {
          flag = (cntFar3 & 7) == 0;
          cntFar3++;
        }

        if (flag) {
          *pBuffer++ = t0;
          *pBuffer++ = t1;
          *pBuffer++ = t2;
          blockPos++;
        }
      }

      a0 += a1 + a2;  // replaces int local_d
      if (a0 < min_d)
        min_d = a0;
    } // end loop on window

    blockEnds.push_back(blockPos);
    Dist += plus1log2shifted4(min_d);  // 1/0.0625 = 16 times log
  } // loop on points of block

  double d = jumpBlock * Dist * 0.0625 + param.lambda * motionEntropy.estimateVector(Vec3<int>(0, 0, 0), LUT_MVestimate);

  // set loop search parameters
  double best_d[2] = { d, d };
  std::set<int32_t> list_tested = { ((0 + 128) << 16) + ((0 + 128) << 8) + 0 + 128 };

  const int searchPattern[3 * 18] = { 1,0,0,   0,1,0,  0,-1,0,  -1,0,0,  0,0,1, 0,0,-1,   1,1,0,     1,-1,0,  -1,1,0,  -1,-1,0,   0,1,1,    0,-1,1,  0,1,-1,  0,-1,-1,   -1,0,1,  1,0,-1,  -1,0,-1,  1,0,1 };
  //const int searchPattern[3 * 6] = {   1,0,0,   0,1,0,  0,-1,0,    -1,0,0,    0,0,1,   0,0,-1};
  //const int searchPattern[3 * 26] = { 1,0,0,   0,1,0,  0,-1,0,   -1,0,0,   0,0,1,   0,0,-1,  1,1,0,  1,-1,0,  -1,1,0,   -1,-1,0,  1,0,1,  0,1,1,  0,-1,1,  -1,0,1,  1,0,-1,  0,1,-1,  0,-1,-1,  -1,0,-1,  1,1,1,  1,1,-1,  1,-1,1,  1,-1,-1,  -1,1,1,  -1,1,-1,  -1,-1,1,  -1,-1,-1 };

  // loop MV search
  int see = 0;
  int Amotion = param.Amotion0;
  while (Amotion >= 1) {

    // pick one MV best seeds
    Vec3<int> Vs = see == 0 ? bestV_NoSplit : bestV_NoSplit_2nd;
    see = 1;

    // loop on searchPattern
    const int* pSearch = searchPattern;
    bool flagBetter = false;
    for (int t = 0; t < 18; t++, pSearch += 3) {

      Vec3<int> V = Vs + Vec3<int>(pSearch[0] * Amotion, pSearch[1] * Amotion, pSearch[2] * Amotion);
      int32_t V1D = ((V[0] + 128) << 16) + ((V[1] + 128) << 8) + V[2] + 128;
      if (list_tested.find(V1D) != list_tested.end())
        continue;
      list_tested.insert(V1D);

      int Vx = V[0];
      int Vy = V[1];
      int Vz = V[2];
      if (std::abs(Vx) > wSize || std::abs(Vy) > wSize || std::abs(Vz) > wSize)
        continue;  // ensure MV does not go out of the window

      Dist = 0;
      int index = 0;
      pBuffer = bufferPoints;
      for (int Nb = 1; Nb <= Nsamples; Nb++) {
        int min_d = max_distance;
        while (index < blockEnds[Nb]) {
          int local_d = std::abs(pBuffer[0] + Vx) + std::abs(pBuffer[1] + Vy) + std::abs(pBuffer[2] + Vz);
          if (local_d < min_d)
            min_d = local_d;
          index++;
          pBuffer += 3;
        }
        Dist += plus1log2shifted4(min_d);  // 1/0.0625 = 16 times log
      }
      d = jumpBlock * Dist * 0.0625 + param.lambda * motionEntropy.estimateVector(V, LUT_MVestimate);

      // keep 2 best MV
      if (d < best_d[0]) {
        best_d[1] = best_d[0];
        bestV_NoSplit_2nd = bestV_NoSplit;
        best_d[0] = d;
        bestV_NoSplit = V;
        see = 0;
        flagBetter = true;
        break;
      }
      else if (d < best_d[1]) {
        best_d[1] = d;
        bestV_NoSplit_2nd = V;
      }

    }  // end loop on searchPattern

    // log reduction of search range
    if (!flagBetter && see == 1) {
      Amotion >>= 1;
      if (!bestV_NoSplit[0] && !bestV_NoSplit[1] && !bestV_NoSplit[2]) {
        bool flag = Amotion > 1;
        Amotion >>= 1;
        if (flag)
          Amotion = std::max(Amotion, 1);
      }
    }

  }  // end loop MV search

  cost_NoSplit = best_d[0];


  // ---------------------------- test split --------------------
  double cost_Split = DBL_MAX;
  PUtree* Split_PU_tree = new PUtree;  // local split tree

  if (local_size > param.motion_min_pu_size && Block0.size() >= 8) {
    // condition on number of points for search acceleration
    int local_size1 = local_size >> 1;

    Vec3<int> xyz0 = Vec3<int>(x0, y0, z0);

    std::array<Vec3<int>, 8> list_xyz = {
      Vec3<int>(0, 0, 0) + xyz0,
      Vec3<int>(0, 0, local_size1) + xyz0,
      Vec3<int>(0, local_size1, 0) + xyz0,
      Vec3<int>(0, local_size1, local_size1) + xyz0,
      Vec3<int>(local_size1, 0, 0) + xyz0,
      Vec3<int>(local_size1, 0, local_size1) + xyz0,
      Vec3<int>(local_size1, local_size1, 0) + xyz0,
      Vec3<int>(local_size1, local_size1, local_size1) + xyz0
    };

    // loop on 8 child PU
    cost_Split = 0.;

    std::vector<Vec3<int>> Block1;
    Block1.reserve(Block0.size());
    std::vector<Vec3<int>> Window1;
    Window1.reserve(Window.size());

    for (int t = 0; t < 8; t++) {
      // child PU coordinates
      int x1 = list_xyz[t][0];
      int y1 = list_xyz[t][1];
      int z1 = list_xyz[t][2];

      // block for child PU
      int xhigh = x1 + local_size1;
      int yhigh = y1 + local_size1;
      int zhigh = z1 + local_size1;

      Block1.resize(0);
      for (const auto& b : Block0) {
        if (b[0] >= x1 && b[0] < xhigh && b[1] >= y1 && b[1] < yhigh && b[2] >= z1 && b[2] < zhigh)
          Block1.push_back(b);
      }

      cost_Split += 1.0;  // the cost due to not coding the occupancy with inter pred

      if (!Block1.size()) {  // empty PU
        Split_PU_tree->popul_flags.push_back(0);
        continue;
      }
      Split_PU_tree->popul_flags.push_back(1);

      // window for child PU
      wSize = param.motion_window_size;
      xhigh += wSize;
      yhigh += wSize;
      zhigh += wSize;
      int xx1 = x1 - wSize;
      int yy1 = y1 - wSize;
      int zz1 = z1 - wSize;

      Window1.resize(0);
      for (const auto& b : Window) {
        if (b[0] >= xx1 && b[0] < xhigh && b[1] >= yy1 && b[1] < yhigh && b[2] >= zz1 && b[2] < zhigh)
          Window1.push_back(b);
      }
      cost_Split += find_motion(param, motionEntropy, Window1, Block1, x1, y1, z1, local_size1, Split_PU_tree, bufferPoints, LUT_MVestimate);
    }
  }

  // ---------------------------- choose split vs no split --------------------
  if (local_size > param.motion_min_pu_size) {
    cost_NoSplit +=
      param.lambda * motionEntropy.hSplitPu[0];  // cost no split flag
  }
  cost_Split += param.lambda * motionEntropy.hSplitPu[1];  // cost split flag

  if (local_size <= param.motion_min_pu_size|| cost_NoSplit <= cost_Split) {  // no split
    // push non split flag, only if size>size_min
    if (local_size > param.motion_min_pu_size) {
      local_PU_tree->split_flags.push_back(0);
    }
    // push MV
    local_PU_tree->MVs.push_back(bestV_NoSplit);

    delete Split_PU_tree;
    return cost_NoSplit;
  }
  else {
    // split
    local_PU_tree->split_flags.push_back(1);  // push split PU flag

    // append Split_PU_tree to  local_PU_tree
    for (const auto& f : Split_PU_tree->popul_flags)
      local_PU_tree->popul_flags.push_back(f);
    for (const auto& f : Split_PU_tree->split_flags)
      local_PU_tree->split_flags.push_back(f);
    for (const auto& v : Split_PU_tree->MVs)
      local_PU_tree->MVs.push_back(v);

    delete Split_PU_tree;
    return cost_Split;
  }
}

//----------------------------------------------------------------------------
void
extracPUsubtree(
  const GeometryParameterSet::Motion& param,
  PUtree* local_PU_tree,
  int block_size,
  int& pos_fs,
  int& pos_fp,
  int& pos_MV,
  PUtree* destination_tree)
{
  // non-split terminal case
  if (
    block_size <= param.motion_min_pu_size
    || !local_PU_tree->split_flags[pos_fs]) {
    if (block_size > param.motion_min_pu_size) {
      destination_tree->split_flags.push_back(0);
      pos_fs++;
    }

    destination_tree->MVs.push_back(local_PU_tree->MVs[pos_MV]);
    pos_MV++;
    return;
  }

  // split case
  destination_tree->split_flags.push_back(1);
  pos_fs++;

  // loop on 8 children
  for (int s = 0; s < 8; s++) {
    if (local_PU_tree->popul_flags[pos_fp]) {  // child populated
      destination_tree->popul_flags.push_back(1);
      pos_fp++;

      extracPUsubtree(
        param, local_PU_tree, block_size >> 1, pos_fs, pos_fp, pos_MV,
        destination_tree);
    } else {  // child not pouplated
      destination_tree->popul_flags.push_back(0);
      pos_fp++;
    }
  }
}


//----------------------------------------------------------------------------
void
buildActiveWindowAndBoundToBB(
  std::vector<std::vector<Vec3<int>>>& lpuActiveWindow,
  int& LPUnumInAxis,
  const int maxBB,
  PCCPointSet3& predPointCloud,
  int motion_window_size,
  const int log2MotionBlockSize,
  Vec3<int> lvlNodeSizeLog2,
  point_t BBorig)
{
  // prepare displacement list
  std::vector<int> th_dists;
  for (int i = 0; i <= motion_window_size >> log2MotionBlockSize; ++i) {
    th_dists.push_back(i << log2MotionBlockSize);
    if (i) th_dists.push_back(-(i << log2MotionBlockSize));
  }
  if (motion_window_size != (motion_window_size >> log2MotionBlockSize) << log2MotionBlockSize) {
    th_dists.push_back(motion_window_size);
    if (motion_window_size)
      th_dists.push_back(-motion_window_size);
  }
  std::sort(th_dists.begin(), th_dists.end());

  // initialize search windows for LPUs
  LPUnumInAxis = (maxBB) >> log2MotionBlockSize;
  if ((LPUnumInAxis << log2MotionBlockSize) != maxBB)
    LPUnumInAxis++;
  lpuActiveWindow.resize(LPUnumInAxis * LPUnumInAxis * LPUnumInAxis);

  // pile up points in Windows for LPU; at the same time remove points outside the box
  int boundx = 1 << lvlNodeSizeLog2[0];
  int boundy = 1 << lvlNodeSizeLog2[1];
  int boundz = 1 << lvlNodeSizeLog2[2];

  int numPoints = 0;
  for (size_t i = 0; i < predPointCloud.getPointCount(); ++i) {
    const auto point = predPointCloud[i] - BBorig;

    // keep only pred poitns in BB
    if (point[0] >= 0 && point[1] >= 0 && point[2] >= 0 && point[0] < boundx && point[1] < boundy && point[2] < boundz) {
      predPointCloud[numPoints++] = point;
    }

    // build up window search for each LPU
    int oldx = INT32_MIN;
    for (size_t m = 0; m < th_dists.size(); m++) {
      int xidx = (point[0] + th_dists[m]) >> log2MotionBlockSize;
      if (oldx == xidx || xidx < 0 || xidx >= LPUnumInAxis)
        continue;

      int oldy = INT32_MIN;
      for (size_t n = 0; n < th_dists.size(); n++) {
        int yidx = (point[1] + th_dists[n]) >> log2MotionBlockSize;
        if (oldy == yidx || yidx < 0 || yidx >= LPUnumInAxis)
          continue;

        int oldz = INT32_MIN;
        for (size_t k = 0; k < th_dists.size(); k++) {
          int zidx = (point[2] + th_dists[k]) >> log2MotionBlockSize;

          if (oldz != zidx && zidx >= 0 && zidx < LPUnumInAxis) {
            int idx = (xidx * LPUnumInAxis + yidx) * LPUnumInAxis + zidx;
            lpuActiveWindow[idx].push_back(point);
          }

          oldz = zidx;
        }
        oldy = yidx;
      }
      oldx = xidx;
    }
  } // end loop on points

  predPointCloud.resize(numPoints);
}
//----------------------------------------------------------------------------
bool
motionSearchForNode(
  const PCCPointSet3& pointCloud,
  const PCCOctree3Node* node0,
  const GeometryParameterSet::Motion& param,
  int nodeSizeLog2,
  EntropyEncoder* arithmeticEncoder,
  int8_t* bufferPoints,
  PUtree* local_PU_tree,
  const std::vector<std::vector<Vec3<int>>>& lpuActiveWindow,
  int lpuIdx)
{
  // if window is empty, no compensation
  if (!lpuActiveWindow[lpuIdx].size())
    return false;

  MotionEntropyEncoder motionEncoder(arithmeticEncoder);

  std::vector<Vec3<int>> Block0;
  Block0.reserve(node0->end - node0->start);
  for (size_t i = node0->start; i < node0->end; ++i)
    Block0.push_back(pointCloud[i]);

  // entropy estimates
  MotionEntropyEstimate mcEstimate(motionEncoder);
  double LUT_MVestimate[64];
  int wSize = param.motion_window_size;
  mcEstimate.prepareEstimateVector(LUT_MVestimate, wSize, param.motion_max_prefix_bits, param.motion_max_suffix_bits);

  // motion search
  Vec3<int32_t> pos = node0->pos << nodeSizeLog2;
  //PUtree local_PU_tree;
  int x0 = pos[0];
  int y0 = pos[1];
  int z0 = pos[2];

  // MV search
  find_motion(
    param, mcEstimate, lpuActiveWindow[lpuIdx], Block0, x0, y0, z0,
    (1 << nodeSizeLog2), local_PU_tree, bufferPoints, LUT_MVestimate);

  return true;
}

void
encode_splitPU_MV_MC(
  PCCOctree3Node* node0,
  PUtree* local_PU_tree,
  const GeometryParameterSet::Motion& param,
  Vec3<int> nodeSizeLog2,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3* compensatedPointCloud,
  std::vector<std::vector<Vec3<int>>>& lpuActiveWindow,
  int numLPUPerLine,
  int log2MotionBlkSize)
{
  const int node_size = 1 << nodeSizeLog2[0];
  MotionEntropyEncoder motionEncoder(arithmeticEncoder);

  // --------------  non-split / terminal case  ----------------
  if (node_size <= param.motion_min_pu_size || !local_PU_tree->split_flags[0]) {
    if (node_size > param.motion_min_pu_size) {
      motionEncoder.encodeSplitPu(0);
    }

    // encode MV
    Vec3<int> MV = local_PU_tree->MVs[0];
    motionEncoder.encodeVector(MV);
    Vec3<int> MVd = MV;

    // find search window
    std::vector<Vec3<int>> Window;
    const Vec3<int32_t> pos = node0->pos << nodeSizeLog2;
    const int lpuX = pos[0] >> log2MotionBlkSize;
    const int lpuY = pos[1] >> log2MotionBlkSize;
    const int lpuZ = pos[2] >> log2MotionBlkSize;
    const int lpuIdx = (lpuX * numLPUPerLine + lpuY) * numLPUPerLine + lpuZ;

    //  use only one ref
    Window.reserve(lpuActiveWindow[lpuIdx].size());
    for (size_t i = 0; i < lpuActiveWindow[lpuIdx].size(); ++i) {
      Window.push_back(lpuActiveWindow[lpuIdx][i]);
    }

    std::vector<Vec3<int>> pointPredictorMC;
    // create the compensated points
    int xlow = pos[0];
    int xhigh = pos[0] + (1 << nodeSizeLog2[0]);
    int ylow = pos[1];
    int yhigh = pos[1] + (1 << nodeSizeLog2[1]);
    int zlow = pos[2];
    int zhigh = pos[2] + (1 << nodeSizeLog2[2]);
    for (auto& w : Window) {
      // apply best motion
      const Vec3<int> wV = w - MVd;
      if (
        wV[0] >= xlow && wV[0] < xhigh && wV[1] >= ylow && wV[1] < yhigh
        && wV[2] >= zlow && wV[2] < zhigh)
        pointPredictorMC.push_back(wV);
    }

    //and make node0 point to them
    node0->predStart = compensatedPointCloud->getPointCount();
    compensatedPointCloud->resize(
      compensatedPointCloud->getPointCount() + pointPredictorMC.size());
    size_t counter = node0->predStart;
    for (const auto& p : pointPredictorMC) {
      auto& predPoint = (*compensatedPointCloud)[counter++];
      predPoint[0] = p[0];
      predPoint[1] = p[1];
      predPoint[2] = p[2];
    }
    node0->predEnd = compensatedPointCloud->getPointCount();
    node0->isCompensated = true;
    return;
  }

  // --------------- split case ----------------------
  motionEncoder.encodeSplitPu(1);
}

//----------------------------------------------------------------------------

void
decode_splitPU_MV_MC(
  PCCOctree3Node* node0,
  const GeometryParameterSet::Motion& param,
  Vec3<int> nodeSizeLog2,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3* compensatedPointCloud,
  std::vector<std::vector<Vec3<int>>>& lpuActiveWindow,
  int numLPUPerLine,
  int log2MotionBlkSize)
{
  int node_size = 1 << nodeSizeLog2[0];
  MotionEntropyDecoder motionDecoder(arithmeticDecoder);

  // decode split flag
  bool split = false;
  if (node_size > param.motion_min_pu_size)
    split = motionDecoder.decodeSplitPu();

  if (!split) {  // not split
                 // decode MV
    Vec3<int> MV = 0;
    Vec3<int> MVd = 0.;
    motionDecoder.decodeVector(&MV);
    MVd = MV;

    // find search window
    std::vector<Vec3<int>> Window;
    Vec3<int32_t> pos = node0->pos << nodeSizeLog2;

    const int lpuX = pos[0] >> log2MotionBlkSize;
    const int lpuY = pos[1] >> log2MotionBlkSize;
    const int lpuZ = pos[2] >> log2MotionBlkSize;
    const int lpuIdx = (lpuX * numLPUPerLine + lpuY) * numLPUPerLine + lpuZ;

    // use only one ref
    Window.reserve(lpuActiveWindow[lpuIdx].size());
    for (size_t i = 0; i < lpuActiveWindow[lpuIdx].size(); ++i) {
      Window.push_back(lpuActiveWindow[lpuIdx][i]);
    }

    std::vector<Vec3<int>> pointPredictorMC;
    // create the compensated points
    const int xlow = pos[0];
    const int xhigh = pos[0] + (1 << nodeSizeLog2[0]);
    const int ylow = pos[1];
    const int yhigh = pos[1] + (1 << nodeSizeLog2[1]);
    const int zlow = pos[2];
    const int zhigh = pos[2] + (1 << nodeSizeLog2[2]);

    for (const auto& w : Window) {
      // apply best motion
      Vec3<int> wV = w - MVd;
      if (
        wV[0] >= xlow && wV[0] < xhigh && wV[1] >= ylow && wV[1] < yhigh
        && wV[2] >= zlow && wV[2] < zhigh)
        pointPredictorMC.push_back(wV);
    }

    //and make node0 point to them
    node0->predStart = compensatedPointCloud->getPointCount();
    compensatedPointCloud->resize(
      compensatedPointCloud->getPointCount() + pointPredictorMC.size());
    size_t counter = node0->predStart;
    for (const auto& p : pointPredictorMC) {
      auto& predPoint = (*compensatedPointCloud)[counter++];
      predPoint[0] = p[0];
      predPoint[1] = p[1];
      predPoint[2] = p[2];
    }
    node0->predEnd = compensatedPointCloud->getPointCount();
    node0->isCompensated = true;
    return;
  }

  // split; nothing to do
}


//----------------------------------------------------------------------------
}  // namespace pcc
