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

  // global
  void encodeIsWorld(bool hasMotion);

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

  //global
  bool decodeIsWorld();

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

//----------------------------------------------------------------------------
// Encode the presence of an horizontal segment
inline void
MotionEntropyEncoder::encodeIsWorld(bool present)
{
  _arithmeticEncoder->encode(present, _ctxIsWorld);
}

//----------------------------------------------------------------------------
// Decode horizontal segment presence flag
inline bool
MotionEntropyDecoder::decodeIsWorld()
{
  return _arithmeticDecoder->decode(_ctxIsWorld);
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
static double
calcCostOfGlobalMotion(
  std::vector<Vec3<int>>& Window,
  std::vector<Vec3<int>>& Block0,
  int* bufferPoints,
  int wSize)
{
  if (!Window.size())
    return DBL_MAX;

  double cost = DBL_MAX;
  const int samples = 4;
  const int decimate = 6;

  if (Window.size() > samples * std::max(int(Block0.size()), 16))
    wSize >>= 1;
  int maxDistance = wSize << 1;

  int dist = 0;
  int* pBuffer = bufferPoints;

  auto* itB = Block0.data();
  int jumpBlock = 1
    + (Block0.size()
       >> decimate);  // (kind of) random sampling of the original block to code
  int nSamples = 0;
  for (int Nb = 0; Nb < Block0.size();
       Nb += jumpBlock, itB += jumpBlock, nSamples++) {
    Vec3<int> b = *itB;
    int min_d = maxDistance;
    for (const auto w : Window) {
      int a[3];
      for (auto i = 0; i < 3; i++) {
        pBuffer[i] = b[i] - w[i];
        a[i] = std::abs(pBuffer[i]);
      }
      if (a[0] <= wSize && a[1] <= wSize && a[2] <= wSize)
        pBuffer += 3;

      a[0] += a[1] + a[2];  // replaces int local_d
      if (a[0] < min_d)
        min_d = a[0];
    }
    dist += plus1log2shifted4(min_d);  // 1/0.0625 = 16 times log
  }
  cost = jumpBlock * dist;
  return cost;
}

//============================================================================

void
populateWindowList(
  VecVecVec3& windowList,
  PCCPointSet3& pointPred,
  std::vector<int> motion_block_size,
  int LPUnumInAxis[3],
  const Vec3<int> min)
{
  for (size_t i = 0; i < pointPred.getPointCount(); ++i) {
    const Vec3<int> point = pointPred[i];
    const int xidx =
      motion_block_size[0] ? (point[0] - min.x()) / motion_block_size[0] : 0;
    if (xidx < 0 || xidx >= LPUnumInAxis[0])
      continue;
    const int yidx =
      motion_block_size[1] ? (point[1] - min.y()) / motion_block_size[1] : 0;
    if (yidx < 0 || yidx >= LPUnumInAxis[1])
      continue;
    const int zidx =
      motion_block_size[2] ? (point[2] - min.z()) / motion_block_size[2] : 0;
    if (zidx < 0 || zidx >= LPUnumInAxis[2])
      continue;
    const int idx = (xidx * LPUnumInAxis[1] + yidx) * LPUnumInAxis[2] + zidx;
    windowList[idx].push_back(point);
  }
}

void
compensateCuboidGlobalMotion(
  PCCPointSet3* compensatedPointCloud,
  PCCPointSet3& pointPredictor,
  PCCPointSet3& pointPredictorWorld,
  std::vector<bool>& isWorldList,
  std::vector<int> motion_block_size,
  int LPUnumInAxis[3],
  const Vec3<int> min)
{
  VecVecVec3 Window_W_List, Window_V_List;
  const int blocksize = LPUnumInAxis[0] * LPUnumInAxis[1] * LPUnumInAxis[2];
  Window_W_List.resize(blocksize);
  Window_V_List.resize(blocksize);
  populateWindowList(
    Window_W_List, pointPredictorWorld, motion_block_size, LPUnumInAxis, min);
  populateWindowList(
    Window_V_List, pointPredictor, motion_block_size, LPUnumInAxis, min);

  std::vector<Vec3<double>> pointPredictorMC;
  for (int idx = 0; idx < blocksize; idx++) {
    const auto& windowList =
      isWorldList[idx] ? Window_W_List[idx] : Window_V_List[idx];
    for (const auto& point : windowList)
      pointPredictorMC.push_back(point);
  }

  compensatedPointCloud->resize(pointPredictorMC.size());
  int counter = 0;
  for (const auto& p : pointPredictorMC) {
    auto& predPoint = (*compensatedPointCloud)[counter++];
    predPoint[0] = p[0];
    predPoint[1] = p[1];
    predPoint[2] = p[2];
  }
}
void
populateCuboidBlocks(
  VecVecVec3& windowList,
  const PCCPointSet3& pointCloud,
  std::vector<int> motion_block_size,
  std::vector<int> th_dists,
  Box3<int32_t> bbox,
  int LPUnumInAxis[3])
{
  for (size_t i = 0; i < pointCloud.getPointCount(); ++i) {
    std::unordered_map<int, bool> lpuToAdd;
    const Vec3<int> point = pointCloud[i];
    for (size_t m = 0; m < th_dists.size(); m++) {
      const int xidx = motion_block_size[0]
        ? (point[0] + th_dists[m] - bbox.min.x()) / motion_block_size[0]
        : 0;
      if (xidx < 0 || xidx >= LPUnumInAxis[0])
        continue;
      for (size_t n = 0; n < th_dists.size(); n++) {
        const int yidx = motion_block_size[1]
          ? (point[1] + th_dists[n] - bbox.min.y()) / motion_block_size[1]
          : 0;
        if (yidx < 0 || yidx >= LPUnumInAxis[1])
          continue;
        for (size_t k = 0; k < th_dists.size(); k++) {
          const int zidx = motion_block_size[2]
            ? (point[2] + th_dists[k] - bbox.min.z()) / motion_block_size[2]
            : 0;
          if (zidx < 0 || zidx >= LPUnumInAxis[2])
            continue;
          const int idx =
            (xidx * LPUnumInAxis[1] + yidx) * LPUnumInAxis[2] + zidx;
          lpuToAdd[idx] = true;
        }
      }
    }
    for (auto idx : lpuToAdd)
      windowList[idx.first].push_back(point);
  }
}
void
encodeCuboidGlobalMotion(
  const PCCPointSet3& pointCloud,
  std::vector<int> motion_block_size,
  int motion_window_size,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3& pointPredictor,
  int* bufferPoints,
  PCCPointSet3* compensatedPointCloud,
  PCCPointSet3& pointPredictorWorld)
{
  MotionEntropyEncoder motionEncoder(arithmeticEncoder);

  const auto bbox = pointPredictor.computeBoundingBox();

  int LPUnumInAxis[3] = {1, 1, 1};
  for (int i = 0; i < 3; i++) {
    if (motion_block_size[i])
      LPUnumInAxis[i] = (bbox.max[i] - bbox.min[i] + motion_block_size[i] - 1)
        / motion_block_size[i];
  }
  const int extendedWindow = motion_window_size;
  std::vector<int> th_dists;
  th_dists.push_back(extendedWindow);
  if (motion_window_size)
    th_dists.push_back(-extendedWindow);

  const int blockSize = LPUnumInAxis[0] * LPUnumInAxis[1] * LPUnumInAxis[2];

  VecVecVec3 Block0_List, Window_W_List, Window_V_List;
  Block0_List.resize(blockSize);
  Window_W_List.resize(blockSize);
  Window_V_List.resize(blockSize);
  populateCuboidBlocks(
    Block0_List, pointCloud, motion_block_size, th_dists, bbox, LPUnumInAxis);
  populateCuboidBlocks(
    Window_W_List, pointPredictorWorld, motion_block_size, th_dists, bbox,
    LPUnumInAxis);
  populateCuboidBlocks(
    Window_V_List, pointPredictor, motion_block_size, th_dists, bbox,
    LPUnumInAxis);

  std::vector<bool> isWorldList(blockSize, true);
  for (int ith = 0; ith < blockSize; ith++) {
    if (
      !Block0_List[ith].size()
      || (!Window_W_List[ith].size() && !Window_V_List[ith].size())) {
      continue;
    }

    double costWorld = calcCostOfGlobalMotion(
      Window_W_List[ith], Block0_List[ith], bufferPoints, motion_window_size);

    double costVehicle = calcCostOfGlobalMotion(
      Window_V_List[ith], Block0_List[ith], bufferPoints, motion_window_size);

    if (!Window_W_List[ith].size() || costWorld >= costVehicle)
      isWorldList[ith] = false;
  }

  Window_V_List.clear();
  Window_W_List.clear();
  Block0_List.clear();

  // encoding
  for (int idx = 0; idx < blockSize; idx++)
    motionEncoder.encodeIsWorld(isWorldList[idx]);

  // compensation
  compensateCuboidGlobalMotion(
    compensatedPointCloud, pointPredictor, pointPredictorWorld, isWorldList,
    motion_block_size, LPUnumInAxis, bbox.min);
}

//----------------------------------------------------------------------------
void
decodeCuboidGlobalMotion(
  std::vector<int> motion_block_size,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3& pointPredictor,
  PCCPointSet3* compensatedPointCloud,
  PCCPointSet3& pointPredictorWorld)
{
  MotionEntropyDecoder motionDecoder(arithmeticDecoder);

  const auto bbox = pointPredictor.computeBoundingBox();
  int LPUnumInAxis[3] = {1, 1, 1};
  for (int i = 0; i < 3; i++) {
    if (motion_block_size[i])
      LPUnumInAxis[i] = (bbox.max[i] - bbox.min[i] + motion_block_size[i] - 1)
        / motion_block_size[i];
  }

  const int blockSize = LPUnumInAxis[0] * LPUnumInAxis[1] * LPUnumInAxis[2];
  std::vector<bool> isWorldList(false, blockSize);

  // decoding
  for (int idx = 0; idx < blockSize; idx++)
    isWorldList[idx] = motionDecoder.decodeIsWorld();

  // compensation
  compensateCuboidGlobalMotion(
    compensatedPointCloud, pointPredictor, pointPredictorWorld, isWorldList,
    motion_block_size, LPUnumInAxis, bbox.min);
}

//============================================================================
void
quantizeGlobalMotion(double Mat_GM[4][3], int32_t Mat_GM_Q[4][3])
{
  double scale = motionParamScale;
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++) {
      if (l == c)
        Mat_GM_Q[l][c] =
          roundIntegerHalfInf((Mat_GM[l][c] - 1.) * scale) + scale;
      else if (l < 3)
        Mat_GM_Q[l][c] = roundIntegerHalfInf(Mat_GM[l][c] * scale);
      else
        Mat_GM_Q[l][c] = roundIntegerHalfInf(Mat_GM[l][c]);
    }
}

//----------------------------------------------------------------------------
void
applyGlobalMotion(std::vector<Vec3<int>>& listPoints, double Mat_GM[4][3])
{
  for (auto& b : listPoints) {
    Vec3<int> point;
    for (auto i = 0; i < 3; i++)
      point[i] = roundIntegerHalfInf(
        Mat_GM[0][i] * b[0] + Mat_GM[1][i] * b[1] + Mat_GM[2][i] * b[2]
        + Mat_GM[3][i]);
    b = point;
  }
}

//----------------------------------------------------------------------------
void
applyGlobalMotion(
  PCCPointSet3& PC,
  const int32_t Mat_GM_Q[4][3],
  Vec3<double> vehicle_position,
  const int32_t global_thresholds[2])
{
  // unquantize
  double Mat_GM[4][3];
  double scale = motionParamScale;
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++) {
      if (l == c)
        Mat_GM[l][c] = double(Mat_GM_Q[l][c]) / scale + 1.;
      else if (l < 3)
        Mat_GM[l][c] = double(Mat_GM_Q[l][c]) / scale;
      else
        Mat_GM[l][c] = double(Mat_GM_Q[l][c]);
    }

  int topZ = global_thresholds[0];
  int bottomZ = global_thresholds[1];

  // apply motion
  for (int n = 0; n < PC.getPointCount(); n++) {
    Vec3<double> b = PC[n] + vehicle_position;

    if (b[2] < bottomZ || b[2] > topZ) {
      for (auto i = 0; i < 3; i++)
        PC[n][i] = roundIntegerHalfInf(
          Mat_GM[0][i] * b[0] + Mat_GM[1][i] * b[1] + Mat_GM[2][i] * b[2]
          + Mat_GM[3][i] - vehicle_position[i]);
    }
  }
}

//----------------------------------------------------------------------------
int
roundIntegerHalfInf(const double x)
{
  return (x >= 0) ? int(x + 0.5) : -int(-x + 0.5);
}

//----------------------------------------------------------------------------
double
map_reference(
  VecVec3& pcWorldTarget,
  const VecVec3& pointPredictorCentered,
  VecVec3& pcWorldRef)
{
  std::vector<int> accu_m;
  int meanM = 0;
  for (const auto& b : pcWorldTarget) {
    int dmin = 1 << 30;
    Vec3<int> closest;
    for (const auto& w : pointPredictorCentered) {
      const int L =
        std::abs(w[0] - b[0]) + std::abs(w[1] - b[1]) + std::abs(w[2] - b[2]);
      if (L < dmin) {
        dmin = L;
        closest = w;
      }
    }
    pcWorldRef.push_back(closest);

    accu_m.push_back(dmin);
    meanM += dmin;
  }

  double err = double(meanM) / pcWorldTarget.size() / 3.;

  // eliminate outliers
  int count = 0;
  auto it_accu = accu_m.begin();
  auto it_target = pcWorldTarget.begin();
  auto it_ref = pcWorldRef.begin();
  auto it_target_fill = pcWorldTarget.begin();
  auto it_ref_fill = pcWorldRef.begin();
  for (; it_accu != accu_m.end(); it_accu++, it_target++, it_ref++) {
    if (*it_accu * accu_m.size() <= 2 * meanM) {
      *it_target_fill++ = *it_target;
      *it_ref_fill++ = *it_ref;
      count++;
    }
  }

  pcWorldTarget.resize(count);
  pcWorldRef.resize(count);

  return err;
}

//----------------------------------------------------------------------------
void
LMS3D(
  std::vector<Vec3<int>>& P1,
  std::vector<Vec3<int>>& P2,
  std::vector<Vec3<int>>& pointPredictor_centered,
  uint32_t maxBB,
  double Mat_GM[4][3])
{
  // determine correlation matrix M in (X,Y,Z,MV_unity)
  const int MV_unity = maxBB >> 4;  //  // for better matrix conditioning
  double M[4][4] = {
    {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}};

  for (auto it1 = P1.begin(); it1 != P1.end(); it1++) {
    Vec3<double> Pref = *it1;
    //{(*it1)[0], (*it1)[1], (*it1)[2]};

    // X
    M[0][0] += Pref[0] * Pref[0];
    M[0][1] += Pref[0] * Pref[1];
    M[0][2] += Pref[0] * Pref[2];
    M[0][3] += Pref[0] * MV_unity;
    // Y
    M[1][1] += Pref[1] * Pref[1];
    M[1][2] += Pref[1] * Pref[2];
    M[1][3] += Pref[1] * MV_unity;
    // Z
    M[2][2] += Pref[2] * Pref[2];
    M[2][3] += Pref[2] * MV_unity;
    // 1
    M[3][3] += MV_unity * MV_unity;
  }
  M[1][0] = M[0][1];
  M[2][0] = M[0][2];
  M[2][1] = M[1][2];
  M[3][0] = M[0][3];
  M[3][1] = M[1][3];
  M[3][2] = M[2][3];

  // inverse M by Gauss pivoting
  double invM[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  for (int pivot = 0; pivot < 3; pivot++)  // descente
  {
    double value_pivot = M[pivot][pivot];
    for (int l = pivot + 1; l < 4; l++) {
      double factor = -M[l][pivot] / value_pivot;
      for (int c = 0; c < 4; c++) {
        M[l][c] += M[pivot][c] * factor;
        invM[l][c] += invM[pivot][c] * factor;
      }
    }
  }

  for (int pivot = 3; pivot > 0; pivot--)  // mont�e
  {
    double value_pivot = M[pivot][pivot];
    for (int l = pivot - 1; l >= 0; l--) {
      double factor = -M[l][pivot] / value_pivot;
      for (int c = 0; c < 4; c++) {
        M[l][c] += M[pivot][c] * factor;
        invM[l][c] += invM[pivot][c] * factor;
      }
    }
  }

  for (int pivot = 0; pivot < 4; pivot++)  // normalisation
  {
    double factor = 1 / M[pivot][pivot];
    for (int c = 0; c < 4; c++)
      invM[pivot][c] *= factor;
  }

  // determine rhs matrix R
  double R[4][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  for (auto it1 = P1.begin(), it2 = P2.begin(); it1 != P1.end();
       it1++, it2++) {
    Vec3<double> Pref = *it1;
    Vec3<double> Ptarget = *it2;

    // X
    R[0][0] += Ptarget[0] * Pref[0];
    R[1][0] += Ptarget[0] * Pref[1];
    R[2][0] += Ptarget[0] * Pref[2];
    R[3][0] += Ptarget[0] * MV_unity;
    // Y
    R[0][1] += Ptarget[1] * Pref[0];
    R[1][1] += Ptarget[1] * Pref[1];
    R[2][1] += Ptarget[1] * Pref[2];
    R[3][1] += Ptarget[1] * MV_unity;
    // Z
    R[0][2] += Ptarget[2] * Pref[0];
    R[1][2] += Ptarget[2] * Pref[1];
    R[2][2] += Ptarget[2] * Pref[2];
    R[3][2] += Ptarget[2] * MV_unity;
  }

  // apply inv M to R to get the transformation matrix T
  double T[4][3];
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++)
      T[l][c] = invM[l][0] * R[0][c] + invM[l][1] * R[1][c]
        + invM[l][2] * R[2][c] + invM[l][3] * R[3][c];

  // deconditioning of 1 <-> MV_unity
  for (int c = 0; c < 3; c++)
    T[3][c] *= double(MV_unity);

  // penalization
  double lambda = 1.0;
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++)
      T[l][c] *= lambda;
  T[0][0] += 1 - lambda;
  T[1][1] += 1 - lambda;
  T[2][2] += 1 - lambda;

  // apply T to global motion matrix  Mat_GM
  double Mat_GM1[4][3];  //copy old GM matrix
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++)
      Mat_GM1[l][c] = Mat_GM[l][c];

  for (int l = 0; l < 3; l++)  // deformation part
    for (int c = 0; c < 3; c++)
      Mat_GM[l][c] = Mat_GM1[l][0] * T[0][c] + Mat_GM1[l][1] * T[1][c]
        + Mat_GM1[l][2] * T[2][c];

  for (int c = 0; c < 3; c++)  // translation part
    Mat_GM[3][c] = Mat_GM1[3][0] * T[0][c] + Mat_GM1[3][1] * T[1][c]
      + Mat_GM1[3][2] * T[2][c] + T[3][c];
}

//----------------------------------------------------------------------------
void
PopulatePCLikelyWorld(
  const int blocknumInAxis,
  const int th_dist,
  const int bsize,
  const int top_z_boundary,
  const int bottom_z_boundary,
  const bool useCuboidalRegionsInGMEstimation,
  const PCCPointSet3& pointCloud,
  const PCCPointSet3& pointPredictor,
  std::vector<Vec3<int>>& pc_likely_world)
{
  int th_dists[2] = {th_dist, -th_dist};
  std::vector<bool> ref_planar_region;

  if (!useCuboidalRegionsInGMEstimation) {
    ref_planar_region.resize(
      blocknumInAxis * blocknumInAxis * blocknumInAxis, false);

    // scan each point in reference frame
    // use cubic block as the process unit
    for (size_t i = 0; i < pointPredictor.getPointCount(); i++) {
      const Vec3<double> point = pointPredictor[i];
      double x = point[0], y = point[1], z = point[2];
      for (size_t m = 0; m < 2; m++) {
        int xidx = (x + th_dists[m]) / bsize;
        if (xidx < 0 || xidx >= blocknumInAxis)
          continue;
        for (size_t n = 0; n < 2; n++) {
          int yidx = (y + th_dists[n]) / bsize;
          if (yidx < 0 || yidx >= blocknumInAxis)
            continue;
          for (size_t k = 0; k < 2; k++) {
            int zidx = (z + th_dists[k]) / bsize;
            if (zidx < 0 || zidx >= blocknumInAxis)
              continue;
            int idx = (xidx * blocknumInAxis + yidx) * blocknumInAxis + zidx;
            ref_planar_region[idx] = true;
          }
        }
      }
    }
    //scan each point in current frame
    for (size_t i = 0; i < pointCloud.getPointCount(); i++) {
      const Vec3<double> point = pointCloud[i];
      int xidx = point[0] / bsize;
      int yidx = point[1] / bsize;
      int zidx = point[2] / bsize;
      int idx = (xidx * blocknumInAxis + yidx) * blocknumInAxis + zidx;
      if (idx >= ref_planar_region.size() || !ref_planar_region[idx])
        continue;
      const Vec3<int> b = {int(point[0]), int(point[1]), int(point[2])};
      if ((b[2] < bottom_z_boundary || b[2] > top_z_boundary)) {
        pc_likely_world.push_back(b);
      }
    }
  } else {  // Use cuboidal regions
    ref_planar_region.resize(blocknumInAxis * blocknumInAxis, false);

    // scan each point in reference frame
    // use cubic block as the process unit
    for (size_t i = 0; i < pointPredictor.getPointCount(); i++) {
      const Vec3<double> point = pointPredictor[i];
      double x = point[0], y = point[1], z = point[2];
      for (size_t m = 0; m < 2; m++) {
        int xidx = (x + th_dists[m]) / bsize;
        if (xidx < 0 || xidx >= blocknumInAxis)
          continue;
        for (size_t n = 0; n < 2; n++) {
          int yidx = (y + th_dists[n]) / bsize;
          if (yidx < 0 || yidx >= blocknumInAxis)
            continue;
          int idx = xidx * blocknumInAxis + yidx;
          ref_planar_region[idx] = true;
        }
      }
    }
    //scan each point in current frame
    for (size_t i = 0; i < pointCloud.getPointCount(); i++) {
      const Vec3<double> point = pointCloud[i];
      int xidx = point[0] / bsize;
      int yidx = point[1] / bsize;
      int idx = xidx * blocknumInAxis + yidx;
      if (idx >= ref_planar_region.size() || !ref_planar_region[idx])
        continue;
      const Vec3<int> b = {int(point[0]), int(point[1]), int(point[2])};
      if ((b[2] < bottom_z_boundary || b[2] > top_z_boundary)) {
        pc_likely_world.push_back(b);
      }
    }
  }
}

//----------------------------------------------------------------------------
void
SearchGlobalMotion(
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointPredictor,
  double QS,
  int bsize,
  int th_dist,
  uint32_t maxBB,
  const bool useCuboidalRegionsInGMEstimation,
  std::vector<int>& gm_matrix,
  Vec3<int>& gm_trans,
  const std::pair<int, int> thresh)
{
  // ------------- first pass: find world-referential-likely LPU ----
  // number of LCU
  uint32_t maxBB_Scalled = (maxBB);

  // loop on LCU
  std::vector<Vec3<int>> pcLikelyWorld;

  int topZ = thresh.first;
  int bottomZ = thresh.second;
  int blocknumInAxis = (maxBB_Scalled % bsize) ? (maxBB_Scalled / bsize + 1)
                                               : (maxBB_Scalled / bsize);

  PopulatePCLikelyWorld(
    blocknumInAxis, th_dist, bsize, topZ, bottomZ,
    useCuboidalRegionsInGMEstimation, pointCloud, pointPredictor,
    pcLikelyWorld);

  Vec3<double> minPositions = {0., 0., 0.};
  Vec3<double> minPositions_local = minPositions * QS;
  Vec3<int> minPositions_int = Vec3<int>(
    int(minPositions_local[0]), int(minPositions_local[1]),
    int(minPositions_local[2]));
  for (auto& b : pcLikelyWorld)
    b += minPositions_int;

  std::vector<Vec3<int>> pointPredictorCentered;
  for (size_t i = 0; i < pointPredictor.getPointCount(); ++i) {
    const Vec3<double> point = pointPredictor[i] + minPositions_local;
    ;
    pointPredictorCentered.push_back(
      Vec3<int>(int(point[0]), int(point[1]), int(point[2])));
  }
  std::vector<Vec3<int>> pointPredictorCentered0 = pointPredictorCentered;

  // global motion found iteratively using LMS

  int NLMS = 1;
  int nb_points = 100;

  int jump = 1 + (pcLikelyWorld.size() / nb_points);

  double Mat_GM[4][3] = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1},
    {0, 0, 0}};  // global motion is identity as a start

  //std::cout << "(points/err) = ";

  for (int i = 0; i < NLMS; i++) {
    // sample pc_likely_world
    std::vector<Vec3<int>> pcWorldTarget;
    for (int N = i; N < pcLikelyWorld.size(); N += jump) {
      auto it = pcLikelyWorld.begin() + N;
      pcWorldTarget.push_back(*it);
    }

    // map reference to pc_world
    std::vector<Vec3<int>> pcWorldRef;
    map_reference(pcWorldTarget, pointPredictorCentered, pcWorldRef);
    // Least Mean Square 3D
    LMS3D(
      pcWorldRef, pcWorldTarget, pointPredictorCentered, maxBB_Scalled,
      Mat_GM);

    // apply global motion
    if (NLMS > 1) {  // Unnecessary when NLMS = 1
      pointPredictorCentered = pointPredictorCentered0;
      applyGlobalMotion(pointPredictorCentered, Mat_GM);
    }
  }
  //std::cout << std::endl;
  int Mat_GM_Q[4][3] = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1},
    {0, 0, 0}};  // global motion is identity as a start
  quantizeGlobalMotion(Mat_GM, Mat_GM_Q);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gm_matrix[3 * i + j] = Mat_GM_Q[j][i];
    }
    gm_trans[i] = Mat_GM_Q[3][i];
  }
}

//----------------------------------------------------------------------------
void
SearchGlobalMotionPerTile(
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointPredictor,
  double QS,
  GeometryBrickHeader& gbh,
  int th_dist,
  const bool useCuboidalRegionsInGMEstimation)
{
  int maxBB = (1 << gbh.maxRootNodeDimLog2) - 1;

  SearchGlobalMotion(
    pointCloud, pointPredictor, QS, gbh.motion_block_size[2], th_dist, maxBB,
    useCuboidalRegionsInGMEstimation, gbh.gm_matrix, gbh.gm_trans,
    gbh.gm_thresh);
}

//----------------------------------------------------------------------------
void
applyGlobalMotion_with_shift(
  PCCPointSet3& pointPredictorWorld,
  const std::vector<int> gm_matrix,
  const Vec3<int> gm_trans,
  const Vec3<int> minimum_position)
{
  static const unsigned int motionParamPrec = 16;
  static const unsigned int motionParamOffset = 1 << (motionParamPrec - 1);

  // apply
  Vec3<int> b = 0;
  for (int n = 0; n < pointPredictorWorld.getPointCount(); n++) {
    b = pointPredictorWorld[n];

    b[0] = b[0] + minimum_position[0];
    b[1] = b[1] + minimum_position[1];
    b[2] = b[2] + minimum_position[2];

    for (int i = 0; i < 3; i++) {
      pointPredictorWorld[n][i] =
        divExp2RoundHalfInfPositiveShift(
          gm_matrix[3 * i] * b[0] + gm_matrix[3 * i + 1] * b[1]
            + gm_matrix[3 * i + 2] * b[2],
          motionParamPrec, motionParamOffset)
        + gm_trans[i] - minimum_position[i];  // translation and offset
    }
  }
}

//----------------------------------------------------------------------------
void
compensateWithRoadObjClassfication(
  PCCPointSet3& pointPredictorWorld,
  const std::vector<int> gm_matrix,
  const Vec3<int> gm_trans,
  const std::pair<int, int> thresh,
  const Vec3<int> minimum_position)
{
  static const unsigned int motionParamPrec = 16;
  static const unsigned int motionParamOffset = 1 << (motionParamPrec - 1);

  // apply
  Vec3<int> b = 0;
  for (int n = 0; n < pointPredictorWorld.getPointCount(); n++) {
    b = pointPredictorWorld[n];

    b[0] = b[0] + minimum_position[0];
    b[1] = b[1] + minimum_position[1];
    b[2] = b[2] + minimum_position[2];

    if ((b[2] < thresh.second) || (b[2] > thresh.first)) {
      for (int i = 0; i < 3; i++) {
        pointPredictorWorld[n][i] =
          divExp2RoundHalfInfPositiveShift(
            gm_matrix[3 * i] * b[0] + gm_matrix[3 * i + 1] * b[1]
              + gm_matrix[3 * i + 2] * b[2],
            motionParamPrec, motionParamOffset)
          + gm_trans[i] - minimum_position[i];  // translation and offset
      }
    }
  }
}

//----------------------------------------------------------------------------
// output: pointPredictorWorld : predPointCloud + applied global motion
//         compensatedPointCloud : predPointCloud + Octree-based partition
void
compensateWithCuboidPartition(
  PCCPointSet3& pointCloud,
  PCCPointSet3& predPointCloud,
  PCCPointSet3& pointPredictorWorld,
  const GeometryBrickHeader& gbh,
  int motion_window_size,
  const Vec3<int> minimum_position,
  EntropyEncoder* arithmeticEncoder)
{
  applyGlobalMotion_with_shift(
    pointPredictorWorld, gbh.gm_matrix, gbh.gm_trans, minimum_position);

  std::unique_ptr<int> bufferPoints;
  bufferPoints.reset(new int[3 * 3000 * 10000]);
  PCCPointSet3 compensatedPointCloud;

  encodeCuboidGlobalMotion(
    pointCloud, gbh.motion_block_size, motion_window_size, arithmeticEncoder,
    predPointCloud, bufferPoints.get(), &compensatedPointCloud,
    pointPredictorWorld);

  bufferPoints.reset();
  pointPredictorWorld.clear();
  pointPredictorWorld = compensatedPointCloud;
}

//----------------------------------------------------------------------------
void
decodeCompensateWithCuboidPartition(
  PCCPointSet3& predPointCloud,
  PCCPointSet3& pointPredictorWorld,
  const GeometryBrickHeader& gbh,
  const Vec3<int> minimum_position,
  EntropyDecoder* arithmeticDecoder)
{
  applyGlobalMotion_with_shift(
    pointPredictorWorld, gbh.gm_matrix, gbh.gm_trans, minimum_position);

  PCCPointSet3 compensatedPointCloud;

  decodeCuboidGlobalMotion(
    gbh.motion_block_size, arithmeticDecoder, predPointCloud,
    &compensatedPointCloud, pointPredictorWorld);

  pointPredictorWorld = compensatedPointCloud;
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
  std::vector<Vec3<int>> list_tested = { {0,0,0,} };

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
      if (std::find(list_tested.begin(), list_tested.end(), V) != list_tested.end())
        continue;
      list_tested.push_back(V);

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

std::vector<std::vector<Vec3<int>>>
buildActiveWindow(
  const int LPUnumInAxis,
  const PCCPointSet3& predPointCloud,
  int motion_window_size,
  const int log2MotionBlockSize)
{
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

  std::vector<std::vector<Vec3<int>>> lpuActiveWindow;

  lpuActiveWindow.resize(LPUnumInAxis * LPUnumInAxis * LPUnumInAxis);
  for (size_t i = 0; i < predPointCloud.getPointCount(); ++i) {
    std::unordered_map<int, bool> lpuToAdd;
    const Vec3<int> point = predPointCloud[i];
    for (size_t m = 0; m < th_dists.size(); m++) {
      int xidx = (point[0] + th_dists[m]) >> log2MotionBlockSize;
      if (xidx < 0 || xidx >= LPUnumInAxis)
        continue;
      for (size_t n = 0; n < th_dists.size(); n++) {
        int yidx = (point[1] + th_dists[n]) >> log2MotionBlockSize;
        if (yidx < 0 || yidx >= LPUnumInAxis)
          continue;
        for (size_t k = 0; k < th_dists.size(); k++) {
          int zidx = (point[2] + th_dists[k]) >> log2MotionBlockSize;
          if (zidx < 0 || zidx >= LPUnumInAxis)
            continue;
          int idx = (xidx * LPUnumInAxis + yidx) * LPUnumInAxis + zidx;
          lpuToAdd[idx] = true;
        }
      }
    }
    for (auto idx : lpuToAdd) {
      lpuActiveWindow[idx.first].push_back(point);
    }
  }

  return lpuActiveWindow;
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
