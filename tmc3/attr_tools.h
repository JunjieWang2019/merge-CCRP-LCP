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

#include "FixedPoint.h"
#include "quantization.h"
#include "hls.h"
#include <vector>
#include "TMC3.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "PCCMisc.h"

using VecAttr = std::vector<std::array<pcc::FixedPoint, 8>>;

namespace pcc {
namespace attr {

enum Mode
{
  Null,
  Intra,
};

inline bool isNull(Mode mode) { return mode == Mode::Null; }
inline bool isIntra(Mode mode) { return mode == Mode::Intra; }

} // namespace attr
} /* namespace pcc */

namespace pcc {

struct MotionVector {
  Vec3<int32_t> position;
  int nodeSize;
  Vec3<int> motionVector;
};

//============================================================================

struct UrahtNode {
  int64_t pos;
  int weight;
  Qps qp;
  bool decoded;
  int64_t offset;
  attr::Mode mode;
  attr::Mode mc;

  uint8_t occupancy;
  std::vector<UrahtNode>::iterator firstChild;
  std::vector<UrahtNode>::iterator lastChild;
};

class RahtKernel {
public:
  RahtKernel(int weightLeft, int weightRight);
  void fwdTransform(FixedPoint& lf, FixedPoint& hf);
  void invTransform(FixedPoint& left, FixedPoint& right);

  int64_t getW0() { return _a.val; }
  int64_t getW1() { return _b.val; }

private:
  FixedPoint _a, _b;
};

//============================================================================
// In-place transform a set of sparse 2x2x2 blocks each using the same weights

void fwdTransformBlock222(
  const int numBufs, VecAttr::iterator buf, const int64_t weights[]);

void invTransformBlock222(
  const int numBufs, VecAttr::iterator buf, const int64_t weights[]);
} /* namespace pcc */
