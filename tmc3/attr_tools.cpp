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

#include "attr_tools.h"
#include "RAHT.h"

#include <cassert>
#include <cinttypes>
#include <climits>
#include <cstddef>
#include <utility>
#include <vector>
#include <stdio.h>

#include "PCCTMC3Common.h"
#include "PCCMisc.h"

namespace pcc {

//============================================================================
// Encapsulation of a RAHT transform stage.

RahtKernel::RahtKernel(int weightLeft, int weightRight)
{
  uint64_t w = weightLeft + weightRight;
  uint64_t isqrtW = irsqrt(w);
  _a.val = (isqrt(uint64_t(weightLeft) << (2 * _a.kFracBits)) * isqrtW) >> 40;
  _b.val = (isqrt(uint64_t(weightRight) << (2 * _b.kFracBits)) * isqrtW) >> 40;
}

void
RahtKernel::fwdTransform(
  FixedPoint& lf, FixedPoint& hf)
{
  // lf = right * b + left * a
  // hf = right * a - left * b
  auto tmp = lf.val * _b.val;
  lf.val *= _a.val;
  lf.val += hf.val * _b.val;
  hf.val *= _a.val;
  hf.val -= tmp;

  lf.fixAfterMultiplication();
  hf.fixAfterMultiplication();
}

void
RahtKernel::invTransform(
  FixedPoint& left, FixedPoint& right)
{
  // left = lf * a - hf * b
  // right = lf * b + hf * a
  auto tmp = right.val * _b.val;
  right.val *= _a.val;
  right.val += left.val * _b.val;
  left.val *= _a.val;
  left.val -= tmp;

  right.fixAfterMultiplication();
  left.fixAfterMultiplication();
}

//============================================================================
// In-place transform a set of sparse 2x2x2 blocks each using the same weights

void
fwdTransformBlock222(
  const int numBufs, VecAttr::iterator buf,
  const int64_t weights[])
{
  static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};
  static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
  for (int i = 0, iw = 0; i < 12; i++, iw += 2) {
    int i0 = a[i];
    int i1 = b[i];
    int64_t w0 = weights[iw + 32];
    int64_t w1 = weights[iw + 33];

    if (!w0 && !w1)
      continue;

    // only one occupied, propagate to next level
    if (!w0 || !w1) {
      if (!w0) {
        for (int k = 0; k < numBufs; k++)
          std::swap(buf[k][i0], buf[k][i1]);
      }
      continue;
    }

    // actual transform
    for (int k = 0; k < numBufs; k++) {
      auto& lf = buf[k][i0];
      auto& hf = buf[k][i1];

      auto tmp = lf.val * w1;
      lf.val = lf.val * w0 + hf.val * w1;
      hf.val = hf.val * w0 - tmp;

      lf.fixAfterMultiplication();
      hf.fixAfterMultiplication();
    }
  }
}

void
invTransformBlock222(
  const int numBufs, VecAttr::iterator buf,
  const int64_t weights[])
{
  static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};
  static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
  for (int i = 11, iw = 22; i >= 0; i--, iw -= 2) {
    int i0 = a[i];
    int i1 = b[i];
    int64_t w0 = weights[iw + 32];
    int64_t w1 = weights[iw + 33];

    if (!a && !b)
      continue;

    // only one occupied, propagate to next level
    if (!w0 || !w1) {
      if (!w0) {
        for (int k = 0; k < numBufs; k++)
          std::swap(buf[k][i0], buf[k][i1]);
      }
      continue;
    }

    // actual transform
    for (int k = 0; k < numBufs; k++) {
      auto& left = buf[k][i0];
      auto& right = buf[k][i1];

      auto tmp = right.val * w1;
      right.val = right.val * w0 + left.val * w1;
      left.val = left.val * w0 - tmp;

      right.fixAfterMultiplication();
      left.fixAfterMultiplication();
    }
  }
}
}  // namespace pcc
