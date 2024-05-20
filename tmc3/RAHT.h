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
#include "attr_tools.h"
#include "ply.h"
#include <vector>
#include "pointset_processing.h"
#include "ringbuf.h"

namespace pcc {

void regionAdaptiveHierarchicalTransform(
  const RahtPredictionParams& rahtPredParams,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  const int attribCount,
  const int voxelCount,
  int64_t* mortonCode,
  int* attributes,
  const int voxelCount_mc,
  int64_t* mortonCode_mc,
  int* attributes_mc,
  int* coefficients,
  attr::ModeEncoder& encoder);

void regionAdaptiveHierarchicalInverseTransform(
  const RahtPredictionParams& rahtPredParams,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  const int attribCount,
  const int voxelCount,
  int64_t* mortonCode,
  int* attributes,
  const int voxelCount_mc,
  int64_t* mortonCode_mc,
  int* attributes_mc,
  int* coefficients,
  attr::ModeDecoder& decoder);

struct PCCRAHTACCoefficientEntropyEstimate {
  PCCRAHTACCoefficientEntropyEstimate()
  {
    init();
  }

  PCCRAHTACCoefficientEntropyEstimate(
    const PCCRAHTACCoefficientEntropyEstimate& other) = default;

  PCCRAHTACCoefficientEntropyEstimate&
    operator=(const PCCRAHTACCoefficientEntropyEstimate&) = default;

  void resStatUpdate(int32_t values, int k);
  void init();
  void updateCostBits(int32_t values, int k);
  double costBits() { return sumCostBits; }
  void resetCostBits() { sumCostBits = 0.; }

private:
  // Encoder side residual cost calculation
  static constexpr unsigned scaleRes = 1 << 20;
  static constexpr unsigned windowLog2 = 6;
  int probResGt0[3];  //prob of residuals larger than 0: 1 for each component
  int probResGt1[3];  //prob of residuals larger than 1: 1 for each component
  double sumCostBits;
};

struct PCCRAHTComputeCCCP {
  int8_t computeCrossChromaComponentPredictionCoeff(int m, int64_t coeffs[][3]);

  void reset()
  {
    sum = {0, 0};
    window.clear();
  }

private:
  struct Elt {
    int64_t k1k2;
    int64_t k1k1;
  };

  Elt sum {0, 0};
  ringbuf<Elt> window = ringbuf<Elt>(128);
};
} /* namespace pcc */
