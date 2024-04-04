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

#include "AttributeEncoder.h"

#include "attribute_raw.h"
#include "constants.h"
#include "entropy.h"
#include "io_hls.h"
#include "quantization.h"
#include "RAHT.h"
#include "FixedPoint.h"

#include <algorithm>

// todo(df): promote to per-attribute encoder parameter
//static const double kAttrPredLambdaR = 0.01;
//static const double kAttrPredLambdaC = 0.14;

namespace pcc {
//============================================================================
// An encapsulation of the entropy coding methods used in attribute coding

class PCCResidualsEncoder : protected AttributeContexts {
public:
  PCCResidualsEncoder(
    const AttributeParameterSet& aps,
    const AttributeBrickHeader& abh,
    const AttributeContexts& ctxtMem);

  EntropyEncoder arithmeticEncoder;

  const AttributeContexts& getCtx() const { return *this; }

  void start(const SequenceParameterSet& sps, int numPoints);
  int stop();

  void encodeRunLength(int runLength);
  void encodeInterPredMode(const int& predMode);
  void encodeSymbol(uint32_t value, int k1, int k2, int k3);
  void encode(int32_t value0, int32_t value1, int32_t value2);
  void encode(int32_t value);

  // Encoder side residual cost calculation
  const int scaleRes = 1 << 20;
  const int windowLog2 = 6;
  int probResGt0[3];  //prob of residuals larger than 0: 1 for each component
  int probResGt1[3];  //prob of residuals larger than 1: 1 for each component
  void resStatUpdateColor(Vec3<int32_t> values);
  void resStatUpdateRefl(int32_t values);
  void resStatReset();
};

//----------------------------------------------------------------------------

PCCResidualsEncoder::PCCResidualsEncoder(
  const AttributeParameterSet& aps,
  const AttributeBrickHeader& abh,
  const AttributeContexts& ctxtMem)
  : AttributeContexts(ctxtMem)
{
  resStatReset();
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::start(const SequenceParameterSet& sps, int pointCount)
{
  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = pointCount * 3 * 2 + 1024;
  arithmeticEncoder.setBuffer(maxAcBufLen, nullptr);
  arithmeticEncoder.enableBypassStream(sps.cabac_bypass_stream_enabled_flag);
  arithmeticEncoder.setBypassBinCodingWithoutProbUpdate(
    sps.bypass_bin_coding_without_prob_update);
  arithmeticEncoder.start();
}

//----------------------------------------------------------------------------

int
PCCResidualsEncoder::stop()
{
  return arithmeticEncoder.stop();
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::resStatReset()
{
  for (int k = 0; k < 3; k++)
    probResGt0[k] = probResGt1[k] = (scaleRes >> 1);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::resStatUpdateColor(Vec3<int32_t> value)
{
  for (int k = 0; k < 3; k++) {
    probResGt0[k] += value[k] ? (scaleRes - probResGt0[k]) >> windowLog2
                              : -((probResGt0[k]) >> windowLog2);
    if (value[k])
      probResGt1[k] += abs(value[k]) > 1
        ? (scaleRes - probResGt1[k]) >> windowLog2
        : -((probResGt1[k]) >> windowLog2);
  }
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::resStatUpdateRefl(int32_t value)
{
  probResGt0[0] += value ? (scaleRes - probResGt0[0]) >> windowLog2
                         : -(probResGt0[0] >> windowLog2);
  if (value)
    probResGt1[0] += abs(value) > 1 ? (scaleRes - probResGt1[0]) >> windowLog2
                                    : -(probResGt1[0] >> windowLog2);
}

//----------------------------------------------------------------------------

void PCCResidualsEncoder::encodeInterPredMode(const int& predMode) {
  const bool& isLayerMode = predMode >= 1;
  arithmeticEncoder.encode(isLayerMode, ctxLayerPred);
  if (isLayerMode)
    arithmeticEncoder.encode(predMode == 1, ctxInterLayerPred);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeRunLength(int runLength)
{
  auto* ctx = ctxRunLen;
  for (int i = 0; i < std::min(3, runLength); i++, ctx++)
    arithmeticEncoder.encode(1, *ctx);

  if (runLength < 3) {
    arithmeticEncoder.encode(0, *ctx);
    return;
  }
  runLength -= 3;

  auto prefix = runLength >> 1;
  for (int i = 0; i < std::min(4, prefix); i++)
    arithmeticEncoder.encode(1, *ctx);

  if (runLength < 8) {
    arithmeticEncoder.encode(0, *ctx);
    arithmeticEncoder.encode(runLength & 1);
    return;
  }
  runLength -= 8;

  arithmeticEncoder.encodeExpGolomb(runLength, 2, *++ctx);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeSymbol(uint32_t value, int k1, int k2, int k3)
{
  arithmeticEncoder.encode(value > 0, ctxCoeffGtN[0][k1]);
  if (!value)
    return;

  arithmeticEncoder.encode(--value > 0, ctxCoeffGtN[1][k2]);
  if (!value)
    return;

  arithmeticEncoder.encodeExpGolomb(
    --value, 1, ctxCoeffRemPrefix[k3], ctxCoeffRemSuffix[k3]);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encode(int32_t value0, int32_t value1, int32_t value2)
{
  int mag0 = abs(value0);
  int mag1 = abs(value1);
  int mag2 = abs(value2);

  int b0 = (mag1 == 0);
  int b1 = (mag1 <= 1);
  int b2 = (mag2 == 0);
  int b3 = (mag2 <= 1);
  encodeSymbol(mag1, 0, 0, 1);
  encodeSymbol(mag2, 1 + b0, 1 + b1, 1);

  auto mag0minusX = b0 && b2 ? mag0 - 1 : mag0;
  assert(mag0minusX >= 0);
  encodeSymbol(mag0minusX, 3 + (b0 << 1) + b2, 3 + (b1 << 1) + b3, 0);

  if (mag0)
    arithmeticEncoder.encode(value0 < 0);
  if (mag1)
    arithmeticEncoder.encode(value1 < 0);
  if (mag2)
    arithmeticEncoder.encode(value2 < 0);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encode(int32_t value)
{
  int mag = abs(value) - 1;
  encodeSymbol(mag, 0, 0, 0);
  arithmeticEncoder.encode(value < 0);
}

//============================================================================
// An encapsulation of the entropy coding methods used in attribute coding

struct PCCResidualsEntropyEstimator {
  size_t freq0[kAttributeResidualAlphabetSize + 1];
  size_t freq1[kAttributeResidualAlphabetSize + 1];
  size_t symbolCount0;
  size_t symbolCount1;
  size_t isZero0Count;
  size_t isZero1Count;
  PCCResidualsEntropyEstimator() { init(); }
  void init();
  double bitsDetail(
    const uint32_t detail,
    const size_t symbolCount,
    const size_t* const freq) const;
  double bits(const uint32_t value0) const;
  void update(const uint32_t value0);
  double bits(
    const uint32_t value0, const uint32_t value1, const uint32_t value2) const;
  void
  update(const uint32_t value0, const uint32_t value1, const uint32_t value2);
};

//----------------------------------------------------------------------------

void
PCCResidualsEntropyEstimator::init()
{
  for (size_t i = 0; i <= kAttributeResidualAlphabetSize; ++i) {
    freq0[i] = 1;
    freq1[i] = 1;
  }
  symbolCount0 = kAttributeResidualAlphabetSize + 1;
  symbolCount1 = kAttributeResidualAlphabetSize + 1;
  isZero1Count = isZero0Count = symbolCount0 / 2;
}

//----------------------------------------------------------------------------

double
PCCResidualsEntropyEstimator::bitsDetail(
  const uint32_t detail,
  const size_t symbolCount,
  const size_t* const freq) const
{
  const uint32_t detailClipped =
    std::min(detail, uint32_t(kAttributeResidualAlphabetSize));
  const double pDetail =
    PCCClip(double(freq[detailClipped]) / symbolCount, 0.001, 0.999);
  double bits = -log2(pDetail);
  if (detail >= kAttributeResidualAlphabetSize) {
    const double x = double(detail) - double(kAttributeResidualAlphabetSize);
    bits += 2.0 * std::floor(log2(x + 1.0)) + 1.0;
  }
  return bits;
}

//----------------------------------------------------------------------------

double
PCCResidualsEntropyEstimator::bits(const uint32_t value0) const
{
  const bool isZero0 = value0 == 0;
  const double pIsZero0 = isZero0
    ? double(isZero0Count) / symbolCount0
    : double(symbolCount0 - isZero0Count) / symbolCount0;
  double bits = -log2(PCCClip(pIsZero0, 0.001, 0.999));
  if (!isZero0) {
    bits += bitsDetail(value0 - 1, symbolCount0, freq0);
  }
  return bits;
}

//----------------------------------------------------------------------------

void
PCCResidualsEntropyEstimator::update(const uint32_t value0)
{
  const bool isZero0 = value0 == 0;
  ++symbolCount0;
  if (!isZero0) {
    ++freq0[std::min(value0 - 1, uint32_t(kAttributeResidualAlphabetSize))];
  } else {
    ++isZero0Count;
  }
}

//----------------------------------------------------------------------------

double
PCCResidualsEntropyEstimator::bits(
  const uint32_t value0, const uint32_t value1, const uint32_t value2) const
{
  const bool isZero0 = value0 == 0;
  const double pIsZero0 = isZero0
    ? double(isZero0Count) / symbolCount0
    : double(symbolCount0 - isZero0Count) / symbolCount0;
  double bits = -log2(PCCClip(pIsZero0, 0.001, 0.999));
  if (!isZero0) {
    bits += bitsDetail(value0 - 1, symbolCount0, freq0);
  }

  const bool isZero1 = value1 == 0 && value2 == 0;
  const double pIsZero1 = isZero1
    ? double(isZero1Count) / symbolCount0
    : double(symbolCount0 - isZero1Count) / symbolCount0;
  bits -= log2(PCCClip(pIsZero1, 0.001, 0.999));
  if (!isZero1) {
    bits += bitsDetail(value1, symbolCount1, freq1);
    bits += bitsDetail(value2, symbolCount1, freq1);
  }
  return bits;
}

//----------------------------------------------------------------------------

void
PCCResidualsEntropyEstimator::update(
  const uint32_t value0, const uint32_t value1, const uint32_t value2)
{
  const bool isZero0 = value0 == 0;
  ++symbolCount0;
  if (!isZero0) {
    ++freq0[std::min(value0 - 1, uint32_t(kAttributeResidualAlphabetSize))];
  } else {
    ++isZero0Count;
  }

  const bool isZero1 = value1 == 0 && value2 == 0;
  symbolCount1 += 2;
  if (!isZero1) {
    ++freq1[std::min(value1, uint32_t(kAttributeResidualAlphabetSize))];
    ++freq1[std::min(value2, uint32_t(kAttributeResidualAlphabetSize))];
  } else {
    ++isZero1Count;
  }
}

//============================================================================
// AttributeEncoderIntf

AttributeEncoderIntf::~AttributeEncoderIntf() = default;

//============================================================================
// AttributeEncoder factory

std::unique_ptr<AttributeEncoderIntf>
makeAttributeEncoder()
{
  return std::unique_ptr<AttributeEncoder>(new AttributeEncoder());
}

//============================================================================
// AttributeEncoder Members

void
AttributeEncoder::encode(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const AttributeDescription& desc,
  const AttributeParameterSet& attr_aps,
  AttributeBrickHeader& abh,
  AttributeContexts& ctxtMem,
  PCCPointSet3& pointCloud,
  PayloadBuffer* payload,
  AttributeInterPredParams& attrInterPredParams,
  attr::ModeEncoder& predEncoder)
{
  if (attr_aps.attr_encoding == AttributeEncoding::kRaw) {
    AttrRawEncoder::encode(sps, desc, attr_aps, abh, pointCloud, payload);
    return;
  }

  // Encoders are able to modify the slice header:
  _abh = &abh;

  QpSet qpSet = deriveQpSet(desc, attr_aps, abh);

  PCCResidualsEncoder encoder(attr_aps, abh, ctxtMem);
  encoder.start(sps, int(pointCloud.getPointCount()));
  abh.attr_layer_code_mode.clear();
  if (desc.attr_num_dimensions_minus1 == 0) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      encodeReflectancesTransformRaht(
        desc, attr_aps, abh, qpSet, pointCloud, encoder, predEncoder,
        attrInterPredParams);
      break;

    case AttributeEncoding::kRaw:
      // Already handled
      break;
    default:
      std::runtime_error("Not supported attributes encoding");
    }
  } else if (desc.attr_num_dimensions_minus1 == 2) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      if (attrInterPredParams.enableAttrInterPred && attr_aps.dual_motion_field_flag
          && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
        if (attr_aps.mcap_to_rec_geom_flag)
          attrInterPredParams.compensatedPointCloud = pointCloud;
        attrInterPredParams.encodeMotionAndBuildCompensated(gps, encoder.arithmeticEncoder, attr_aps.mcap_to_rec_geom_flag);
      }
      encodeColorsTransformRaht(
        desc, attr_aps, abh, qpSet, pointCloud, encoder, predEncoder,
        attrInterPredParams);
      break;

    case AttributeEncoding::kRaw:
      // Already handled
      break;
    default:
      std::runtime_error("Not supported attributes encoding");
    }
  } else {
    assert(
      desc.attr_num_dimensions_minus1 == 0
      || desc.attr_num_dimensions_minus1 == 2);
  }

  uint32_t acDataLen = encoder.stop();

  // write abh
  write(sps, attr_aps, abh, payload);
  _abh = nullptr;

  std::copy_n(
    encoder.arithmeticEncoder.buffer(), acDataLen,
    std::back_inserter(*payload));

  // save the context state for re-use by a future slice if required
  ctxtMem = encoder.getCtx();
}

//----------------------------------------------------------------------------

template<const int attribCount>
inline void
encodeRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder,
  attr::ModeEncoder& predEncoder,
  const AttributeInterPredParams& attrInterPredParams)
{
  const int voxelCount = pointCloud.getPointCount();

  // Allocate arrays.
  std::vector<int64_t> mortonCode;
  std::vector<int> attributes;
  std::vector<Qps> pointQpOffsets;
  std::vector<int> coefficients(attribCount * voxelCount);

  // Populate input arrays.
  auto indexOrd =
    sortedPointCloud(attribCount, pointCloud, mortonCode, attributes);
  pointQpOffsets.reserve(voxelCount);
  for (auto index : indexOrd) {
    pointQpOffsets.push_back(qpSet.regionQpOffset(pointCloud[index]));
  }

  if (attrInterPredParams.hasLocalMotion()) {
    predEncoder.set(&encoder.arithmeticEncoder);
    const int voxelCount_mc =
      int(attrInterPredParams.compensatedPointCloud.getPointCount());

    // Allocate arrays.
    std::vector<int64_t> mortonCode_mc;
    std::vector<int> attributes_mc;
    sortedPointCloud(
      attribCount, attrInterPredParams.compensatedPointCloud, mortonCode_mc,
      attributes_mc);
    uint64_t maxMortonCode = mortonCode.back();
    int depth = (ilog2(maxMortonCode|7) + 2) / 3;
    abh.attr_layer_code_mode.resize(depth, 0);
    // Transform.
    regionAdaptiveHierarchicalTransform(
      aps.rahtPredParams, abh, qpSet, pointQpOffsets.data(), attribCount,
      voxelCount, mortonCode.data(), attributes.data(), voxelCount_mc,
      mortonCode_mc.data(), attributes_mc.data(), coefficients.data(),
      predEncoder);
  } else {
    predEncoder.reset();
    predEncoder.set(&encoder.arithmeticEncoder);

    // Transform.
    regionAdaptiveHierarchicalTransform(
      aps.rahtPredParams, abh, qpSet, pointQpOffsets.data(), attribCount,
      voxelCount, mortonCode.data(), attributes.data(), 0, nullptr, nullptr,
      coefficients.data(), predEncoder);
  }

  // Entropy encode.
  int zeroRun = 0;
  if (attribCount == 3) {
    int values[attribCount];
    for (int n = 0; n < voxelCount; ++n) {
      for (int d = 0; d < attribCount; ++d) {
        values[d] = coefficients[voxelCount * d + n];
      }
      if (!values[0] && !values[1] && !values[2])
        ++zeroRun;
      else {
        encoder.encodeRunLength(zeroRun);
        encoder.encode(values[0], values[1], values[2]);
        zeroRun = 0;
      }
    }
  } else if (attribCount == 1) {
    for (int n = 0; n < voxelCount; ++n) {
      auto value = coefficients[n];
      if (!value)
        ++zeroRun;
      else {
        encoder.encodeRunLength(zeroRun);
        encoder.encode(value);
        zeroRun = 0;
      }
    }
  }
  if (zeroRun)
    encoder.encodeRunLength(zeroRun);

  if (aps.rahtPredParams.raht_enable_inter_intra_layer_RDO) {
    int codeModeSize = abh.attr_layer_code_mode.size();
    for (int layerIdx = 0; layerIdx < codeModeSize; ++layerIdx) {
      int predMode = abh.attr_layer_code_mode[layerIdx];
      encoder.encodeInterPredMode(predMode);
    }
  }

  predEncoder.flush();

  int clipMax = (1 << desc.bitdepth) - 1;
  auto attribute = attributes.begin();
  if (attribCount == 3) {
    for (auto index : indexOrd) {
      auto& color = pointCloud.getColor(index);
      color[0] = attr_t(PCCClip(*attribute++, 0, clipMax));
      color[1] = attr_t(PCCClip(*attribute++, 0, clipMax));
      color[2] = attr_t(PCCClip(*attribute++, 0, clipMax));
    }
  } else if (attribCount == 1) {
    for (auto index : indexOrd) {
      auto& refl = pointCloud.getReflectance(index);
      refl = attr_t(PCCClip(*attribute++, 0, clipMax));
    }
  }
}

void
AttributeEncoder::encodeReflectancesTransformRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder,
  attr::ModeEncoder& predEncoder,
  const AttributeInterPredParams& attrInterPredParams)
{
  encodeRaht<1>(
    desc, aps,abh, qpSet, pointCloud, encoder, predEncoder, attrInterPredParams);
}

void
AttributeEncoder::encodeColorsTransformRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder,
  attr::ModeEncoder& predEncoder,
  const AttributeInterPredParams& attrInterPredParams)
{
  encodeRaht<3>(
    desc, aps, abh, qpSet, pointCloud, encoder, predEncoder, attrInterPredParams);
}

//============================================================================

} /* namespace pcc */
