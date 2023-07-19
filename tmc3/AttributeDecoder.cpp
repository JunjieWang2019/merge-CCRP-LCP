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

#include "AttributeDecoder.h"

#include "AttributeCommon.h"
#include "attribute_raw.h"
#include "constants.h"
#include "entropy.h"
#include "hls.h"
#include "io_hls.h"
#include "RAHT.h"
#include "FixedPoint.h"

namespace pcc {

//============================================================================
// An encapsulation of the entropy decoding methods used in attribute coding

class PCCResidualsDecoder : protected AttributeContexts {
public:
  PCCResidualsDecoder(
    const AttributeBrickHeader& abh, const AttributeContexts& ctxtMem);

  EntropyDecoder arithmeticDecoder;

  const AttributeContexts& getCtx() const { return *this; }

  void start(const SequenceParameterSet& sps, const char* buf, int buf_len);
  void stop();

  int decodeRunLength();
  int decodeSymbol(int k1, int k2, int k3);
  void decode(int32_t values[3]);
  int32_t decode();
};

//----------------------------------------------------------------------------

PCCResidualsDecoder::PCCResidualsDecoder(
  const AttributeBrickHeader& abh, const AttributeContexts& ctxtMem)
  : AttributeContexts(ctxtMem)
{}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::start(
  const SequenceParameterSet& sps, const char* buf, int buf_len)
{
  arithmeticDecoder.setBuffer(buf_len, buf);
  arithmeticDecoder.enableBypassStream(sps.cabac_bypass_stream_enabled_flag);
  arithmeticDecoder.setBypassBinCodingWithoutProbUpdate(sps.bypass_bin_coding_without_prob_update);
  arithmeticDecoder.start();
}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::stop()
{
  arithmeticDecoder.stop();
}

//----------------------------------------------------------------------------

int
PCCResidualsDecoder::decodeRunLength()
{
  int runLength = 0;
  auto* ctx = ctxRunLen;
  for (; runLength < 3; runLength++, ctx++) {
    int bin = arithmeticDecoder.decode(*ctx);
    if (!bin)
      return runLength;
  }

  for (int i = 0; i < 4; i++) {
    int bin = arithmeticDecoder.decode(*ctx);
    if (!bin) {
      runLength += arithmeticDecoder.decode();
      return runLength;
    }
    runLength += 2;
  }

  runLength += arithmeticDecoder.decodeExpGolomb(2, *++ctx);
  return runLength;
}

//----------------------------------------------------------------------------

int
PCCResidualsDecoder::decodeSymbol(int k1, int k2, int k3)
{
  if (!arithmeticDecoder.decode(ctxCoeffGtN[0][k1]))
    return 0;

  if (!arithmeticDecoder.decode(ctxCoeffGtN[1][k2]))
    return 1;

  int coeff_abs_minus2 = arithmeticDecoder.decodeExpGolomb(
    1, ctxCoeffRemPrefix[k3], ctxCoeffRemSuffix[k3]);

  return coeff_abs_minus2 + 2;
}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::decode(int32_t value[3])
{
  value[1] = decodeSymbol(0, 0, 1);
  int b0 = value[1] == 0;
  int b1 = value[1] <= 1;
  value[2] = decodeSymbol(1 + b0, 1 + b1, 1);
  int b2 = value[2] == 0;
  int b3 = value[2] <= 1;
  value[0] = decodeSymbol(3 + (b0 << 1) + b2, 3 + (b1 << 1) + b3, 0);

  if (b0 && b2)
    value[0] += 1;

  if (value[0] && arithmeticDecoder.decode())
    value[0] = -value[0];
  if (value[1] && arithmeticDecoder.decode())
    value[1] = -value[1];
  if (value[2] && arithmeticDecoder.decode())
    value[2] = -value[2];
}

//----------------------------------------------------------------------------

int32_t
PCCResidualsDecoder::decode()
{
  auto mag = decodeSymbol(0, 0, 0) + 1;
  bool sign = arithmeticDecoder.decode();
  return sign ? -mag : mag;
}

//============================================================================
// AttributeDecoderIntf

AttributeDecoderIntf::~AttributeDecoderIntf() = default;

//============================================================================
// AttributeDecoder factory

std::unique_ptr<AttributeDecoderIntf>
makeAttributeDecoder()
{
  return std::unique_ptr<AttributeDecoder>(new AttributeDecoder());
}

//============================================================================
// AttributeDecoder Members

void
AttributeDecoder::decode(
  const SequenceParameterSet& sps,
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  const char* payload,
  size_t payloadLen,
  AttributeContexts& ctxtMem,
  PCCPointSet3& pointCloud  ,
   const AttributeInterPredParams& attrInterPredParams
  )
{
  if (attr_aps.attr_encoding == AttributeEncoding::kRaw) {
    AttrRawDecoder::decode(
      attr_desc, attr_aps, abh, payload, payloadLen, pointCloud);
    return;
  }

  QpSet qpSet = deriveQpSet(attr_desc, attr_aps, abh);

  PCCResidualsDecoder decoder(abh, ctxtMem);
  decoder.start(sps, payload, payloadLen);

  if (attr_desc.attr_num_dimensions_minus1 == 0) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeReflectancesRaht(attr_desc, attr_aps, qpSet, decoder, pointCloud);
      break;

    case AttributeEncoding::kRaw:
      // Already handled
      break;
    }
  } else if (attr_desc.attr_num_dimensions_minus1 == 2) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeColorsRaht(attr_desc, attr_aps, qpSet, decoder, pointCloud);
      break;

    case AttributeEncoding::kRaw:
      // Already handled
      break;
    }
  } else {
    assert(
      attr_desc.attr_num_dimensions_minus1 == 0
      || attr_desc.attr_num_dimensions_minus1 == 2);
  }

  decoder.stop();

  // save the context state for re-use by a future slice if required
  ctxtMem = decoder.getCtx();
}

//----------------------------------------------------------------------------

template<const int attribCount>
inline void
decodeRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const int voxelCount = pointCloud.getPointCount();

  // Morton codes
  std::vector<int64_t> mortonCode;
  std::vector<int> attributes;
  auto indexOrd =
    sortedPointCloud(attribCount, pointCloud, mortonCode, attributes);
  attributes.resize(voxelCount * attribCount);

  // Entropy decode
  std::vector<int> coefficients(attribCount * voxelCount, 0);
  std::vector<Qps> pointQpOffsets;
  pointQpOffsets.reserve(voxelCount);

  int n = 0;
  for (auto index : indexOrd) {
    pointQpOffsets.push_back(qpSet.regionQpOffset(pointCloud[index]));
  }

  // Decode coefficients
  if (attribCount == 3) {
    int32_t values[3];

    for (int n = 0; n < voxelCount; n++) {
      int zeroRun = decoder.decodeRunLength();
      if (zeroRun) {
        n += zeroRun;
        if (n >= voxelCount)
          break;
      }

      decoder.decode(values);
      for (int d = 0; d < 3; d++)
        coefficients[n + voxelCount * d] = values[d];
    }
  } else if (attribCount == 1) {
    for (int n = 0; n < voxelCount; n++) {
      int zeroRun = decoder.decodeRunLength();
      if (zeroRun) {
        n += zeroRun;
        if (n >= voxelCount)
          break;
      }
      coefficients[n] = decoder.decode();
    }
  }

  regionAdaptiveHierarchicalInverseTransform(
    aps.rahtPredParams, qpSet,
    pointQpOffsets.data(), mortonCode.data(), attributes.data(), attribCount,
    voxelCount, coefficients.data());

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
AttributeDecoder::decodeReflectancesRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  decodeRaht<1>(desc, aps, qpSet, decoder, pointCloud);
}

void
AttributeDecoder::decodeColorsRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  decodeRaht<3>(desc, aps, qpSet, decoder, pointCloud);
}

//============================================================================

} /* namespace pcc */
