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
  int decodeInterPredMode(
    const RahtPredictionParams& rahtPredParams,
    int layerIndex,
    int layerDepth);
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
PCCResidualsDecoder::decodeInterPredMode(
  const RahtPredictionParams& rahtPredParams,
  int layerIndex,
  int layerDepth)
{
  int distanceToRoot = layerIndex + 1;
  int layerD = layerDepth - distanceToRoot;
  int predCtxLevel = 0;
  if (rahtPredParams.enable_inter_prediction) {
    predCtxLevel = layerD - rahtPredParams.mode_level;
    if (predCtxLevel >= NUMBER_OF_LEVELS_MODE)
      predCtxLevel = NUMBER_OF_LEVELS_MODE - 1;
  } else if (rahtPredParams.prediction_enabled_flag) {
    predCtxLevel = layerD - rahtPredParams.intra_mode_level;
    if (predCtxLevel >= NUMBER_OF_LEVELS_MODE)
      predCtxLevel = NUMBER_OF_LEVELS_MODE - 1;
  }

  bool upperInferMode =
    distanceToRoot < rahtPredParams.upper_mode_level
    && distanceToRoot < layerDepth - rahtPredParams.mode_level + 1;

  bool realInferInLowerLevel =
    rahtPredParams.enable_average_prediction
    ? (distanceToRoot >= (layerDepth - rahtPredParams.mode_level
        + rahtPredParams.lower_mode_level_for_average_prediction + 1)
      && !upperInferMode)
    : predCtxLevel < 0 && !upperInferMode;

  bool isLayerMode;
  isLayerMode = arithmeticDecoder.decode(ctxLayerPred);
  if (!isLayerMode)
    return 0;

  if (upperInferMode)
    return 2;

  if (realInferInLowerLevel)
    return 1;

  bool isInterLayerMode = arithmeticDecoder.decode(ctxInterLayerPred);
  if (isInterLayerMode)
    return 1;
  return 2;
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
  const GeometryParameterSet& gps,
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  AttributeBrickHeader& abh,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  const char* payload,
  size_t payloadLen,
  PCCPointSet3& pointCloud,
  AttributeInterPredParams& attrInterPredParams,
  attr::ModeDecoder& predDecoder)
{
  if (attr_aps.attr_encoding == AttributeEncoding::kRaw) {
    AttrRawDecoder::decode(
      attr_desc, attr_aps, abh, payload, payloadLen, pointCloud);
    return;
  }

  auto& decoder = *_pDecoder.get();

  if (attr_desc.attr_num_dimensions_minus1 == 0) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeReflectancesRaht(
        attr_desc, attr_aps, abh, _qpSet, decoder, pointCloud, predDecoder,
        attrInterPredParams);
      break;

    case AttributeEncoding::kRaw:
      // Already handled
      break;
    }
  } else if (attr_desc.attr_num_dimensions_minus1 == 2) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      if (attrInterPredParams.enableAttrInterPred && attr_aps.dual_motion_field_flag
          && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
        if (attr_aps.mcap_to_rec_geom_flag)
          attrInterPredParams.compensatedPointCloud = pointCloud;
        attrInterPredParams.decodeMotionAndBuildCompensated(
          attr_aps.motion, decoder.arithmeticDecoder,
          attr_aps.mcap_to_rec_geom_flag);
      }
      decodeColorsRaht(
        attr_desc, attr_aps, abh, _qpSet, decoder, pointCloud, predDecoder,
        attrInterPredParams);
      break;

    case AttributeEncoding::kRAHTperBlock:
      if (
        attrInterPredParams.enableAttrInterPred
        && attr_aps.dual_motion_field_flag
        && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
        if (attr_aps.mcap_to_rec_geom_flag)
          attrInterPredParams.compensatedPointCloud = pointCloud;
        attrInterPredParams.decodeMotionAndBuildCompensated(
          attr_aps.motion, decoder.arithmeticDecoder,
          attr_aps.mcap_to_rec_geom_flag);
      }
      decodeRAHTperBlock(
        attr_desc, attr_aps, abh, _qpSet, decoder, pointCloud, predDecoder,
        attrInterPredParams);
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
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeSlab(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  AttributeBrickHeader& abh,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  const char* payload,
  size_t payloadLen,
  PCCPointSet3& slabPointCloud,
  AttributeInterPredParams& attrInterPredParams,
  attr::ModeDecoder& predDecoder)
{
  if (attr_aps.attr_encoding == AttributeEncoding::kRaw) {
    throw std::runtime_error("Not Supported yet");
    return;
  }

  PCCPointSet3 tmp;
  if (attrInterPredParams.enableAttrInterPred
      && attr_aps.dual_motion_field_flag
      && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
    // compensatedPointCloud is needed by geometry
    tmp.swap(attrInterPredParams.compensatedPointCloud);
    attrInterPredParams.mortonCode_mc.clear();
    attrInterPredParams.attributes_mc.clear();
    if (attr_aps.mcap_to_rec_geom_flag)
      attrInterPredParams.compensatedPointCloud = slabPointCloud;
    attrInterPredParams.decodeMotionAndBuildCompensated(
      attr_aps.motion, _pDecoder->arithmeticDecoder, attr_aps.mcap_to_rec_geom_flag);
  }
  if (attr_desc.attr_num_dimensions_minus1 == 0) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeReflectancesRaht(
        attr_desc, attr_aps, abh, _qpSet, *_pDecoder, slabPointCloud, predDecoder,
        attrInterPredParams);
      break;

    case AttributeEncoding::kRaw:
      // Already handled
      break;
    }
  } else if (attr_desc.attr_num_dimensions_minus1 == 2) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeColorsRaht(
        attr_desc, attr_aps, abh, _qpSet, *_pDecoder, slabPointCloud, predDecoder,
        attrInterPredParams);
      break;

    case AttributeEncoding::kRAHTperBlock:
      decodeRAHTperBlock(
        attr_desc, attr_aps, abh, _qpSet, *_pDecoder, slabPointCloud, predDecoder,
        attrInterPredParams);
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
  if (attrInterPredParams.enableAttrInterPred
      && attr_aps.dual_motion_field_flag
      && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
    tmp.swap(attrInterPredParams.compensatedPointCloud);
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::startDecode(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  const char* payload,
  size_t payloadLen,
  const AttributeContexts& ctxtMem)
{
  if (attr_aps.attr_encoding == AttributeEncoding::kRaw) {
    return;
  }

  _qpSet = deriveQpSet(attr_desc, attr_aps, abh);

  _pDecoder.reset(new PCCResidualsDecoder(abh, ctxtMem));
  _pDecoder->start(sps, payload, payloadLen);
}

//----------------------------------------------------------------------------

void
AttributeDecoder::finishDecode(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  AttributeContexts& ctxtMem)
{
  if (attr_aps.attr_encoding == AttributeEncoding::kRaw) {
    return;
  }

  _pDecoder->stop();

  // save the context state for re-use by a future slice if required
  ctxtMem = _pDecoder->getCtx();
}

//----------------------------------------------------------------------------

template<const int attribCount>
void
rahtEntropyDecoder(
  const int voxelCount, int* coefficients, PCCResidualsDecoder& decoder)
{
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
}

//----------------------------------------------------------------------------

template<const int attribCount>
inline void
decodeRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud,
  attr::ModeDecoder& predDecoder,
  const AttributeInterPredParams& attrInterPredParams)
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

  for (auto index : indexOrd) {
    pointQpOffsets.push_back(qpSet.regionQpOffset(pointCloud[index]));
  }
  abh.attr_layer_code_mode.clear();

  rahtEntropyDecoder<attribCount>(voxelCount, coefficients.data(), decoder);

  if (attrInterPredParams.hasLocalMotion()) {
    predDecoder.set(&decoder.arithmeticDecoder);
    const int voxelCount_mc =
      int(attrInterPredParams.compensatedPointCloud.getPointCount());
    uint64_t maxMortonCode = mortonCode.back() | 7;
    int depth = 0;
    while (maxMortonCode) {
      maxMortonCode >>= 3;
      depth++;
    }
    abh.attr_layer_code_mode.resize(depth, 0);
    int codeModeSize = abh.attr_layer_code_mode.size();
    for (int layerIdx = 0; layerIdx < codeModeSize; ++layerIdx) {
      int& predMode = abh.attr_layer_code_mode[layerIdx];
      predMode = decoder.decodeInterPredMode(
        aps.rahtPredParams, layerIdx, codeModeSize);
    }

    auto& mortonCode_mc = attrInterPredParams.mortonCode_mc;
    auto& attributes_mc = attrInterPredParams.attributes_mc;
    if (mortonCode_mc.empty())
      sortedPointCloud(
        attribCount, attrInterPredParams.compensatedPointCloud, mortonCode_mc,
        attributes_mc);

    regionAdaptiveHierarchicalInverseTransform(
      aps.rahtPredParams, abh, qpSet, pointQpOffsets.data(), attribCount,
      voxelCount, mortonCode.data(), attributes.data(), voxelCount_mc,
      mortonCode_mc.data(), attributes_mc.data(), coefficients.data(),
      predDecoder);
  } else {
    predDecoder.reset();
    predDecoder.set(&decoder.arithmeticDecoder);

    regionAdaptiveHierarchicalInverseTransform(
      aps.rahtPredParams, abh, qpSet, pointQpOffsets.data(), attribCount,
      voxelCount, mortonCode.data(), attributes.data(), 0, nullptr, nullptr,
      coefficients.data(), predDecoder);
  }

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

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeRAHTperBlock(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud,
  attr::ModeDecoder& predDecoder,
  const AttributeInterPredParams& attrInterPredParams)
{
  const int attribCount = 3;
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

  for (auto index : indexOrd) {
    pointQpOffsets.push_back(qpSet.regionQpOffset(pointCloud[index]));
  }

  const int prefix_shift = 3 * aps.block_size_log2;
  abh.attr_layer_code_mode.clear();
  if (attrInterPredParams.hasLocalMotion()) {
    predDecoder.set(&decoder.arithmeticDecoder);
    const int voxelCount_mc =
      int(attrInterPredParams.compensatedPointCloud.getPointCount());
    //std::cout << "Using inter MC for prediction" << std::endl;

    auto& mortonCode_mc = attrInterPredParams.mortonCode_mc;
    auto& attributes_mc = attrInterPredParams.attributes_mc;
    if (mortonCode_mc.empty())
      sortedPointCloud(
        attribCount, attrInterPredParams.compensatedPointCloud, mortonCode_mc,
        attributes_mc);

    abh.attr_layer_code_mode.resize(aps.block_size_log2);

    int block_pc_begin = 0;
    int block_mc_begin = 0;
    while (block_pc_begin < voxelCount) {
      int64_t prefix = mortonCode[block_pc_begin] >> prefix_shift;

      int block_pc_end = block_pc_begin + 1;
      while (block_pc_end < voxelCount
             && (mortonCode[block_pc_end] >> prefix_shift) == prefix)
        block_pc_end++;

      while (block_mc_begin < voxelCount_mc
             && (mortonCode_mc[block_mc_begin] >> prefix_shift) < prefix)
        block_mc_begin++;
      int block_mc_end = block_mc_begin;
      while (block_mc_end < voxelCount_mc
             && (mortonCode_mc[block_mc_end] >> prefix_shift) == prefix)
        block_mc_end++;

      rahtEntropyDecoder<3>(
        block_pc_end - block_pc_begin,
        coefficients.data() + attribCount * block_pc_begin, decoder);

      int codeModeSize = abh.attr_layer_code_mode.size();
      for (int layerIdx = 0; layerIdx < codeModeSize; ++layerIdx) {
        int& predMode = abh.attr_layer_code_mode[layerIdx];
        predMode = decoder.decodeInterPredMode(
          aps.rahtPredParams, layerIdx, codeModeSize);
      }

      // Transform.
      regionAdaptiveHierarchicalInverseTransform(
        aps.rahtPredParams, abh, qpSet, pointQpOffsets.data() + block_pc_begin,
        attribCount, block_pc_end - block_pc_begin,
        mortonCode.data() + block_pc_begin,
        attributes.data() + attribCount * block_pc_begin,
        block_mc_end - block_mc_begin, mortonCode_mc.data() + block_mc_begin,
        attributes_mc.data() + attribCount * block_mc_begin,
        coefficients.data() + attribCount * block_pc_begin, predDecoder);

      block_pc_begin = block_pc_end;
      block_mc_begin = block_mc_end;
    }
  } else {
    predDecoder.reset();
    predDecoder.set(&decoder.arithmeticDecoder);

    int block_pc_begin = 0;
    while (block_pc_begin < voxelCount) {
      int64_t prefix = mortonCode[block_pc_begin] >> prefix_shift;

      int block_pc_end = block_pc_begin + 1;
      while (block_pc_end < voxelCount
             && (mortonCode[block_pc_end] >> prefix_shift) == prefix)
        block_pc_end++;

      rahtEntropyDecoder<3>(
        block_pc_end - block_pc_begin,
        coefficients.data() + attribCount * block_pc_begin, decoder);

      // Transform.
      regionAdaptiveHierarchicalInverseTransform(
        aps.rahtPredParams, abh, qpSet, pointQpOffsets.data() + block_pc_begin,
        attribCount, block_pc_end - block_pc_begin,
        mortonCode.data() + block_pc_begin,
        attributes.data() + attribCount * block_pc_begin, 0, nullptr, nullptr,
        coefficients.data() + attribCount * block_pc_begin, predDecoder);

      block_pc_begin = block_pc_end;
    }
  }

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
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud,
  attr::ModeDecoder& predDecoder,
  const AttributeInterPredParams& attrInterPredParams)
{
  decodeRaht<1>(
    desc, aps, abh, qpSet, decoder, pointCloud, predDecoder, attrInterPredParams);
}

void
AttributeDecoder::decodeColorsRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud,
  attr::ModeDecoder& predDecoder,
  const AttributeInterPredParams& attrInterPredParams)
{
  decodeRaht<3>(
    desc, aps, abh,qpSet, decoder, pointCloud, predDecoder, attrInterPredParams);
}

//============================================================================

} /* namespace pcc */
