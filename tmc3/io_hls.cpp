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

#include "BitReader.h"
#include "BitWriter.h"
#include "PCCMisc.h"
#include "hls.h"
#include "io_hls.h"
#include "motionWip.h"

#include <iomanip>
#include <iterator>
#include <algorithm>
#include <sstream>

namespace pcc {

//============================================================================

Oid::operator std::string() const
{
  std::stringstream ss;
  int subidentifier = 0;
  int firstSubidentifier = 1;
  for (auto byte : this->contents) {
    if (byte & 0x80) {
      subidentifier = (subidentifier << 7) | (byte & 0x7f);
      continue;
    }

    // end of subidentifer.
    // NB: the first subidentifier encodes two oid components
    if (firstSubidentifier) {
      firstSubidentifier = 0;
      if (subidentifier < 40)
        ss << '0';
      else if (subidentifier < 80) {
        ss << '1';
        subidentifier -= 40;
      } else {
        ss << '2';
        subidentifier -= 80;
      }
    }

    ss << '.' << std::to_string(subidentifier);
    subidentifier = 0;
  }

  return ss.str();
}

//----------------------------------------------------------------------------

bool
operator==(const Oid& lhs, const Oid& rhs)
{
  // NB: there is a unique encoding for each OID.  Equality may be determined
  // comparing just the content octets
  return lhs.contents == rhs.contents;
}

//----------------------------------------------------------------------------

template<typename T>
void
writeOid(T& bs, const Oid& oid)
{
  // write out the length according to the BER definite shoft form.
  // NB: G-PCC limits the length to 127 octets.
  constexpr int oid_reserved_zero_bit = 0;
  int oid_len = oid.contents.size();
  bs.writeUn(1, oid_reserved_zero_bit);
  bs.writeUn(7, oid_len);

  const auto& oid_contents = oid.contents;
  for (int i = 0; i < oid_len; i++)
    bs.writeUn(8, oid_contents[i]);
}

//----------------------------------------------------------------------------

template<typename T>
void
readOid(T& bs, Oid* oid)
{
  int oid_reserved_zero_bit;
  int oid_len;
  bs.readUn(1, &oid_reserved_zero_bit);
  bs.readUn(7, &oid_len);
  oid->contents.resize(oid_len);
  for (int i = 0; i < oid_len; i++)
    bs.readUn(8, &oid->contents[i]);
}

//----------------------------------------------------------------------------

int
lengthOid(const Oid& oid)
{
  return 1 + oid.contents.size();
}

//============================================================================

std::ostream&
operator<<(std::ostream& os, const AttributeLabel& label)
{
  switch (label.known_attribute_label) {
  case KnownAttributeLabel::kColour: os << "color"; break;
  case KnownAttributeLabel::kReflectance: os << "reflectance"; break;
  case KnownAttributeLabel::kFrameIndex: os << "frame index"; break;
  case KnownAttributeLabel::kFrameNumber: os << "frame number"; break;
  case KnownAttributeLabel::kMaterialId: os << "material id"; break;
  case KnownAttributeLabel::kOpacity: os << "opacity"; break;
  case KnownAttributeLabel::kNormal: os << "normal"; break;
  case KnownAttributeLabel::kOid: os << std::string(label.oid);
  default:
    // An unknown known attribute
    auto iosFlags = os.flags(std::ios::hex);
    os << std::setw(8) << label.known_attribute_label;
    os.flags(iosFlags);
  }

  return os;
}

//============================================================================

template<typename Bs>
void
writeAttrParamCicp(Bs& bs, const AttributeParameters& param)
{
  bs.writeUe(param.cicp_colour_primaries_idx);
  bs.writeUe(param.cicp_transfer_characteristics_idx);
  bs.writeUe(param.cicp_matrix_coefficients_idx);
  bs.write(param.cicp_video_full_range_flag);
  bs.byteAlign();
}

//----------------------------------------------------------------------------

template<typename Bs>
void
parseAttrParamCicp(Bs& bs, AttributeParameters* param)
{
  bs.readUe(&param->cicp_colour_primaries_idx);
  bs.readUe(&param->cicp_transfer_characteristics_idx);
  bs.readUe(&param->cicp_matrix_coefficients_idx);
  bs.read(&param->cicp_video_full_range_flag);
  param->cicpParametersPresent = 1;
  bs.byteAlign();
}

//============================================================================

template<typename Bs>
void
writeAttrParamScaling(Bs& bs, const AttributeParameters& param)
{
  int attr_offset_bits = numBits(std::abs(param.attr_offset));
  int attr_scale_bits = numBits(param.attr_scale_minus1);
  bs.writeUe(attr_offset_bits);
  bs.writeSn(attr_offset_bits, param.attr_offset);
  bs.writeUe(attr_scale_bits);
  bs.writeUn(attr_scale_bits, param.attr_scale_minus1);
  bs.writeUe(param.attr_frac_bits);
  bs.byteAlign();
}

//----------------------------------------------------------------------------

template<typename Bs>
void
parseAttrParamScaling(Bs& bs, AttributeParameters* param)
{
  int attr_offset_bits, attr_scale_bits;
  bs.readUe(&attr_offset_bits);
  bs.readSn(attr_offset_bits, &param->attr_offset);
  bs.readUe(&attr_scale_bits);
  bs.readUn(attr_scale_bits, &param->attr_scale_minus1);
  bs.readUe(&param->attr_frac_bits);
  param->scalingParametersPresent = 1;
  bs.byteAlign();
}

//============================================================================

template<typename Bs>
void
writeAttrParamDefaultValue(
  const AttributeDescription& desc, Bs& bs, const AttributeParameters& param)
{
  bs.writeUn(desc.bitdepth, param.attr_default_value[0]);
  for (int k = 1; k <= desc.attr_num_dimensions_minus1; k++)
    bs.writeUn(desc.bitdepth, param.attr_default_value[k]);
  bs.byteAlign();
}

//----------------------------------------------------------------------------

template<typename Bs>
void
parseAttrParamDefaultValue(
  const AttributeDescription& desc, Bs& bs, AttributeParameters* param)
{
  param->attr_default_value.resize(desc.attr_num_dimensions_minus1 + 1);

  bs.readUn(desc.bitdepth, &param->attr_default_value[0]);
  for (int k = 1; k <= desc.attr_num_dimensions_minus1; k++)
    bs.readUn(desc.bitdepth, &param->attr_default_value[k]);
  bs.byteAlign();
}

//============================================================================

template<typename Bs>
void
writeAttrParamOpaque(Bs& bs, const OpaqueAttributeParameter& param)
{
  if (param.attr_param_type == AttributeParameterType::kItuT35) {
    bs.writeUn(8, param.attr_param_itu_t_t35_country_code);
    if (param.attr_param_itu_t_t35_country_code == 0xff)
      bs.writeUn(8, param.attr_param_itu_t_t35_country_code_extension);
  } else if (param.attr_param_type == AttributeParameterType::kOid)
    writeOid(bs, param.attr_param_oid);

  for (auto attr_param_byte : param.attr_param_byte)
    bs.writeUn(8, attr_param_byte);

  bs.byteAlign();
}

//----------------------------------------------------------------------------

template<typename Bs>
OpaqueAttributeParameter
parseAttrParamOpaque(
  Bs& bs, AttributeParameterType attr_param_type, int attrParamLen)
{
  bs.byteAlign();

  OpaqueAttributeParameter param;
  param.attr_param_type = attr_param_type;

  if (param.attr_param_type == AttributeParameterType::kItuT35) {
    bs.readUn(8, &param.attr_param_itu_t_t35_country_code);
    attrParamLen--;
    if (param.attr_param_itu_t_t35_country_code == 0xff) {
      bs.readUn(8, &param.attr_param_itu_t_t35_country_code_extension);
      attrParamLen--;
    }
  } else if (param.attr_param_type == AttributeParameterType::kOid) {
    readOid(bs, &param.attr_param_oid);
    attrParamLen -= lengthOid(param.attr_param_oid);
  }

  if (attrParamLen > 0) {
    param.attr_param_byte.resize(attrParamLen);
    for (int i = 0; i < attrParamLen; i++)
      bs.readUn(8, &param.attr_param_byte[i]);
  }

  return param;
}

//============================================================================
// NB: this writes all present parameters, whereas parse parses only one
// The encoder works in a fixed order.  However, this is non-normative.

template<typename T>
void
writeAttributeParameters(
  const AttributeDescription& attr, T& bs, const AttributeParameters& params)
{
  int num_attr_parameters = params.numParams();
  bs.writeUe(num_attr_parameters);
  bs.byteAlign();

  if (!params.attr_default_value.empty()) {
    int attr_param_len = 0;
    auto bsCounter = makeBitWriter(InsertionCounter(&attr_param_len));
    writeAttrParamDefaultValue(attr, bsCounter, params);

    auto attr_param_type = AttributeParameterType::kDefaultValue;
    bs.writeUn(8, attr_param_type);
    bs.writeUn(8, attr_param_len);
    writeAttrParamDefaultValue(attr, bs, params);
  }

  if (params.cicpParametersPresent) {
    int attr_param_len = 0;
    auto bsCounter = makeBitWriter(InsertionCounter(&attr_param_len));
    writeAttrParamCicp(bsCounter, params);

    auto attr_param_type = AttributeParameterType::kCicp;
    bs.writeUn(8, attr_param_type);
    bs.writeUn(8, attr_param_len);
    writeAttrParamCicp(bs, params);
  }

  if (params.scalingParametersPresent) {
    int attr_param_len = 0;
    auto bsCounter = makeBitWriter(InsertionCounter(&attr_param_len));
    writeAttrParamScaling(bsCounter, params);

    auto attr_param_type = AttributeParameterType::kScaling;
    bs.writeUn(8, attr_param_type);
    bs.writeUn(8, attr_param_len);
    writeAttrParamScaling(bs, params);
  }

  for (const auto& param : params.opaqueParameters) {
    int attr_param_len = 0;
    auto bsCounter = makeBitWriter(InsertionCounter(&attr_param_len));
    writeAttrParamOpaque(bsCounter, param);

    bs.writeUn(8, param.attr_param_type);
    bs.writeUn(8, attr_param_len);
    writeAttrParamOpaque(bs, param);
  }
}

//----------------------------------------------------------------------------

template<typename T>
void
parseAttributeParameter(
  const AttributeDescription& attr, T& bs, AttributeParameters& params)
{
  AttributeParameterType attr_param_type;
  int attr_param_len;
  bs.readUn(8, &attr_param_type);
  bs.readUn(8, &attr_param_len);
  // todo(df): check that all attr_param_len bytes are consumed
  switch (attr_param_type) {
    using Type = AttributeParameterType;
  case Type::kCicp: parseAttrParamCicp(bs, &params); break;
  case Type::kScaling: parseAttrParamScaling(bs, &params); break;
  case Type::kDefaultValue:
    parseAttrParamDefaultValue(attr, bs, &params);
    break;

  case Type::kItuT35:
  case Type::kOid:
  default:
    params.opaqueParameters.emplace_back(
      parseAttrParamOpaque(bs, attr_param_type, attr_param_len));
  }
}

//============================================================================

PayloadBuffer
write(const SequenceParameterSet& sps)
{
  PayloadBuffer buf(PayloadType::kSequenceParameterSet);
  auto bs = makeBitWriter(std::back_inserter(buf));

  // NB: if taking bits from reserved_profile_compatibility_XXbits, be sure
  // not to change the bitstream position of the other constraint flags.
  bs.writeUn(1, sps.profile.main_profile_compatibility_flag);
  bs.writeUn(21, sps.profile.reserved_profile_compatibility_21bits);
  bs.writeUn(1, sps.profile.slice_reordering_constraint_flag);
  bs.writeUn(1, sps.profile.unique_point_positions_constraint_flag);

  bs.writeUn(8, sps.level);
  bs.writeUn(4, sps.sps_seq_parameter_set_id);

  bs.writeUn(5, sps.frame_ctr_bits);
  bs.writeUn(5, sps.slice_tag_bits);

  bs.writeUe(sps.sps_bounding_box_offset_bits);
  if (auto bbOriginBits = sps.sps_bounding_box_offset_bits) {
    auto sps_bounding_box_offset_xyz =
      toXyz(sps.geometry_axis_order, sps.seqBoundingBoxOrigin);

    bs.writeSn(bbOriginBits, sps_bounding_box_offset_xyz.x());
    bs.writeSn(bbOriginBits, sps_bounding_box_offset_xyz.y());
    bs.writeSn(bbOriginBits, sps_bounding_box_offset_xyz.z());

    int seq_bounding_box_offset_log2_scale = 0;
    bs.writeUe(seq_bounding_box_offset_log2_scale);
  }

  bs.writeUe(sps.sps_bounding_box_size_bits);
  if (auto bbSizeBits = sps.sps_bounding_box_size_bits) {
    auto seq_bounding_box_size_minus1 =
      toXyz(sps.geometry_axis_order, sps.seqBoundingBoxSize - 1);

    bs.writeUn(bbSizeBits, seq_bounding_box_size_minus1.x());
    bs.writeUn(bbSizeBits, seq_bounding_box_size_minus1.y());
    bs.writeUn(bbSizeBits, seq_bounding_box_size_minus1.z());
  }

  int seq_unit_numerator_minus1 = sps.seqGeomScale.numerator - 1;
  int seq_unit_denominator_minus1 = sps.seqGeomScale.denominator - 1;
  bs.writeUe(seq_unit_numerator_minus1);
  bs.writeUe(seq_unit_denominator_minus1);
  bs.writeUn(1, sps.seq_geom_scale_unit_flag);

  bs.writeUe(sps.global_scale_mul_log2());
  bs.writeUe(sps.global_scale_fp_bits());
  bs.writeUn(sps.global_scale_fp_bits(), sps.global_scale_rem());

  int num_attribute_sets = int(sps.attributeSets.size());
  bs.writeUe(num_attribute_sets);
  for (const auto& attr : sps.attributeSets) {
    bs.writeUe(attr.attr_num_dimensions_minus1);
    bs.writeUe(attr.attr_instance_id);

    int attr_bitdepth_minus1 = attr.bitdepth - 1;
    bs.writeUe(attr_bitdepth_minus1);

    const auto& label = attr.attributeLabel;
    bs.write(label.known_attribute_label_flag());
    if (label.known_attribute_label_flag())
      bs.writeUe(label.known_attribute_label);
    else
      writeOid(bs, label.oid);

    writeAttributeParameters(attr, bs, attr.params);
  }

  bs.writeUn(3, sps.geometry_axis_order);
  bs.write(sps.cabac_bypass_stream_enabled_flag);
  bs.write(sps.entropy_continuation_enabled_flag);

  bool sps_extension_flag = true;
  bs.write(sps_extension_flag);

  bs.write(sps.inter_frame_prediction_enabled_flag);
  bs.write(sps.bypass_bin_coding_without_prob_update);

  if (sps.inter_frame_prediction_enabled_flag)
    bs.write(sps.inter_frame_trisoup_enabled_flag);

  if (sps.inter_frame_trisoup_enabled_flag)
    bs.write(sps.inter_frame_trisoup_align_slices_flag);

  if(sps.inter_frame_trisoup_align_slices_flag)
    bs.writeUe(sps.inter_frame_trisoup_align_slices_step);

  bs.byteAlign();

  return buf;
}

//----------------------------------------------------------------------------

SequenceParameterSet
parseSps(const PayloadBuffer& buf)
{
  SequenceParameterSet sps;
  assert(buf.type == PayloadType::kSequenceParameterSet);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUn(1, &sps.profile.main_profile_compatibility_flag);
  bs.readUn(21, &sps.profile.reserved_profile_compatibility_21bits);
  bs.readUn(1, &sps.profile.slice_reordering_constraint_flag);
  bs.readUn(1, &sps.profile.unique_point_positions_constraint_flag);

  bs.readUn(8, &sps.level);
  bs.readUn(4, &sps.sps_seq_parameter_set_id);

  bs.readUn(5, &sps.frame_ctr_bits);
  bs.readUn(5, &sps.slice_tag_bits);

  sps.seqBoundingBoxOrigin = 0;
  bs.readUe(&sps.sps_bounding_box_offset_bits);
  if (auto bbOriginBits = sps.sps_bounding_box_offset_bits) {
    Vec3<int> seq_bounding_box_offset;
    bs.readSn(bbOriginBits, &seq_bounding_box_offset.x());
    bs.readSn(bbOriginBits, &seq_bounding_box_offset.y());
    bs.readSn(bbOriginBits, &seq_bounding_box_offset.z());

    int seq_bounding_box_offset_log2_scale;
    bs.readUe(&seq_bounding_box_offset_log2_scale);
    seq_bounding_box_offset *= 1 << seq_bounding_box_offset_log2_scale;

    // NB: this is in XYZ axis order until the SPS is converted to STV
    sps.seqBoundingBoxOrigin = seq_bounding_box_offset;
  }

  sps.seqBoundingBoxSize = 0;
  bs.readUe(&sps.sps_bounding_box_size_bits);
  if (auto bbSizeBits = sps.sps_bounding_box_size_bits) {
    Vec3<int> seq_bounding_box_whd_minus1;
    bs.readUn(bbSizeBits, &seq_bounding_box_whd_minus1.x());
    bs.readUn(bbSizeBits, &seq_bounding_box_whd_minus1.y());
    bs.readUn(bbSizeBits, &seq_bounding_box_whd_minus1.z());

    // NB: this is in XYZ axis order until the SPS is converted to STV
    sps.seqBoundingBoxSize = seq_bounding_box_whd_minus1 + 1;
  }

  int seq_unit_numerator_minus1;
  int seq_unit_denominator_minus1;
  bs.readUe(&seq_unit_numerator_minus1);
  bs.readUe(&seq_unit_denominator_minus1);
  bs.readUn(1, &sps.seq_geom_scale_unit_flag);

  sps.seqGeomScale.numerator = seq_unit_numerator_minus1 + 1;
  sps.seqGeomScale.denominator = seq_unit_denominator_minus1 + 1;

  bs.readUe(&sps.global_scale_mul_log2());
  bs.readUe(&sps.global_scale_fp_bits());
  bs.readUn(sps.global_scale_fp_bits(), &sps.global_scale_rem());

  int num_attribute_sets = int(bs.readUe());
  for (int i = 0; i < num_attribute_sets; i++) {
    sps.attributeSets.emplace_back();
    auto& attr = sps.attributeSets.back();
    bs.readUe(&attr.attr_num_dimensions_minus1);
    bs.readUe(&attr.attr_instance_id);

    int attr_bitdepth_minus1;
    bs.readUe(&attr_bitdepth_minus1);
    attr.bitdepth = attr_bitdepth_minus1 + 1;

    auto& label = attr.attributeLabel;
    bool known_attribute_label_flag = bs.read();
    if (known_attribute_label_flag)
      bs.readUe(&label.known_attribute_label);
    else {
      label.known_attribute_label = KnownAttributeLabel::kOid;
      readOid(bs, &label.oid);
    }

    int num_attribute_parameters;
    bs.readUe(&num_attribute_parameters);
    bs.byteAlign();
    for (int i = 0; i < num_attribute_parameters; i++)
      parseAttributeParameter(attr, bs, attr.params);
  }

  bs.readUn(3, &sps.geometry_axis_order);
  bs.read(&sps.cabac_bypass_stream_enabled_flag);
  bs.read(&sps.entropy_continuation_enabled_flag);

  // conformance check: reordering constraint must be set with continuation
  if (sps.entropy_continuation_enabled_flag)
    assert(sps.profile.slice_reordering_constraint_flag);

  sps.inter_frame_prediction_enabled_flag = false;
  sps.bypass_bin_coding_without_prob_update = false;
  sps.inter_frame_trisoup_enabled_flag = false;
  sps.inter_frame_trisoup_align_slices_flag = false;
  sps.inter_frame_trisoup_align_slices_step = 1;
  bool sps_extension_flag = bs.read();
  if (sps_extension_flag) {
    bs.read(&sps.inter_frame_prediction_enabled_flag);
    bs.read(&sps.bypass_bin_coding_without_prob_update);

    if (sps.inter_frame_prediction_enabled_flag)
      bs.read(&sps.inter_frame_trisoup_enabled_flag);

    if (sps.inter_frame_trisoup_enabled_flag)
      bs.read(&sps.inter_frame_trisoup_align_slices_flag);

    if(sps.inter_frame_trisoup_align_slices_flag)
      bs.readUe(&sps.inter_frame_trisoup_align_slices_step);
  }
  bs.byteAlign();

  return sps;
}

//----------------------------------------------------------------------------

void
convertXyzToStv(SequenceParameterSet* sps)
{
  // permute the bounding box from xyz to internal stv order
  sps->seqBoundingBoxOrigin =
    fromXyz(sps->geometry_axis_order, sps->seqBoundingBoxOrigin);

  sps->seqBoundingBoxSize =
    fromXyz(sps->geometry_axis_order, sps->seqBoundingBoxSize);
}

//============================================================================

PayloadBuffer
write(const SequenceParameterSet& sps, const GeometryParameterSet& gps)
{
  PayloadBuffer buf(PayloadType::kGeometryParameterSet);
  auto bs = makeBitWriter(std::back_inserter(buf));

  bs.writeUn(4, gps.gps_geom_parameter_set_id);
  bs.writeUn(4, gps.gps_seq_parameter_set_id);
  bs.write(gps.geom_box_log2_scale_present_flag);
  if (!gps.geom_box_log2_scale_present_flag)
    bs.writeUe(gps.gps_geom_box_log2_scale);

  bs.write(gps.geom_unique_points_flag);

  bs.write(gps.octree_point_count_list_present_flag);

  bs.writeUn(2, /*gps.inferred_direct_coding_mode*/0);
  /*if (gps.inferred_direct_coding_mode)
    bs.write(gps.joint_2pt_idcm_enabled_flag);*/

  bs.write(gps.qtbt_enabled_flag);
  bs.writeUn(3, gps.neighbour_avail_boundary_log2_minus1);

  bs.write(gps.geom_scaling_enabled_flag);
  if (gps.geom_scaling_enabled_flag) {
    bs.writeUe(gps.geom_base_qp);
    bs.writeUn(2, gps.geom_qp_multiplier_log2);
    /*if (gps.inferred_direct_coding_mode)
      bs.writeSe(gps.geom_idcm_qp_offset);*/
  }

  // NB: bitstreams conforming to the first edition must set
  // gps_extension_flag equal to 0.
  bool gps_extension_flag = sps.profile.isDraftProfile();
  bs.write(gps_extension_flag);
  if (gps_extension_flag) {
    bs.write(gps.trisoup_enabled_flag);
    if (gps.trisoup_enabled_flag){
      bs.write(gps.non_cubic_node_start_edge);
      bs.write(gps.non_cubic_node_end_edge);
    }

    // inter
    bs.write(gps.interPredictionEnabledFlag);
    if (gps.interPredictionEnabledFlag) {

      bs.writeUe(gps.motion.motion_block_size);
      bs.writeUe(gps.motion.motion_min_pu_size);

      bs.write(gps.gof_geom_entropy_continuation_enabled_flag);

      if (gps.trisoup_enabled_flag) {
        bs.write(gps.trisoup_skip_mode_enabled_flag);
      }
    }
  }
  bs.byteAlign();

  return buf;
}

//----------------------------------------------------------------------------

GeometryParameterSet
parseGps(const PayloadBuffer& buf)
{
  GeometryParameterSet gps;
  assert(buf.type == PayloadType::kGeometryParameterSet);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUn(4, &gps.gps_geom_parameter_set_id);
  bs.readUn(4, &gps.gps_seq_parameter_set_id);
  bs.read(&gps.geom_box_log2_scale_present_flag);
  if (!gps.geom_box_log2_scale_present_flag)
    bs.readUe(&gps.gps_geom_box_log2_scale);

  bs.read(&gps.geom_unique_points_flag);

  gps.octree_point_count_list_present_flag = false;

  bs.read(&gps.octree_point_count_list_present_flag);

  bs.readUn(2, &gps.inferred_direct_coding_mode);
  /*if (gps.inferred_direct_coding_mode)
    bs.read(&gps.joint_2pt_idcm_enabled_flag);*/

  bs.read(&gps.qtbt_enabled_flag);
  bs.readUn(3, &gps.neighbour_avail_boundary_log2_minus1);

  gps.geom_base_qp = 0;
  gps.geom_qp_multiplier_log2 = 0;
  gps.geom_idcm_qp_offset = 0;
  bs.read(&gps.geom_scaling_enabled_flag);
  if (gps.geom_scaling_enabled_flag) {
    bs.readUe(&gps.geom_base_qp);
    bs.readUn(2, &gps.geom_qp_multiplier_log2);
  }

  gps.trisoup_enabled_flag = false;
  gps.interPredictionEnabledFlag = false;
  gps.gof_geom_entropy_continuation_enabled_flag = false;
  gps.trisoup_skip_mode_enabled_flag = false;
  bool gps_extension_flag = bs.read();
  if (gps_extension_flag) {

    bs.read(&gps.trisoup_enabled_flag);
    if (gps.trisoup_enabled_flag) {
      bs.read(&gps.non_cubic_node_start_edge);
      bs.read(&gps.non_cubic_node_end_edge);
    }

    // inter
    bs.read(&gps.interPredictionEnabledFlag);
    if (gps.interPredictionEnabledFlag) {

      bs.readUe(&gps.motion.motion_block_size);
      bs.readUe(&gps.motion.motion_min_pu_size);

      bs.read(&gps.gof_geom_entropy_continuation_enabled_flag);

      if (gps.trisoup_enabled_flag) {
        bs.read(&gps.trisoup_skip_mode_enabled_flag);
      }
    }
  }
  bs.byteAlign();

  return gps;
}

//----------------------------------------------------------------------------

void
convertXyzToStv(const SequenceParameterSet& sps, GeometryParameterSet* gps)
{
}

//============================================================================

PayloadBuffer
write(const SequenceParameterSet& sps, const AttributeParameterSet& aps)
{
  PayloadBuffer buf(PayloadType::kAttributeParameterSet);
  auto bs = makeBitWriter(std::back_inserter(buf));

  bs.writeUn(4, aps.aps_attr_parameter_set_id);
  bs.writeUn(4, aps.aps_seq_parameter_set_id);
  bs.writeUe(aps.attr_encoding);

  bs.writeUe(aps.init_qp_minus4);
  bs.writeSe(aps.aps_chroma_qp_offset);
  bs.write(aps.aps_slice_qp_deltas_present_flag);

  if (aps.attr_encoding == AttributeEncoding::kRAHTransform) {
    bs.write(aps.rahtPredParams.prediction_enabled_flag);
    if (aps.rahtPredParams.prediction_enabled_flag) {
      bs.writeUe(aps.rahtPredParams.prediction_threshold0);
      bs.writeUe(aps.rahtPredParams.prediction_threshold1);
      bs.write(aps.rahtPredParams.enable_inter_prediction);
      bs.write(aps.rahtPredParams.raht_enable_inter_intra_layer_RDO);
      if(aps.rahtPredParams.enable_inter_prediction) {
        bs.writeUe(aps.rahtPredParams.mode_level);
        bs.writeUe(aps.rahtPredParams.upper_mode_level);
      } else {
        bs.writeUe(aps.rahtPredParams.intra_mode_level);
      }
    }
    bs.write(aps.rahtPredParams.integer_haar_enable_flag);
    bs.write(aps.rahtPredParams.cross_chroma_component_prediction_flag);
  }

  if (aps.attr_encoding == AttributeEncoding::kRaw)
    bs.write(aps.raw_attr_variable_len_flag);

  // NB: bitstreams conforming to the first edition must set
  // aps_extension_flag equal to 0.
  bool aps_extension_flag = sps.profile.isDraftProfile();
  bs.write(aps_extension_flag);
  if (aps_extension_flag) {
    bs.write(aps.attrInterPredictionEnabled);
    if (aps.attrInterPredictionEnabled) {
      bs.write(aps.dual_motion_field_flag);
      if (aps.dual_motion_field_flag) {
        bs.write(aps.mcap_to_rec_geom_flag);
        bs.writeUe(aps.motion.motion_block_size);
        bs.writeUe(aps.motion.motion_min_pu_size);
      }
    }

    if (
      aps.attr_encoding == AttributeEncoding::kRAHTransform
      && aps.rahtPredParams.prediction_enabled_flag) {
      bs.write(aps.rahtPredParams.prediction_skip1_flag);
      bs.write(aps.rahtPredParams.subnode_prediction_enabled_flag);
      if (aps.rahtPredParams.enable_inter_prediction) {
        bs.write(aps.rahtPredParams.enable_average_prediction);
        if (aps.rahtPredParams.enable_average_prediction) {
          bs.writeUe(aps.rahtPredParams.upper_mode_level_for_average_prediction);
          bs.writeUe(aps.rahtPredParams.lower_mode_level_for_average_prediction);
        }
      }
      if (aps.rahtPredParams.subnode_prediction_enabled_flag) 
        for (int i = 0; i < 5; i++) 
          bs.writeUe(aps.rahtPredParams.prediction_weights[i]);
    }
  }

  bs.byteAlign();

  return buf;
}

//----------------------------------------------------------------------------

AttributeParameterSet
parseAps(const PayloadBuffer& buf)
{
  AttributeParameterSet aps;
  assert(buf.type == PayloadType::kAttributeParameterSet);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUn(4, &aps.aps_attr_parameter_set_id);
  bs.readUn(4, &aps.aps_seq_parameter_set_id);
  bs.readUe(&aps.attr_encoding);

  bs.readUe(&aps.init_qp_minus4);
  bs.readSe(&aps.aps_chroma_qp_offset);
  bs.read(&aps.aps_slice_qp_deltas_present_flag);

  aps.rahtPredParams.raht_enable_inter_intra_layer_RDO = false;
  if (aps.attr_encoding == AttributeEncoding::kRAHTransform) {
    aps.rahtPredParams.intra_mode_level = 0;
    aps.rahtPredParams.mode_level = 0;
    aps.rahtPredParams.upper_mode_level = 0;
    bs.read(&aps.rahtPredParams.prediction_enabled_flag);
    if (aps.rahtPredParams.prediction_enabled_flag) {
      bs.readUe(&aps.rahtPredParams.prediction_threshold0);
      bs.readUe(&aps.rahtPredParams.prediction_threshold1);
      bs.read(&aps.rahtPredParams.enable_inter_prediction);
      bs.read(&aps.rahtPredParams.raht_enable_inter_intra_layer_RDO);
      if(aps.rahtPredParams.enable_inter_prediction) {
        bs.readUe(&aps.rahtPredParams.mode_level);
        bs.readUe(&aps.rahtPredParams.upper_mode_level);
      } else {
        bs.readUe(&aps.rahtPredParams.intra_mode_level);
      }
    }
    bs.read(&aps.rahtPredParams.integer_haar_enable_flag);
    bs.read(&aps.rahtPredParams.cross_chroma_component_prediction_flag);
  }

  if (aps.attr_encoding == AttributeEncoding::kRaw)
    bs.read(&aps.raw_attr_variable_len_flag);

  bool aps_extension_flag = bs.read();
  aps.rahtPredParams.prediction_skip1_flag = false;
  aps.rahtPredParams.subnode_prediction_enabled_flag = false;
  aps.attrInterPredictionEnabled = false;
  aps.dual_motion_field_flag = false;
  aps.mcap_to_rec_geom_flag = false;
  aps.rahtPredParams.enable_average_prediction = false;
  aps.rahtPredParams.upper_mode_level_for_average_prediction = 0;
  aps.rahtPredParams.lower_mode_level_for_average_prediction = 0;

  if (aps_extension_flag) {
    bs.read(&aps.attrInterPredictionEnabled);
    if (aps.attrInterPredictionEnabled) {
      bs.read(&aps.dual_motion_field_flag);
      if (aps.dual_motion_field_flag) {
        bs.read(&aps.mcap_to_rec_geom_flag);
        bs.readUe(&aps.motion.motion_block_size);
        bs.readUe(&aps.motion.motion_min_pu_size);
      }
    }

    if (
      aps.attr_encoding == AttributeEncoding::kRAHTransform
      && aps.rahtPredParams.prediction_enabled_flag) {
      bs.read(&aps.rahtPredParams.prediction_skip1_flag);
      bs.read(&aps.rahtPredParams.subnode_prediction_enabled_flag);
      if (aps.rahtPredParams.enable_inter_prediction) {
        bs.read(&aps.rahtPredParams.enable_average_prediction);
        if (aps.rahtPredParams.enable_average_prediction) {
          bs.readUe(&aps.rahtPredParams.upper_mode_level_for_average_prediction);
          bs.readUe(&aps.rahtPredParams.lower_mode_level_for_average_prediction);
        }
      }
      if (aps.rahtPredParams.subnode_prediction_enabled_flag) {
        aps.rahtPredParams.prediction_weights.resize(5);
        for (int i = 0; i < 5; i++)
          bs.readUe(&aps.rahtPredParams.prediction_weights[i]);
        aps.rahtPredParams.setPredictionWeights();
      }
    }
  }

  bs.byteAlign();

  return aps;
}

//----------------------------------------------------------------------------

void
convertXyzToStv(const SequenceParameterSet& sps, AttributeParameterSet* aps)
{
}

//============================================================================

void
write(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PayloadBuffer* buf)
{
  assert(buf->type == PayloadType::kGeometryBrick);
  auto bs = makeBitWriter(std::back_inserter(*buf));

  int gbh_reserved_zero_3bits = 0;
  bs.writeUn(4, gbh.geom_geom_parameter_set_id);
  bs.writeUn(3, gbh_reserved_zero_3bits);
  bs.writeUe(gbh.geom_slice_id);
  bs.writeUn(sps.slice_tag_bits, gbh.slice_tag);
  bs.writeUn(sps.frame_ctr_bits, gbh.frame_ctr_lsb);

  if (sps.entropy_continuation_enabled_flag) {
    bs.write(gbh.entropy_continuation_flag);
    if (gbh.entropy_continuation_flag)
      bs.writeUe(gbh.prev_slice_id);
  }

  int geomBoxLog2Scale = gbh.geomBoxLog2Scale(gps);
  auto geom_box_origin = toXyz(sps.geometry_axis_order, gbh.geomBoxOrigin);
  geom_box_origin.x() >>= geomBoxLog2Scale;
  geom_box_origin.y() >>= geomBoxLog2Scale;
  geom_box_origin.z() >>= geomBoxLog2Scale;

  if (gps.geom_box_log2_scale_present_flag)
    bs.writeUe(gbh.geom_box_log2_scale);

  bs.writeUe(gbh.geom_box_origin_bits_minus1);
  if (auto originBits = gbh.geom_box_origin_bits_minus1 + 1) {
    bs.writeUn(originBits, geom_box_origin.x());
    bs.writeUn(originBits, geom_box_origin.y());
    bs.writeUn(originBits, geom_box_origin.z());
  }

  int tree_depth_minus1 = gbh.tree_depth_minus1();
  bs.writeUe(tree_depth_minus1);
  if (gps.qtbt_enabled_flag)
    for (int i = 0; i <= tree_depth_minus1; i++)
      bs.writeUn(3, gbh.tree_lvl_coded_axis_list[i]);

  bs.writeUe(gbh.geom_stream_cnt_minus1);

  if (gps.geom_scaling_enabled_flag) {
    bs.writeSe(gbh.geom_slice_qp_offset);
  }

  if (gps.trisoup_enabled_flag) {
    bs.writeUe(gbh.trisoup_node_size);
    bs.writeUe(gbh.trisoup_QP);
    bs.write(gbh.trisoup_centroid_vertex_residual_flag);
    if( gbh.trisoup_centroid_vertex_residual_flag ){
      bs.write(gbh.trisoup_face_vertex_flag);
    }
    bs.write(gbh.trisoup_halo_flag);
    bs.write(gbh.trisoup_vertex_merge_flag);
    bs.writeUe(gbh.trisoup_thickness);

    int _write_bits = 0;
    if( gps.non_cubic_node_start_edge ){
      bs.writeUe( gbh.slice_bb_pos_bits );
      if( gbh.slice_bb_pos_bits > 0 ) {
        bs.writeUe( gbh.slice_bb_pos_log2_scale );
        _write_bits = gbh.slice_bb_pos_bits;
        bs.writeUn( _write_bits, gbh.slice_bb_pos[0] );
        bs.writeUn( _write_bits, gbh.slice_bb_pos[1] );
        bs.writeUn( _write_bits, gbh.slice_bb_pos[2] );
      }
    }
    if( gps.non_cubic_node_end_edge ){
      bs.writeUe( gbh.slice_bb_width_bits );
      if( gbh.slice_bb_width_bits > 0 ) {
        bs.writeUe( gbh.slice_bb_width_log2_scale );
        _write_bits = gbh.slice_bb_width_bits;
        bs.writeUn( _write_bits, gbh.slice_bb_width[0] );
        bs.writeUn( _write_bits, gbh.slice_bb_width[1] );
        bs.writeUn( _write_bits, gbh.slice_bb_width[2] );
      }
    }
  }

  if (gps.interPredictionEnabledFlag)
    bs.write(gbh.interPredictionEnabledFlag);
  else
    assert(!gbh.interPredictionEnabledFlag);

  bs.byteAlign();
}

//----------------------------------------------------------------------------

GeometryBrickHeader
parseGbh(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const PayloadBuffer& buf,
  int* bytesReadHead,
  int* bytesReadFoot)
{
  GeometryBrickHeader gbh;
  assert(buf.type == PayloadType::kGeometryBrick);
  auto bs = makeBitReader(buf.begin(), buf.end());

  int gbh_reserved_zero_3bits;
  bs.readUn(4, &gbh.geom_geom_parameter_set_id);
  bs.readUn(3, &gbh_reserved_zero_3bits);
  bs.readUe(&gbh.geom_slice_id);
  bs.readUn(sps.slice_tag_bits, &gbh.slice_tag);
  bs.readUn(sps.frame_ctr_bits, &gbh.frame_ctr_lsb);

  gbh.entropy_continuation_flag = false;
  if (sps.entropy_continuation_enabled_flag) {
    bs.read(&gbh.entropy_continuation_flag);
    if (gbh.entropy_continuation_flag)
      bs.readUe(&gbh.prev_slice_id);
  }

  if (gps.geom_box_log2_scale_present_flag)
    bs.readUe(&gbh.geom_box_log2_scale);

  Vec3<int> geom_box_origin;
  bs.readUe(&gbh.geom_box_origin_bits_minus1);
  if (auto originBits = gbh.geom_box_origin_bits_minus1 + 1) {
    bs.readUn(originBits, &geom_box_origin.x());
    bs.readUn(originBits, &geom_box_origin.y());
    bs.readUn(originBits, &geom_box_origin.z());
  }
  gbh.geomBoxOrigin = fromXyz(sps.geometry_axis_order, geom_box_origin);
  gbh.geomBoxOrigin *= 1 << gbh.geomBoxLog2Scale(gps);

  gbh.geom_stream_cnt_minus1 = 0;

  int tree_depth_minus1;
  bs.readUe(&tree_depth_minus1);

  gbh.tree_lvl_coded_axis_list.resize(tree_depth_minus1 + 1, 7);
  if (gps.qtbt_enabled_flag)
    for (int i = 0; i <= tree_depth_minus1; i++)
      bs.readUn(3, &gbh.tree_lvl_coded_axis_list[i]);

  bs.readUe(&gbh.geom_stream_cnt_minus1);

  gbh.geom_slice_qp_offset = 0;
  if (gps.geom_scaling_enabled_flag) {
    bs.readSe(&gbh.geom_slice_qp_offset);
  }

  if (gps.trisoup_enabled_flag) {
    bs.readUe(&gbh.trisoup_node_size);
    bs.readUe(&gbh.trisoup_QP);
    bs.read(&gbh.trisoup_centroid_vertex_residual_flag);
    if( gbh.trisoup_centroid_vertex_residual_flag ){
      bs.read(&gbh.trisoup_face_vertex_flag);
    }
    bs.read(&gbh.trisoup_halo_flag);
    bs.read(&gbh.trisoup_vertex_merge_flag);
    bs.readUe(&gbh.trisoup_thickness);

    gbh.slice_bb_pos_bits = 0;
    gbh.slice_bb_pos_log2_scale = 0;
    gbh.slice_bb_pos = 0;
    gbh.slice_bb_width_bits = 0;
    gbh.slice_bb_width_log2_scale = 0;
    gbh.slice_bb_width = 0;
    int _read_bits = 0;
    if( gps.non_cubic_node_start_edge ){
      bs.readUe( &gbh.slice_bb_pos_bits );
      if( gbh.slice_bb_pos_bits > 0 ) {
        bs.readUe( &gbh.slice_bb_pos_log2_scale );
        _read_bits = gbh.slice_bb_pos_bits;
        bs.readUn( _read_bits, &gbh.slice_bb_pos[0] );
        bs.readUn( _read_bits, &gbh.slice_bb_pos[1] );
        bs.readUn( _read_bits, &gbh.slice_bb_pos[2] );
      }
    }
    if( gps.non_cubic_node_end_edge ){
      bs.readUe( &gbh.slice_bb_width_bits );
      if( gbh.slice_bb_width_bits > 0 ) {
        bs.readUe( &gbh.slice_bb_width_log2_scale );
        _read_bits = gbh.slice_bb_width_bits;
        bs.readUn( _read_bits, &gbh.slice_bb_width[0] );
        bs.readUn( _read_bits, &gbh.slice_bb_width[1] );
        bs.readUn( _read_bits, &gbh.slice_bb_width[2] );
      }
    }
  }

  if (gps.interPredictionEnabledFlag)
    bs.read(&gbh.interPredictionEnabledFlag);
  else
    gbh.interPredictionEnabledFlag = false;

  bs.byteAlign();

  if (bytesReadHead)
    *bytesReadHead = int(std::distance(buf.begin(), bs.pos()));

  // To avoid having to make separate calls, the footer is parsed here
  gbh.footer = parseGbf(gps, gbh, buf, bytesReadFoot);

  return gbh;
}

//----------------------------------------------------------------------------

GeometryBrickHeader
parseGbhIds(const PayloadBuffer& buf)
{
  GeometryBrickHeader gbh;
  assert(buf.type == PayloadType::kGeometryBrick);
  auto bs = makeBitReader(buf.begin(), buf.end());

  int gbh_reserved_zero_3bits;
  bs.readUn(4, &gbh.geom_geom_parameter_set_id);
  bs.readUn(3, &gbh_reserved_zero_3bits);
  bs.readUe(&gbh.geom_slice_id);
  // NB: to decode slice_tag requires sps activation

  /* NB: this function only decodes ids at the start of the header. */
  /* NB: do not attempt to parse any further */

  return gbh;
}

//============================================================================

void
write(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& /* gbh */,
  const GeometryBrickFooter& gbf,
  PayloadBuffer* buf)
{
  assert(buf->type == PayloadType::kGeometryBrick);
  auto bs = makeBitWriter(std::back_inserter(*buf));

  // NB: if modifying this footer, it is essential that the decoder can
  // either decode backwards, or seek to the start.

  if (gps.octree_point_count_list_present_flag) {
    for (int i = 0; i < gbf.octree_lvl_num_points_minus1.size(); i++)
      bs.writeUn(24, gbf.octree_lvl_num_points_minus1[i]);
  }

  bs.writeUn(24, gbf.geom_num_points_minus1);
}

//----------------------------------------------------------------------------

GeometryBrickFooter
parseGbf(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const PayloadBuffer& buf,
  int* bytesRead)
{
  GeometryBrickFooter gbf;
  assert(buf.type == PayloadType::kGeometryBrick);

  // todo(df): this would be simpler to parse if it were written in reverse
  auto bufStart = buf.end() - 3; /* geom_num_points_minus1 */
  if (gps.octree_point_count_list_present_flag)
    bufStart -= gbh.tree_depth_minus1() * 3;

  if (bytesRead)
    *bytesRead = int(std::distance(bufStart, buf.end()));

  auto bs = makeBitReader(bufStart, buf.end());

  if (gps.octree_point_count_list_present_flag) {
    auto tree_depth_minus1 = gbh.tree_depth_minus1();
    gbf.octree_lvl_num_points_minus1.resize(tree_depth_minus1);
    for (int i = 0; i < tree_depth_minus1; i++)
      bs.readUn(24, &gbf.octree_lvl_num_points_minus1[i]);
  }

  bs.readUn(24, &gbf.geom_num_points_minus1);

  return gbf;
}

//============================================================================

void
write(
  const SequenceParameterSet& sps,
  const AttributeParameterSet& aps,
  const AttributeBrickHeader& abh,
  PayloadBuffer* buf)
{
  assert(buf->type == PayloadType::kAttributeBrick);
  auto bs = makeBitWriter(std::back_inserter(*buf));

  int abh_reserved_zero_3bits = 0;
  bs.writeUn(4, abh.attr_attr_parameter_set_id);
  bs.writeUn(3, abh_reserved_zero_3bits);
  bs.writeUe(abh.attr_sps_attr_idx);
  bs.writeUe(abh.attr_geom_slice_id);

  assert(abh.attr_sps_attr_idx < sps.attributeSets.size());

  if (aps.aps_slice_qp_deltas_present_flag) {
    bs.writeSe(abh.attr_qp_delta_luma);
    bs.writeSe(abh.attr_qp_delta_chroma);
  }

  bool attr_layer_qp_present_flag = !abh.attr_layer_qp_delta_luma.empty();
  bs.write(attr_layer_qp_present_flag);
  if (attr_layer_qp_present_flag) {
    int attr_num_qp_layers_minus1 = abh.attr_num_qp_layers_minus1();
    bs.writeUe(attr_num_qp_layers_minus1);
    for (int i = 0; i <= attr_num_qp_layers_minus1; i++) {
      bs.writeSe(abh.attr_layer_qp_delta_luma[i]);
      bs.writeSe(abh.attr_layer_qp_delta_chroma[i]);
    }
  }

  int attr_num_regions = abh.qpRegions.size();
  bs.writeUe(attr_num_regions);

  if (attr_num_regions)
    bs.writeUe(abh.attr_region_bits_minus1);

  for (int i = 0; i < attr_num_regions; i++) {
    // NB: only one region is currently permitted.
    auto& region = abh.qpRegions[i];

    auto attr_region_origin =
      toXyz(sps.geometry_axis_order, region.regionOrigin);

    auto attr_region_whd_minus1 =
      toXyz(sps.geometry_axis_order, region.regionSize - 1);

    auto regionBits = abh.attr_region_bits_minus1 + 1;
    bs.writeUn(regionBits, attr_region_origin.x());
    bs.writeUn(regionBits, attr_region_origin.y());
    bs.writeUn(regionBits, attr_region_origin.z());
    bs.writeUn(regionBits, attr_region_whd_minus1.x());
    bs.writeUn(regionBits, attr_region_whd_minus1.y());
    bs.writeUn(regionBits, attr_region_whd_minus1.z());
    bs.writeSe(region.attr_region_qp_offset[0]);
    if (sps.attributeSets[abh.attr_sps_attr_idx].attr_num_dimensions_minus1)
      bs.writeSe(region.attr_region_qp_offset[1]);
  }

  bs.write(abh.disableAttrInterPred);

  bs.byteAlign();
}

//----------------------------------------------------------------------------

AttributeBrickHeader
parseAbhIds(const PayloadBuffer& buf)
{
  AttributeBrickHeader abh;
  assert(buf.type == PayloadType::kAttributeBrick);
  auto bs = makeBitReader(buf.begin(), buf.end());

  int abh_reserved_zero_3bits;
  bs.readUn(4, &abh.attr_attr_parameter_set_id);
  bs.readUn(3, &abh_reserved_zero_3bits);
  bs.readUe(&abh.attr_sps_attr_idx);
  bs.readUe(&abh.attr_geom_slice_id);

  /* NB: this function only decodes ids at the start of the header. */
  /* NB: do not attempt to parse any further */

  return abh;
}

//----------------------------------------------------------------------------

AttributeBrickHeader
parseAbh(
  const SequenceParameterSet& sps,
  const AttributeParameterSet& aps,
  const PayloadBuffer& buf,
  int* bytesRead)
{
  AttributeBrickHeader abh;
  assert(buf.type == PayloadType::kAttributeBrick);
  auto bs = makeBitReader(buf.begin(), buf.end());

  int abh_reserved_zero_3bits;
  bs.readUn(4, &abh.attr_attr_parameter_set_id);
  bs.readUn(3, &abh_reserved_zero_3bits);
  bs.readUe(&abh.attr_sps_attr_idx);
  bs.readUe(&abh.attr_geom_slice_id);

  assert(abh.attr_sps_attr_idx < sps.attributeSets.size());

  if (aps.aps_slice_qp_deltas_present_flag) {
    bs.readSe(&abh.attr_qp_delta_luma);
    bs.readSe(&abh.attr_qp_delta_chroma);
  }

  bool attr_layer_qp_present_flag;
  bs.read(&attr_layer_qp_present_flag);
  if (attr_layer_qp_present_flag) {
    int attr_num_qp_layers_minus1;
    bs.readUe(&attr_num_qp_layers_minus1);
    abh.attr_layer_qp_delta_luma.resize(attr_num_qp_layers_minus1 + 1);
    abh.attr_layer_qp_delta_chroma.resize(attr_num_qp_layers_minus1 + 1);
    for (int i = 0; i <= attr_num_qp_layers_minus1; i++) {
      bs.readSe(&abh.attr_layer_qp_delta_luma[i]);
      bs.readSe(&abh.attr_layer_qp_delta_chroma[i]);
    }
  }

  // NB: Number of regions is restricted in this version of specification.
  int attr_num_regions;
  bs.readUe(&attr_num_regions);
  assert(attr_num_regions <= 1);

  if (attr_num_regions)
    bs.readUe(&abh.attr_region_bits_minus1);

  abh.qpRegions.resize(attr_num_regions);
  for (int i = 0; i < attr_num_regions; i++) {
    auto regionBits = abh.attr_region_bits_minus1 + 1;
    auto& region = abh.qpRegions[i];

    Vec3<int> attr_region_origin;
    bs.readUn(regionBits, &attr_region_origin.x());
    bs.readUn(regionBits, &attr_region_origin.y());
    bs.readUn(regionBits, &attr_region_origin.z());
    region.regionOrigin = fromXyz(sps.geometry_axis_order, attr_region_origin);

    Vec3<int> attr_region_whd_minus1;
    bs.readUn(regionBits, &attr_region_whd_minus1.x());
    bs.readUn(regionBits, &attr_region_whd_minus1.y());
    bs.readUn(regionBits, &attr_region_whd_minus1.z());
    region.regionSize =
      fromXyz(sps.geometry_axis_order, attr_region_whd_minus1 + 1);

    bs.readSe(&region.attr_region_qp_offset[0]);
    if (sps.attributeSets[abh.attr_sps_attr_idx].attr_num_dimensions_minus1)
      bs.readSe(&region.attr_region_qp_offset[1]);
  }

  bs.read(&abh.disableAttrInterPred);

  bs.byteAlign();

  if (bytesRead)
    *bytesRead = int(std::distance(buf.begin(), bs.pos()));

  return abh;
}

//============================================================================

ConstantAttributeDataUnit
parseConstantAttribute(
  const SequenceParameterSet& sps, const PayloadBuffer& buf)
{
  ConstantAttributeDataUnit cadu;
  assert(buf.type == PayloadType::kConstantAttribute);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUn(4, &cadu.constattr_attr_parameter_set_id);
  bs.readUe(&cadu.constattr_sps_attr_idx);
  bs.readUe(&cadu.constattr_geom_slice_id);

  // todo(df): check bounds
  const auto& attrDesc = sps.attributeSets[cadu.constattr_sps_attr_idx];

  cadu.constattr_default_value.resize(attrDesc.attr_num_dimensions_minus1 + 1);
  bs.readUn(attrDesc.bitdepth, &cadu.constattr_default_value[0]);
  for (int k = 1; k <= attrDesc.attr_num_dimensions_minus1; k++)
    bs.readUn(attrDesc.bitdepth, &cadu.constattr_default_value[k]);

  return cadu;
}

//============================================================================

void
write(
  const SequenceParameterSet& sps,
  const FrameBoundaryMarker& fbm,
  PayloadBuffer* buf)
{
  assert(buf->type == PayloadType::kFrameBoundaryMarker);
  auto bs = makeBitWriter(std::back_inserter(*buf));

  int fbdu_frame_ctr_lsb_bits = sps.frame_ctr_bits;
  bs.writeUn(5, fbdu_frame_ctr_lsb_bits);
  bs.writeUn(fbdu_frame_ctr_lsb_bits, fbm.fbdu_frame_ctr_lsb);
  bs.byteAlign();
}

//----------------------------------------------------------------------------

FrameBoundaryMarker
parseFrameBoundaryMarker(const PayloadBuffer& buf)
{
  FrameBoundaryMarker fbm;
  assert(buf.type == PayloadType::kFrameBoundaryMarker);
  auto bs = makeBitReader(buf.begin(), buf.end());

  int fbdu_frame_ctr_lsb_bits;
  bs.readUn(5, &fbdu_frame_ctr_lsb_bits);
  bs.readUn(fbdu_frame_ctr_lsb_bits, &fbm.fbdu_frame_ctr_lsb);

  return fbm;
}

//============================================================================

PayloadBuffer
write(const SequenceParameterSet& sps, const TileInventory& inventory)
{
  PayloadBuffer buf(PayloadType::kTileInventory);
  auto bs = makeBitWriter(std::back_inserter(buf));

  bs.writeUn(4, inventory.ti_seq_parameter_set_id);

  bs.writeUn(5, inventory.ti_frame_ctr_bits);
  bs.writeUn(inventory.ti_frame_ctr_bits, inventory.ti_frame_ctr);

  int num_tiles = inventory.tiles.size();
  bs.writeUn(16, num_tiles);
  if (!num_tiles) {
    bs.byteAlign();
    return buf;
  }

  bs.writeUn(5, inventory.tile_id_bits);
  bs.writeUn(8, inventory.tile_origin_bits_minus1);
  bs.writeUn(8, inventory.tile_size_bits_minus1);

  for (const auto& entry : inventory.tiles) {
    bs.writeUn(inventory.tile_id_bits, entry.tile_id);

    auto tile_origin = toXyz(sps.geometry_axis_order, entry.tileOrigin);
    if (auto tileOriginBits = inventory.tile_origin_bits_minus1 + 1) {
      bs.writeSn(tileOriginBits, tile_origin.x());
      bs.writeSn(tileOriginBits, tile_origin.y());
      bs.writeSn(tileOriginBits, tile_origin.z());
    }

    auto tile_size_minus1 = toXyz(sps.geometry_axis_order, entry.tileSize) - 1;
    if (auto tileSizeBits = inventory.tile_size_bits_minus1 + 1) {
      bs.writeUn(tileSizeBits, tile_size_minus1.x());
      bs.writeUn(tileSizeBits, tile_size_minus1.y());
      bs.writeUn(tileSizeBits, tile_size_minus1.z());
    }
  }

  // NB: this is at the end of the inventory to aid fixed-width parsing
  bs.writeUe(inventory.ti_origin_bits_minus1);
  if (auto tiOriginBits = inventory.ti_origin_bits_minus1 + 1) {
    auto ti_origin_xyz = toXyz(sps.geometry_axis_order, inventory.origin);
    bs.writeSn(tiOriginBits, ti_origin_xyz.x());
    bs.writeSn(tiOriginBits, ti_origin_xyz.y());
    bs.writeSn(tiOriginBits, ti_origin_xyz.z());
  }

  int ti_origin_log2_scale = 0;
  bs.writeUe(ti_origin_log2_scale);

  bs.byteAlign();

  return buf;
}

//----------------------------------------------------------------------------

TileInventory
parseTileInventory(const PayloadBuffer& buf)
{
  TileInventory inventory;
  assert(buf.type == PayloadType::kTileInventory);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUn(4, &inventory.ti_seq_parameter_set_id);

  bs.readUn(5, &inventory.ti_frame_ctr_bits);
  bs.readUn(inventory.ti_frame_ctr_bits, &inventory.ti_frame_ctr);

  int num_tiles;
  bs.readUn(16, &num_tiles);
  if (!num_tiles) {
    bs.byteAlign();
    return inventory;
  }

  bs.readUn(5, &inventory.tile_id_bits);
  bs.readUn(8, &inventory.tile_origin_bits_minus1);
  bs.readUn(8, &inventory.tile_size_bits_minus1);

  for (int i = 0; i < num_tiles; i++) {
    int tile_id = i;
    bs.readUn(inventory.tile_id_bits, &tile_id);

    Vec3<int> tile_origin;
    if (auto tileOriginBits = inventory.tile_origin_bits_minus1 + 1) {
      bs.readSn(tileOriginBits, &tile_origin.x());
      bs.readSn(tileOriginBits, &tile_origin.y());
      bs.readSn(tileOriginBits, &tile_origin.z());
    }

    Vec3<int> tile_size_minus1;
    if (auto tileSizeBits = inventory.tile_size_bits_minus1 + 1) {
      bs.readUn(tileSizeBits, &tile_size_minus1.x());
      bs.readUn(tileSizeBits, &tile_size_minus1.y());
      bs.readUn(tileSizeBits, &tile_size_minus1.z());
    }
    // NB: this is in XYZ axis order until the inventory is converted to STV
    TileInventory::Entry entry;
    entry.tile_id = tile_id;
    entry.tileOrigin = tile_origin;
    entry.tileSize = tile_size_minus1 + 1;
    inventory.tiles.push_back(entry);
  }

  Vec3<int> ti_origin_xyz;
  bs.readUe(&inventory.ti_origin_bits_minus1);
  if (auto tiOriginBits = inventory.ti_origin_bits_minus1 + 1) {
    bs.readSn(tiOriginBits, &ti_origin_xyz.x());
    bs.readSn(tiOriginBits, &ti_origin_xyz.y());
    bs.readSn(tiOriginBits, &ti_origin_xyz.z());
  }

  int ti_origin_log2_scale;
  bs.readUe(&ti_origin_log2_scale);
  ti_origin_xyz *= 1 << ti_origin_log2_scale;

  // NB: this is in XYZ axis order until converted to STV
  inventory.origin = ti_origin_xyz;

  bs.byteAlign();

  return inventory;
}

//----------------------------------------------------------------------------

void
convertXyzToStv(const SequenceParameterSet& sps, TileInventory* inventory)
{
  for (auto& tile : inventory->tiles) {
    tile.tileOrigin = fromXyz(sps.geometry_axis_order, tile.tileOrigin);
    tile.tileSize = fromXyz(sps.geometry_axis_order, tile.tileSize);
  }
}

//============================================================================

PayloadBuffer
write(
  const SequenceParameterSet& sps,
  const AttributeParamInventoryHdr& inv,
  const AttributeParameters& params)
{
  PayloadBuffer buf(PayloadType::kGeneralizedAttrParamInventory);
  auto bs = makeBitWriter(std::back_inserter(buf));

  assert(inv.attr_param_seq_parameter_set_id == sps.sps_seq_parameter_set_id);
  bs.writeUn(4, inv.attr_param_seq_parameter_set_id);
  int attr_param_frame_ctr_lsb_bits = sps.frame_ctr_bits;
  bs.writeUn(5, attr_param_frame_ctr_lsb_bits);
  bs.writeUn(attr_param_frame_ctr_lsb_bits, inv.attr_param_frame_ctr_lsb);
  bs.writeUe(inv.attr_param_sps_attr_idx);

  assert(inv.attr_param_sps_attr_idx < sps.attributeSets.size());
  auto& desc = sps.attributeSets[inv.attr_param_sps_attr_idx];
  writeAttributeParameters(desc, bs, params);

  return buf;
}

//----------------------------------------------------------------------------

AttributeParamInventoryHdr
parseAttrParamInventoryHdr(const PayloadBuffer& buf)
{
  AttributeParamInventoryHdr inv;
  assert(buf.type == PayloadType::kGeneralizedAttrParamInventory);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUn(4, &inv.attr_param_seq_parameter_set_id);
  int attr_param_frame_ctr_lsb_bits;
  bs.readUn(5, &attr_param_frame_ctr_lsb_bits);
  bs.readUn(attr_param_frame_ctr_lsb_bits, &inv.attr_param_frame_ctr_lsb);
  bs.readUe(&inv.attr_param_sps_attr_idx);

  return inv;
}

//----------------------------------------------------------------------------

AttributeParameters&
parseAttrParamInventory(
  const AttributeDescription& attr,
  const PayloadBuffer& buf,
  AttributeParameters& params)
{
  AttributeParamInventoryHdr inv;
  assert(buf.type == PayloadType::kGeneralizedAttrParamInventory);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUn(4, &inv.attr_param_seq_parameter_set_id);
  int attr_param_frame_ctr_lsb_bits;
  bs.readUn(5, &attr_param_frame_ctr_lsb_bits);
  bs.readUn(attr_param_frame_ctr_lsb_bits, &inv.attr_param_frame_ctr_lsb);
  bs.readUe(&inv.attr_param_sps_attr_idx);

  int num_attr_parameters;
  bs.readUe(&num_attr_parameters);
  bs.byteAlign();
  for (auto i = 0; i < num_attr_parameters; i++)
    parseAttributeParameter(attr, bs, params);

  return params;
}

//============================================================================

PayloadBuffer
write(const UserData& ud)
{
  PayloadBuffer buf(PayloadType::kUserData);
  auto bs = makeBitWriter(std::back_inserter(buf));

  writeOid(bs, ud.user_data_oid);
  // todo(df): write a blob of userdata

  return buf;
}

//----------------------------------------------------------------------------

UserData
parseUserData(const PayloadBuffer& buf)
{
  UserData ud;
  assert(buf.type == PayloadType::kUserData);
  auto bs = makeBitReader(buf.begin(), buf.end());

  readOid(bs, &ud.user_data_oid);
  // todo(df): read the blob of userdata

  return ud;
}

//============================================================================
// Helpers for Global scaling

SequenceParameterSet::GlobalScale::GlobalScale(Rational x)
{
  // x must be representable using fixed point
  if (popcntGt1(uint32_t(x.denominator)))
    throw std::runtime_error("cannot convert Rational to GlobalScale");

  // Direct conversion
  denominatorLog2 = ilog2(uint32_t(x.numerator));
  numeratorModDenominator = x.numerator - (1 << denominatorLog2);
  numeratorMulLog2 = denominatorLog2 - ilog2(uint32_t(x.denominator));

  // Simplify: remove powers of two from numeratorModDenominator
  while (!(numeratorModDenominator & 1) && denominatorLog2) {
    numeratorModDenominator >>= 1;
    denominatorLog2--;
  }
}

//----------------------------------------------------------------------------

SequenceParameterSet::GlobalScale::operator Rational() const
{
  int numeratorPreMul = ((1 << denominatorLog2) + numeratorModDenominator);

  // Simplify 2^numeratorMulLog2 / 2^denominatorLog2
  int numeratorS = std::max(0, numeratorMulLog2 - denominatorLog2);
  int denominatorS = denominatorLog2 - (numeratorMulLog2 - numeratorS);

  // Simplify numeratorPreMul / 2^denominatorLog2
  while (!(numeratorPreMul & 1) && denominatorS) {
    numeratorPreMul >>= 1;
    denominatorS--;
  }

  return Rational(numeratorPreMul << numeratorS, 1 << denominatorS);
}

//============================================================================

}  // namespace pcc
