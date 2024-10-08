/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2019, ISO/IEC
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

#include <memory>

#include "hls.h"
#include "PayloadBuffer.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "attr_tools.h"

namespace pcc {

//============================================================================

class AttributeContexts;
struct AttributeInterPredParams;

//============================================================================

class AttributeDecoderIntf {
public:
  virtual ~AttributeDecoderIntf();

  virtual void decode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    AttributeBrickHeader& abh,
    int geom_num_points_minus1,
    int minGeomNodeSizeLog2,
    const char* payload,
    size_t payloadLen,
    PCCPointSet3& pointCloud,
    AttributeInterPredParams& attrInterPredParams,
    attr::ModeDecoder& predDecoder
    ) = 0;

  virtual void decodeSlab(
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
    attr::ModeDecoder& predDecoder
  ) = 0;

  virtual void startDecode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    const AttributeBrickHeader& abh,
    const char* payload,
    size_t payloadLen,
    const AttributeContexts& ctxtMem
  ) = 0;

  virtual void finishDecode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    const AttributeBrickHeader& abh,
    AttributeContexts& ctxtMem
  ) = 0;
};

//----------------------------------------------------------------------------

std::unique_ptr<AttributeDecoderIntf> makeAttributeDecoder();

//============================================================================

class AttributeEncoderIntf {
public:
  virtual ~AttributeEncoderIntf();

  virtual void encode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& attr_aps,
    PCCPointSet3& pointCloud,
    PayloadBuffer* payload,
    AttributeInterPredParams &attrInterPredParams,
    attr::ModeEncoder& predEncoder
  ) = 0;

  virtual void encodeSlab(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& attr_aps,
    PCCPointSet3& slabPointCloud,
    PayloadBuffer* payload,
    AttributeInterPredParams& attrInterPredParams,
    attr::ModeEncoder& predEncoder
  ) = 0;

  virtual void startEncode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& attr_aps,
    const AttributeBrickHeader& abh,
    const AttributeContexts& ctxtMem,
    uint32_t pointCountInPointCloud
  ) = 0;

  virtual void finishEncode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& attr_aps,
    AttributeContexts& ctxtMem,
    PayloadBuffer* payload
  ) = 0;
};

//----------------------------------------------------------------------------

std::unique_ptr<AttributeEncoderIntf> makeAttributeEncoder();

//============================================================================

}  // namespace pcc
