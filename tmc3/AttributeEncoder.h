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

#include <stdint.h>
#include <vector>

#include "Attribute.h"
#include "AttributeCommon.h"
#include "PayloadBuffer.h"
#include "PCCTMC3Common.h"
#include "hls.h"
#include "quantization.h"
#include "attr_tools.h"

namespace pcc {

//============================================================================
// Opaque definitions (Internal detail)

class PCCResidualsEncoder;
struct PCCResidualsEntropyEstimator;

//============================================================================

class AttributeEncoder : public AttributeEncoderIntf {
public:
  void encode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& attr_aps,
    PCCPointSet3& pointCloud,
    PayloadBuffer* payload,
    AttributeInterPredParams &attrInterPredParams,
    attr::ModeEncoder& predEncoder
  ) override;

  void encodeSlab(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& attr_aps,
    PCCPointSet3& slabPointCloud,
    PayloadBuffer* payload,
    AttributeInterPredParams& attrInterPredParams,
    attr::ModeEncoder& predEncoder
  ) override;

  void startEncode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& attr_aps,
    const AttributeBrickHeader& abh,
    const AttributeContexts& ctxtMem,
    uint32_t pointCountInPointCloud
  ) override;

  void finishEncode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& attr_aps,
    AttributeContexts& ctxtMem,
    PayloadBuffer* payload
  ) override;

protected:
  // todo(df): consider alternative encapsulation

  void encodeReflectancesTransformRaht(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    AttributeBrickHeader& abh,
    const QpSet& qpSet,
    PCCPointSet3& pointCloud,
    PCCResidualsEncoder& encoder,
    attr::ModeEncoder& predEncoder,
    const AttributeInterPredParams& attrInterPredParams);

  void encodeColorsTransformRaht(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    AttributeBrickHeader& abh,
    const QpSet& qpSet,
    PCCPointSet3& pointCloud,
    PCCResidualsEncoder& encoder,
    attr::ModeEncoder& predEncoder,
    const AttributeInterPredParams& attrInterPredParams);

  void encodeRAHTperBlock(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    AttributeBrickHeader& abh,
    const QpSet& qpSet,
    PCCPointSet3& pointCloud,
    PCCResidualsEncoder& encoder,
    attr::ModeEncoder& predEncoder,
    const AttributeInterPredParams& attrInterPredParams);

private:
  // The current attribute slice header
  AttributeBrickHeader _abh;

  // for local attributes
  std::unique_ptr<PCCResidualsEncoder> _pEncoder;
  QpSet _qpSet;
};

//============================================================================

} /* namespace pcc */
