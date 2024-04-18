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

#include <stdint.h>
#include <vector>

#include "entropy.h"
#include "hls.h"
#include "PCCTMC3Common.h"
#include "motionWip.h"

namespace pcc {

//============================================================================

std::vector<int> sortedPointCloud(
  const int attribCount,
  const PCCPointSet3& pointCloud,
  std::vector<int64_t>& mortonCode,
  std::vector<int>& attributes);

//---------------------------------------------------------------------------

void sortedPointCloud(
  const int attribCount,
  const PCCPointSet3& pointCloud,
  const std::vector<int>& indexOrd,
  std::vector<int>& attributes);

//---------------------------------------------------------------------------

struct MSOctree;
struct EncoderParams;

struct AttributeInterPredParams: InterPredParams {
  bool enableAttrInterPred;
#if 0
  int getPointCount() const { return referencePointCloud.getPointCount(); }
  void clear() { referencePointCloud.clear(); }
#endif
  bool hasLocalMotion() const { return compensatedPointCloud.getPointCount() > 0; }

  void findMotion(
    const EncoderParams* params,
    const EncodeMotionSearchParams& msParams,
    const ParameterSetMotion& mvPS,
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    PCCPointSet3& pointCloud
  );

  void encodeMotionAndBuildCompensated(
    const ParameterSetMotion& mvPS,
    EntropyEncoder& arithmeticEncoder,
    bool mcap_to_rec_geom_flag
  );

  void buildCompensatedSlab(
    const ParameterSetMotion& mvPS,
    bool mcap_to_rec_geom_flag
  );

  void prepareEncodeMotion(
    const ParameterSetMotion& mvPS,
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    PCCPointSet3& pointCloud
  );

  void prepareDecodeMotion(
    const ParameterSetMotion& mvPS,
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    PCCPointSet3& pointCloud
  );

  void decodeMotionAndBuildCompensated(
    const ParameterSetMotion& mvPS,
    EntropyDecoder& arithmeticDecoder,
    bool mcap_to_rec_geom_flag
  );

  void extractMotionForSlab(int startX, int thickness)
  {
    slabStartX = startX;
    slabThickness = thickness;
    dualMotion = MVField(
        mvField,
        point_t(startX, 0, 0),
        point_t(startX + thickness, INT32_MAX, INT32_MAX));
  }

  void copyMotion()
  {
    dualMotion = mvField;
  }
protected:
  int slabStartX = -1;
  int slabThickness = -1;
  MVField  dualMotion;
  MSOctree mSOctreeCurr;
  std::vector<std::pair<int/*dualMotion's puNodeIdx*/, int/*mSOctreeCurr's nodeIdx*/> > motionPUTrees;
};

//----------------------------------------------------------------------------

class AttributeContexts {
public:
  void reset();

protected:
  AdaptiveBitModel ctxRunLen[5];
  AdaptiveBitModel ctxLayerPred;
  AdaptiveBitModel ctxInterLayerPred;
  AdaptiveBitModel ctxCoeffGtN[2][7];
  AdaptiveBitModel ctxCoeffRemPrefix[2][12];
  AdaptiveBitModel ctxCoeffRemSuffix[2][12];
};

//----------------------------------------------------------------------------

inline void
AttributeContexts::reset()
{
  this->~AttributeContexts();
  new (this) AttributeContexts;
}

//============================================================================

}  // namespace pcc
