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

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "Attribute.h"
#include "PayloadBuffer.h"
#include "PCCMath.h"
#include "PCCPointSet.h"
#include "frame.h"
#include "geometry.h"
#include "geometry_params.h"
#include "hls.h"
#include "partitioning.h"
#include "pointset_processing.h"
#include "TMC3.h"
#include "attr_tools.h"
#include "AttributeCommon.h"


namespace pcc {

//============================================================================

struct EncodeMotionSearchParams {
  // window_size, max_prefix_bits and max_suffix_bits
  // are used for historical reasons, and kept for motion entropy estimation
  // TODO: Some cleanups should be made to avoid using them
  int window_size = 0;
  int max_prefix_bits = 0;
  int max_suffix_bits = 0;

  // Controls the use of an approximate or exact nearest neighbour search
  bool approximate_nn = false;

  // Number of iterations for the mean Nearest Neighbour intialization
  // of the starting motion vector
  int K = 1;

  // Starting search distance around starting motion vector
  int Amotion0 = 0;

  // Provide lambda value for Lagrangian Rate distorsion optimization
  double lambda = 0.;

  // Factor applied to the color difference of a neighbour point in the
  // computation of its distance to a current point in distorsion estimation
  double dgeom_color_factor = 0.;

  // A decimation factor on the number of points being considered during the
  // distorsion estimation is obtained by: MotionBlockSize >> decimate
  int decimate = 0;
};

//----------------------------------------------------------------------------

struct EncoderAttributeParams {
  // NB: this only makes sense for setting configurable parameters
  AttributeBrickHeader abh;

  // parameters for attributes motion search
  EncodeMotionSearchParams motion;
};

//----------------------------------------------------------------------------

struct EncoderParams {
  SequenceParameterSet sps;
  GeometryParameterSet gps;
  GeometryBrickHeader gbh;

  // NB: information about attributes is split between the SPS and the APS.
  //  => The SPS enumerates the attributes, the APS controls coding params.
  std::vector<AttributeParameterSet> aps;

  // Encoder specific parameters for attributes
  std::vector<EncoderAttributeParams> attr;

  // todo(df): this should go away
  std::map<std::string, int> attributeIdxMap;

  // Determine the sequence bonuding box using the first input frame
  bool autoSeqBbox;

  // Length of the source point cloud unit vectors.
  double srcUnitLength;

  // Scale factor used to define the coordinate system used for coding.
  // This is the coordinate system where slicing is performed.
  //  P_cod = P_src * codedGeomScale
  double codedGeomScale;

  // Scale factor used to define the sequence coordinate system.
  //  P_seq = P_src * seqGeomScale
  double seqGeomScale;

  // Scale factor used to define the external coordinate system.
  //  P_ext = P_src * extGeomScale
  double extGeomScale;

  // Number of fractional bits used in output position representation.
  int outputFpBits;

  // Encoder specific parameters for geometry
  OctreeEncOpts geom;

  // Encoder specific parameters for trisoup
  TrisoupEncOpts trisoup;

  // Parameters that control partitioning
  PartitionParams partition;

  std::vector<point_t> fixedSliceOrigin;

  // attribute recolouring parameters
  RecolourParams recolour;

  // per-slice trisoup node sizes
  std::vector<int> trisoupNodeSizes;

  // Enable enforcement of level limits (encoder will abort if exceeded)
  bool enforceLevelLimits;

  // Qp used for IDCM quantisation (used to derive HLS values)
  int idcmQp;

  // Period of random access points (managed by SequenceEncoder)
  int randomAccessPeriod;

  // local motion
  int motionPreset;

  // parameters for geometry motion search
  EncodeMotionSearchParams motion;

  // Encoder will use localized attributes even during the encoding
  bool localized_attributes_encoding;

  // localized attributes slab thickness
  int localized_attributes_slab_thickness;
};

//============================================================================

class PCCTMC3Encoder3 {
public:
  class Callbacks;

  PCCTMC3Encoder3();
  PCCTMC3Encoder3(const PCCTMC3Encoder3&) = delete;
  PCCTMC3Encoder3(PCCTMC3Encoder3&&) = default;
  PCCTMC3Encoder3& operator=(const PCCTMC3Encoder3& rhs) = delete;
  PCCTMC3Encoder3& operator=(PCCTMC3Encoder3&& rhs) = default;
  ~PCCTMC3Encoder3();

  int compress(
    const PCCPointSet3& inputPointCloud,
    EncoderParams* params,
    Callbacks*,
    CloudFrame* reconstructedCloud = nullptr);

  void compressPartition(
    const PCCPointSet3& inputPointCloud,
    const PCCPointSet3& originPartCloud,
    EncoderParams* params,
    Callbacks*,
    CloudFrame* reconstructedCloud = nullptr);

  static void deriveParameterSets(EncoderParams* params);
  static void fixupParameterSets(EncoderParams* params);
  void setInterForCurrPic(bool x) { _codeCurrFrameAsInter = x; }
  void deriveMotionParams(EncoderParams* params);

  void processNextSlabAttributes(PCCPointSet3& slabPointCloud, uint32_t xStartSlab, bool isLast);

private:
  void appendSlice(PCCPointSet3& cloud);

  void startEncodeGeometryBrick(const EncoderParams*);
  void encodeGeometryBrick(const EncoderParams*, PayloadBuffer* buf, AttributeInterPredParams& attrInterPredParams);

  SrcMappedPointSet quantization(const PCCPointSet3& src);

private:
  attr::ModeEncoder predCoder;
  PCCPointSet3 pointCloud;

  // Point positions in spherical coordinates of the current slice
  std::vector<point_t> _posSph;

  // Scale factor used to decimate the input point cloud.
  // Decimation is performed as if the input were scaled by
  //   Round(P_src * inputDecimationScale)
  // and duplicate points removed.
  // todo: expose this parameter?
  double _inputDecimationScale;

  // Scale factor that defines coding coordinate system
  double _srcToCodingScale;

  // Sequence origin in terms of coding coordinate system
  Vec3<int> _originInCodingCoords;

  // Position of the slice in the translated+scaled co-ordinate system.
  Vec3<int> _sliceOrigin;

  // Size of the current slice
  Vec3<int> _sliceBoxWhd;

  // The active parameter sets
  const SequenceParameterSet* _sps;
  const GeometryParameterSet* _gps;
  std::vector<const AttributeParameterSet*> _aps;

  // Cached copy of the curent _gbh (after encoding geometry)
  GeometryBrickHeader _gbh;

  // Indicates that this is the start of a new frame
  bool _firstSliceInFrame;

  // Current identifier of payloads with the same geometry
  int _sliceId;

  // Identifies the previous slice in bistream order
  int _prevSliceId;

  // Identifies the current tile
  int _tileId;

  // Current frame number.
  // NB: only the log2_max_frame_ctr LSBs are sampled for frame_ctr
  int _frameCounter;

  // Memorized context buffers
  std::unique_ptr<GeometryOctreeContexts> _ctxtMemOctreeGeom;
  std::vector<AttributeContexts> _ctxtMemAttrs;
  std::vector<int> _ctxtMemAttrSliceIds;
  // Code current picture as inter prediction
  bool _codeCurrFrameAsInter;

  AttributeInterPredParams attrInterPredParams;

  pcc::point_t minPos_ref;

  CloudFrame _refFrame;

  // For local encoding
  PayloadBuffer
    payload_geom;

  std::vector<PayloadBuffer>
    payload_attr;

  pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock>
    clock_user_geom;

  std::vector<pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock>>
    clock_user_attr;

  std::vector<AttributeBrickHeader>
    _abh;

  std::vector<decltype(makeAttributeEncoder())>
    attrEncoder;

  // LocalizedAttributes

  enum class Type {
    kNone = -1, // global decoder
    kLocalEncoder = 1,
    kGlobalEncoder = 2
  };

  Type type;

  uint32_t slabThickness;

  bool isInter;

  int currSlabIdx;

  // For encoder
  const PCCPointSet3* originPartCloud;
  point_t origin;
  double targetToSourceScaleFactor;
  Box3<int32_t> bBoxOrigin;
  const EncoderParams* params;

  // This is for global encoder
  std::vector<std::pair<int/*startx*/,int/*numPoints*/>> startXAndNumPointsPerSlab;
};

//----------------------------------------------------------------------------

class PCCTMC3Encoder3::Callbacks {
public:
  virtual void onOutputBuffer(const PayloadBuffer&) = 0;
  virtual void onPostRecolour(const PCCPointSet3&) = 0;
};

//============================================================================

}  // namespace pcc
