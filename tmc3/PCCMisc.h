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

#include <array>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <iterator>
#include <limits>
#include <utility>
#include <string>
#include <vector>

#if _MSC_VER
#  define DEPRECATED_MSVC __declspec(deprecated)
#  define DEPRECATED
#else
#  define DEPRECATED_MSVC
#  define DEPRECATED __attribute__((deprecated))
#endif

#if _MSC_VER && !defined(__attribute__)
#  define __attribute__(...)
#endif

namespace pcc {
const uint32_t PCC_UNDEFINED_INDEX = -1;

enum PCCEndianness
{
  PCC_BIG_ENDIAN = 0,
  PCC_LITTLE_ENDIAN = 1
};

inline PCCEndianness
PCCSystemEndianness()
{
  uint32_t num = 1;
  return (*(reinterpret_cast<char*>(&num)) == 1) ? PCC_LITTLE_ENDIAN
                                                 : PCC_BIG_ENDIAN;
}

//---------------------------------------------------------------------------
// Replace any occurence of %d with formatted number.  The %d format
// specifier may use the formatting conventions of snprintf().
std::string expandNum(const std::string& src, int num);

//---------------------------------------------------------------------------
// Population count -- return the number of bits set in @x.
//
inline int
popcnt(uint32_t x)
{
  x = x - ((x >> 1) & 0x55555555u);
  x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
  return ((x + (x >> 4) & 0xF0F0F0Fu) * 0x1010101u) >> 24;
}

//---------------------------------------------------------------------------
// Population count -- return the number of bits set in @x.
//
inline int
popcnt(uint8_t x)
{
  uint32_t val = x * 0x08040201u;
  val >>= 3;
  val &= 0x11111111u;
  val *= 0x11111111u;
  return val >> 28;
}

//---------------------------------------------------------------------------
// Test if population count is greater than 1.
// Returns non-zero if true.
//
inline uint32_t
popcntGt1(uint32_t x)
{
  return x & (x - 1);
}

//---------------------------------------------------------------------------
// Round @x up to next power of two.
//
inline uint32_t
ceilpow2(uint32_t x)
{
  x--;
  x = x | (x >> 1);
  x = x | (x >> 2);
  x = x | (x >> 4);
  x = x | (x >> 8);
  x = x | (x >> 16);
  return x + 1;
}

//---------------------------------------------------------------------------
// Round @x up to next power of two.
//
inline uint64_t
ceilpow2(uint64_t x)
{
  x--;
  x = x | (x >> 1);
  x = x | (x >> 2);
  x = x | (x >> 4);
  x = x | (x >> 8);
  x = x | (x >> 16);
  x = x | (x >> 32);
  return x + 1;
}

//---------------------------------------------------------------------------
// Compute \left\floor \text{log}_2(x) \right\floor.
// NB: ilog2(0) = -1.

inline int
ilog2(uint32_t x)
{
#if defined __GNUC__
  return x ? 31 - __builtin_clz(x) : -1;
#elif defined _MSC_VER
  unsigned long r = -1;
  if (x)
    _BitScanReverse(&r, x);
  return r;
#else
  x = ceilpow2(x + 1) - 1;
  return popcnt(x) - 1;
#endif
}

//---------------------------------------------------------------------------
// Compute \left\floor \text{log}_2(x) \right\floor.
// NB: ilog2(0) = -1.

inline int
ilog2(uint64_t x)
{
#if defined __GNUC__
  return x ? 63 - __builtin_clzll(x) : -1;
#elif defined _MSC_VER
  unsigned long r = -1;
  if (x)
    _BitScanReverse64(&r, x);
  return r;
#else
  x = ceilpow2(x + 1) - 1;
  return popcnt(uint32_t(x >> 32)) + popcnt(uint32_t(x)) - 1;
#endif
}

//---------------------------------------------------------------------------
// Compute \left\ceil \text{log}_2(x) \right\ceil.
// NB: ceillog2(0) = 32.

inline int
ceillog2(uint32_t x)
{
  return ilog2(x - 1) + 1;
}

//---------------------------------------------------------------------------
// The number of bits required to represent x.
// NB: x must be >= 0.
// NB: numBits(0) = 1.

inline int
numBits(int x)
{
  return std::max(0, ilog2(uint32_t(x))) + 1;
}

//---------------------------------------------------------------------------
// The number of bits required to represent constexpr X.
// NB: X must be >= 0.
// NB: NumBits<0>::val = 0.
template <uint64_t X>
struct NumBits {
  enum { val = NumBits<(X >> 1)>::val + 1 };
};
// terminal case
template <>
struct NumBits<1ULL> {
  enum { val = 1 };
};
// special case
template <>
struct NumBits<0ULL> {
  enum { val = 1 };
};

//---------------------------------------------------------------------------
// rotate left by n bits.
// If n is negative, the effect is to rotate right.
// NB: Signed types are treated as unsigned.

template<typename T>
T
rotateLeft(T val, int n)
{
  n &= sizeof(T) * 8 - 1;
  using unsigned_t = typename std::make_unsigned<T>::type;
  return (unsigned_t(val) << n) | (unsigned_t(val) >> (8 * sizeof(val) - n));
}

//---------------------------------------------------------------------------
// rotate right by n bits.
// If n is negative, the effect is to rotate left.
// NB: Signed types are treated as unsigned.

template<typename T>
T
rotateRight(T val, int n)
{
  n &= sizeof(T) * 8 - 1;
  using unsigned_t = typename std::make_unsigned<T>::type;
  return (unsigned_t(val) >> n) | (unsigned_t(val) << (8 * sizeof(val) - n));
}

//---------------------------------------------------------------------------
// Compute an approximation of \left\floor \sqrt{x} \right\floor

uint32_t isqrt(uint64_t x) __attribute__((const));
uint32_t fastIsqrt(uint64_t x) __attribute__((const));

//---------------------------------------------------------------------------
// Compute an approximation of reciprocal sqrt

int64_t irsqrt(uint64_t a64) __attribute__((const));
int64_t fastIrsqrt(int64_t a64) __attribute__((const));

//---------------------------------------------------------------------------
// Compute an approximation of atan2

int iatan2(int y, int x);

//---------------------------------------------------------------------------
// Reduce fixed point by one level of fixed point precision (N bits)
// may be used for restoring original fixed point precision after a
// multiplication by fixed point value;
// or for rounding fixed point value to integer value.
template <int N, typename T>
inline T fpReduce(T val)
{
  return (val >= 0) * (((1 << N - 1) + val) >> N)
    - (val < 0) * (((1 << N - 1) - val) >> N);
}

//---------------------------------------------------------------------------
// Expand by one level of fixed point precision (N bits)
// may be used for expanding precision before a division by a fixed point
// value;
// or for converting an integer value to a fixed point value.
template <int N, typename T>
inline T fpExpand(T val)
{
  return (val >= 0) * (val << N) - (val < 0) * ((-val) << N);
}

#if 0 // Not used currently
//---------------------------------------------------------------------------
// Estimate log2 function (with decimal part) for fixed points numbers
//
// InFPP: fixed point precision of the input
// OutFPP: fixed point precision of the output
// K: number of iterations for the approximation
//
// In general (i.e. without taking into account possible range of input values),
// OutFPP should not be more than 63-6 (63-5 should always work if InFPP = 32).
// InFPP may have any value in range [0..64] and similar approximation error
//   should be obtained for any value.
//
// Method found in a discussion on StackExchange
//   https://math.stackexchange.com/questions/1706939/approximation-log-2x
//
template <int InFPP, int OutFPP=InFPP, unsigned K=2, bool Faster=true>
int64_t fpLog2(uint64_t x) {
  if (!x) return std::numeric_limits<int64_t>::min();

  static_assert(!Faster || K < 5, "invk[] does not support K >= 5 currently");

  const int64_t log2 = ilog2(x);
  const int64_t shift = log2 - InFPP;

  int64_t res = shift << OutFPP;
  // compute decimal part
  if (K > 0) {
    // try taking maximum fixed point precision
    int64_t xD;
    if (log2 >= 32)
      xD = x >> log2 - 32;
    else
      xD = x << 32 - log2;
    // highest bit is lost, x is set equal to {(x offset by shift) - 1}
    // with 64 bits FP precision
    x <<= 64 - log2;
    constexpr uint64_t oneD = 1ULL << 32;
    uint64_t t = x / (xD + oneD);

    // Note: accum could get more precision than 32 bits if needed
    uint64_t accum = t;
    uint64_t tPowered = t;

    static constexpr uint64_t invk[10] = {
      0, 0, 0, 341, 0, 205, 0, 146, 0, 114 };

    for (unsigned k = 3; k < 2*K+1; k+=2) {
      // Note: we could also increase tPowered precision along iterations
      tPowered *= t;
      tPowered >>= 32;
      tPowered *= t;
      if (Faster) {
        tPowered >>= 10;
        accum += tPowered * invk[k] >> 32;
        tPowered >>= 22;
      } else {
        accum += tPowered / k >> 32;
        tPowered >>= 32;
      }
    }
    // 2 / log(2) ~= 2.88539 @ 32 fpp bits
    constexpr uint64_t twoOverLogTwo = 0x2e2a8eca5ULL;
    // Note: in theory it should not overflow
    //assert(accum <= 0x58b90bfbULL);
    res += (twoOverLogTwo * accum) >> 64 - OutFPP;
  }
  return res;
}
#endif

//---------------------------------------------------------------------------

#if 0 // Not used currently
template <unsigned FPP, typename T>
double fpToDouble(T x)
{
  if (FPP < 64)
    return double(x) / double(T(1) << FPP);
  else
    return double(x) / std::pow(2., FPP);
}
#endif

#if 0 // Not used currently
template <unsigned FPP, typename T>
T fpFromDouble(double x)
{
  if (FPP < 64)
    return T(x * double(T(1) << FPP) + 0.5);
  else
    return T(x * std::pow(2., FPP) + 0.5);
}
#endif

//---------------------------------------------------------------------------
// Luted log2
static const int64_t log2P_LUT[130] = {
  917504, 458752, 393216, 354880, 327680, 306582, 289344, 274769, 262144,
  251008, 241046, 232035, 223808, 216240, 209233, 202710, 196608, 190876,
  185472, 180360, 175510, 170897, 166499, 162296, 158272, 154412, 150704,
  147136, 143697, 140379, 137174, 134074, 131072, 128163, 125340, 122599,
  119936, 117345, 114824, 112368, 109974, 107639, 105361, 103136, 100963,
  98838,  96760,  94726,  92736,  90786,  88876,  87004,  85168,  83367,
  81600,  79865,  78161,  76488,  74843,  73227,  71638,  70075,  68538,
  67025,  65536,  64070,  62627,  61205,  59804,  58424,  57063,  55722,
  54400,  53096,  51809,  50540,  49288,  48052,  46832,  45627,  44438,
  43264,  42103,  40957,  39825,  38706,  37600,  36507,  35427,  34358,
  33302,  32257,  31224,  30202,  29190,  28190,  27200,  26220,  25250,
  24290,  23340,  22399,  21468,  20546,  19632,  18727,  17831,  16943,
  16064,  15192,  14329,  13473,  12625,  11785,  10952,  10126,  9307,
  8496,   7691,   6893,   6102,   5317,   4539,   3767,   3002,   2242,
  1489,   742,    0,      0};


// entry one is 1 << InFPP; output log(0.5) is 1 << OutFPP
template <int InFPP = 16, int OutFPP = InFPP>
int64_t fpEntropyProbaLUT(int64_t x) {
  constexpr int shiftIn = 9 + InFPP - 16;

  constexpr int64_t one = 1 << shiftIn;
  int64_t idx = x >> shiftIn;
  int64_t lambda = x - (idx << shiftIn);
  int64_t out = (one - lambda) * log2P_LUT[idx] + lambda * log2P_LUT[idx + 1];

  constexpr int shiftOut = OutFPP - 16 - shiftIn;
  if (shiftOut >= 0)
    return out << shiftOut;
  else
    return out + (1 << -shiftOut - 1) >> -shiftOut;
}

//---------------------------------------------------------------------------
// Decrement the @axis-th dimension of 3D morton code @x.
//
inline int64_t
morton3dAxisDec(int64_t val, int axis)
{
  const int64_t mask0 = 0x9249249249249249llu << axis;
  return ((val & mask0) - 1 & mask0) | (val & ~mask0);
}

//---------------------------------------------------------------------------
// add the three dimentional addresses @a + @b;
//
inline uint64_t
morton3dAdd(uint64_t a, uint64_t b)
{
  constexpr uint64_t mask1 = 0x9249249249249249llu;
  constexpr uint64_t mask2 = mask1 << 1;
  constexpr uint64_t mask3 = mask1 << 2;

  return (a | ~mask1) + (b & mask1) & mask1
    | (a | ~mask2) + (b & mask2) & mask2
    | (a | ~mask3) + (b & mask3) & mask3;
}

//---------------------------------------------------------------------------
// Sort the elements in range [@first, @last) using a counting sort.
//
// The value of each element is determined by the function
// @value_of(RandomIt::value_type) and must be in the range [0, Radix).
//
// A supplied output array of @counts represents the histogram of values,
// and may be used to calculate the output position of each value span.
//
// NB: This is an in-place implementation and is not a stable sort.

template<class RandomIt, class ValueOp, std::size_t Radix>
void
countingSort(
  RandomIt first,
  RandomIt last,
  std::array<int, Radix>& counts,
  ValueOp value_of)
{
  // step 1: count each radix
  for (auto it = first; it != last; ++it) {
    counts[value_of(*it)]++;
  }

  // step 2: determine the output offsets
  std::array<RandomIt, Radix> ptrs = {{first}};
  for (int i = 1; i < Radix; i++) {
    ptrs[i] = std::next(ptrs[i - 1], counts[i - 1]);
  }

  // step 3: re-order, completing each radix in turn.
  RandomIt ptr_orig_last = first;
  for (int i = 0; i < Radix; i++) {
    std::advance(ptr_orig_last, counts[i]);
    while (ptrs[i] != ptr_orig_last) {
      int radix = value_of(*ptrs[i]);
      std::iter_swap(ptrs[i], ptrs[radix]);
      ++ptrs[radix];
    }
  }
}

#if 0 // Currently unused
//---------------------------------------------------------------------------

struct NoOp {
  template<typename... Args>
  void operator()(Args...)
  {}
};

#endif
#if 0 // Currently unused
//---------------------------------------------------------------------------

template<typename It, typename ValueOp, typename AccumOp>
void
radixSort8WithAccum(int maxValLog2, It begin, It end, ValueOp op, AccumOp acc)
{
  std::array<int, 8> counts = {};
  countingSort(begin, end, counts, [=](decltype(*begin)& it) {
    return op(maxValLog2, it);
  });

  acc(maxValLog2, counts);

  if (--maxValLog2 < 0)
    return;

  auto childBegin = begin;
  for (int i = 0; i < counts.size(); i++) {
    if (!counts[i])
      continue;
    auto childEnd = std::next(childBegin, counts[i]);
    radixSort8WithAccum(maxValLog2, childBegin, childEnd, op, acc);
    childBegin = childEnd;
  }
}

//---------------------------------------------------------------------------

template<typename It, typename ValueOp>
void
radixSort8(int maxValLog2, It begin, It end, ValueOp op)
{
  radixSort8WithAccum(maxValLog2, begin, end, op, NoOp());
}

#endif
//============================================================================
// A wrapper to reverse the iteration order of a range based for loop

template<typename T>
struct fwd_is_reverse_iterator {
  T& obj;
};

template<typename T>
auto
begin(fwd_is_reverse_iterator<T> w) -> decltype(w.obj.rbegin())
{
  return w.obj.rbegin();
}

template<typename T>
auto
end(fwd_is_reverse_iterator<T> w) -> decltype(w.obj.rend())
{
  return w.obj.rend();
}

template<typename T>
fwd_is_reverse_iterator<T>
inReverse(T&& obj)
{
  return {obj};
}

//============================================================================
// numeric tools

template<typename T>
T
gcd(T a, T b)
{
  while (b != 0) {
    T r = a % b;
    a = b;
    b = r;
  }
  return a;
}

//----------------------------------------------------------------------------

template<typename T>
T
lcm(const T& a, const T& b)
{
  return std::abs(a * b) / gcd(a, b);
}

//----------------------------------------------------------------------------

template<typename Iterator>
typename std::iterator_traits<Iterator>::value_type
lcm_all(Iterator it, Iterator it_end)
{
  typedef typename std::iterator_traits<Iterator>::value_type T;

  T res = T();

  if (it != it_end)
    res = *it++;

  while (it != it_end)
    res = lcm(res, *it++);

  return res;
}

//----------------------------------------------------------------------------

// apply value of differently typed parameters depending on a condition

template<bool cond>
struct conditional_value {};

template<>
struct conditional_value<true> {
  template<typename TA, typename TB>
  static const TA& from(const TA& a, const TB& b)
  { return a; }

  template<typename TA, typename TB>
  static TA& from(TA& a, TB& b)
  { return a; }

  template<typename TA, typename TB>
  static TA from(TA&& a, TB&& b)
  { return a; }
};

template<>
struct conditional_value<false> {
  template<typename TA, typename TB>
  static const TB& from(const TA& a, const TB& b)
  { return b; }

  template<typename TA, typename TB>
  static TB& from(TA& a, TB& b)
  { return b; }

  template<typename TA, typename TB>
  static TB from(TA&& a, TB&& b)
  { return b; }
};

//----------------------------------------------------------------------------

// conditional_value may require to compute both typed alternative values to
// pass them as arguments of conditional_value::value() function.
// conditional_call relies on INVOKE expression (e.g. lambda) to only compute
// the required value by conditionally selecting the INVOKE expression
// and using it.

template<bool cond>
struct conditional_call {
  /*
  template<typename A, typename B>
  static typename std::conditional<
    cond,
    typename std::result_of<A()>::type,
    typename std::result_of<B()>::type
  >::type from(A,B);
  */
};

template<>
struct conditional_call<true> {
  template<typename A, typename B>
  static typename std::result_of<A()>::type
  from(const A& a, const B& b)
  { return a(); }
};


template<>
struct conditional_call<false> {
  template<typename A, typename B>
  static typename std::result_of<B()>::type
  from(const A& a, const B& b)
  { return b(); }
};

//----------------------------------------------------------------------------

}  // namespace pcc
