// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W.  Doak (doakpw@ornl.gov)
//

#ifndef DCA_UTIL_TYPE_HELP_HPP
#define DCA_UTIL_TYPE_HELP_HPP

#include <type_traits>
#include <complex>
#include <tuple>

// namespace dca {
// namespace util {

// #ifndef DCA_HAVE_GPU
// template <typename Real>
// struct Real2CudaComplex {
//   using type = void;
// };
// #endif
// }  // namespace util
// }  // namespace dca

#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu_complex.h"
#include "dca/linalg/util/gpu_type_mapping.hpp"
#endif

namespace dca {
namespace util {

#ifdef DCA_HAVE_GPU
template <typename Real>
using CudaComplex = typename Real2CudaComplex<Real>::type;
#endif

template <typename T>
struct IsComplex_t : public std::false_type {};
template <typename T>
struct IsComplex_t<std::complex<T>> : public std::true_type {};

// template <typename T>
// using IsComplex_t< CudaComplex<T>> : public std::true_type {};

template <typename T>
using IsComplex = std::enable_if_t<IsComplex_t<T>::value, bool>;

template <typename T>
using IsReal = std::enable_if_t<std::is_floating_point<T>::value, bool>;

template <typename T, typename = bool>
struct RealAlias_impl {};

template <typename T>
struct RealAlias_impl<T, IsReal<T>> {
  using value_type = T;
};

#ifdef DCA_HAVE_GPU
template <>
struct RealAlias_impl<cuComplex, bool> {
  using value_type = float;
};

template <>
struct RealAlias_impl<cuDoubleComplex, bool> {
  using value_type = double;
};
#endif

template <typename T>
struct RealAlias_impl<T, IsComplex<T>> {
  using value_type = typename T::value_type;
};

/** If you have a function templated on a value that can be real or complex
 *   and you need to get the base Real type if its complex or just the real.
 *
 *  If you try to do this on anything but a fp or a std::complex<fp> you will
 *  get a compilation error.
 */
template <typename T>
using RealAlias = typename RealAlias_impl<T>::value_type;

template <typename T, typename = bool>
struct ComplexAlias_impl {};

template <typename T>
struct ComplexAlias_impl<T, IsReal<T>> {
  using value_type = std::complex<T>;
};

template <typename T>
struct ComplexAlias_impl<T, IsComplex<T>> {
  using value_type = T;
};

#ifdef DCA_HAVE_GPU
template <>
struct ComplexAlias_impl<double2> {
  using value_type = cuDoubleComplex;
};

template <>
struct ComplexAlias_impl<float2> {
  using value_type = cuComplex;
};

template <typename T>
struct ComplexAlias_impl<T, IsCudaComplex<T>> {
  using value_type = CudaComplex<std::remove_pointer<T*>>;
};

template <typename T, typename = bool>
struct CudaScalar_impl {};

template <typename T>
struct CudaScalar_impl<T, IsReal<T>> {
  using value_type = T;
};

template <typename T>
struct CudaScalar_impl<T, IsComplex<T>> {
  using value_type = dca::util::CudaComplex<RealAlias<T>>;
};

template <typename T>
using CudaScalar = typename CudaScalar_impl<T>::value_type;

template <typename T>
struct CudaScalarStruct {
  typename CudaScalar_impl<T>::value_type value;
};
#endif

template <typename T>
using ComplexAlias = typename ComplexAlias_impl<T>::value_type;

template <typename REAL, bool complex>
struct ScalarSelect {
  using type = REAL;
};

template <typename REAL>
struct ScalarSelect<REAL, false> {
  using type = REAL;
};

template <typename REAL>
struct ScalarSelect<REAL, true> {
  using type = std::complex<REAL>;
};

template <typename T, typename = bool>
struct TheOne;
template <typename T, typename = bool>
struct TheZero;

template <typename T>
struct TheOne<T, IsReal<T>> {
  static constexpr T value = 1.0;
};

template <typename T>
struct TheZero<T, IsReal<T>> {
  static constexpr T value = 0.0;
};


#ifdef DCA_HAVE_GPU
template <typename T>
struct TheOne<T, IsCudaComplex<T>> {
  static constexpr T value{1.0, 0.0};
};

template <typename T>
struct TheZero<T, IsCudaComplex<T>> {
  static constexpr T value{0.0, 0.0};
};

template <typename T>
std::enable_if_t<IsCudaComplex_t<T>::value, void> makeOne(T& one) {
  one = T{1.0, 0.0};
}

template <typename T>
std::enable_if_t<IsCudaComplex_t<T>::value, void> makeZero(T& zero) {
  zero = T{0.0, 0.0};
}

#endif

template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, void> makeOne(T& one) {
  one = 1.0;
}

template <typename T>
std::enable_if_t<IsComplex_t<T>::value, void> makeOne(T& one) {
  one = T{1.0, 0.0};
}

template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, void> makeZero(T& zero) {
  zero = 1.0;
}

template <typename T>
std::enable_if_t<IsComplex_t<T>::value, void> makeZero(T& zero) {
  zero = T{1.0, 0.0};
}

template <typename T>
struct TheOne<T, IsComplex<T>> {
  static constexpr T value{1.0, 0.0};
};

  
template <typename T>
struct TheZero<T, IsComplex<T>> {
  static constexpr T value = {0.0, 0.0};
};

template <typename T, typename T2>
auto makeMaybe(
    T2 t2, typename std::enable_if_t<IsComplex_t<T>::value || std::is_floating_point<T>::value>* = 0) {
  return T(t2);
}

#ifdef DCA_HAVE_GPU
/** to handle making double2 and float2 values from Real in a generic way.
 *  static cast required to deal with possibility of narrowing conversion from literal expressed as double.
 */
template <typename T, typename T2>
auto makeMaybe(T2 t2, typename std::enable_if_t<IsCudaComplex_t<T>::value>* = 0) {
  using Real = RealAlias<T>;
  return T{static_cast<Real>(t2), 0.0};
}

template <typename T>
inline auto GPUTypeConversion(T var, typename std::enable_if_t<IsCudaComplex_t<T>::value>* = 0) {
  return HOSTTypeMap<T>{var.x, var.y};
}

template <typename T>
inline auto GPUTypeConversion(
    T var, typename std::enable_if_t<IsComplex_t<T>::value && (!IsCudaComplex_t<T>::value)>* = 0) {
  return CUDATypeMap<T>{var.real(), var.imag()};
}

template <typename T>
inline auto GPUTypeConversion(T var,
                              typename std::enable_if_t<std::is_floating_point<T>::value>* = 0) {
  return var;
}

template <typename T>
inline auto GPUTypeConversion(T var, typename std::enable_if_t<std::is_integral<T>::value>* = 0) {
  return var;
}
#endif

template <typename ARRAY, std::size_t... SIZE>
auto Array2Tuple_impl(const ARRAY& a, std::index_sequence<SIZE...>) {
  return std::make_tuple(a[SIZE]...);
}

template <typename T, std::size_t N, typename Indices = std::make_index_sequence<N>>
auto Array2Tuple(const std::array<T, N>& a) {
  return Array2Tuple_impl(a, Indices{});
}

template <typename T, typename = bool>
struct SignType_impl {};

template <typename T>
struct SignType_impl<T, IsComplex<T>> {
  using type = T;
};

template <typename T>
struct SignType_impl<T, IsReal<T>> {
  using type = std::int8_t;
};

template <typename T>
using SignType = typename SignType_impl<T>::type;

}  // namespace util
}  // namespace dca

#endif
