// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: John Biddiscombe (john.biddiscombe@cscs.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides utility functions that operate on types.

#ifndef DCA_UTIL_TYPE_UTILS_HPP
#define DCA_UTIL_TYPE_UTILS_HPP

#include <complex>
#ifdef DCA_HAVE_CUDA
#include <cuComplex.h>
#endif

//#include "dca/util/type_list.hpp"

namespace dca {
namespace util {
// dca::util::

// Throws an assertion if the two types do not match.
// Extends std::is_same<> by forcing the compiler to print the types in the error message which
// helps with debugging.
template <typename T1, typename T2>
struct assert_same {
  assert_same() {
    static_assert(std::is_same<T1, T2>::value, "Types must be equal.");
  }
  static_assert(std::is_same<T1, T2>::value, "Types must be equal.");
};

// Determine if a type is complex or not.
template <class T>
struct IsComplex : std::false_type {};
template <class T>
struct IsComplex<std::complex<T>> : std::true_type {};

namespace {
template <class T>
struct RealImpl {
  using type = T;
};
template <class T>
struct RealImpl<std::complex<T>> {
  using type = T;
};
#ifdef DCA_HAVE_CUDA
template <>
struct RealImpl<cuComplex> {
  using type = float;
};
template <>
struct RealImpl<cuDoubleComplex> {
  using type = double;
};
#endif  // DCA_HAVE_CUDA
}  // namespace

template <class T>
using Real = typename RealImpl<T>::type;

namespace {
template <class T>
struct ComplexImpl {
  using type = std::complex<T>;
};
template <class T>
struct ComplexImpl<std::complex<T>> {
  using type = std::complex<T>;
};
}  // namespace

template <class T>
using Complex = typename ComplexImpl<T>::type;

namespace {
template <bool single_precision, bool complex>
struct ScalarImpl;
template <>
struct ScalarImpl<true, true> {
  using type = std::complex<float>;
};
template <>
struct ScalarImpl<false, true> {
  using type = std::complex<double>;
};
template <>
struct ScalarImpl<true, false> {
  using type = float;
};
template <>
struct ScalarImpl<false, false> {
  using type = double;
};
}  // namespace

template <bool single_precision, bool complex>
using Scalar = typename ScalarImpl<single_precision, complex>::type;

}  // namespace util
}  // namespace dca

#endif  // DCA_UTIL_TYPE_UTILS_HPP
