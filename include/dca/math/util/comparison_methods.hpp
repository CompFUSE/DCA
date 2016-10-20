// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// The file provides comparing methods used e.g. for sorting containers.

#ifndef DCA_MATH_UTIL_COMPARISON_METHODS_HPP
#define DCA_MATH_UTIL_COMPARISON_METHODS_HPP

#include <cmath>
#include <complex>
#include <utility>

namespace dca {
namespace math {
namespace util {
// dca::math::util::

template <typename T1, typename T2>
bool pairLess(const std::pair<T1, T2>& x, const std::pair<T1, T2>& y) {
  return x.first < y.first;
}

template <typename T1, typename T2>
bool ComplexPairLess(const std::pair<std::complex<T1>, T2>& x,
                     const std::pair<std::complex<T1>, T2>& y) {
  return x.first.real() < y.first.real();
}

template <typename T1, typename T2>
bool pairAbsGreater(const std::pair<T1, T2>& x, const std::pair<T1, T2>& y) {
  return std::abs(x.first) > std::abs(y.first);
}

template <typename T1, typename T2>
bool susceptibilityPairGreater(const std::pair<std::complex<T1>, T2>& x,
                               const std::pair<std::complex<T1>, T2>& y) {
  return std::abs(x.first - 1.) > std::abs(y.first - 1.);
}

}  // util
}  // math
}  // dca

#endif  // DCA_MATH_UTIL_COMPARISON_METHODS_HPP
