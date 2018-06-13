// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
bool complexPairLess(const std::pair<std::complex<T1>, T2>& x,
                     const std::pair<std::complex<T1>, T2>& y) {
  return x.first.real() < y.first.real();
}

template <typename T1, typename T2>
bool pairAbsGreater(const std::pair<T1, T2>& x, const std::pair<T1, T2>& y) {
  return std::abs(x.first) > std::abs(y.first);
}

// Returns true, if x.first is closer to 1 than y.first, otherwise returns false.
// Used in BseLatticeSolver to sort the eigenvalues.
template <typename T1, typename T2>
bool susceptibilityPairLess(const std::pair<T1, T2>& x, const std::pair<T1, T2>& y) {
  return std::abs(x.first - 1.) < std::abs(y.first - 1.);
}

}  // util
}  // math
}  // dca

#endif  // DCA_MATH_UTIL_COMPARISON_METHODS_HPP
