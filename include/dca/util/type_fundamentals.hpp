// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W.  Doak (doakpw@ornl.gov)
//

#ifndef DCA_UTIL_TYPE_FUNDAMENTALS_HPP
#define DCA_UTIL_TYPE_FUNDAMENTALS_HPP

#include <complex>

namespace dca {
namespace util {
template <typename T, typename = bool>
struct TheOne;
template <typename T, typename = bool>
struct TheZero;

template <typename T>
using IsReal = std::enable_if_t<std::is_floating_point<T>::value, bool>;

template <typename T>
struct IsComplex_t : public std::false_type {};
template <typename T>
struct IsComplex_t<std::complex<T>> : public std::true_type {};

// template <typename T>
// using IsComplex_t< CudaComplex<T>> : public std::true_type {};

template <typename T>
using IsComplex = std::enable_if_t<IsComplex_t<T>::value, bool>;

template <typename T, typename = bool>
struct RealAlias_impl {};

template <typename T>
struct RealAlias_impl<T, IsReal<T>> {
  using value_type = T;
};

  
}  // namespace util
}  // namespace dca

#endif
