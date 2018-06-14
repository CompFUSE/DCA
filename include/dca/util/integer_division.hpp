// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides a function that performs rounded-up integer division.

#ifndef DCA_UTIL_INTEGER_DIVISION_HPP
#define DCA_UTIL_INTEGER_DIVISION_HPP

#include <type_traits>

namespace dca {
namespace util {
// dca::util::

// Preconditions: a >= 0, and b > 0.
template <typename IntegerType>
constexpr typename std::enable_if<std::is_integral<IntegerType>::value, IntegerType>::type ceilDiv(
    IntegerType a, IntegerType b) {
  return (a + b - 1) / b;
}

}  // util
}  // dca

#endif  // DCA_UTIL_INTEGER_DIVISION_HPP
