// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
