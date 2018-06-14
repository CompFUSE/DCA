// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides a function to determine the sign of number.
// Reference:
// http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c

#ifndef DCA_MATH_UTIL_SGN_HPP
#define DCA_MATH_UTIL_SGN_HPP

namespace dca {
namespace math {
namespace util {
// dca::math::util::

template <typename T>
int sgn(const T val) {
  return (T(0) < val) - (val < T(0));
}

}  // util
}  // math
}  // dca

#endif  // DCA_MATH_UTIL_SGN_HPP
