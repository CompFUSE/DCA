// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides utilities for vector manipulation.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_VECTOR_OPERATIONS_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_VECTOR_OPERATIONS_HPP

#include <vector>

namespace dca {
namespace math {
namespace transform {
namespace details {
// dca::math::transform::details::

template <typename Scalar>
Scalar dot_prod(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
  assert(x.size() == y.size());
  Scalar result = 0;
  for (std::size_t l = 0; l < x.size(); l++)
    result += x[l] * y[l];
  return result;
}

template <typename Scalar>
Scalar dot_prod(Scalar x, Scalar y) {
  return x * y;
}

template <typename Scalar>
std::vector<Scalar> sum(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
  assert(x.size() == y.size());
  std::vector<Scalar> result(x.size());
  for (std::size_t l = 0; l < x.size(); l++)
    result[l] = x[l] + y[l];
  return result;
}

template <typename Scalar>
std::vector<Scalar> diff(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
  assert(x.size() == y.size());
  std::vector<Scalar> result(x.size());
  for (std::size_t l = 0; l < x.size(); l++)
    result[l] = x[l] - y[l];
  return result;
}

}  // details
}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_VECTOR_OPERATIONS_HPP
