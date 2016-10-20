// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the different types of basis expansions.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_BASIS_EXPANSIONS_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_BASIS_EXPANSIONS_HPP

#include <string>

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

enum BASIS_EXPANSIONS {
  KRONECKER_DELTA,

  HARMONICS,
  COSINE,
  SINE,

  HERMITE_LINEAR_SPLINE,
  HERMITE_CUBIC_SPLINE,

  LEGENDRE_P,
  LEGENDRE_Q,
  LEGENDRE_LM
};

std::string to_str(BASIS_EXPANSIONS BS);

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_BASIS_EXPANSIONS_HPP
