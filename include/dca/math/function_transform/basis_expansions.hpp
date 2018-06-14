// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
