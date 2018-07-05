// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements basis_expansions.hpp.

#include "dca/math/function_transform/basis_expansions.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

std::string to_str(BASIS_EXPANSIONS BS) {
  switch (BS) {
    case KRONECKER_DELTA:
      return "KRONECKER_DELTA";

    case HARMONICS:
      return "HARMONICS";

    case SINE:
      return "SINE";

    case COSINE:
      return "COSINE";

    case HERMITE_LINEAR_SPLINE:
      return "HERMITE_LINEAR_SPLINE";

    case HERMITE_CUBIC_SPLINE:
      return "HERMITE_CUBIC_SPLINE";

    case LEGENDRE_P:
      return "LEGENDRE_P";

    case LEGENDRE_Q:
      return "LEGENDRE_Q";

    case LEGENDRE_LM:
      return "LEGENDRE_LM";

    default:
      return "NOT DEFINED";
  }
}

}  // transform
}  // math
}  // dca
