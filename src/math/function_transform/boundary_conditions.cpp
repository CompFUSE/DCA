// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements boundary_conditions.hpp.

#include "dca/math/function_transform/boundary_conditions.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

std::string to_str(BOUNDARY_CONDITIONS BC) {
  switch (BC) {
    case INTERVAL:
      return "INTERVAL";

    case PERIODIC:
      return "PERIODIC";

    case ANTIPERIODIC:
      return "ANTIPERIODIC";

    default:
      return "NOT DEFINED";
  }
}

}  // transform
}  // math
}  // dca
