// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
