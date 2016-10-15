// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements domain_representations.hpp.

#include "dca/math/function_transform/domain_representations.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

std::string to_str(DOMAIN_REPRESENTATIONS DR) {
  switch (DR) {
    case DISCRETE:
      return "DISCRETE";

    case CONTINUOUS:
      return "CONTINUOUS";

    case EXPANSION:
      return "EXPANSION";

    default:
      return "NOT DEFINED";
  }
}

}  // transform
}  // math
}  // dca
