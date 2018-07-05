// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
