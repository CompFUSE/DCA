// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements element_spacings.hpp.

#include "dca/math/function_transform/element_spacings.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

std::string to_str(ELEMENT_SPACINGS IS) {
  switch (IS) {
    case EQUIDISTANT:
      return "EQUIDISTANT";

    case NONEQUIDISTANT:
      return "NONEQUIDISTANT";

    default:
      return "NOT DEFINED";
  }
}

}  // transform
}  // math
}  // dca
