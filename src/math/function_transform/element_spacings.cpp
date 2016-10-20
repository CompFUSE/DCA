// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
