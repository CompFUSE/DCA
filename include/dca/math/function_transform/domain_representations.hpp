// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the different types of domain representations.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_DOMAIN_REPRESENTATIONS_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_DOMAIN_REPRESENTATIONS_HPP

#include <string>

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

enum DOMAIN_REPRESENTATIONS { DISCRETE, CONTINUOUS, EXPANSION };

std::string to_str(DOMAIN_REPRESENTATIONS DR);

}  // transform
}  // math
}  // dca

#endif  //  DCA_MATH_FUNCTION_TRANSFORM_DOMAIN_REPRESENTATIONS_HPP
