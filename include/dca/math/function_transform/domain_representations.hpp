// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
