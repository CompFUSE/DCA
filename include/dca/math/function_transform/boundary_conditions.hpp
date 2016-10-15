// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the different types of boundary conditions.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_BOUNDARY_CONDTIONS_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_BOUNDARY_CONDTIONS_HPP

#include <string>

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

enum BOUNDARY_CONDITIONS { INTERVAL, PERIODIC, ANTIPERIODIC };

std::string to_str(BOUNDARY_CONDITIONS BC);

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_BOUNDARY_CONDTIONS_HPP
