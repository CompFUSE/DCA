// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
