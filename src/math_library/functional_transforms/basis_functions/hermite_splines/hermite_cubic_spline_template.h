// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_HERMITE_SPLINES_HERMITE_CUBIC_SPLINE_TEMPLATE_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_HERMITE_SPLINES_HERMITE_CUBIC_SPLINE_TEMPLATE_H

#include "math_library/functional_transforms/typedefs.hpp"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

template <typename lh_dmn_type, typename rh_dmn_type, ELEMENT_SPACINGS ES, BOUNDARY_CONDITIONS BC,
          int DIMENSION>
struct hermite_cubic_spline {};

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_HERMITE_SPLINES_HERMITE_CUBIC_SPLINE_TEMPLATE_H
