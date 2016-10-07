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

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_TEMPLATE_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_TEMPLATE_H

#include "dca/function/function.hpp"
#include "math_library/functional_transforms/typedefs.hpp"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

template <typename type_input, DOMAIN_REPRESENTATIONS DMN_REP_LHS, typename type_output,
          DOMAIN_REPRESENTATIONS DMN_REP_RHS, int DMN_INDEX>
struct TRANSFORM_DOMAIN {
  const static bool VERBOSE = true;

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output);
};

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_TEMPLATE_H
