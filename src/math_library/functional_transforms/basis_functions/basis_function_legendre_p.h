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

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_BASIS_FUNCTION_LEGENDRE_P_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_BASIS_FUNCTION_LEGENDRE_P_H

#include <cassert>
#include "math_library/functional_transforms/basis_functions/basis_function_template.h"
#include "math_library/functional_transforms/typedefs.hpp"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

template <typename lh_dmn_type, BASIS_EXPANSIONS BS_LHS, typename rh_dmn_type>
class basis_function<lh_dmn_type, BS_LHS, rh_dmn_type, LEGENDRE_P> {
public:
  typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
  typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;

  typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
  typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;

  typedef typename lh_spec_dmn_type::element_type lh_element_type;
  typedef typename rh_spec_dmn_type::element_type rh_element_type;

  typedef lh_scalar_type f_scalar_type;

public:
  static f_scalar_type execute(int i, int j) {
    return execute(lh_dmn_type::get_elements()[i], rh_dmn_type::get_elements()[j]);
  }

  static f_scalar_type execute(lh_element_type& lh_elem, rh_element_type& rh_elem) {
    lh_scalar_type delta = lh_dmn_type::get_volume();
    lh_scalar_type min = lh_dmn_type::get_min()[0];

    return compute(2 * (lh_elem - min) / delta - 1, rh_elem);
  }

private:
  inline static f_scalar_type compute(lh_scalar_type t, int l) {
    assert(t > -1. - 1.e-6 and t < 1. + 1.e-6);

    f_scalar_type f_t = 0;
    lh_scalar_type l_v = l;

    switch (l) {
      case 0:
        f_t = 1;
        break;

      case 1:
        f_t = t;
        break;

      default:
        f_t += ((2 * l_v - 1) * t * compute(t, l - 1) - (l_v - 1) * compute(t, l - 2)) / l_v;
    }

    return f_t;
  }
};

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_BASIS_FUNCTION_LEGENDRE_P_H
