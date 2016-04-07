//-*-C++-*-
#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_BASIS_FUNCTION_HERMITE_SPLINE_CUBLIC_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_BASIS_FUNCTION_HERMITE_SPLINE_CUBLIC_H

#include "math_library/typedefs.hpp"
#include "math_library/functional_transforms/basis_functions/basis_function_template.h"
#include "math_library/functional_transforms/basis_functions/hermite_splines/hermite_spline.h"
#include "math_library/functional_transforms/basis_functions/hermite_splines/hermite_cubic_spline_equidistant_periodic_nD.h"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

template <typename lhs_dmn_type, BASIS_EXPANSIONS BS_LHS, typename rhs_dmn_type>
class basis_function<lhs_dmn_type, BS_LHS, rhs_dmn_type, HERMITE_CUBIC_SPLINE> {
public:
  typedef typename lhs_dmn_type::dmn_specifications_type lhs_spec_dmn_type;
  typedef typename rhs_dmn_type::dmn_specifications_type rhs_spec_dmn_type;

  typedef typename lhs_spec_dmn_type::scalar_type lhs_scalar_type;
  typedef typename rhs_spec_dmn_type::scalar_type rhs_scalar_type;

  typedef typename lhs_spec_dmn_type::element_type lhs_element_type;
  typedef typename rhs_spec_dmn_type::element_type rhs_element_type;

  const static ELEMENT_SPACINGS ELEMENT_SPACING = rhs_spec_dmn_type::ELEMENT_SPACING;
  const static BOUNDARY_CONDITIONS BOUNDARY_CONDITION = rhs_spec_dmn_type::BOUNDARY_CONDITION;

  typedef rhs_scalar_type f_scalar_type;

public:
  static rhs_scalar_type execute(int i, int j) {
    return hermite_cubic_spline<lhs_dmn_type, rhs_dmn_type, ELEMENT_SPACING, BOUNDARY_CONDITION,
                                rhs_dmn_type::DIMENSION>::execute(i, j);
  }

  static rhs_scalar_type execute(lhs_element_type& x, rhs_element_type& y) {
    const static lhs_scalar_type a = -0.5;

    lhs_scalar_type d = rhs_dmn_type::get_elements()[1] - rhs_dmn_type::get_elements()[0];

    return hermite_spline::cubic(x, y, d, a);
  }
};

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_BASIS_FUNCTION_HERMITE_SPLINE_CUBLIC_H
