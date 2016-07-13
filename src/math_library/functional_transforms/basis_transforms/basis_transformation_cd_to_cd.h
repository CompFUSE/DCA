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

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_BASIS_TRANSFORMATION_CD_TO_CD_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_BASIS_TRANSFORMATION_CD_TO_CD_H

#include <cassert>
#include <string>

#include "comp_library/linalg/linalg.hpp"
#include "math_library/functional_transforms/basis_functions/basis_functions.hpp"
#include "math_library/functional_transforms/basis_transforms/basis_transformation_template.h"
#include "math_library/functional_transforms/basis_transforms/inner_product_domain.h"
#include "math_library/functional_transforms/typedefs.hpp"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

template <typename input_type, typename output_type>
class basis_transformation<input_type, CONTINUOUS, output_type, CONTINUOUS> {
public:
  typedef input_type rh_dmn_type;
  typedef output_type lh_dmn_type;

  typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
  typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;

  typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
  typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;

  typedef typename lh_spec_dmn_type::element_type lh_element_type;
  typedef typename rh_spec_dmn_type::element_type rh_element_type;

  typedef basis_function<lh_dmn_type, lh_spec_dmn_type::BASIS_EXPANSION, rh_dmn_type,
                         rh_spec_dmn_type::BASIS_EXPANSION>
      basis_function_type;

  typedef typename basis_function_type::f_scalar_type f_scalar_type;

  typedef LIN_ALG::matrix<f_scalar_type, LIN_ALG::CPU> matrix_type;

  typedef inner_product_domain<output_type> inner_product_dmn;

  typedef basis_transformation<input_type, CONTINUOUS, inner_product_dmn, CONTINUOUS> input_transform_type;
  typedef basis_transformation<output_type, CONTINUOUS, inner_product_dmn, CONTINUOUS> output_transform_type;

public:
  static bool& is_initialized() {
    static bool initialized = false;
    return initialized;
  }

  static std::string& get_name() {
    static std::string name = "basis-transformation";
    return name;
  }

  static matrix_type& get_transformation_matrix() {
    static matrix_type T;

    if (not is_initialized())
      initialize_transformation_matrix();

    return T;
  }

  static void initialize_transformation_matrix() {
    is_initialized() = true;

    int M = lh_dmn_type::get_size();
    int N = rh_dmn_type::get_size();

    assert(M > 0 and N > 0);

    matrix_type& T = get_transformation_matrix();

    T.reserve(std::pair<int, int>(M, N));

    for (int j = 0; j < N; j++)
      for (int i = 0; i < M; i++)
        T(i, j) = basis_function_type::execute(i, j);
  }
};

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_BASIS_TRANSFORMATION_CD_TO_CD_H
