// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_BASIS_TRANSFORMATION_TEMPLATE_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_BASIS_TRANSFORMATION_TEMPLATE_H

#include <cassert>
#include <string>

#include "comp_library/linalg/linalg.hpp"
#include "math_library/functional_transforms/basis_functions/basis_functions.hpp"
#include "math_library/typedefs.hpp"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

template <typename input_type, DOMAIN_REPRESENTATIONS DMN_REP_INPUT, typename output_type,
          DOMAIN_REPRESENTATIONS DMN_REP_OUTPUT>
class basis_transformation {
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

public:
  static bool& is_initialized() {
    static bool initialized = false;
    return initialized;
  }

  static std::string& get_name() {
    static std::string name = "basis-transformation";
    return name;
  }

  static LIN_ALG::matrix<f_scalar_type, LIN_ALG::CPU>& get_transformation_matrix() {
    static LIN_ALG::matrix<f_scalar_type, LIN_ALG::CPU> T;

    if (not is_initialized())
      initialize_transformation_matrix();

    return T;
  }

  static void initialize_transformation_matrix() {
    is_initialized() = true;

    int M = lh_dmn_type::get_size();
    int N = rh_dmn_type::get_size();

    assert(M > 0 and N > 0);

    LIN_ALG::matrix<f_scalar_type, LIN_ALG::CPU>& T = get_transformation_matrix();

    T.reserve(std::pair<int, int>(M, N));

    for (int j = 0; j < N; j++)
      for (int i = 0; i < M; i++)
        T(i, j) = basis_function_type::execute(i, j);
  }
};

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_BASIS_TRANSFORMATION_TEMPLATE_H
