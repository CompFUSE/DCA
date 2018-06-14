// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Helper class for basis_transform.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_BASIS_TRANSFORM_BASIS_TRANSFORMATION_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_BASIS_TRANSFORM_BASIS_TRANSFORMATION_HPP

#include <cassert>
#include <string>

#include "dca/linalg/linalg.hpp"
#include "dca/math/function_transform/basis_function.hpp"
#include "dca/math/function_transform/basis_transform/inner_product_domain.hpp"
#include "dca/math/function_transform/domain_representations.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

// General template
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

  typedef dca::linalg::Matrix<f_scalar_type, dca::linalg::CPU> matrix_type;

public:
  static bool& is_initialized() {
    static bool initialized = false;
    return initialized;
  }

  static std::string& get_name() {
    static std::string name = "basis-transformation";
    return name;
  }

  static dca::linalg::Matrix<f_scalar_type, dca::linalg::CPU>& get_transformation_matrix() {
    static dca::linalg::Matrix<f_scalar_type, dca::linalg::CPU> T;

    if (not is_initialized())
      initialize_transformation_matrix();

    return T;
  }

  static void initialize_transformation_matrix() {
    is_initialized() = true;

    int M = lh_dmn_type::get_size();
    int N = rh_dmn_type::get_size();

    assert(M > 0 and N > 0);

    dca::linalg::Matrix<f_scalar_type, dca::linalg::CPU>& T = get_transformation_matrix();

    T.reserve(std::pair<int, int>(M, N));

    for (int j = 0; j < N; j++)
      for (int i = 0; i < M; i++)
        T(i, j) = basis_function_type::execute(i, j);
  }
};

// Specialization for continous --> continous
#include "basis_transformation_cd_to_cd.inc"

// Specialization for continous --> discrete
#include "basis_transformation_cd_to_dd.inc"

// Specialization for continous --> expansion
#include "basis_transformation_cd_to_ed.inc"

// Specialization for discrete --> expansion
#include "basis_transformation_dd_to_ed.inc"

// Specialization for expansion --> continous
#include "basis_transformation_ed_to_cd.inc"

// Specialization for expansion --> discrete
#include "basis_transformation_ed_to_dd.inc"

// Specialization for expansion --> expansion
#include "basis_transformation_ed_to_ed.inc"

}  // transform
}  // math
}  // dca

#endif  //  DCA_MATH_FUNCTION_TRANSFORM_BASIS_TRANSFORM_BASIS_TRANSFORMATION_HPP
