// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class transforms bases.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_BASIS_TRANSFORM_BASIS_TRANSFORM_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_BASIS_TRANSFORM_BASIS_TRANSFORM_HPP

#include <string>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/math/function_transform/basis_transform/basis_transformation.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

template <typename input_type, typename output_type>
class basis_transform {
public:
  typedef typename input_type::dmn_specifications_type input_spec_dmn_type;
  typedef typename output_type::dmn_specifications_type output_spec_dmn_type;

  typedef basis_transformation<input_type, input_spec_dmn_type::DOMAIN_REPRESENTATION, output_type,
                               output_spec_dmn_type::DOMAIN_REPRESENTATION>
      basis_transformation_type;

  typedef typename basis_transformation_type::matrix_type matrix_type;

public:
  static bool& is_initialized() {
    return basis_transformation_type::is_initialized();
  }

  static std::string& get_name() {
    return basis_transformation_type::get_name();
  }

  static matrix_type& get_transformation_matrix() {
    return basis_transformation_type::get_transformation_matrix();
  }
};

// Specialization for func::dmn_0 domains.
template <typename input_type, typename output_type>
class basis_transform<func::dmn_0<input_type>, func::dmn_0<output_type>> {
public:
  typedef typename input_type::dmn_specifications_type input_spec_dmn_type;
  typedef typename output_type::dmn_specifications_type output_spec_dmn_type;

  typedef basis_transformation<input_type, input_spec_dmn_type::DOMAIN_REPRESENTATION, output_type,
                               output_spec_dmn_type::DOMAIN_REPRESENTATION>
      basis_transformation_type;

  typedef typename basis_transformation_type::matrix_type matrix_type;

public:
  static bool& is_initialized() {
    return basis_transformation_type::is_initialized();
  }

  static std::string& get_name() {
    return basis_transformation_type::get_name();
  }

  static matrix_type& get_transformation_matrix() {
    return basis_transformation_type::get_transformation_matrix();
  }
};

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_BASIS_TRANSFORM_BASIS_TRANSFORM_HPP
