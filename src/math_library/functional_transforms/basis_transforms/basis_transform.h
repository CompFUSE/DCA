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

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_BASIS_TRANSFORM_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_BASIS_TRANSFORM_H

#include <string>

#include "math_library/functional_transforms/basis_transforms/basis_transformation_cd_to_cd.h"
#include "math_library/functional_transforms/basis_transforms/basis_transformation_cd_to_dd.h"
#include "math_library/functional_transforms/basis_transforms/basis_transformation_cd_to_ed.h"
#include "math_library/functional_transforms/basis_transforms/basis_transformation_dd_to_ed.h"
#include "math_library/functional_transforms/basis_transforms/basis_transformation_ed_to_cd.h"
#include "math_library/functional_transforms/basis_transforms/basis_transformation_ed_to_dd.h"
#include "math_library/functional_transforms/basis_transforms/basis_transformation_ed_to_ed.h"
#include "math_library/functional_transforms/basis_transforms/basis_transformation_template.h"
#include "comp_library/function_library/domains/special_domains/dmn_0.h"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::
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

template <typename input_type, typename output_type>
class basis_transform<dmn_0<input_type>, dmn_0<output_type>> {
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

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_BASIS_TRANSFORM_H
