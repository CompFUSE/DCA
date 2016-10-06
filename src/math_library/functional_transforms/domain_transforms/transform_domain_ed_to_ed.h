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

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_ED_TO_ED_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_ED_TO_ED_H

#include <iostream>

#include "dca/function/function.hpp"
#include "comp_library/linalg/linalg.hpp"
#include "math_library/functional_transforms/basis_transforms/basis_transforms.hpp"
#include "math_library/functional_transforms/domain_transforms/transform_domain_template.h"
#include "math_library/functional_transforms/domain_transforms/transformation_characteristics.h"
#include "math_library/functional_transforms/typedefs.hpp"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

template <typename type_input, typename type_output, int DMN_INDEX>
class TRANSFORM_DOMAIN<type_input, EXPANSION, type_output, EXPANSION, DMN_INDEX>
    : public TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX> {
private:
  const static bool VERBOSE = false;

  typedef basis_transformation<type_input, EXPANSION, type_output, EXPANSION> basis_transformation_type;
  typedef typename basis_transformation_type::matrix_type matrix_type;

public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output) {
    default_execute(f_input, f_output);
  }

private:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void default_execute(func::function<scalartype_input, domain_input>& f_input,
                              func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      std::cout << "\n\t default-transform (continuous -> expansion) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

    matrix_type& T = basis_transformation_type::get_transformation_matrix();

    transform(f_input, f_output, T);
  }
};

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_ED_TO_ED_H
