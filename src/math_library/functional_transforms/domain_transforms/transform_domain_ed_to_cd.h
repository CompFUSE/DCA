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

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_ED_TO_CD_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_ED_TO_CD_H

#include <iostream>

#include "comp_library/linalg/linalg.hpp"
#include "comp_library/function_library/include_function_library.h"
#include "math_library/functional_transforms/basis_transforms/basis_transforms.hpp"
#include "math_library/functional_transforms/domain_transforms/transform_domain_template.h"
#include "math_library/functional_transforms/domain_transforms/transformation_characteristics.h"
#include "math_library/functional_transforms/typedefs.hpp"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

template <typename type_input, typename type_output, int DMN_INDEX>
class TRANSFORM_DOMAIN<type_input, EXPANSION, type_output, CONTINUOUS, DMN_INDEX>
    : public TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX> {
private:
  const static bool VERBOSE = false;

  typedef basis_transformation<type_input, EXPANSION, type_output, CONTINUOUS> basis_transformation_type;
  typedef typename basis_transformation_type::matrix_type matrix_type;

public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(FUNC_LIB::function<scalartype_input, domain_input>& f_input,
                      FUNC_LIB::function<scalartype_output, domain_output>& f_output) {
    default_execute(f_input, f_output);
  }

private:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void default_execute(FUNC_LIB::function<scalartype_input, domain_input>& f_input,
                              FUNC_LIB::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      std::cout << "\n\t default-transform (expansion -> continuous) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

    matrix_type& T = basis_transformation_type::get_transformation_matrix();

    TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(f_input, f_output, T);
  }
};

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_ED_TO_CD_H
