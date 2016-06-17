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

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_FUNCTION_TRANSFORMS_TRANSFORM_FUNCTION_DOMAINWISE_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_FUNCTION_TRANSFORMS_TRANSFORM_FUNCTION_DOMAINWISE_H

#include <iostream>

#include "dca/util/type_list.hpp"

#include "comp_library/linalg/linalg.hpp"
#include "comp_library/function_library/include_function_library.h"
#include "math_library/functional_transforms/typedefs.hpp"
#include "math_library/functional_transforms/domain_transforms/domain_transforms.hpp"

namespace math_algorithms {
namespace functional_transforms {

template <typename domain_input, typename domain_output, typename type_input, typename type_output>
struct TRANSFORM_DOMAINWISE {
  const static bool VERBOSE = false;

  using TRANSFORMED_DOMAIN =
      typename dca::util::SWAP_FIRST<domain_input, type_input, type_output>::Result;

  const static int CURR_DMN_INDEX =
      dca::util::IndexOf<type_input, typename domain_input::this_type>::value;
  const static int NEXT_DMN_INDEX =
      dca::util::IndexOf<type_input, typename TRANSFORMED_DOMAIN::this_type>::value;

  typedef typename type_input::dmn_specifications_type input_specs_type;
  typedef typename type_output::dmn_specifications_type output_specs_type;

  const static DOMAIN_REPRESENTATIONS DMN_REP_LHS = input_specs_type::DOMAIN_REPRESENTATION;
  const static DOMAIN_REPRESENTATIONS DMN_REP_RHS = output_specs_type::DOMAIN_REPRESENTATION;

  template <typename scalartype_input, typename scalartype_output>
  static void execute_on_first(FUNC_LIB::function<scalartype_input, domain_input>& f_input,
                               FUNC_LIB::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      std::cout << "\n\n\t" << __FUNCTION__ << "\t" << f_input.get_name() << " --> "
                << f_output.get_name() << "\n\n";

    assert(CURR_DMN_INDEX > -1 and CURR_DMN_INDEX < f_input.signature());

    math_algorithms::functional_transforms::TRANSFORM_DOMAIN<
        type_input, DMN_REP_LHS, type_output, DMN_REP_RHS, CURR_DMN_INDEX>::execute(f_input,
                                                                                    f_output);
  }

  template <typename scalartype_input, typename scalartype_output, typename scalartype_T>
  static void execute_on_first(FUNC_LIB::function<scalartype_input, domain_input>& f_input,
                               FUNC_LIB::function<scalartype_output, domain_output>& f_output,
                               LIN_ALG::matrix<scalartype_T, LIN_ALG::CPU>& T) {
    if (VERBOSE)
      std::cout << "\n\n\t" << __FUNCTION__ << "\t" << f_input.get_name() << " --> "
                << f_output.get_name() << "\n\n";

    assert(CURR_DMN_INDEX > -1 and CURR_DMN_INDEX < f_input.signature());

    TRANSFORM_DOMAIN_PROCEDURE<CURR_DMN_INDEX>::transform(f_input, f_output, T);
  }

  template <typename scalartype_input, typename scalartype_output>
  static void execute_on_all(FUNC_LIB::function<scalartype_input, domain_input>& f_input,
                             FUNC_LIB::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE) {
      std::cout << "\n\n\t" << __FUNCTION__ << "\t" << f_input.get_name() << " --> "
                << f_output.get_name() << "\n\n";

      f_input.print_fingerprint();
    }

    FUNC_LIB::function<scalartype_output, TRANSFORMED_DOMAIN> f_output_new("f_output_new");

    TRANSFORM_DOMAIN<type_input, DMN_REP_LHS, type_output, DMN_REP_RHS, CURR_DMN_INDEX>::execute(
        f_input, f_output_new);

    if (NEXT_DMN_INDEX == -1) {
      if (VERBOSE)
        std::cout << "\n\n\tSTOP\t  start the copy " << std::endl;

      assert(f_output_new.size() == f_output.size());

      for (int l = 0; l < f_output.size(); l++)
        f_output(l) = f_output_new(l);

      if (VERBOSE) {
        f_output.print_fingerprint();

        std::cout << "\n\n\tSTOP\t finished the copy\n\n\n" << std::endl;
      }
    }
    else {
      if (VERBOSE) {
        std::cout << "\n\n\tSTART\t" << __FUNCTION__ << std::endl;

        f_output_new.print_fingerprint();
      }

      TRANSFORM_DOMAINWISE<TRANSFORMED_DOMAIN, domain_output, type_input, type_output>::execute_on_all(
          f_output_new, f_output);
    }
  }

  template <typename scalartype_input, typename scalartype_output, typename scalartype_T>
  static void execute_on_all(FUNC_LIB::function<scalartype_input, domain_input>& f_input,
                             FUNC_LIB::function<scalartype_output, domain_output>& f_output,
                             LIN_ALG::matrix<scalartype_T, LIN_ALG::CPU>& T) {
    if (VERBOSE) {
      std::cout << "\n\n\t" << __FUNCTION__ << "\t" << f_input.get_name() << " --> "
                << f_output.get_name() << "\n\n";

      f_input.print_fingerprint();
    }

    FUNC_LIB::function<scalartype_output, TRANSFORMED_DOMAIN> f_output_new("f_output_new");

    TRANSFORM_DOMAIN_PROCEDURE<CURR_DMN_INDEX>::transform(f_input, f_output_new, T);

    if (NEXT_DMN_INDEX == -1) {
      if (VERBOSE)
        std::cout << "\n\n\tSTOP\t  start the copy " << std::endl;

      assert(f_output_new.size() == f_output.size());

      for (int l = 0; l < f_output.size(); l++)
        f_output(l) = f_output_new(l);

      if (VERBOSE) {
        f_output.print_fingerprint();

        std::cout << "\n\n\tSTOP\t finished the copy\n\n\n" << std::endl;
      }
    }
    else {
      if (VERBOSE) {
        std::cout << "\n\n\tSTART\t" << __FUNCTION__ << std::endl;

        f_output_new.print_fingerprint();
      }

      TRANSFORM_DOMAINWISE<TRANSFORMED_DOMAIN, domain_output, type_input, type_output>::execute_on_all(
          f_output_new, f_output, T);
    }
  }
};

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_FUNCTION_TRANSFORMS_TRANSFORM_FUNCTION_DOMAINWISE_H
