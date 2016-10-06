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

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_FUNCTION_TRANSFORMS_TRANSFORM_FUNCTION_DMN_0_TO_DMN_0_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_FUNCTION_TRANSFORMS_TRANSFORM_FUNCTION_DMN_0_TO_DMN_0_H

#include <iostream>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/util/type_list.hpp"

#include "comp_library/linalg/linalg.hpp"
#include "math_library/functional_transforms/function_transforms/transform_function_domainwise.h"
#include "math_library/functional_transforms/function_transforms/transform_function_template.h"

namespace math_algorithms {
namespace functional_transforms {

template <typename type_input, typename type_output>
class TRANSFORM<func::dmn_0<type_input>, func::dmn_0<type_output>> {
  const static bool VERBOSE = false;

public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      print_types(f_input, f_output);

    using TRANSFORMED_DOMAIN =
        typename dca::func::SWAP_FIRST<domain_input, type_input, type_output>::Result;
    dca::util::assert_same<TRANSFORMED_DOMAIN, domain_output>();

    TRANSFORM_DOMAINWISE<domain_input, domain_output, type_input, type_output>::execute_on_first(
        f_input, f_output);
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output,
            class domain_output, typename scalartype_T>
  static void execute(func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output,
                      dca::linalg::Matrix<scalartype_T, dca::linalg::CPU>& T) {
    if (VERBOSE)
      print_types(f_input, f_output);

    using TRANSFORMED_DOMAIN =
        typename dca::func::SWAP_FIRST<domain_input, type_input, type_output>::Result;
    dca::util::assert_same<TRANSFORMED_DOMAIN, domain_output>();

    TRANSFORM_DOMAINWISE<domain_input, domain_output, type_input, type_output>::execute_on_first(
        f_input, f_output, T);
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute_on_all(func::function<scalartype_input, domain_input>& f_input,
                             func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      print_types(f_input, f_output);

    using TRANSFORMED_DOMAIN =
        typename dca::func::SWAP_ALL<domain_input, type_input, type_output>::Result;

    dca::util::assert_same<TRANSFORMED_DOMAIN, domain_output>();

    TRANSFORM_DOMAINWISE<domain_input, domain_output, type_input, type_output>::execute_on_all(
        f_input, f_output);
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output,
            class domain_output, typename scalartype_T>
  static void execute_on_all(func::function<scalartype_input, domain_input>& f_input,
                             func::function<scalartype_output, domain_output>& f_output,
                             dca::linalg::Matrix<scalartype_T, dca::linalg::CPU>& T) {
    if (VERBOSE)
      print_types(f_input, f_output);

    using TRANSFORMED_DOMAIN =
        typename dca::func::SWAP_ALL<domain_input, type_input, type_output>::Result;

    dca::util::assert_same<TRANSFORMED_DOMAIN, domain_output>();

    TRANSFORM_DOMAINWISE<domain_input, domain_output, type_input, type_output>::execute_on_all(
        f_input, f_output, T);
  }

private:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void print_types(func::function<scalartype_input, domain_input>& f_input,
                          func::function<scalartype_output, domain_output>& f_output,
                          bool do_all_domains = false) {
    typedef typename domain_input::this_type type_list_input;
    typedef typename domain_output::this_type type_list_output;

    dca::util::print_type<type_input>::to_JSON(std::cout);
    std::cout << "\n\n";
    dca::util::print_type<type_output>::to_JSON(std::cout);
    std::cout << "\n\n";

    if (do_all_domains) {
      dca::util::print_type<type_list_input>::to_JSON(std::cout);
      std::cout << "\n\n";

      using TRANSFORMED_DOMAIN =
          typename dca::func::SWAP_ALL<domain_input, type_input, type_output>::Result;

      dca::util::print_type<typename TRANSFORMED_DOMAIN::this_type>::to_JSON(std::cout);
      std::cout << "\n\n";

      dca::util::print_type<type_list_output>::to_JSON(std::cout);
      std::cout << "\n\n";

      func::function<scalartype_output, TRANSFORMED_DOMAIN> T;
      T.print_fingerprint();
      f_output.print_fingerprint();
    }
    else {
      dca::util::print_type<type_list_input>::to_JSON(std::cout);
      std::cout << "\n\n";

      using TRANSFORMED_DOMAIN =
          typename dca::func::SWAP_FIRST<domain_input, type_input, type_output>::Result;

      dca::util::print_type<typename TRANSFORMED_DOMAIN::this_type>::to_JSON(std::cout);
      std::cout << "\n\n";

      dca::util::print_type<type_list_output>::to_JSON(std::cout);
      std::cout << "\n\n";

      func::function<scalartype_output, TRANSFORMED_DOMAIN> T(
          "func::function<scalartype_output, TRANSFORMED_DOMAIN>");

      f_input.print_fingerprint();
      T.print_fingerprint();
      f_output.print_fingerprint();
    }
  }
};

}  // functional_transforms
}  // math_algorithms

#endif
