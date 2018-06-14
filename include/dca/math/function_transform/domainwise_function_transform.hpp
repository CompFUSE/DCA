// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Helper class for FunctionTransform.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_DOMAINWISE_FUNCTION_TRANSFORM_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_DOMAINWISE_FUNCTION_TRANSFORM_HPP

#include <iostream>

#include "dca/function/domains/domain_type_operations.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/math/function_transform/domain_representations.hpp"
#include "dca/math/function_transform/transform_domain.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

template <typename domain_input, typename domain_output, typename type_input, typename type_output>
struct DomainwiseFunctionTransform {
  const static bool VERBOSE = false;

  using TRANSFORMED_DOMAIN =
      typename dca::func::SWAP_FIRST<domain_input, type_input, type_output>::Result;

  const static int CURR_DMN_INDEX =
      dca::util::IndexOf<type_input, typename domain_input::this_type>::value;
  const static int NEXT_DMN_INDEX =
      dca::util::IndexOf<type_input, typename TRANSFORMED_DOMAIN::this_type>::value;

  typedef typename type_input::dmn_specifications_type input_specs_type;
  typedef typename type_output::dmn_specifications_type output_specs_type;

  const static DOMAIN_REPRESENTATIONS DMN_REP_LHS = input_specs_type::DOMAIN_REPRESENTATION;
  const static DOMAIN_REPRESENTATIONS DMN_REP_RHS = output_specs_type::DOMAIN_REPRESENTATION;

  template <typename scalartype_input, typename scalartype_output>
  static void execute_on_first(const func::function<scalartype_input, domain_input>& f_input,
                               func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      std::cout << "\n\n\t" << __FUNCTION__ << "\t" << f_input.get_name() << " --> "
                << f_output.get_name() << "\n\n";

    assert(CURR_DMN_INDEX > -1 and CURR_DMN_INDEX < f_input.signature());

    TRANSFORM_DOMAIN<type_input, DMN_REP_LHS, type_output, DMN_REP_RHS, CURR_DMN_INDEX>::execute(
        f_input, f_output);
  }

  template <typename scalartype_input, typename scalartype_output, typename scalartype_T>
  static void execute_on_first(const func::function<scalartype_input, domain_input>& f_input,
                               func::function<scalartype_output, domain_output>& f_output,
                               const linalg::Matrix<scalartype_T, linalg::CPU>& T) {
    if (VERBOSE)
      std::cout << "\n\n\t" << __FUNCTION__ << "\t" << f_input.get_name() << " --> "
                << f_output.get_name() << "\n\n";

    assert(CURR_DMN_INDEX > -1 and CURR_DMN_INDEX < f_input.signature());

    TRANSFORM_DOMAIN_PROCEDURE<CURR_DMN_INDEX>::transform(f_input, f_output, T);
  }

  template <typename scalartype_input, typename scalartype_output>
  static void execute_on_all(const func::function<scalartype_input, domain_input>& f_input,
                             func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE) {
      std::cout << "\n\n\t" << __FUNCTION__ << "\t" << f_input.get_name() << " --> "
                << f_output.get_name() << "\n\n";

      f_input.print_fingerprint();
    }

    func::function<scalartype_output, TRANSFORMED_DOMAIN> f_output_new("f_output_new");

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

      DomainwiseFunctionTransform<TRANSFORMED_DOMAIN, domain_output, type_input, type_output>::execute_on_all(
          f_output_new, f_output);
    }
  }

  template <typename scalartype_input, typename scalartype_output, typename scalartype_T>
  static void execute_on_all(const func::function<scalartype_input, domain_input>& f_input,
                             func::function<scalartype_output, domain_output>& f_output,
                             const linalg::Matrix<scalartype_T, linalg::CPU>& T) {
    if (VERBOSE) {
      std::cout << "\n\n\t" << __FUNCTION__ << "\t" << f_input.get_name() << " --> "
                << f_output.get_name() << "\n\n";

      f_input.print_fingerprint();
    }

    func::function<scalartype_output, TRANSFORMED_DOMAIN> f_output_new("f_output_new");

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

      DomainwiseFunctionTransform<TRANSFORMED_DOMAIN, domain_output, type_input, type_output>::execute_on_all(
          f_output_new, f_output, T);
    }
  }
};

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_DOMAINWISE_FUNCTION_TRANSFORM_HPP
