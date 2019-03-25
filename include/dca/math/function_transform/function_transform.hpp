// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class transforms dca::func::function objects.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_FUNCTION_TRANSFORM_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_FUNCTION_TRANSFORM_HPP

#include <iostream>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/math/function_transform/domainwise_function_transform.hpp"
#include "dca/math/function_transform/special_transforms/momentum_to_space.hpp"
#include "dca/math/function_transform/special_transforms/space_to_momentum.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

// Empty class template.
template <typename type_input, typename type_output>
class FunctionTransform {};

// Specialization for dmn_0 to dmn_0.
template <typename type_input, typename type_output>
class FunctionTransform<func::dmn_0<type_input>, func::dmn_0<type_output>> {
  const static bool VERBOSE = false;

public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(const func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      print_types(f_input, f_output);

    using TRANSFORMED_DOMAIN =
        typename dca::func::SWAP_FIRST<domain_input, type_input, type_output>::Result;
    dca::util::assert_same<TRANSFORMED_DOMAIN, domain_output>();

    DomainwiseFunctionTransform<domain_input, domain_output, type_input, type_output>::execute_on_first(
        f_input, f_output);
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output,
            class domain_output, typename scalartype_T>
  static void execute(const func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output,
                      const linalg::Matrix<scalartype_T, linalg::CPU>& T) {
    if (VERBOSE)
      print_types(f_input, f_output);

    using TRANSFORMED_DOMAIN =
        typename dca::func::SWAP_FIRST<domain_input, type_input, type_output>::Result;
    dca::util::assert_same<TRANSFORMED_DOMAIN, domain_output>();

    DomainwiseFunctionTransform<domain_input, domain_output, type_input, type_output>::execute_on_first(
        f_input, f_output, T);
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute_on_all(const func::function<scalartype_input, domain_input>& f_input,
                             func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      print_types(f_input, f_output);

    using TRANSFORMED_DOMAIN =
        typename dca::func::SWAP_ALL<domain_input, type_input, type_output>::Result;

    dca::util::assert_same<TRANSFORMED_DOMAIN, domain_output>();

    DomainwiseFunctionTransform<domain_input, domain_output, type_input, type_output>::execute_on_all(
        f_input, f_output);
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output,
            class domain_output, typename scalartype_T>
  static void execute_on_all(const func::function<scalartype_input, domain_input>& f_input,
                             func::function<scalartype_output, domain_output>& f_output,
                             const linalg::Matrix<scalartype_T, linalg::CPU>& T) {
    if (VERBOSE)
      print_types(f_input, f_output);

    using TRANSFORMED_DOMAIN =
        typename dca::func::SWAP_ALL<domain_input, type_input, type_output>::Result;

    dca::util::assert_same<TRANSFORMED_DOMAIN, domain_output>();

    DomainwiseFunctionTransform<domain_input, domain_output, type_input, type_output>::execute_on_all(
        f_input, f_output, T);
  }

private:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void print_types(const func::function<scalartype_input, domain_input>& f_input,
                          const func::function<scalartype_output, domain_output>& f_output,
                          const bool do_all_domains = false) {
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

// Specialization for momentum to space transform.
template <typename DmnScalar, int dim, phys::domains::CLUSTER_NAMES c_name,
          phys::domains::CLUSTER_SHAPE c_shape>
class FunctionTransform<
    func::dmn_0<phys::domains::cluster_domain<DmnScalar, dim, c_name, phys::domains::MOMENTUM_SPACE, c_shape>>,
    func::dmn_0<phys::domains::cluster_domain<DmnScalar, dim, c_name, phys::domains::REAL_SPACE, c_shape>>> {
  using KDmn = func::dmn_0<
      phys::domains::cluster_domain<DmnScalar, dim, c_name, phys::domains::MOMENTUM_SPACE, c_shape>>;
  using RDmn =
      func::dmn_0<phys::domains::cluster_domain<DmnScalar, dim, c_name, phys::domains::REAL_SPACE, c_shape>>;

public:
  template <typename ScalarInp, typename ScalarOut, class DomainInput, class DomainOutput>
  static void execute(const func::function<ScalarInp, DomainInput>& f_input,
                      func::function<ScalarOut, DomainOutput>& f_output) {
    MomentumToSpaceTransform<KDmn, RDmn>::execute(f_input, f_output);
  }
};

// Specialization for space to momentum transformation.
template <typename DmnScalar, int dim, phys::domains::CLUSTER_NAMES c_name,
          phys::domains::CLUSTER_SHAPE c_shape>
class FunctionTransform<
    func::dmn_0<phys::domains::cluster_domain<DmnScalar, dim, c_name, phys::domains::REAL_SPACE, c_shape>>,
    func::dmn_0<phys::domains::cluster_domain<DmnScalar, dim, c_name, phys::domains::MOMENTUM_SPACE, c_shape>>> {
  using KDmn = func::dmn_0<
      phys::domains::cluster_domain<DmnScalar, dim, c_name, phys::domains::MOMENTUM_SPACE, c_shape>>;
  using RDmn =
      func::dmn_0<phys::domains::cluster_domain<DmnScalar, dim, c_name, phys::domains::REAL_SPACE, c_shape>>;

public:
  template <typename ScalarInp, typename ScalarOut, class DomainInput, class DomainOutput>
  static void execute(const func::function<ScalarInp, DomainInput>& f_input,
                      func::function<ScalarOut, DomainOutput>& f_output) {
    SpaceToMomentumTransform<RDmn, KDmn>::execute(f_input, f_output);
  }
};

}  // namespace transform
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_FUNCTION_TRANSFORM_HPP
