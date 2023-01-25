// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file refactors code shared between the many function template class tests
//

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"

#ifndef DCA_FUNCTION_TESTING_HPP
#define DCA_FUNCTION_TESTING_HPP

using dca::func::dmn;
using dca::func::dmn_0;
using dca::func::dmn_variadic;

namespace dca {
namespace testing {

// A selection of domain types we can use for testing.
using Domain0a = dmn_0<dmn<1, double>>;
using Domain0b = dmn_0<dmn<2, double>>;
using Domain0c = dmn_0<dmn<4, double>>;
using Domain0d = dmn_0<dmn<8, double>>;
using Domain1d = dmn_variadic<Domain0d>;
using Domain2a = dmn_variadic<Domain0a, Domain0b>;
using Domain2c = dmn_variadic<Domain0c, Domain0d>;
using Domain4a = dmn_variadic<Domain2a, Domain2c>;
using Domain16 = dmn_variadic<Domain4a, Domain4a, Domain2c>;
using Domain0v = dmn_variadic<Domain0d>;
using Domain16v = dmn_variadic<Domain4a, Domain4a, Domain2c>;
using Domain2c0c0c = dmn_variadic<Domain2c, Domain0c, Domain0c>;

template <typename T1>
struct function_test {};

template <int N, typename Dmn>
struct function_test<dca::func::function<double, dmn<N, Dmn>>> {
  using fType = dca::func::function<double, dmn<N, Dmn>>;

  function_test(fType& func) : f(func) {}

  template <typename Arg>
  bool check_1(Arg /*arg*/) {
    /*
    std::cout << "Sub branch size " << std::endl;
    for (int i = 0; i < f.get_Nb_branch_domains(); i++) {
      std::cout << f.get_branch_size(i) << "\t";
    }
    std::cout << std::endl;

    s
    td::cout << "Sub branch steps " << std::endl;
    for (int i = 0; i < f.get_Nb_branch_domains(); i++) {
      std::cout << f.get_branch_domain_steps()[i] << "\t";
    }
    std::cout << std::endl;
    */
    return true;
  }

  fType& f;
};


  
template <typename DMN>
struct function_test<dca::func::function<double, DMN>> {
  using FuncType = dca::func::function<double, DMN>;
  using ScalarType = typename FuncType::this_scalar_type;
  using Domain = DMN;

  function_test(FuncType& func) : f(func) {}

  int signature() {
    return f.signature();
  }
  int size() {
    return f.size();
  }

  void fill_sequence(int start = 0) {
    int N = f.size();
    for (int i = 0; i < N; ++i) {
      f(i) = i + start;
      // if (i<1024) std::cout << i << ",";
    }
  }

  void check_sequence() {
    int N = f.size();
    for (int i = 0; i < N; ++i) {
      if (f(i) != i)
        throw(std::runtime_error("fault"));
    }
    // std::cout << "Ntypes " << Ntypes << " signature " << signature() << " size " << size()
    //           << std::endl;
  }

  template <typename Arg>
  bool check_1(Arg /*arg*/) {
    std::cout << "Sub branch size " << std::endl;
    for (int i = 0; i < f.get_domain().get_Nb_branch_domains(); i++) {
      std::cout << f.get_domain().get_branch_size(i) << "\t";
    }
    std::cout << std::endl;

    std::cout << "Sub branch steps " << std::endl;
    for (int i = 0; i < f.get_domain().get_Nb_branch_domains(); i++) {
      std::cout << f.get_domain().get_branch_domain_steps()[i] << "\t";
    }
    std::cout << std::endl;

    // IsSame<mp_append<test_domain_0a::this_type, test_domain_0b::this_type>::type, double>();
    // IsSame<test_domain_0a::this_type, test_domain_0b::this_type>();
    // IsSame<test_domain_0a, test_domain_0b>();
    // IsSame<test_domain_2a, test_domain_2c>();

    std::cout << std::endl;

    // std::cout << "\nTesting Typelist count "
    //           << dca::util::type_name<dca::util::TypeAt<2, test_list>>().c_str() << std::endl;
    // std::cout << "\nTesting Typelist count "
    //           << dca::util::type_name<dca::util::TypeAt<2, test_list>>().c_str() << std::endl;

    // typedef typename TypeAt<typename Domain::domain_typelist_0, 0>::Result dom_0;
    // std::cout << "Getting first subdomain "
    //           << "Type Id is " << typeid(dom_0).name() << std::endl;
    // dca::func::function<double, dom_0> sub_function;
    // function_test<decltype(sub_function)> sub_domain(sub_function);
    // sub_domain.check_1(1);

    return true;
  }

  template <typename... Args>
  ScalarType expand(Args... /*args*/) {
    return ScalarType(0);
  }

  template <typename... Args>
  bool check_value(Args... args) {
    // if (f(args...) == arg1 * offset<f, 1> + arg2 * offset<f, 2> +) {
    // }
    return f.operator()(args...) == f(args...);
    // return check_value(args...);
  }
  FuncType& f;
};

}  // namespace testing
}  // namespace dca

#endif
