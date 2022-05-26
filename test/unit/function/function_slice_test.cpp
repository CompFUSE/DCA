// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests the slicing functionality of the function template class
//

#include <numeric>
#include "function_testing.hpp"
#include "gtest/gtest.h"
#include "dca/util/to_string.hpp"
using dca::func::dmn;
using dca::func::dmn_0;
using dca::func::dmn_variadic;

namespace dca {
namespace testing {}
}  // namespace dca

TEST(Function, Slice) {
  using namespace dca::testing;
  dca::func::function<double, Domain2c> f2c;
  function_test<decltype(f2c)> f2c_test(f2c);
  f2c_test.fill_sequence();
  dca::func::function<double, dmn_variadic<Domain0c>> new_func;
  int subind[]{0, 3};
  f2c.slice(0, subind, static_cast<double*>(new_func.values()));
  // std::cout << dca::vectorToString(f2c.getValues()) << '\n';
  // std::cout << dca::vectorToString(new_func.getValues()) << '\n';
  dca::func::function<double, dmn_variadic<Domain0c>> expected_func{12, 13, 14, 15};
  EXPECT_EQ(new_func, expected_func);

  // Produce a slice in the Domain0d
  // The slice in Domain0c == 2
  int subind02[]{2, 0};
  // function to receive slice
  dca::func::function<double, dmn_variadic<Domain0d>> f02;
  f2c.slice(1, subind02, static_cast<double*>(f02.values()));
  dca::func::function<double, dmn_variadic<Domain0d>> expected_f02{2, 6, 10, 14, 18, 22, 26, 30};
  EXPECT_EQ(f02, expected_f02);

  // Two index slicing
  dca::func::function<double, Domain2c0c0c> f2c0c0c;
  function_test<decltype(f2c0c0c)> f2c0c0c_test(f2c0c0c);
  f2c0c0c_test.fill_sequence();
  dca::func::function<double, Domain2c> expected_f2c;
  dca::func::function<double, Domain2c> result_f2c;
  std::vector<int> subind03(3);
  for (int iw = 0; iw < Domain0c::dmn_size(); ++iw)
    for (int ik = 0; ik < Domain0c::dmn_size(); ++ik) {
      for (int ileading = 0; ileading < Domain2c::dmn_size(); ++ileading)
        expected_f2c(ileading) = f2c0c0c(ileading, ik, iw);
      subind03 = {0,ik,iw};
      f2c0c0c.slice(0, subind03, result_f2c.values());
      EXPECT_EQ(result_f2c, expected_f2c) << ik << ':' << iw << '\n'
                                          << dca::vectorToString(result_f2c.getValues()) << '\n'
                                          << dca::vectorToString(expected_f2c.getValues());
    }
}

TEST(Function, Distribute) {
  using namespace dca::testing;
  dca::func::function<double, Domain2c> f2c;
  function_test<decltype(f2c)> f2c_test(f2c);
  f2c_test.fill_sequence();
  double pre_sum = std::accumulate(f2c.begin(), f2c.end(), 0.0, std::plus<double>());
  dca::func::function<double, dmn_variadic<Domain0c>> new_func;
  int subind[]{0, 3};
  f2c.slice(0, subind, static_cast<double*>(new_func.values()));
  // std::cout << dca::vectorToString(f2c.getValues()) << '\n';
  // std::cout << dca::vectorToString(new_func.getValues()) << '\n';
  dca::func::function<double, dmn_variadic<Domain0c>> expected_func{12, 13, 14, 15};
  EXPECT_EQ(new_func, expected_func);
  new_func += 1;
  f2c.distribute(0, subind, static_cast<double*>(new_func.values()));
  double post_sum = std::accumulate(f2c.begin(), f2c.end(), 0.0, std::plus<double>());
  EXPECT_EQ(post_sum, pre_sum + Domain0c::dmn_size());

  // // Produce a slice in the Domain0d
  // // The slice in Domain0c == 2
  // int subind02[]{2,0};
  // // function to receive slice
  // dca::func::function<double, dmn_variadic<Domain0d>> f02;
  // f2c.slice(1,subind02, static_cast<double*>(f02.values()));
  // dca::func::function<double, dmn_variadic<Domain0d>> expected_f02{2,6,10,14,18,22,26,30};
  // EXPECT_EQ(f02, expected_f02);;

  // // Two index slicing
  // dca::func::function<double, Domain2c0c0c> f2c0c0c;
  // function_test<decltype(f2c0c0c)> f2c0c0c_test(f2c0c0c);
  // f2c0c0c_test.fill_sequence();
  // int subind03[]{2,2,0};
  // dca::func::function<double, Domain2c> result_f2c;
  // f2c0c0c.slice(1,2,subind03, static_cast<double*>(result_f2c.values()));
  // // std::cout << dca::vectorToString(result_f2c.getValues()) << '\n';
  // dca::func::function<double, Domain2c>
  // expected_f2c{2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82,86,90,94,98,102,106,110,114,118,122,126};
  // EXPECT_EQ(result_f2c, expected_f2c);
}
