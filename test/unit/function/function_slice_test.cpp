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
#include "dca/function/functionToString.hpp"

using dca::func::dmn;
using dca::func::dmn_0;
using dca::func::dmn_variadic;

namespace dca {
namespace testing {}
}  // namespace dca

TEST(Function, SingleIndexSlice) {
  using namespace dca::testing;
  // Domain2c -> <Domain0c, Domain0d> ->  <dmn_0<4, double>, dmn_0<8, double>

  dca::func::function<double, Domain2c> f2c;
  function_test<decltype(f2c)> f2c_test(f2c);
  f2c_test.fill_sequence();
  /**  0  1  2  3
       4  5  6  7
       8  9  10 11
       12 13 14 15
       16 17 18 19
       20 21 22 23
       24 25 26 27
       28 29 30 21
  */
  dca::func::function<double, dmn_variadic<Domain0c>> new_func;

  // tests of the deprecated ptr to int interface
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

  //  prefered reference to std::vector<int> inteface
  std::vector<int> subind_2 = {3, 0};
  dca::func::function<double, dmn_variadic<Domain0d>> new_func_0d;
  dca::func::function<double, dmn_variadic<Domain0d>> expected_func_0d{3,  7,  11, 15,
                                                                       19, 23, 27, 31};
  f2c.slice(1, subind_2, static_cast<double*>(new_func_0d.values()));
  EXPECT_TRUE(new_func_0d == expected_func_0d) << functionToString(new_func_0d);

  // multi level domains single slice
  dca::func::function<double, dmn_variadic<dmn_variadic<Domain0c, Domain0d>, Domain0b>> ml_func;
  function_test<decltype(ml_func)> mlf_test(ml_func);
  mlf_test.fill_sequence();
  /* 0, 1, 2, 3, 4, 5, 6, 7 ,8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 ,25...
     32, 33, 34 ...
  */
  dca::func::function<double, Domain2c> new_func_2c;
  dca::func::function<double, Domain2c> expected_func_2c{32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
                                                         43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
                                                         54, 55, 56, 57, 58, 59, 60, 61, 62, 63};
  ml_func.slice(0, {0, 1}, new_func_2c.values());
  EXPECT_EQ(new_func_2c, expected_func_2c);

  // doubly nested
  dca::func::function<double, dmn_variadic<dmn_variadic<Domain2c, Domain0c, Domain0d>, Domain0b>>
      ml_func_double;
  function_test<decltype(ml_func_double)> mlf_test_double(ml_func_double);
  mlf_test_double.fill_sequence();
  // exercise for the reader what the layout looks like.
  // but the leftmost Domain

  dca::func::function<double, dmn_variadic<Domain2c, Domain0c, Domain0d>> first_func;
  ml_func_double.slice(0, {0, 1}, first_func.values());
  first_func.slice(0, {0, 0, 1}, new_func_2c.values());

  function_test<decltype(expected_func_2c)> expected_double(expected_func_2c);
  expected_double.fill_sequence(1152);
  EXPECT_EQ(new_func_2c, expected_func_2c);
}

TEST(Function, TwoIndexSlice) {
  using namespace dca::testing;

  // Two index slicing
  dca::func::function<double, Domain2c0c0c> f2c0c0c;
  function_test<decltype(f2c0c0c)> f2c0c0c_test(f2c0c0c);
  f2c0c0c_test.fill_sequence();
  /**
   */
  dca::func::function<double, Domain2c> expected_f2c;
  dca::func::function<double, Domain2c> result_f2c;
  std::vector<int> subind03(3);
  for (int iw = 0; iw < Domain0c::dmn_size(); ++iw)
    for (int ik = 0; ik < Domain0c::dmn_size(); ++ik) {
      for (int ileading = 0; ileading < Domain2c::dmn_size(); ++ileading)
        expected_f2c(ileading) = f2c0c0c(ileading, ik, iw);
      subind03 = {0, ik, iw};
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
