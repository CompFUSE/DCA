// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file tests the helper methods of function_utils.hpp.

#include "dca/function/function_utils.hpp"

#include "gtest/gtest.h"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"

template <int n_elem, typename Scalar = double>
using Function = dca::func::function<Scalar, dca::func::dmn_0<dca::func::dmn<n_elem>>>;

TEST(FunctionUtilsTest, Difference) {
  Function<4> f1, f2;
  f1 = 1.;

  f2(0) = 1.3;
  f2(1) = 0.7;
  f2(2) = 2;
  f2(3) = 1.4;

  dca::func::utils::Difference diff = dca::func::utils::difference(f1, f2);

  EXPECT_NEAR(0.5, diff.l1, 3e-7);
  EXPECT_NEAR(0.5787918, diff.l2, 3e-7);
  EXPECT_NEAR(1, diff.l_inf, 3e-7);

  Function<6> f3;
  f3 = 1.;

  EXPECT_THROW(dca::func::utils::difference(f1, f3), std::logic_error);

  Function<4, float> f_float;
  for (int i = 0; i < f_float.size(); ++i)
    f_float(i) = f2(i);
  auto diff2 = dca::func::utils::difference(f2, f_float);
  EXPECT_NEAR(0., diff2.l2, 1e-7);

  // Test with different complex types.
  Function<4, std::complex<float>> f_c;
  Function<4, std::complex<double>> f_z;
  for (int i = 0; i < f_c.size(); ++i) {
    f_c(i) = std::complex<float>(i, i * i);
    f_z(i) = std::complex<double>(i, i * i);
  }
  auto diff3 = dca::func::utils::difference(f_c, f_z);
  EXPECT_NEAR(0., diff3.l2, 5e-7);
}

TEST(FunctionUtilsTest, Complex) {
  Function<2> f;
  f(0) = 0.5, f(1) = 1.;

  auto f_cmplx = dca::func::utils::complex(f);

  EXPECT_DOUBLE_EQ(0.5, f_cmplx(0).real());
  EXPECT_DOUBLE_EQ(0., f_cmplx(0).imag());
  EXPECT_DOUBLE_EQ(1., f_cmplx(1).real());
  EXPECT_DOUBLE_EQ(0., f_cmplx(1).imag());
}

TEST(FunctionUtilsTest, Real) {
  Function<2, std::complex<double>> f;
  f(0) = std::complex<double>(0.5, 0);
  f(1) = std::complex<double>(1, 1);

  auto f_real = dca::func::utils::real(f, false);

  EXPECT_DOUBLE_EQ(0.5, f_real(0));
  EXPECT_DOUBLE_EQ(1., f_real(1));

  // If the second argument is set to true, the function's imaginary part must be zero.
  EXPECT_THROW(dca::func::utils::real(f, true), std::logic_error);

  f(1).imag(0.);
  f_real = dca::func::utils::real(f, true);
  EXPECT_DOUBLE_EQ(1., f_real(1));
}
