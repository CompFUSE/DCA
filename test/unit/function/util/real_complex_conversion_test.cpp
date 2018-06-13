// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file tests real_complex_conversion.hpp.

#include "dca/function/util/real_complex_conversion.hpp"

#include "gtest/gtest.h"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"

template <int num_elements, typename Scalar = double>
using Function = dca::func::function<Scalar, dca::func::dmn_0<dca::func::dmn<num_elements>>>;

TEST(RealComplexConversionTest, RealToComplex) {
  Function<2> f;
  f(0) = 0.5;
  f(1) = 1.;

  auto f_complex = dca::func::util::complex(f);

  EXPECT_DOUBLE_EQ(0.5, f_complex(0).real());
  EXPECT_DOUBLE_EQ(0., f_complex(0).imag());
  EXPECT_DOUBLE_EQ(1., f_complex(1).real());
  EXPECT_DOUBLE_EQ(0., f_complex(1).imag());
}

TEST(RealComplexConversionTest, ComplexToReal) {
  Function<2, std::complex<double>> f;
  f(0) = std::complex<double>(0.5, 0.);
  f(1) = std::complex<double>(1., 1.);

  auto f_real = dca::func::util::real(f, false);

  EXPECT_DOUBLE_EQ(0.5, f_real(0));
  EXPECT_DOUBLE_EQ(1., f_real(1));

  // If the second argument is true, the function's imaginary part must be zero.
  EXPECT_THROW(dca::func::util::real(f, true), std::logic_error);

  f(1).imag(0.);
  f_real = dca::func::util::real(f, true);

  EXPECT_DOUBLE_EQ(0.5, f_real(0));
  EXPECT_DOUBLE_EQ(1., f_real(1));
}
