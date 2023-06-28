// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests the complex operators on cuda complex types

#include "dca/platform/dca_gpu.h"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include <stdexcept>
#include "dca/testing/gtest_h_w_warning_blocking.h"

TEST(ComplexOpCuda, Assign) {
  double2 d1_a{0.0, 0.0};
  double2 d2_a{1.0,2.0};
  dca::linalg::assign(d1_a,d2_a);
  EXPECT_EQ(d1_a.x, d2_a.x);
  EXPECT_EQ(d1_a.y, d2_a.y);

  std::complex<double> c1{1.3,2.4};
  dca::linalg::assign(d1_a,c1);
  EXPECT_EQ(d1_a.x, c1.real());
  EXPECT_EQ(d1_a.y, c1.imag());

  std::complex<double> c2{0.0,0.5};
  dca::linalg::assign(c2, d1_a);
  EXPECT_EQ(c2.real(), c1.real());
  EXPECT_EQ(c2.imag(),c1.imag());

  std::int8_t i81 = 1;
  dca::linalg::assign<double>(d1_a, i81);
}
