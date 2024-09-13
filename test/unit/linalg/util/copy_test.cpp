// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests generic dca::util::memorycopy functions

#include "dca/linalg/util/copy.hpp"

#include <complex>
#include <stdexcept>

#include "dca/testing/type_testing.hpp"
#include "dca/util/type_help.hpp"
#include "dca/platform/dca_gpu.h"
#include "dca/linalg/util/util_gpublas.hpp"

#include "dca/testing/gtest_h_w_warning_blocking.h"

TEST(DCAMemoryCopy, copy) {
  //  dca::linalg::util::initializeMagma();
  std::vector<std::complex<double>> some_complex(10, {1, 1});
  std::vector<double[2]> some_double2(10);
  std::for_each(some_double2.begin(), some_double2.end(), [](auto& d_arr) { d_arr[0] = 0.0; d_arr[1] = 2.0; });
  EXPECT_EQ(some_double2[0][1], 2);
  std::complex<double> expected{0,2};
  EXPECT_EQ(sizeof(std::remove_pointer<decltype(some_double2.data())>),sizeof(std::remove_pointer<decltype(some_complex.data())>));
  dca::linalg::util::memoryCopyCpu(some_complex.data(), some_double2.data(), 10);
  EXPECT_EQ(some_complex[9].imag(), 2);
}

