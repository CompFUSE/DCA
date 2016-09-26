// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests serial_collective_sum.hpp.

#include "dca/parallel/no_concurrency/serial_collective_sum.hpp"
#include "gtest/gtest.h"
#include "comp_library/function_library/domains/special_domains/dmn.h"
#include "comp_library/function_library/domains/special_domains/dmn_0.h"

class SerialCollectiveSumTest : public ::testing::Test {
protected:
  using Domain = dmn_0<dmn<3, int>>;
  dca::parallel::SerialCollectiveSum sum_interface_;
};

TEST_F(SerialCollectiveSumTest, Sum) {
  // Scalar
  double d = 3.14;
  sum_interface_.sum(d);
  EXPECT_EQ(3.14, d);

  // std::vector
  const std::vector<int> v_check{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<int> v(v_check);
  sum_interface_.sum(v);
  EXPECT_EQ(v_check, v);

  // std::map<std::string, std::vector>

  // function<Scalar, Domain>
  FUNC_LIB::function<double, Domain> f;
  f(0) = 3.14;
  f(1) = 2.72;
  f(2) = 42.;
  sum_interface_.sum(f);

  EXPECT_EQ(3.14, f(0));
  EXPECT_EQ(2.72, f(1));
  EXPECT_EQ(42., f(2));

  // function<Scalar, Domain> (separate f_in and f_out)

  // function<std::vector, Domain>

  // dca::linalg::Vector

  // dca::linalg::Matrix
}

TEST_F(SerialCollectiveSumTest, SumAndAverage) {
  const int measurements = 2000;

  // Scalar
  double d = 3.14;
  sum_interface_.sum_and_average(d, measurements);
  EXPECT_EQ(0.00157, d);

  // function<Scalar>
  FUNC_LIB::function<double, Domain> f;
  f(0) = 3.14;
  f(1) = 2.72;
  f(2) = 42.;

  sum_interface_.sum_and_average(f, measurements);

  EXPECT_EQ(0.00157, f(0));
  EXPECT_EQ(0.00136, f(1));
  EXPECT_EQ(0.021, f(2));
}
