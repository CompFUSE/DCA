// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests set_to_zero.hpp.

#include "dca/function/set_to_zero.hpp"
#include <complex>
#include "gtest/gtest.h"

TEST(SetToZeroTest, FundamentalTypes) {
  int i = 42;
  dca::func::setToZero(i);
  EXPECT_EQ(0, i);

  double d = 3.14;
  dca::func::setToZero(d);
  EXPECT_EQ(0., d);
}

TEST(SetToZeroTest, Complex) {
  std::complex<double> c(1.2, 3.4);
  dca::func::setToZero(c);
  EXPECT_EQ(std::complex<double>(0., 0.), c);
}

TEST(SetToZeroTest, CustomTypes) {
  struct MyStruct {
    int a = -42;
    double b = 9.9;
  } my_struct;

  dca::func::setToZero(my_struct);
  EXPECT_EQ(-42, my_struct.a);
  EXPECT_EQ(9.9, my_struct.b);
}
