// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests dmn.hpp.

#include "dca/function/domains/dmn.hpp"
#include <vector>
#include "gtest/gtest.h"

TEST(DmnTest, WithoutElements) {
  // Default element type
  using Dmn1 = dca::func::dmn<3>;
  EXPECT_EQ(3, Dmn1::get_size());
  EXPECT_EQ(Dmn1::dmn_size(), Dmn1::get_size());
  EXPECT_EQ("dca::func::dmn<3, int>", Dmn1::get_name());
  EXPECT_THROW(Dmn1::get_elements(), std::logic_error);

  // Specify element type
  using Dmn2 = dca::func::dmn<42, double>;
  EXPECT_EQ(42, Dmn2::get_size());
  EXPECT_EQ("dca::func::dmn<42, double>", Dmn2::get_name());
  EXPECT_THROW(Dmn2::get_elements(), std::logic_error);
}

TEST(DmnTest, WithElements) {
  const std::vector<int> elements_1{1, 0, -3, 42};
  using Dmn1 = dca::func::dmn<4, int>;
  EXPECT_EQ(4, Dmn1::get_size());
  EXPECT_EQ("dca::func::dmn<4, int>", Dmn1::get_name());

  EXPECT_THROW(Dmn1::get_elements(), std::logic_error);
  Dmn1::set_elements(elements_1);
  EXPECT_NO_THROW(Dmn1::get_elements());
  EXPECT_EQ(elements_1, Dmn1::get_elements());

  // Try to set the elements of the domain with a vector of wrong size.
  const std::vector<float> elements_2{3.14, 1.23f, 42};
  using Dmn2 = dca::func::dmn<2, float>;
  EXPECT_THROW(Dmn2::set_elements(elements_2), std::logic_error);

  const std::vector<long> elements_3{123};
  using Dmn3 = dca::func::dmn<2, long>;
  EXPECT_THROW(Dmn3::set_elements(elements_3), std::logic_error);
}
