// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests get_bounds.hpp.

#include "dca/parallel/util/get_bounds.hpp"

#include "gtest/gtest.h"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"

TEST(GetBoundsTest, Domain) {
  const int num_threads = 3;

  // number of threads < domain-size
  using Domain1 = dca::func::dmn_0<dca::func::dmn<11, int>>;
  Domain1 dmn1;

  std::pair<int, int> bounds_0 = dca::parallel::util::getBounds(0, num_threads, dmn1);
  EXPECT_EQ(0, bounds_0.first);
  EXPECT_EQ(3, bounds_0.second);

  std::pair<int, int> bounds_1 = dca::parallel::util::getBounds(1, num_threads, dmn1);
  EXPECT_EQ(3, bounds_1.first);
  EXPECT_EQ(7, bounds_1.second);

  std::pair<int, int> bounds_2 = dca::parallel::util::getBounds(2, num_threads, dmn1);
  EXPECT_EQ(7, bounds_2.first);
  EXPECT_EQ(11, bounds_2.second);

  // number of threads = domain-size
  using Domain2 = dca::func::dmn_0<dca::func::dmn<3, int>>;
  Domain2 dmn2;

  bounds_0 = dca::parallel::util::getBounds(0, num_threads, dmn2);
  EXPECT_EQ(0, bounds_0.first);
  EXPECT_EQ(1, bounds_0.second);

  bounds_1 = dca::parallel::util::getBounds(1, num_threads, dmn2);
  EXPECT_EQ(1, bounds_1.first);
  EXPECT_EQ(2, bounds_1.second);

  bounds_2 = dca::parallel::util::getBounds(2, num_threads, dmn2);
  EXPECT_EQ(2, bounds_2.first);
  EXPECT_EQ(3, bounds_2.second);

  // number of threads > domain-size
  using Domain3 = dca::func::dmn_0<dca::func::dmn<2, int>>;
  Domain3 dmn3;

  bounds_0 = dca::parallel::util::getBounds(0, num_threads, dmn3);
  EXPECT_EQ(0, bounds_0.first);
  EXPECT_EQ(1, bounds_0.second);

  bounds_1 = dca::parallel::util::getBounds(1, num_threads, dmn3);
  EXPECT_EQ(1, bounds_1.first);
  EXPECT_EQ(2, bounds_1.second);

  bounds_2 = dca::parallel::util::getBounds(2, num_threads, dmn3);
  EXPECT_EQ(-1, bounds_2.first);
  EXPECT_EQ(-1, bounds_2.second);
}

TEST(GetBoundsTest, Bounds) {
  const int num_threads = 3;

  // number of threads < current range
  const std::pair<int, int> current_bounds_1(3, 14);

  std::pair<int, int> bounds_0 = dca::parallel::util::getBounds(0, num_threads, current_bounds_1);
  EXPECT_EQ(3, bounds_0.first);
  EXPECT_EQ(6, bounds_0.second);

  std::pair<int, int> bounds_1 = dca::parallel::util::getBounds(1, num_threads, current_bounds_1);
  EXPECT_EQ(6, bounds_1.first);
  EXPECT_EQ(10, bounds_1.second);

  std::pair<int, int> bounds_2 = dca::parallel::util::getBounds(2, num_threads, current_bounds_1);
  EXPECT_EQ(10, bounds_2.first);
  EXPECT_EQ(14, bounds_2.second);

  // number of threads = current range
  const std::pair<int, int> current_bounds_2(-1, 2);

  bounds_0 = dca::parallel::util::getBounds(0, num_threads, current_bounds_2);
  EXPECT_EQ(-1, bounds_0.first);
  EXPECT_EQ(0, bounds_0.second);

  bounds_1 = dca::parallel::util::getBounds(1, num_threads, current_bounds_2);
  EXPECT_EQ(0, bounds_1.first);
  EXPECT_EQ(1, bounds_1.second);

  bounds_2 = dca::parallel::util::getBounds(2, num_threads, current_bounds_2);
  EXPECT_EQ(1, bounds_2.first);
  EXPECT_EQ(2, bounds_2.second);

  // number of threads > current range
  const std::pair<int, int> current_bounds_3(0, 2);

  bounds_0 = dca::parallel::util::getBounds(0, num_threads, current_bounds_3);
  EXPECT_EQ(0, bounds_0.first);
  EXPECT_EQ(1, bounds_0.second);

  bounds_1 = dca::parallel::util::getBounds(1, num_threads, current_bounds_3);
  EXPECT_EQ(1, bounds_1.first);
  EXPECT_EQ(2, bounds_1.second);

  bounds_2 = dca::parallel::util::getBounds(2, num_threads, current_bounds_3);
  EXPECT_EQ(-1, bounds_2.first);
  EXPECT_EQ(-1, bounds_2.second);
}
