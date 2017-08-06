// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests time_domain_left_oriented.hpp.

#include "dca/phys/domains/time_and_frequency/time_domain_left_oriented.hpp"
#include <vector>
#include <gtest/gtest.h>
#include "dca/phys/domains/time_and_frequency/time_domain_left_oriented.hpp"

using namespace dca::phys::domains;

TEST(TimeDomainLeftOriented, Basic) {
  EXPECT_EQ("time-domain-left-oriented", time_domain_left_oriented::get_name());
  EXPECT_FALSE(time_domain_left_oriented::is_initialized());
  EXPECT_THROW(time_domain_left_oriented::initialize(), std::logic_error);

  // First we need to initialize time_domain.
  const double beta = 4.;
  const int time_slices = 2;
  const double eps = 1.e-10;
  time_domain::initialize(beta, time_slices, eps);

  // Check initialization of time_domain_left_oriented.
  time_domain_left_oriented::initialize();

  EXPECT_TRUE(time_domain_left_oriented::is_initialized());
  EXPECT_EQ(time_domain::get_size() - 2, time_domain_left_oriented::get_size());

  const std::vector<double> elements_check{-beta + eps, -beta / time_slices, 0, beta / time_slices};
  EXPECT_EQ(elements_check, time_domain_left_oriented::get_elements());

  EXPECT_THROW(time_domain_left_oriented::initialize(), std::logic_error);
}
