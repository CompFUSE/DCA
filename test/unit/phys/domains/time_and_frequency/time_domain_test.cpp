// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests time_domain.hpp.

#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include <vector>
#include <gtest/gtest.h>

using namespace dca::phys::domains;

TEST(TimeDomainTest, Basic) {
  const double beta = 4.;
  const int time_slices = 2;
  const double eps = 1.e-10;

  const int level = 1;  // 2^level = 2 elements per time slice.
  std::vector<double> weights;
  std::vector<double> nodes;

  EXPECT_EQ("time-domain", time_domain::get_name());
  EXPECT_FALSE(time_domain::is_initialized());
  EXPECT_THROW(time_domain::initialize_integration_domain(level, weights, nodes), std::logic_error);

  // Check initialization.
  time_domain::initialize(beta, time_slices, eps);

  EXPECT_TRUE(time_domain::is_initialized());
  EXPECT_EQ(6, time_domain::get_size());  // 6 = 2*(time_slices+1).

  const std::vector<double> elements_check{-beta + eps, -beta / time_slices, -eps,
                                           eps,         beta / time_slices,  beta - eps};
  EXPECT_EQ(elements_check, time_domain::get_elements());

  EXPECT_THROW(time_domain::initialize(beta, time_slices, eps), std::logic_error);

  // Check initialization of integration domain.
  time_domain::initialize_integration_domain(level, weights, nodes);

  EXPECT_EQ(weights.size(), nodes.size());
  EXPECT_EQ(8, nodes.size());  // 8 = (elements.size()-2) * 2^level.

  const std::vector<double> weights_check(8, beta);
  const std::vector<double> nodes_check{-beta + eps,
                                        (-beta + eps - beta / time_slices) / 2.,
                                        -beta / time_slices,
                                        (-beta / time_slices - eps) / 2.,
                                        eps,
                                        (eps + beta / time_slices) / 2.,
                                        beta / time_slices,
                                        (beta / time_slices + beta - eps) / 2.};

  EXPECT_EQ(weights_check, weights);
  EXPECT_EQ(nodes_check, nodes);
}
