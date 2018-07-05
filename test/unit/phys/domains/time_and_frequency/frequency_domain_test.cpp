// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests frequency_domain.hpp.

#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

#include <cmath>  // M_PI
#include <stdexcept>
#include <vector>

#include "gtest/gtest.h"

using dca::phys::domains::frequency_domain;

TEST(FrequencyDomainTest, Full) {
  const double beta = 2.4;
  const int num_freqs = 3;

  const std::vector<int> indices_expected{-5, -3, -1, 1, 3, 5};
  const std::vector<double> elements_expected{-5 * M_PI / beta, -3 * M_PI / beta, -M_PI / beta,
                                              M_PI / beta,      3 * M_PI / beta,  5 * M_PI / beta};

  EXPECT_FALSE(frequency_domain::is_initialized());
  EXPECT_EQ("frequency-domain", frequency_domain::get_name());
  EXPECT_DEBUG_DEATH(frequency_domain::get_size(), "initialized_");
  EXPECT_DEBUG_DEATH(frequency_domain::get_elements(), "initialized_");
  EXPECT_DEBUG_DEATH(frequency_domain::get_indices(), "initialized_");

  frequency_domain::initialize(beta, num_freqs);

  EXPECT_TRUE(frequency_domain::is_initialized());
  EXPECT_EQ("frequency-domain", frequency_domain::get_name());
  EXPECT_EQ(2 * num_freqs, frequency_domain::get_size());
  EXPECT_EQ(elements_expected, frequency_domain::get_elements());
  EXPECT_EQ(indices_expected, frequency_domain::get_indices());

  // The domain can only be initialized once.
  EXPECT_THROW(frequency_domain::initialize(0.9, 10), std::logic_error);
}
