// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests frequency_domain.hpp.

#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

#include "gtest/gtest.h"

TEST(FrequencyDomainTest, Initialize) {
  using dca::phys::domains::frequency_domain;
  EXPECT_FALSE(frequency_domain::is_initialized());

  const double beta = 2.4;
  const int n_positive_frequencies = 3;
  frequency_domain::initialize(beta, n_positive_frequencies);
  const std::vector<double> expected_frequencies{-5 * M_PI / beta, -3 * M_PI / beta,
                                                 -M_PI / beta,     M_PI / beta,
                                                 3 * M_PI / beta,  5 * M_PI / beta};
  EXPECT_EQ(expected_frequencies, frequency_domain::get_elements());
  EXPECT_EQ(2 * n_positive_frequencies, frequency_domain::get_size());

  // The frequency domain can be initialized only once.
  EXPECT_THROW(frequency_domain::initialize(1, 1), std::logic_error);
}
