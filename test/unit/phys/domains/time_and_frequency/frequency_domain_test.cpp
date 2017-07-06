// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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

TEST(FrequencyDomainTest, Basic) {
  const double beta = 2.4;
  const int num_freqs = 3;

  const std::vector<double> elements_expected{-5 * M_PI / beta, -3 * M_PI / beta, -M_PI / beta,
                                              M_PI / beta,      3 * M_PI / beta,  5 * M_PI / beta};
  const std::vector<int> integer_wave_vectors_expected{-5, -3, -1, 1, 3, 5};

  EXPECT_FALSE(frequency_domain::is_initialized());

  frequency_domain::initialize(beta, num_freqs);

  EXPECT_TRUE(frequency_domain::is_initialized());
  EXPECT_EQ(2 * num_freqs, frequency_domain::get_size());
  EXPECT_EQ("frequency-domain", frequency_domain::get_name());

  EXPECT_EQ(elements_expected, frequency_domain::get_elements());
  EXPECT_EQ(integer_wave_vectors_expected, frequency_domain::get_integer_wave_vectors());

  EXPECT_EQ(2.*M_PI/beta, frequency_domain::get_basis()[0]);
  EXPECT_EQ(beta/(2.*M_PI), frequency_domain::get_inverse_basis()[0]);

  // The frequency domain can only be initialized once.
  EXPECT_THROW(frequency_domain::initialize(1, 1), std::logic_error);
}
