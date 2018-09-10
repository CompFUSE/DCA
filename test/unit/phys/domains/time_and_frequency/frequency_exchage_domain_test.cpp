// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests frequency_exchange_domain.hpp and momentum_exchange_domain.hpp.

#include "dca/phys/domains/time_and_frequency/frequency_exchange_domain.hpp"

#include <gtest/gtest.h>

#include "dca/math/util/vector_operations.hpp"

struct MockParameters {
  bool compute_all_transfers() const {
    return compute_all_transfers_;
  }
  int get_four_point_frequency_transfer() const {
    return fp_freq_exchange_;
  }

  bool compute_all_transfers_;
  int fp_freq_exchange_;
};

using dca::phys::domains::FrequencyExchangeDomain;
using dca::math::util::isSameVector;

TEST(FrequencyExchangeDomainTest, SingleExchange) {
  const int exchange_id = -2;
  MockParameters parameters{false, exchange_id};

  FrequencyExchangeDomain::initialize(parameters);
  EXPECT_TRUE(FrequencyExchangeDomain::isInitialized());

  EXPECT_EQ(1, FrequencyExchangeDomain::get_size());
  EXPECT_TRUE(isSameVector(std::vector<int>{-2}, FrequencyExchangeDomain::get_elements()));
  EXPECT_EQ(2, FrequencyExchangeDomain::get_extension_size());
}

TEST(FrequencyExchangeDomainTest, MultipleExchanges) {
  const int max_exchange_id = 3;
  MockParameters parameters{true, max_exchange_id};

  FrequencyExchangeDomain::initialize(parameters);
  EXPECT_TRUE(FrequencyExchangeDomain::isInitialized());

  EXPECT_EQ(max_exchange_id + 1, FrequencyExchangeDomain::get_size());
  EXPECT_TRUE(isSameVector(std::vector<int>{0, 1, 2, 3}, FrequencyExchangeDomain::get_elements()));
  EXPECT_EQ(3, FrequencyExchangeDomain::get_extension_size());

  // The maximum exchange frequency must be non-negative.
  parameters.fp_freq_exchange_ = -1;
  EXPECT_THROW(FrequencyExchangeDomain::initialize(parameters), std::logic_error);
}
