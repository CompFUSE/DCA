// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests momentum_exchange_domain.hpp.

#include "dca/phys/domains/cluster/momentum_exchange_domain.hpp"

#include <numeric>
#include <gtest/gtest.h>

#include "dca/math/util/vector_operations.hpp"

constexpr int cluster_size = 6;

struct MockParameters {
  struct KClusterDmn {
    static int dmn_size() {
      return cluster_size;
    }
  };

  bool compute_all_transfers() const {
    return compute_all_transfers_;
  }
  int get_four_point_momentum_transfer_index() const {
    return transfer_index_;
  }

  bool compute_all_transfers_ = false;
  int transfer_index_ = -1;
};

using dca::phys::domains::MomentumExchangeDomain;
using dca::math::util::isSameVector;

TEST(MomentumExchangeDomainTest, SingleExchange) {
  const int exchange_id = -2;
  MockParameters parameters{false, exchange_id};

  MomentumExchangeDomain::initialize(parameters);
  EXPECT_TRUE(MomentumExchangeDomain::isInitialized());

  EXPECT_EQ(1, MomentumExchangeDomain::get_size());
  EXPECT_TRUE(isSameVector(std::vector<int>{-2}, MomentumExchangeDomain::get_elements()));
}

TEST(MomentumExchangeDomainTest, MultipleExchanges) {
  MockParameters parameters{true};

  MomentumExchangeDomain::initialize(parameters);
  EXPECT_TRUE(MomentumExchangeDomain::isInitialized());

  EXPECT_EQ(cluster_size, MomentumExchangeDomain::get_size());

  std::vector<int> expected_elements(cluster_size);
  std::iota(expected_elements.begin(), expected_elements.end(), 0);
  EXPECT_TRUE(isSameVector(expected_elements, MomentumExchangeDomain::get_elements()));
}
