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
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

using dca::phys::domains::MomentumExchangeDomain;
using dca::math::util::isSameVector;

constexpr char input_file[] = DCA_SOURCE_DIR "/test/unit/phys/domains/cluster/input.json";
using MomentumExchangeDomainTest =
    dca::testing::G0Setup<dca::testing::LatticeSquare, dca::phys::solver::CT_AUX, input_file>;


TEST_F(MomentumExchangeDomainTest, MultipleExchanges) {
  MomentumExchangeDomain::initialize(parameters_);
  EXPECT_TRUE(MomentumExchangeDomain::isInitialized());

  // Cluster sites with their id:
  // 12 13 14 15
  // 08 09 10 11
  // 04 05 06 07
  // 00 01 02 03
  //
  // Independent sites due to rotations and reflections: 00, 01, 02, 05, 06, 10.

  EXPECT_EQ(6, MomentumExchangeDomain::get_size());
  const std::vector<int> expected_elements{00, 01, 02, 05, 06, 10};
  EXPECT_TRUE(expected_elements == MomentumExchangeDomain::get_elements());
}
