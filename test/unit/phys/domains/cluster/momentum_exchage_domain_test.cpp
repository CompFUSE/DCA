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

using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "dca/phys/domains/cluster/momentum_exchange_domain.hpp"

#include <numeric>
#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/math/util/vector_operations.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

using dca::phys::domains::MomentumExchangeDomain;
using dca::math::util::isSameVector;

using Scalar = double;

constexpr char input_file[] = DCA_SOURCE_DIR "/test/unit/phys/domains/cluster/input.json";
using MomentumExchangeDomainTest =
  dca::testing::G0Setup<Scalar, dca::testing::LatticeSquare, dca::ClusterSolverId::CT_AUX, input_file>;


TEST_F(MomentumExchangeDomainTest, MultipleExchanges) {
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

