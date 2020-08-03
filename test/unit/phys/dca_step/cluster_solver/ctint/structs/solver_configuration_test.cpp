// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// CT-INT configuration test.

#include "dca/phys/dca_step/cluster_solver/ctint/structs/solver_configuration.hpp"

#include <functional>
#include <random>
#include <vector>

#include "gtest/gtest.h"

using dca::phys::solver::ctint::Vertex;
std::mt19937_64 rgen(0);
std::uniform_real_distribution<double> distro(0, 1);
auto rng = std::bind(distro, rgen);

TEST(SolverConfigurationTest, MoveAndShrink) {
  dca::phys::solver::ctint::InteractionVertices interactions;
  interactions.insertElement({{0, 0, 0, 0}, {0, 0, 1, 1}, 1});  // up-down
  interactions.insertElement({{1, 1, 1, 1}, {0, 0, 1, 1}, 1});
  interactions.insertElement({{2, 2, 2, 2}, {0, 0, 1, 1}, 1});

  using HostVector = dca::linalg::util::HostVector<int>;
  std::array<HostVector, 2> indices_source;
  std::array<HostVector, 2> indices_remove;

  auto execute_test = [&](int conf_size, std::vector<int> remove) {
    dca::phys::solver::ctint::SolverConfiguration config(1, 1, interactions);
    for (int i = 0; i < conf_size; ++i) {
      config.insertRandom(rng);
    }

    EXPECT_EQ(conf_size, config.size());
    EXPECT_TRUE(config.checkConsistency());

    std::set<std::uint64_t> removed_tags(remove.begin(), remove.end());

    config.moveAndShrink(indices_source, indices_remove, remove);
    EXPECT_EQ(config.size(), conf_size - remove.size());

    for (int i = 0; i < config.size(); ++i) {
      const auto tag = config[i].tag;
      EXPECT_LT(tag, conf_size);
      EXPECT_EQ(0, removed_tags.count(tag));
    }

    EXPECT_TRUE(config.checkConsistency());
  };

  execute_test(4, {0, 2});
  execute_test(8, {0, 1, 2, 3, 4});
  execute_test(4, {0, 1, 2});
}

TEST(SolverConfigurationTest, randomRemovalCandidate) {
  dca::phys::solver::ctint::InteractionVertices interactions;
  interactions.insertElement({{0, 0, 0, 0}, {0, 0, 1, 1}, 1});  // up-down

  dca::phys::solver::ctint::SolverConfiguration config(1, 1, interactions);

  std::set<uint64_t> tags;
  for (uint64_t tag = 0; tag < 20; ++tag) {
    config.insertRandom(rng);
    config.commitInsertion(tag);
    tags.insert(tag);
  }

  auto execute_test = [&](int n_attempts) {
    for (int i = 0; i < n_attempts; ++i) {
      const auto candidates = config.randomRemovalCandidate(rng());
      EXPECT_EQ(candidates[1], -1);
      const auto tag = config[candidates[0]].tag;
      EXPECT_EQ(tags.count(tag), 1);
    }
  };

  execute_test(10);

  // Remove some indices from the possible candidates.
  for (int idx : {0, 3, 7, 4, 2, 19}) {
    config.markForRemoval(idx);
    tags.erase(idx);
  }

  execute_test(10);

  // Remove more indices to have less than 10 candidates.
  for (int idx : {1, 5, 6, 8, 9, 10, 11}) {
    config.markForRemoval(idx);
    tags.erase(idx);
  }

  execute_test(10);
}
