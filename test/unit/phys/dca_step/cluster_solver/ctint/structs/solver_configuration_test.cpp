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

#include "test/unit/phys/dca_step/cluster_solver/stub_rng.hpp"

TEST(SolverConfigurationTest, InsertAndSwap) {
  dca::phys::solver::ctint::InteractionVertices interactions;
  interactions.insertElement({{0, 0, 0, 0}, {0, 1, 2, 3}, 1});
  interactions.insertElement({{0, 0, 0, 0}, {0, 0, 2, 2}, 1});
  interactions.insertElement({{0, 0, 0, 0}, {1, 0, 3, 2}, 1});

  dca::phys::solver::ctint::SolverConfiguration config(1, 2, interactions, 1);
  using Vector = std::vector<double>;
  // Select numbers for: first vertex(ndd), tau, aux_spin, double insertion, tau2, aux_spin2.
  dca::testing::StubRng rng(Vector{0.9, 0.66, 0.2, 0, 0.66, 0.2});

  config.insertRandom(rng);
  EXPECT_EQ(2, config.size());
  EXPECT_EQ(2, config.lastInsertionSize());
  EXPECT_EQ(0, config[0].interaction_id);

  // Select density-density
  rng.setNewValues(Vector{0.5, 0, 0});
  config.insertRandom(rng);
  EXPECT_EQ(1, config.lastInsertionSize());
  EXPECT_EQ(3, config.size());
  EXPECT_EQ(1, config[2].interaction_id);

  config.swapVertices(0, 1);
  EXPECT_EQ(2, config[0].interaction_id);
  EXPECT_EQ(0, config[1].interaction_id);

  // Pop is used by the class CtintWalker after the matrix configuration has been reordered
  // accordingly. Using it right after a reorder of the configuration will leave it inconsistent.
  config.pop();
  EXPECT_EQ(2, config.size());
  EXPECT_FALSE(config.checkConsistency());
}

TEST(SolverConfigurationTest, MatrixConfigurationUpdate) {
  dca::phys::solver::ctint::InteractionVertices interactions;
  interactions.insertElement({{0, 0, 0, 0}, {0, 0, 1, 1}, 1});  // up-up
  interactions.insertElement({{0, 0, 0, 0}, {2, 2, 3, 3}, 1});  // down-down
  interactions.insertElement({{0, 0, 0, 0}, {0, 0, 3, 3}, 1});  // up-down
  interactions.insertElement({{0, 0, 0, 0}, {2, 2, 1, 1}, 1});  // down-up

  dca::phys::solver::ctint::SolverConfiguration config(1, 2, interactions);
  using dca::phys::solver::ctint::Vertex;
  std::uint64_t tag = 0;
  config.push_back(Vertex{0, 2, ++tag, 0});  // up-down
  config.push_back(Vertex{0, 3, ++tag, 0});  // down-up
  config.push_back(Vertex{0, 0, ++tag, 0});  // up-up

  EXPECT_EQ(4, config.getSector(0).size());
  EXPECT_EQ(2, config.getSector(1).size());

  config.swapVertices(1, 2);
  config.pop();  // pops a down-up vertex

  EXPECT_EQ(3, config.getSector(0).size());
  EXPECT_EQ(1, config.getSector(1).size());

  // See above test regarding the consistency of the configuration.
  // TODO: always leave config consitent.
  EXPECT_FALSE(config.checkConsistency());
}

TEST(SolverConfigurationTest, MoveAndShrink) {
  dca::phys::solver::ctint::InteractionVertices interactions;
  interactions.insertElement({{0, 0, 0, 0}, {0, 0, 1, 1}, 1});  // up-down
  interactions.insertElement({{1, 1, 1, 1}, {0, 0, 1, 1}, 1});
  interactions.insertElement({{2, 2, 2, 2}, {0, 0, 1, 1}, 1});

  dca::phys::solver::ctint::SolverConfiguration config(1, 1, interactions);
  using dca::phys::solver::ctint::Vertex;
  std::mt19937_64 rgen(0);
  std::uniform_real_distribution<double> distro(0, 1);
  auto rng = std::bind(distro, rgen);

  for (int i = 0; i < 4; ++i)
    config.insertRandom(rng);

  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(config.getTag(i), config.getSector(0).getTag(i));
    EXPECT_EQ(config.getTag(i), config.getSector(1).getTag(i));
  }

  std::vector<int> from{0};
  std::vector<int> to{0, 2};

  std::array<std::vector<int>, 2> from_sectors{from, from};
  std::array<std::vector<int>, 2> to_sectors{to, to};

  config.moveAndShrink(from_sectors, to_sectors, from, to);

  EXPECT_EQ(config.size(), 2);
  for (int i = 0; i < config.size(); ++i) {
    EXPECT_EQ(config.getTag(i), config.getSector(0).getTag(i));
    EXPECT_EQ(config.getTag(i), config.getSector(1).getTag(i));
  }

  EXPECT_TRUE(config.checkConsistency());
}
