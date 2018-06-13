// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#include "gtest/gtest.h"
#include <vector>

#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration.hpp"

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

  config.pop();
  EXPECT_EQ(2, config.size());
}

TEST(SolverConfigurationTest, MatrixConfigurationUpdate){
  dca::phys::solver::ctint::InteractionVertices interactions;
  interactions.insertElement({{0, 0, 0, 0}, {0, 0, 1, 1}, 1}); // up-up
  interactions.insertElement({{0, 0, 0, 0}, {2, 2, 3, 3}, 1}); // down-down
  interactions.insertElement({{0, 0, 0, 0}, {0, 0, 3, 3}, 1}); // up-down
  interactions.insertElement({{0, 0, 0, 0}, {2, 2, 1, 1}, 1}); // down-up

  dca::phys::solver::ctint::SolverConfiguration config(1, 2, interactions);
  using dca::phys::solver::ctint::Vertex;
  config.push_back(Vertex{0, 2, 0}); // up-down
  config.push_back(Vertex{0, 3, 0}); // down-up
  config.push_back(Vertex{0, 0, 0}); // up-up

  EXPECT_EQ(4, config.getSector(0).size());
  EXPECT_EQ(2, config.getSector(1).size());

  config.swapVertices(1,2);
  config.pop(); // pops a down-up vertex

  EXPECT_EQ(3, config.getSector(0).size());
  EXPECT_EQ(1, config.getSector(1).size());
}

