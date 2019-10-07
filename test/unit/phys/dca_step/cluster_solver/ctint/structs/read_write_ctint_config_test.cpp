// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the reading and writing of the CT-INT configuration to a buffer.

#include "dca/phys/dca_step/cluster_solver/ctint/structs/solver_configuration.hpp"

#include <random>
#include <functional>
#include "gtest/gtest.h"

#include "test/unit/phys/dca_step/cluster_solver/stub_rng.hpp"

TEST(ReadWriteConfigTest, All) {
  dca::phys::solver::ctint::InteractionVertices interactions;
  interactions.insertElement({{0, 0, 0, 0}, {0, 1, 2, 3}, 1});
  interactions.insertElement({{0, 0, 0, 0}, {0, 0, 2, 2}, 1});
  interactions.insertElement({{0, 0, 0, 0}, {1, 0, 3, 2}, 1});

  dca::phys::solver::ctint::SolverConfiguration config(1, 2, interactions, false);
  dca::phys::solver::ctint::SolverConfiguration config2(1, 2, interactions, false);

  std::mt19937_64 rng_internal(0);
  std::uniform_real_distribution<double> distro(0., 1.);
  auto rng = std::bind(distro, rng_internal);

  const int n = 15;
  for(int i =0; i < n; ++i)
    config.insertRandom(rng);

  dca::io::Buffer buff;
  buff << config;

  buff >> config2;

  EXPECT_EQ(config, config2);
}
