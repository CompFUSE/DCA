// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests  ct_int_matrix_configuration.hpp.

#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"

#include "gtest/gtest.h"

using dca::phys::solver::ctint::InteractionVertices;
using dca::phys::solver::ctint::InteractionElement;
using dca::phys::solver::ctint::Vertex;

class MatrixConfigurationWrapper : public dca::phys::solver::ctint::MatrixConfiguration {
public:
  using BaseClass = dca::phys::solver::ctint::MatrixConfiguration;

  MatrixConfigurationWrapper(const InteractionVertices* H_int, int bands)
      : BaseClass(H_int, bands) {}
  using BaseClass::addVertex;
  using BaseClass::findIndices;
  using BaseClass::pop;
};

TEST(MatrixConfigurationTest, InsertAndSwap) {
  InteractionVertices vertices;
  vertices.insertElement(InteractionElement{{0, 0, 0, 0}, {0, 0, 2, 2}, 1});  // up-down
  vertices.insertElement(InteractionElement{{0, 0, 0, 0}, {0, 0, 1, 1}, 1});  // up-up
  MatrixConfigurationWrapper config(&vertices, 2);

  const double tau0(0), tau1(0.1);
  std::uint64_t tag1(0), tag2(1);
  Vertex v1{0, 0, tag1, tau0};  // based on first element
  Vertex v2{0, 1, tag2, tau1};  // based on second element

  config.addVertex(v1);
  config.addVertex(v2);
  EXPECT_EQ(3, config.size(0));
  EXPECT_EQ(1, config.size(1));

  auto indices_up = config.findIndices(tag2, 0);
  auto indices_down = config.findIndices(tag2, 1);
  const std::vector<int> expected_up{1, 2};
  EXPECT_EQ(expected_up, indices_up);
  EXPECT_EQ(0, indices_down.size());

  config.swapSectorLabels(0, 2, 0);  // push up-down vertex to the end.
  config.pop(1, 1);
  EXPECT_EQ(2, config.size(0));
  EXPECT_EQ(0, config.size(1));

  indices_up = config.findIndices(tag2, 0);
  const std::vector<int> expected_up2{0, 1};
  EXPECT_EQ(expected_up2, indices_up);
}
