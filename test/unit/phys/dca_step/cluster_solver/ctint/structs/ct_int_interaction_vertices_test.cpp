// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the InteractionElement class in the contest of a density-density interaction
// model.

#include "gtest/gtest.h"

#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"

#include "dca/math/random/random.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

TEST(InteractionVerticesTest, InsertElement) {
  dca::phys::solver::ctint::InteractionVertices H;
  using Element = dca::phys::solver::ctint::InteractionElement;

  H.insertElement(Element{{0, 0, 0, 0}, {0, 1, 2, 3}, 1.});  // pair-hop 1
  H.insertElement(Element{{0, 0, 0, 0}, {0, 0, 2, 2}, 1.});  // density-density
  H.insertElement(Element{{0, 0, 0, 0}, {1, 0, 3, 2}, 1.});  // pair-hop 2

  EXPECT_EQ(2, H[0].partner_id);
  EXPECT_EQ(-1, H[1].partner_id);
  EXPECT_EQ(0, H[2].partner_id);
}

TEST(InteractionVerticesTest, SamplingProb) {
  dca::phys::solver::ctint::InteractionVertices H;
  using dca::phys::solver::ctint::InteractionElement;
  const InteractionElement v1{{1, 0, 0, 0}, {0, 0, 0, 0}, 2};
  const InteractionElement v2{{2, 0, 0, 0}, {0, 0, 0, 0}, 2};
  const InteractionElement v3{{3, 0, 0, 0}, {0, 0, 0, 0}, 6};
  H.insertElement(v1), H.insertElement(v2), H.insertElement(v3);

  // Test if v1 and v2 are sampled with probability 1/5.
  // The interval [0,1) is mapped int he following way: [0,0.6) -> v3, [0.6,0.8) -> v2 and [0.8,1)
  // -> v1
  EXPECT_EQ(2, H.getInsertionIndices(0.3).first);
  EXPECT_EQ(1, H.getInsertionIndices(0.7).first);
  EXPECT_EQ(1, H.getInsertionIndices(0.79).first);
  EXPECT_EQ(0, H.getInsertionIndices(0.8).first);
}

// TEST(InteractionVertices, Statistic) {
//  dca::phys::solver::ctint::InteractionVertices vertices;
//  using dca::phys::solver::ctint::InteractionElement;
//  const InteractionElement zero{{3, 0, 0, 0}, {0, 0, 0, 0}};
//
//  vertices.insertElement(zero, 1);
//  vertices.insertElement(zero, 1);
//  vertices.insertElement(zero, 2);
//
//  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 0);
//  std::array<int, 3> results{0, 0, 0};
//  const int n_samples = 1e4;
//
//  // Count the number of vertices of each type.
//  for (int i = 0; i < n_samples; i++) {
//    const int res = &vertices.getRandomElement(rng()) - &vertices.getElements()[0];
//    results[res]++;
//  }
//
//  const std::array<double, 3> p{0.25, 0.25, 0.5};
//  for (int i = 0; i < p.size(); i++) {
//    const double var = std::sqrt(p[i] * (1 - p[i]) * n_samples);
//    EXPECT_NEAR(p[i], results[i] / n_samples, 3 * var);
//  }
//}

using G0Setup = dca::testing::G0Setup<dca::testing::LatticeHund>;
TEST_F(G0Setup, InitializeFromHamiltonians) {
  //  ****************************
  //  **  Density-density test  **
  //  ****************************
  {
    auto& H_int = G0Setup::data->H_interactions;
    const double U = 1;
    H_int = 0.;
    H_int(0, 0, 0, 1, 0) = U;
    H_int(0, 1, 0, 0, 0) = U;

    dca::phys::solver::ctint::InteractionVertices vertices;
    vertices.initializeFromHamiltonian(H_int, true /*double counting*/);
    EXPECT_EQ(U * RDmn::dmn_size(), vertices.integratedInteraction());

    vertices.reset();
    vertices.initializeFromHamiltonian(H_int, false /*double counting*/);
    EXPECT_EQ(2 * U * RDmn::dmn_size(), vertices.integratedInteraction());
  }

  //  ********************************
  //  **  Non density-density test  **
  //  ********************************
  // Note: H_nd is set as \sum_{i != j} Jh c+_{i, up} c_{j, up} *
  //                                     * (c+_{j, down} c_{i, down} +  c+_{i,down} c_{j, down})
  {
    using Nu = dca::func::dmn_variadic<BDmn, SDmn>;
    using Rdmn = G0Setup::RDmn;
    using Lattice = G0Setup::LatticeType;

    dca::func::function<double, dca::func::dmn_variadic<Nu, Nu, Nu, Nu, Rdmn>> H_nd;
    Lattice::initializeNonDensityInteraction(H_nd, G0Setup::parameters);

    dca::phys::solver::ctint::InteractionVertices vertices;
    vertices.initializeFromNonDensityHamiltonian(H_nd);

    const double Jh = G0Setup::parameters.get_Jh();
    // Double insertion is proposed with w**2 probability.
    EXPECT_EQ(Jh, vertices[0].w);
    EXPECT_NEAR(Jh * 16, vertices.integratedInteraction(), 1e-10);

    const auto& vertex = vertices[0];
    std::array<int, 4> expected{0, 1, 2, 3};  // 0up, 1up, 0down, 1up.
    for (int i = 0; i < expected.size(); i++)
      EXPECT_EQ(expected[i], vertex.nu[i]);
  }
}
