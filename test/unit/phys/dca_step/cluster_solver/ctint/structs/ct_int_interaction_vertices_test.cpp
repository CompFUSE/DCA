// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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

  EXPECT_EQ(std::vector<ushort>{2}, H[0].partners_id);
  EXPECT_EQ(std::vector<ushort>{}, H[1].partners_id);
  EXPECT_EQ(std::vector<ushort>{0}, H[2].partners_id);
}

TEST(InteractionVerticesTest, SamplingProb) {
  dca::phys::solver::ctint::InteractionVertices H;
  using dca::phys::solver::ctint::InteractionElement;
  const InteractionElement v1{{1, 0, 0, 0}, {0, 0, 0, 0}, 2};
  const InteractionElement v2{{2, 0, 0, 0}, {0, 0, 0, 0}, 2};
  const InteractionElement v3{{3, 0, 0, 0}, {0, 0, 0, 0}, 6};
  H.insertElement(v1), H.insertElement(v2), H.insertElement(v3);

  // TODO: re-enable.
  // Test if v1 and v2 are sampled with probability 1/5.
  // The interval [0,1) is mapped int he following way: [0,0.6) -> v3, [0.6,0.8) -> v2 and [0.8,1)
  // -> v1
  //  EXPECT_EQ(2, H.getInsertionIndices(0.3).first);
  //  EXPECT_EQ(1, H.getInsertionIndices(0.7).first);
  //  EXPECT_EQ(1, H.getInsertionIndices(0.79).first);
  //  EXPECT_EQ(0, H.getInsertionIndices(0.8).first);
}

using G0Setup = dca::testing::G0Setup<dca::testing::LatticeHund>;
TEST_F(G0Setup, InitializeFromHamiltonians) {
  //  ****************************
  //  **  Density-density test  **
  //  ****************************
  {
    auto& H_int = G0Setup::data_->H_interactions;
    const double U = 1;
    H_int = 0.;
    H_int(0, 0, 0, 1, 0) = U;
    H_int(0, 1, 0, 0, 0) = U;

    dca::phys::solver::ctint::InteractionVertices vertices;
    vertices.initializeFromHamiltonian(H_int);
    EXPECT_EQ(U * RDmn::dmn_size(), vertices.integratedInteraction());

    // Note: use this check if the double counting is removed.
    // EXPECT_EQ(2 * RDmn::dmn_size(), vertices.integratedInteraction());
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
    Lattice::initializeNonDensityInteraction(H_nd, G0Setup::parameters_);

    dca::phys::solver::ctint::InteractionVertices vertices;
    vertices.initializeFromNonDensityHamiltonian(H_nd);

    const double Jh = G0Setup::parameters_.get_Jh();
    EXPECT_EQ(Jh, vertices[0].w);
    EXPECT_NEAR(Jh * 16, vertices.integratedInteraction(), 1e-10);

    const auto& vertex = vertices[0];
    std::array<int, 4> expected{0, 1, 2, 3};  // 0up, 1up, 0down, 1up.
    for (int i = 0; i < expected.size(); i++)
      EXPECT_EQ(expected[i], vertex.nu[i]);
  }
}
