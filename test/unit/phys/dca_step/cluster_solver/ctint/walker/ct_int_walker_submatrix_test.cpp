// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Jérémie Bouquet   (bouquetj@gmail.com).
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch).
//
// This class tests the CPU walker used by the ctint cluster solver. The fast updated matrix
// are compared with their direct computation.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"
#include "gtest/gtest.h"

#include "walker_wrapper.hpp"
#include "walker_wrapper_submatrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr char input_name[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctint/walker/submatrix_input.json";

using G0Setup =
    typename dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_INT, input_name>;
using namespace dca::phys::solver;
using SubmatrixWalker = testing::phys::solver::ctint::WalkerWrapperSubmatrix<G0Setup::Parameters>;
using Walker = testing::phys::solver::ctint::WalkerWrapper<G0Setup::Parameters>;
using Matrix = Walker::Matrix;
using MatrixPair = std::array<Matrix, 2>;

// Compare the submatrix update with a direct computation of the M matrix, and compare the
// acceptance probability to
// the CTINT walker with no submatrix update.
TEST_F(G0Setup, doSteps) {
  std::vector<double> setup_rngs{0., 0.00, 0.9,  0.5, 0.01, 0,    0.75, 0.02,
                                 0,  0.6,  0.03, 1,   0.99, 0.04, 0.99};
  G0Setup::RngType rng(setup_rngs);

  ctint::G0Interpolation<dca::linalg::CPU> g0(
      dca::phys::solver::ctint::details::shrinkG0(data_->G0_r_t));
  G0Setup::LabelDomain label_dmn;
  Walker::setDMatrixBuilder(g0, RDmn::parameter_type::get_subtract_matrix(),
                            label_dmn.get_branch_domain_steps());
  Walker::setDMatrixAlpha(parameters_.getAlphas(), false);
  Walker::setInteractionVertices(parameters_, *data_);

  // ************************************
  // Test vertex insertion / removal ****
  // ************************************
  // Set rng values.
  //
  // Insertion, vertex_id, tau, aux_spin, acceptance_rng
  // Removal, vertex_id, acceptance_rng
  // ...
  // Note: if acceptance_rng <= 0 the move is always accepted, if it is > 1 the move is always
  // rejected.
  const std::vector<double> rng_vals{
      0, 0,    0.1, 0.8, -1,  // Insertion.
      0, 0.99, 0.2, 0.8, -1,  // Insertion.
      0, 0,    0.3, 0.8, 2,   // Insertion. Rejected.
      1, 0,    -1,            // Remove pre-existing.
      1, 0.99, -1,            // Remove recently inserted.
      1, 0.99, 2,             // Remove recently inserted. Rejected
      1, 0,    2,             // Remove . Rejected
      0, 0.99, 0.4, 0.2, -1,  // Insertion
  };

  // NOTE: the non-submatrix walker reorders the configuration after each removal, therefore we need
  // to change the vertex_id value for some removals, in order to obtain the same results.
  const std::vector<double> rng_nosub_vals{
      0, 0,    0.1, 0.8, -1,  // Insertion.
      0, 0.99, 0.2, 0.8, -1,  // Insertion.
      0, 0,    0.3, 0.8, 2,   // Insertion. Rejected.
      1, 0,    -1,            // Remove tau = 0.00
      1, 0.,   -1,            // Remove tau = 0.2
      1, 0.99, 2,             // Remove tau = 0.1. Rejected
      1, 0.2,  2,             // Remove tau = 0.01 . Rejected
      0, 0.99, 0.4, 0.2, -1,  // Insertion
  };

  for (int steps = 1; steps <= 8; ++steps) {
    rng.setNewValues(setup_rngs);
    SubmatrixWalker walker(parameters_, rng);

    MatrixPair old_M(walker.getM());
    rng.setNewValues(rng_vals);
    walker.doStep(steps);
    MatrixPair new_M(walker.getM());
    auto config = walker.getWalkerConfiguration();

    // Compute directly the new M.
    walker.setMFromConfig();
    MatrixPair direct_M(walker.getM());

    for (int s = 0; s < 2; ++s)
      EXPECT_TRUE(dca::linalg::matrixop::areNear(direct_M[s], new_M[s], 1e-7));

    // Compare with non submatrix walker.
    rng.setNewValues(setup_rngs);
    Walker walker_nosub(parameters_, rng);

    rng.setNewValues(rng_nosub_vals);
    for (int i = 0; i < steps; ++i)
      walker_nosub.doStep();

    EXPECT_NEAR(walker.getAcceptanceProbability(), walker_nosub.getAcceptanceProbability(), 1e-5);

    auto config_nosubm = walker_nosub.getWalkerConfiguration();
    ASSERT_EQ(config.size(), config_nosubm.size());
    // The final configuration is the same up to a permutation.
    for (int i = 0; i < config.size(); ++i) {
      bool found = false;
      for (int j = 0; j < config_nosubm.size(); ++j)
        if (config[i] == config_nosubm[j]) {
          found = true;
          break;
        }
      EXPECT_TRUE(found);
    }
  }
}
