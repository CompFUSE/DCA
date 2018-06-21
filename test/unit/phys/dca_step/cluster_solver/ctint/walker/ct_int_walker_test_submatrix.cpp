// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Jérémie Bouquet (bouquetj@gmail.com).
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

using G0Setup =
    typename dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_INT>;
using namespace dca::phys::solver;
using SubmatrixWalker = testing::phys::solver::ctint::WalkerWrapperSubmatrix<G0Setup::Parameters>;
using Walker = testing::phys::solver::ctint::WalkerWrapper<G0Setup::Parameters>;
using Matrix = Walker::Matrix;
using MatrixPair = std::array<Matrix, 2>;

// Compare the submatrix update with a direct computation of the M matrix, and compare the acceptance probability to
// the CTINT walker with no submatrix update.
TEST_F(G0Setup, doSteps) {
  std::vector<double> setup_rngs(200);
  for (double& x : setup_rngs)
    x = double(std::rand()) / RAND_MAX;
  G0Setup::RngType rng(setup_rngs);

  ctint::G0Interpolation<dca::linalg::CPU> g0(
      dca::phys::solver::ctint::details::shrinkG0(data->G0_r_t));
  G0Setup::LabelDomain label_dmn;
  ctint::DMatrixBuilder<dca::linalg::CPU> builder(g0, RDmn::parameter_type::get_subtract_matrix(),
                                                  label_dmn.get_branch_domain_steps(),
                                                  parameters.getAlphas());

  // *******************************
  // Test vertex insert/removal ****
  // *******************************
  // Set rng values.
  //
  // Insertion, vertex_id, tau, aux_spin.
  // Removal, vertex_id.
  // ...
  std::vector<double> new_vals{
      0.3, 0.01, 0.454, 0.8,  // Insertion.
      1, 0.,                  // Remove pre-existing.
      1, 0.99,                // Remove recently inserted.
      0.1, 0.45, 0.934, 0.2,  //
      1, 0.5                  //
  };

  for (int steps = 1; steps <= 5; ++steps) {
    rng.setNewValues(setup_rngs);
    SubmatrixWalker walker(parameters, rng, G0Setup::interaction_vertices, builder);
    walker.forceAcceptance();

    MatrixPair old_M(walker.getM());
    rng.setNewValues(new_vals);
    walker.doStep(steps);
    MatrixPair new_M(walker.getM());

    // Compute directly the new M.
    walker.setMFromConfig();
    MatrixPair direct_M(walker.getM());

    for (int s = 0; s < 2; ++s)
      EXPECT_TRUE(dca::linalg::matrixop::areNear(direct_M[s], new_M[s], 1e-7));

    // Compare with non submatrix walker.
    rng.setNewValues(setup_rngs);
    Walker walker_nosub(parameters, rng, G0Setup::interaction_vertices, builder);
    walker_nosub.forceAcceptance();

    rng.setNewValues(new_vals);
    for (int i = 0; i < steps; ++i)
      walker_nosub.doStep();

    EXPECT_NEAR(walker.getAcceptanceProbability(), walker_nosub.getAcceptanceProbability(), 1e-5);
  }
}
