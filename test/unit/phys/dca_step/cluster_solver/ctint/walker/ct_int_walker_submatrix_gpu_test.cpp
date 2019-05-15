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
// This class tests the GPU walker used by the ctint cluster solver by comparing it with the CPU
// version.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_gpu_submatrix.hpp"

#include "gtest/gtest.h"

#include "dca/linalg/matrixop.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"
#include "walker_wrapper_submatrix.hpp"

using dca::linalg::CPU;
using dca::linalg::GPU;

constexpr char input_name[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctint/walker/submatrix_input.json";

using G0Setup =
    typename dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_INT, input_name>;
using namespace dca::phys::solver;
template <dca::linalg::DeviceType device>
using SubmatrixWalker =
    testing::phys::solver::ctint::WalkerWrapperSubmatrix<G0Setup::Parameters, device>;
using Matrix = dca::linalg::Matrix<double, CPU>;
using MatrixPair = std::array<Matrix, 2>;

// Compare the submatrix update with a direct computation of the M matrix, and compare the
// acceptance probability to
// the CTINT walker with no submatrix update.
TEST_F(G0Setup, doSteps) {
  std::vector<double> setup_rngs{0., 0.00, 0.9,  0.5, 0.01, 0,    0.75, 0.02,
                                 0,  0.6,  0.03, 1,   0.99, 0.04, 0.99};
  G0Setup::RngType rng(setup_rngs);

  const auto g0_func = dca::phys::solver::ctint::details::shrinkG0(data_->G0_r_t);
  ctint::G0Interpolation<CPU> g0_cpu(g0_func);
  ctint::G0Interpolation<GPU> g0_gpu(g0_func);
  G0Setup::LabelDomain label_dmn;

  // TODO: improve API.
  SubmatrixWalker<CPU>::setDMatrixBuilder(g0_cpu, RDmn::parameter_type::get_subtract_matrix(),
                                          label_dmn.get_branch_domain_steps(),
                                          parameters_.getAlphas());
  SubmatrixWalker<GPU>::setDMatrixBuilder(g0_gpu, RDmn::parameter_type::get_subtract_matrix(),
                                          label_dmn.get_branch_domain_steps(),
                                          parameters_.getAlphas());
  SubmatrixWalker<CPU>::setInteractionVertices(parameters_, *data_);
  SubmatrixWalker<GPU>::setInteractionVertices(parameters_, *data_);


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

  for (const int initial_size : std::array<int, 2>{0, 5}) {
    parameters_.setInitialConfigurationSize(initial_size);

    for (int steps = 1; steps <= 8; ++steps) {
      rng.setNewValues(setup_rngs);
      SubmatrixWalker<CPU> walker_cpu(parameters_, rng);
      rng.setNewValues(setup_rngs);
      SubmatrixWalker<GPU> walker_gpu(parameters_, rng);

      rng.setNewValues(rng_vals);
      walker_cpu.doStep(steps);
      rng.setNewValues(rng_vals);
      walker_gpu.doStep(steps);

      auto M_cpu = walker_cpu.getM();
      auto M_gpu = walker_gpu.getM();
      for (int s = 0; s < 2; ++s)
        EXPECT_TRUE(dca::linalg::matrixop::areNear(M_cpu[s], M_gpu[s], 1e-7));

      // The final configuration is the same.
      const auto& config1 = walker_cpu.getWalkerConfiguration();
      const auto& config2 = walker_gpu.getWalkerConfiguration();
      ASSERT_EQ(config1.size(), config2.size());
      for (int i = 0; i < config1.size(); ++i)
        EXPECT_EQ(config1[i], config2[i]);
    }
  }
}
