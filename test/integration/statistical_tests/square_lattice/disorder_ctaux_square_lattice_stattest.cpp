// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Anirudha Mirmira (anirudha.mirmira@gmail.com)
//
// Disorder stat test, Step 1 (clean limit): drives the CT-AUX r-space disorder pipeline with a
// zero potential (V(R) = 0) and checks that the QMC-sampled local_G_k_w statistically agrees
// with the clean ED reference. With V = 0 the per-configuration r-space Dyson, the translational
// average, and the FT r->k must collapse to the clean answer, so any corruption introduced by the
// disorder machinery (driven here by genuine M sampled at U = 4) shows up as a failed p-value.
// This is the clean-limit regression of the full disorder code path.
//
// Only built when DCA_WITH_DISORDER is ON (see CMakeLists.txt), which is CPU-only, so there is no
// GPU branch here; the clean path is already covered by ctaux_square_lattice_validation_stattest.

#include <string>
#include <iostream>
#include <cmath>

#include "gtest/gtest.h"

#include "dca/math/statistical_testing/statistical_testing.hpp"
#include "dca/config/profiler.hpp"

// The shared setup header uses a global Scalar and includes ctint_cluster_solver.hpp, which needs
// dca::config::McOptions. Provide both before pulling in the setup header (same pattern as the
// rashba stat test). The square lattice Green's function is real, so Scalar is double.
using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "test/integration/statistical_tests/square_lattice/square_lattice_setup.hpp"

dca::testing::DcaMpiTestEnvironment* dca_test_env;

TEST(DisorderCtauxSquareLatticeCleanLimitTest, GreensFunction) {
  using namespace dca::testing;
  const std::string ed_data_name = "data.ed.hdf5";

  const int id = dca_test_env->concurrency.id();
  const int number_of_samples = dca_test_env->concurrency.number_of_processors();

  if (id == dca_test_env->concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  ParametersType<ClusterSolverId::CT_AUX> parameters(dca::util::GitVersion::string(),
                                                     dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(dca_test_env->input_file_name);
  parameters.update_model();
  parameters.update_domains();

  parameters.set_measurements(parameters.get_measurements().back() * number_of_samples / 50);

  DcaData<ClusterSolverId::CT_AUX> data(parameters);
  data.initialize();

  // The stat test drives the solver directly (initialize -> integrate -> local_G_k_w) and never
  // calls DcaLoop::workTheClusters, which is what normally populates the disordered bare
  // propagator data.disordered_G0_r_r_{t,w}_cl_exl that the walker and the post-integration
  // r-space Dyson read. Populate it here with a zero potential (V(R) = 0); the default-constructed
  // DisorderConfiguration is identically zero.
  typename DcaData<ClusterSolverId::CT_AUX>::DisorderConfiguration zero_config;
  data.makeDisorderedG0(zero_config);

  // Do one QMC iteration.
  QuantumClusterSolver<ClusterSolverId::CT_AUX> qmc_solver(parameters, data, nullptr);
  qmc_solver.initialize(0);
  qmc_solver.integrate();

  using dca::func::function;
  function<double, SigmaCutDomain> G_k_w_sample =
      cutFrequency(qmc_solver.local_G_k_w(), n_frequencies);
  G_k_w_sample.set_name("G_k_w");

  // Catastrophic-failure guard: the clean-limit disorder pipeline must produce finite numbers.
  for (int i = 0; i < G_k_w_sample.size(); ++i)
    ASSERT_TRUE(std::isfinite(G_k_w_sample(i)))
        << "Non-finite value in disorder clean-limit local_G_k_w at linear index " << i;

  // Read the expected (clean) ED result.
  function<double, SigmaCutDomain> G_k_w_expected;
  if (id == dca_test_env->concurrency.first()) {
    function<std::complex<double>, SigmaDomain> G_k_w_full("cluster_greens_function_G_k_w");
    dca::io::HDF5Reader reader;
    reader.open_file(ed_data_name);
    reader.open_group("functions");
    reader.execute(G_k_w_full);
    reader.close_group();
    reader.close_file();
    G_k_w_expected = cutFrequency(G_k_w_full, n_frequencies);
  }
  dca_test_env->concurrency.broadcast(G_k_w_expected);

  // Compute covariance and average the QMC result across ranks.
  function<double, CovarianceDomain> G_k_w_covariance("G_k_w_covariance");
  dca_test_env->concurrency.computeCovarianceAndAvg(G_k_w_sample, G_k_w_covariance);

  // Compute the p-value against the clean ED reference.
  if (id == dca_test_env->concurrency.first()) {
    dca::math::StatisticalTesting test(G_k_w_sample, G_k_w_expected, G_k_w_covariance, 1);
    double p_value = test.computePValue(false, number_of_samples);
    test.printInfo("disorder_ctaux_square_testinfo.out", true);
    double p_value_default = 0.05;
    std::cout << "\n***\nThe p-value is " << p_value << "\n***\n";
    EXPECT_LT(p_value_default, p_value);
  }
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca::parallel::MPIConcurrency concurrency(argc, argv);
  dca_test_env = new dca::testing::DcaMpiTestEnvironment(
      concurrency, dca::testing::test_directory + "square_lattice_input.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
