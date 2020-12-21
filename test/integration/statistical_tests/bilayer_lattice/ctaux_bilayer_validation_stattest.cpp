// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Verification test of CT-INT against a reference run

#include <string>
#include <iostream>
#include <cmath>

#include "gtest/gtest.h"

#include "dca/math/statistical_testing/statistical_testing.hpp"
#include "test/integration/statistical_tests/bilayer_lattice/bilayer_lattice_setup.hpp"

dca::testing::DcaMpiTestEnvironment* dca_test_env;

TEST(CtauxBilayerValidationTest, GreensFunction) {
  //  dca::linalg::util::initializeMagma();

  using namespace dca::testing;
  const std::string ed_data_name = dca::testing::test_directory + "/data.ed.hdf5";

  const int id = dca_test_env->concurrency.id();
  const int number_of_samples = dca_test_env->concurrency.number_of_processors();

  if (id == dca_test_env->concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  ParametersType<CT_AUX> parameters(dca::util::GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(dca_test_env->input_file_name);
  parameters.update_model();
  parameters.update_domains();

  parameters.set_measurements(parameters.get_measurements().back() * number_of_samples);

  DcaData<CT_AUX> data(parameters);
  data.initialize();

  // Do one QMC iteration
  QuantumClusterSolver<CT_AUX, CPU> qmc_solver(parameters, data);
  qmc_solver.initialize(0);
  qmc_solver.integrate();

  // stores quantities from integration.
  using dca::func::function;
  function<double, SigmaCutDomain> G_k_w_sample =
      cutFrequency(qmc_solver.local_G_k_w(), n_frequencies);
  G_k_w_sample.set_name("G_k_w");

  // read the expected result
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

  // compute covariance and average ctin result.
  function<double, CovarianceDomain> G_k_w_covariance("G_k_w_covariance");
  dca_test_env->concurrency.computeCovarianceAndAvg(G_k_w_sample, G_k_w_covariance);

  //   compute p-value
  if (id == dca_test_env->concurrency.first()) {
    // read the stored reference data
    dca::math::StatisticalTesting test(G_k_w_sample, G_k_w_expected, G_k_w_covariance, 1);
    double p_value = test.computePValue(false, number_of_samples);
    test.printInfo("ctaux_bilayer_testinfo.out", true);
    double p_value_default = 0.05;
    std::cout << "\n***\nThe p-value is " << p_value << "\n***\n";
    EXPECT_LT(p_value_default, p_value);
  }

  // Uncomment to write integrator output.
  dca::phys::DcaLoopData<ParametersType<CT_AUX>> loop_data;
  qmc_solver.finalize(loop_data);
  if (id == 0) {
    std::cout << "\nProcessor " << id << " is writing data " << std::endl;
    dca::io::HDF5Writer writer;
    writer.open_file("ctint_bilayer_results.hdf5");
    writer.open_group("functions");
    writer.execute(data.G_k_w);
    writer.close_group();
    writer.close_file();
  }
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca::parallel::MPIConcurrency concurrency(argc, argv);
  dca_test_env = new dca::testing::DcaMpiTestEnvironment(
      concurrency, dca::testing::test_directory + "bilayer_lattice_input.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
