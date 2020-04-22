// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Validation test of CT-INT against ED with a non density-density coupling.

#include <string>
#include <iostream>
#include <cmath>

#include "gtest/gtest.h"

#include "dca/math/statistical_testing/statistical_testing.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/math/random/random.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/models/analytic_hamiltonians/hund_lattice.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"

using Model = dca::phys::models::TightBindingModel<
    dca::phys::models::HundLattice<dca::phys::domains::no_symmetry<2>>>;
using RandomNumberGenerator = dca::math::random::StdRandomWrapper<std::mt19937_64>;

using ParametersType =
    dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType,
                                  dca::parallel::stdthread, dca::profiling::NullProfiler, Model,
                                  RandomNumberGenerator, dca::phys::solver::CT_INT>;

using DcaData = dca::phys::DcaData<ParametersType>;
using Solver = dca::phys::solver::CtintClusterSolver<dca::linalg::CPU, ParametersType, true>;

using SigmaCutDomain = dca::math::util::SigmaCutDomain<dca::math::util::details::Kdmn>;
using SigmaDomain = dca::math::util::SigmaDomain<dca::math::util::details::Kdmn>;
using CovarianceDomain = dca::math::util::CovarianceDomain<dca::math::util::details::Kdmn>;
using dca::math::util::cutFrequency;

dca::testing::DcaMpiTestEnvironment* dca_test_env;
const std::string test_directory =
    DCA_SOURCE_DIR "/test/integration/statistical_tests/bilayer_lattice/";
const int n_frequencies = 5;

TEST(CtintHundValidationTest, GreensFunction) {
  using namespace dca::testing;
  const std::string ed_data_name = test_directory + "/hund_data.ed.hdf5";

  const int id = dca_test_env->concurrency.id();
  const int number_of_samples = dca_test_env->concurrency.number_of_processors();

  if (id == dca_test_env->concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  ParametersType parameters(dca::util::GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(test_directory +
                                                           "hund_lattice_input.json");
  parameters.update_model();
  parameters.update_domains();

  parameters.set_measurements(parameters.get_measurements().back() * number_of_samples);

  DcaData data(parameters);
  data.initialize();

  // Do one QMC iteration
  Solver qmc_solver(parameters, data);
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
    test.printInfo("ctint_hund_testinfo.out", true);
    double p_value_default = 0.05;
    std::cout << "\n***\nThe p-value is " << p_value << "\n***\n";
    EXPECT_LT(p_value_default, p_value);
  }

  // Uncomment to write integrator output.
  // dca::phys::DcaLoopData<ParametersType> loop_data;
  // qmc_solver.finalize(loop_data);
  // if (id == 0) {
  //   std::cout << "\nProcessor " << id << " is writing data " << std::endl;
  //   dca::io::HDF5Writer writer;
  //   writer.open_file("ctint_square_results.hdf5");
  //   writer.open_group("functions");
  //   writer.execute(data.G_k_w);
  //   writer.close_group();
  //   writer.close_file();
  // }
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env = new dca::testing::DcaMpiTestEnvironment(argc, argv, "");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
