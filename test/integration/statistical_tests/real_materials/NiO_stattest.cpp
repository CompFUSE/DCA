// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
// 
// Statistical integration test for the CT-HYB solver using  a NiO lattice.
// This is a verification test with reference data taken from an exact diagonalization
// It can run with any number of MPI ranks.

#include "gtest/gtest.h"

#include "NiO_setup.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/models/material_hamiltonians/material_lattice.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

// needed only if solver output is written
//#include "dca/phys/dca_loop/dca_loop_data.hpp"

dca::testing::DcaMpiTestEnvironment* dca_test_env;

const std::string file_coarsegraining =
    DCA_SOURCE_DIR "/test/integration/cluster_solver/ss_ct_hyb/NiO_coarse_grained.hdf5";

using dca::func::dmn_0;
using dca::func::dmn_variadic;

TEST(NiO, ExactDiagonalization) {
  using namespace dca::testing;
  constexpr int dim = 3;

  using Concurrency = std::remove_pointer_t<decltype(dca_test_env)>::ConcurrencyType;
  using SigmaDomain = dca::math::util::SigmaDomain<dca::math::util::details::Kdmn<dim>>;
  using SigmaCutDomain = dca::math::util::SigmaCutDomain<dca::math::util::details::Kdmn<dim>>;
  using CovarianceDomain = dca::math::util::CovarianceDomain<dca::math::util::details::Kdmn<dim>>;
  
  const std::string ed_data_name = "NiO.data.ed.hdf5";
  const int id = dca_test_env->concurrency.id();
  const int number_of_samples = dca_test_env->concurrency.number_of_processors();

  if (id == dca_test_env->concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  dca::testing::TestParameters<dca::parallel::MPIConcurrency, dca::ClusterSolverId::SS_CT_HYB> parameters(dca::util::GitVersion::string(),
                                                          dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(test_directory + "input_NiO.json");
  parameters.set_t_ij_file_name(test_directory + "t_ij_NiO.txt");
  parameters.set_U_ij_file_name(test_directory + "U_ij_NiO_8_lit.txt");
  parameters.update_model();
  parameters.update_domains();

  parameters.set_measurements(parameters.get_measurements().back() * number_of_samples);

  DcaData<dca::ClusterSolverId::SS_CT_HYB, Concurrency> data(parameters);
  data.initializeH0_and_H_i();

  // Read and broadcast the rest of the initialization from full DCA results.
  if (id == 0) {
    dca::io::HDF5Reader reader;
    reader.open_file(file_coarsegraining);
    reader.open_group("functions");
    reader.execute(data.G0_k_w);
    reader.execute(data.G0_r_t);
    reader.execute(data.G_k_w);
    reader.close_group();
    reader.open_group("additional_functions");
    reader.execute(data.Sigma_cluster);
    reader.close_group();
    reader.close_file();
  }
  dca_test_env->concurrency.broadcast(data.G0_k_w);
  dca_test_env->concurrency.broadcast(data.G0_r_t);
  dca_test_env->concurrency.broadcast(data.G_k_w);
  dca_test_env->concurrency.broadcast(data.Sigma_cluster);

  //data.initialize();

  // Do one QMC iteration
  QuantumClusterSolver<dca::ClusterSolverId::SS_CT_HYB, Concurrency, CPU> qmc_solver(parameters, data, nullptr);
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
    reader.open_file(test_directory + ed_data_name);
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
    test.printInfo("NiO_testinfo.out", true);
    double p_value_default = 0.05;
    std::cout << "\n***\nThe p-value is " << p_value << "\n***\n";
    EXPECT_LT(p_value_default, p_value);
  }

  // Uncomment to write integrator output.
  // dca::phys::DcaLoopData<TestParameters<dca::ClusterSolverId::CT_AUX>> loop_data;
  // qmc_solver.finalize(loop_data);
  // if (id == 0) {
  //   std::cout << "\nProcessor " << id << " is writing data " << std::endl;
  //   dca::io::HDF5Writer writer;
  //   writer.open_file("ctint_bilayer_results.hdf5");
  //   writer.open_group("functions");
  //   writer.execute(data.G_k_w);
  //   writer.close_group();
  //   writer.close_file();
  // }
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  dca::parallel::MPIConcurrency concurrency(argc, argv);
  dca_test_env =
    new dca::testing::DcaMpiTestEnvironment(concurrency, dca::testing::test_directory + "input_NiO.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }
  result = RUN_ALL_TESTS();
  return result;
}
