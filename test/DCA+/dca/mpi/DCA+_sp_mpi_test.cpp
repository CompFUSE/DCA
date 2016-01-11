//====================================================================
// Copyright 2015 ETH Zurich.
//
// Description
//
// Author: Urs Haehner (haehneru@itp.phys.ethz.ch), ETH Zurich
//====================================================================

#include <string>
#include <iostream>
#include <cmath>

#include "gitVersion.hpp"
#include "modules.hpp"
#include "include_files.h"
#include "type_definitions.h"
#include "gtest/gtest.h"
#include "../../../common/minimalist_printer.hpp"
#include "../../../common/dca_mpi_test_environment.hpp"

dca_mpi_test_environment* dca_test_env;

TEST(DCA_loop_sp, Self_Energy) {
  using namespace DCA;

  static const LIN_ALG::device_type LIN_ALG_DEVICE = LIN_ALG::CPU;
  static const CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME = CT_AUX_CLUSTER_SOLVER;

  using parameters_type = Parameters<dca_mpi_test_environment::concurrency_type,
                                     model, CLUSTER_SOLVER_NAME>;
  using MOMS_type = DCA_data<parameters_type>;
  using quantum_cluster_solver_type =
    cluster_solver<CLUSTER_SOLVER_NAME, LIN_ALG_DEVICE, parameters_type, MOMS_type>;
  using Monte_Carlo_Integrator_type = quantum_cluster_solver_type;
  using DCA_calculation_type =
    DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>;

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    std::cout << "\nDCA main starting.\n"
              << "MPI-world set up: "
              << dca_test_env->concurrency.number_of_processors()
              << " processes.\n" << std::endl;

    GitVersion::print();
    Modules::print();
  }

  parameters_type parameters(GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast(dca_test_env->input_file);
  parameters.update_model();
  parameters.update_domains();

  MOMS_type MOMS(parameters);
  MOMS.initialize();

  DCA_calculation_type dca_object(parameters, MOMS, dca_test_env->concurrency);
  dca_object.initialize();
  dca_object.execute();
  dca_object.finalize();

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    std::cout << "\nProcessor " << dca_test_env->concurrency.id()
              << " is checking data " << std::endl;

    // Final QMC self energy
    ASSERT_NEAR(-0.19763191899502217,
                     std::real(MOMS.Sigma(0, 0, 0, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.39543084934934963,
                     std::imag(MOMS.Sigma(0, 0, 0, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.1838397635149297,
                     std::real(MOMS.Sigma(0, 0, 1, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.39980886467162741,
                     std::imag(MOMS.Sigma(0, 0, 1, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.1838397635149297,
                     std::real(MOMS.Sigma(0, 0, 2, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.39980886467162741,
                     std::imag(MOMS.Sigma(0, 0, 2, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.16942929468603241,
                     std::real(MOMS.Sigma(0, 0, 3, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.39971634937441847,
                     std::imag(MOMS.Sigma(0, 0, 3, w::dmn_size()/2)), 1.e-14);

    // Final coarsegrained self energy
    ASSERT_NEAR(-0.18606344993297538,
                     std::real(MOMS.Sigma_cluster(0, 0, 0, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.39857826244853495,
                     std::imag(MOMS.Sigma_cluster(0, 0, 0, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.18369581993172768,
                     std::real(MOMS.Sigma_cluster(0, 0, 1, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.39875266456335279,
                     std::imag(MOMS.Sigma_cluster(0, 0, 1, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.18369581993172768,
                     std::real(MOMS.Sigma_cluster(0, 0, 2, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.39875266456335279,
                     std::imag(MOMS.Sigma_cluster(0, 0, 2, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.1813281899304795,
                     std::real(MOMS.Sigma_cluster(0, 0, 3, w::dmn_size()/2)), 1.e-14);
    ASSERT_NEAR(-0.39892706667816974,
                     std::imag(MOMS.Sigma_cluster(0, 0, 3, w::dmn_size()/2)), 1.e-14);
  }

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.last())
    std::cout << "\nDCA main ending.\n" << std::endl;
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env = new dca_mpi_test_environment(
      argc, argv, "input.U=4_d=0.900_Nc=4A_T=1.0.DCA+_sp.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
