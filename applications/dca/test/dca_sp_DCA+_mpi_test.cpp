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
#include "minimalist_printer.hpp"
#include "dca_mpi_test_environment.hpp"

dca_mpi_test_environment* dca_test_env;

TEST(dca_sp_DCAplus_mpi, Self_energy) {
  using namespace DCA;

  using parameters_type = Parameters<dca_mpi_test_environment::concurrency_type,
                                     model, CT_AUX_CLUSTER_SOLVER>;
  using MOMS_type = DCA_data<parameters_type>;
  using Monte_Carlo_Integrator_type =
    cluster_solver<CT_AUX_CLUSTER_SOLVER, LIN_ALG::CPU, parameters_type, MOMS_type>;
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

    // Read self-energy from check_data file.
    FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> > Sigma_check("Self_Energy");
    IO::reader<IO::HDF5> reader;
    reader.open_file("check_data.dca_sp_DCA+_mpi_test.hdf5");
    reader.open_group("functions");
    reader.execute(Sigma_check);
    reader.close_file();

    // Compare the computed self-energy with the expected result.
    for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind) {
      for (int k_ind = 0; k_ind < k_DCA::dmn_size(); ++k_ind) {
        for (int nu_ind_2 = 0; nu_ind_2 < nu::dmn_size(); ++nu_ind_2) {
          for (int nu_ind_1 = 0; nu_ind_1 < nu::dmn_size(); ++nu_ind_1) {
            EXPECT_NEAR(Sigma_check(nu_ind_1, nu_ind_2, k_ind, w_ind).real(),
                        MOMS.Sigma(nu_ind_1, nu_ind_2, k_ind, w_ind).real(), 1.e-12);
            EXPECT_NEAR(Sigma_check(nu_ind_1, nu_ind_2, k_ind, w_ind).imag(),
                        MOMS.Sigma(nu_ind_1, nu_ind_2, k_ind, w_ind).imag(), 1.e-12);
          }
        }
      }
    }
  }
  
  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.last()) {
    std::cout << "\nProcessor " << dca_test_env->concurrency.id()
              << " is writing data " << std::endl;
    dca_object.write();
    
    std::cout << "\nDCA main ending.\n" << std::endl;
  }
 }
  
int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env = new dca_mpi_test_environment(
      argc, argv, "input.dca_sp_DCA+_mpi_test.json");
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
