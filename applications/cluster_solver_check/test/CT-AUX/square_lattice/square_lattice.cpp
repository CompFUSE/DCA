//====================================================================
// Copyright 2016 ETH Zurich.
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

TEST(DCA_loop_sp, Self_Energy) {
  using namespace DCA;

  using parameters_type = Parameters<dca_mpi_test_environment::concurrency_type,
                                     model, CT_AUX_CLUSTER_SOLVER>;
  using MOMS_type = DCA_data<parameters_type>;
  // using ED_solver_type =
  //     cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, LIN_ALG::CPU,
  //                    parameters_type, MOMS_type>;
  using quantum_cluster_solver_type =
    cluster_solver<CT_AUX_CLUSTER_SOLVER, LIN_ALG::CPU, parameters_type, MOMS_type>;
  using QMC_solver_type = quantum_cluster_solver_type;

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

  DCA_calculation_data DCA_info_struct;

  MOMS_type MOMS_imag(parameters);
  MOMS_imag.initialize();

  // Read ED data
  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    IO::reader<IO::HDF5> reader;
    reader.open_file("data.ED.hdf5");
    reader.open_group("functions");
    reader.execute(MOMS_imag.Sigma);
    reader.execute(MOMS_imag.G0_k_w_cluster_excluded);
    reader.execute(MOMS_imag.G0_r_t_cluster_excluded);
    reader.close_file();
  }

  dca_test_env->concurrency.broadcast(MOMS_imag.Sigma);
  dca_test_env->concurrency.broadcast(MOMS_imag.G0_k_w_cluster_excluded);
  dca_test_env->concurrency.broadcast(MOMS_imag.G0_r_t_cluster_excluded);
  
  FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> > Sigma_ED(MOMS_imag.Sigma);
  
  // QMC solver
  // The QMC solver uses the free Greens function G0 computed by the ED solver.
  // It is passed via the MOMS_imag object.
  QMC_solver_type QMC_obj(parameters, MOMS_imag);
  QMC_obj.initialize(1);
  QMC_obj.integrate();
  QMC_obj.finalize(DCA_info_struct);

  FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> > Sigma_QMC(MOMS_imag.Sigma);

  // Check errors
  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); ++k_ind) {
      
      double ED_min_QMC_inf = 0.;
      double ED_inf = 0.;
      for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind) {
        ED_min_QMC_inf = std::max(ED_min_QMC_inf, std::abs(Sigma_ED(0, 0, 0, 0, k_ind, w_ind)
                                                           - Sigma_QMC(0, 0, 0, 0, k_ind, w_ind)));
        ED_inf = std::max(ED_inf, std::abs(Sigma_ED(0, 0, 0, 0, k_ind, w_ind)));
      }
      EXPECT_LT(ED_min_QMC_inf, 5.e-2);
      EXPECT_LT(ED_min_QMC_inf/ED_inf, 5.e-2);
    }
  }

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    std::cout << "\nCluster-solver-check ending.\n" << std::endl;
  }
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env = new dca_mpi_test_environment(
      argc, argv, "input.square_lattice.json");
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
