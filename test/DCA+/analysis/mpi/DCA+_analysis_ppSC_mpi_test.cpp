//====================================================================
// Copyright 2015 ETH Zurich.
//
// Description
//
// Author: Urs Haehner (haehneru@itp.phys.ethz.ch), ETH Zurich
//====================================================================

#include <string>
#include <iostream>
#include <complex>

#include "gitVersion.hpp"
#include "modules.hpp"
#include "include_files.h"
#include "type_definitions.h"
#include "gtest/gtest.h"
#include "../../../common/minimalist_printer.hpp"
#include "../../../common/dca_mpi_test_environment.hpp"

dca_mpi_test_environment* dca_test_env;

TEST(DCA_analysis_ppSC, leading_eigenvalues) {
  using namespace DCA;

  using parameters_type = Parameters<dca_mpi_test_environment::concurrency_type,
                                     model, CT_AUX_CLUSTER_SOLVER>;
  using MOMS_type = DCA_data<parameters_type>;

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    std::cout << "Analysis starting.\n"
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
  MOMS.read(parameters.get_directory() + parameters.get_output_file_name());

  BSE_solver<parameters_type, MOMS_type> analysis_obj(parameters, MOMS);
  analysis_obj.calculate_susceptibilities_2();

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    const static int N_LAMBDAS = 10;
    typedef dmn_0<dmn<N_LAMBDAS, int> > lambda_dmn_type;

    FUNC_LIB::function<std::complex<double>, lambda_dmn_type>& leading_eigenvalues =
        analysis_obj.get_leading_eigenvalues();

    std::cout << "\nProcessor " << dca_test_env->concurrency.id()
              << " is checking data " << std::endl;

    ASSERT_NEAR(0.13607121357030313, std::real(leading_eigenvalues(0)), 1.e-14);
    ASSERT_NEAR(0.034431977492612967, std::real(leading_eigenvalues(1)), 1.e-14);
    ASSERT_NEAR(0.012266372581790334, std::real(leading_eigenvalues(2)), 1.e-14);
    ASSERT_NEAR(0.0054573275829339421, std::real(leading_eigenvalues(3)), 1.e-14);
    ASSERT_NEAR(0.0034328623407512224, std::real(leading_eigenvalues(4)), 1.e-14);
    ASSERT_NEAR(0.0020805637006257741, std::real(leading_eigenvalues(5)), 1.e-14);
    ASSERT_NEAR(0.0014371558441993106, std::real(leading_eigenvalues(6)), 1.e-14);
    ASSERT_NEAR(0.00090401084880077046, std::real(leading_eigenvalues(7)), 1.e-14);
    ASSERT_NEAR(0.0004745991433807814, std::real(leading_eigenvalues(8)), 1.e-14);
    ASSERT_NEAR(0.00041011935715527505, std::real(leading_eigenvalues(9)), 1.e-14);
  }

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.last())
    std::cout << "\nAnalysis ending.\n" << std::endl;
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env = new dca_mpi_test_environment(
      argc, argv, "input.U=4_d=0.900_Nc=4A_T=1.0.DCA+_analysis_ppSC.json");
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
