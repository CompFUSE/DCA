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
#include "gtest/gtest.h"
#include "minimalist_printer.hpp"
#include "dca_mpi_test_environment.hpp"

dca_mpi_test_environment* dca_test_env;

TEST(analysis_DCAplus_mpi, leading_eigenvalues) {
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
    std::cout << "\nProcessor " << dca_test_env->concurrency.id()
              << " is checking data " << std::endl;
    
    const static int N_LAMBDAS = 10;
    typedef dmn_0<dmn<N_LAMBDAS, int> > lambda_dmn_type;

    FUNC_LIB::function<std::complex<double>, lambda_dmn_type>& leading_eigenvalues =
        analysis_obj.get_leading_eigenvalues();

    // Read eigenvalues from check_data file.
    FUNC_LIB::function<std::complex<double>, lambda_dmn_type> leading_eigenvalues_check("leading-eigenvalues");
    IO::reader<IO::HDF5> reader;
    reader.open_file("check_data.analysis_DCA+_mpi_test.hdf5");
    reader.open_group("analysis-functions");
    reader.execute(leading_eigenvalues_check);
    reader.close_file();

    // Compare the computed eigenvalues with the expected result.
    for (int i = 0; i < lambda_dmn_type::dmn_size(); ++i) {
      EXPECT_NEAR(leading_eigenvalues_check(i).real(),
                  leading_eigenvalues(i).real(), 1.e-14);
      EXPECT_NEAR(leading_eigenvalues_check(i).imag(),
                  leading_eigenvalues(i).imag(), 1.e-14);
    }
  }
  
  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.last()) {
    std::cout << "\nProcessor " << dca_test_env->concurrency.id() << " is writing data "
              << std::endl;
    analysis_obj.write();
    
    std::cout << "\nAnalysis ending.\n" << std::endl;
  }
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env = new dca_mpi_test_environment(
      argc, argv, "input.analysis_DCA+_mpi_test.json");
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
