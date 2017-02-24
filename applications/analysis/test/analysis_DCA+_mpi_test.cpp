// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// No-change test for a concurrent (using MPI) DCA+ analysis calculation.
// It runs a simulation of a tight-binding model on 2D square lattice.

#define DCA_WITH_REDUCED_VERTEX_FUNCTION

#include <string>
#include <iostream>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/parallel/pthreading/pthreading.hpp"
#include "dca/phys/dca_analysis/bse_solver/bse_solver.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

dca::testing::DcaMpiTestEnvironment* dca_test_env;

TEST(analysis_DCAplus_mpi, leading_eigenvalues) {
  using DcaPointGroupType = dca::phys::domains::D4;
  using LatticeType = dca::phys::models::square_lattice<DcaPointGroupType>;
  using ModelType = dca::phys::models::TightBindingModel<LatticeType>;
  using Threading = dca::parallel::Pthreading;
  using ParametersType =
      dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType, Threading,
                                    dca::profiling::NullProfiler, ModelType, void /*RngType*/,
                                    dca::phys::solver::CT_AUX>;
  using DcaDataType = dca::phys::DcaData<ParametersType>;

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    std::cout << "Analysis starting.\n"
              << "MPI-world set up: " << dca_test_env->concurrency.number_of_processors()
              << " processes.\n"
              << std::endl;

    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  ParametersType parameters(dca::util::GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(dca_test_env->input_file_name);
  parameters.update_model();
  parameters.update_domains();

  DcaDataType dca_data(parameters);
  dca_data.initialize();
  dca_data.read(parameters.get_directory() + parameters.get_filename_dca());

  dca::phys::analysis::BseSolver<ParametersType, DcaDataType> analysis_obj(parameters, dca_data);
  analysis_obj.calculateSusceptibilities();

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    std::cout << "\nProcessor " << dca_test_env->concurrency.id() << " is checking data "
              << std::endl;

    const static int num_lambdas = 10;
    using lambda_dmn_type = dca::func::dmn_0<dca::func::dmn<num_lambdas, int>>;

    dca::func::function<std::complex<double>, lambda_dmn_type>& leading_eigenvalues =
        analysis_obj.get_leading_eigenvalues();

    // Read eigenvalues from check_data file.
    dca::func::function<std::complex<double>, lambda_dmn_type> leading_eigenvalues_check(
        "leading-eigenvalues");
    dca::io::HDF5Reader reader;
    reader.open_file(DCA_SOURCE_DIR
                     "/applications/analysis/test/check_data.analysis_DCA+_mpi_test.hdf5");
    reader.open_group("analysis-functions");
    reader.execute(leading_eigenvalues_check);
    reader.close_file();

    // Compare the computed eigenvalues with the expected result.
    for (int i = 0; i < lambda_dmn_type::dmn_size(); ++i) {
      EXPECT_NEAR(leading_eigenvalues_check(i).real(), leading_eigenvalues(i).real(), 1.e-14);
      EXPECT_NEAR(leading_eigenvalues_check(i).imag(), leading_eigenvalues(i).imag(), 1.e-14);
    }
  }

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.last()) {
    std::cout << "\nProcessor " << dca_test_env->concurrency.id() << " is writing data " << std::endl;
    analysis_obj.write();

    std::cout << "\nAnalysis ending.\n" << std::endl;
  }
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env = new dca::testing::DcaMpiTestEnvironment(
      argc, argv, DCA_SOURCE_DIR "/applications/analysis/test/input.analysis_DCA+_mpi_test.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
