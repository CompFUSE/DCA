// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// No-change test for a concurrent (using MPI) and CUDA-accelerated DCA+ calculation using the
// CT-AUX cluster solver.
// The test simulates the tight-binding model on a 2D square lattice.

#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/config/cmake_options.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/filesystem.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/linalg/util/info_cuda.hpp"
#include "dca/linalg/util/util_cublas.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
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

TEST(dca_sp_DCAplus_mpi, Self_energy) {
  using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
  using DcaPointGroupType = dca::phys::domains::D4;
  using LatticeType = dca::phys::models::square_lattice<DcaPointGroupType>;
  using ModelType = dca::phys::models::TightBindingModel<LatticeType>;
  using Threading = dca::parallel::NoThreading;
  using ParametersType =
      dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType, Threading,
                                    dca::profiling::NullProfiler, ModelType, RngType,
                                    dca::phys::solver::CT_AUX>;
  using DcaDataType = dca::phys::DcaData<ParametersType>;
  using ClusterSolverType =
      dca::phys::solver::CtauxClusterSolver<dca::linalg::CPU, ParametersType, DcaDataType>;
  using DcaLoopType =
      dca::phys::DcaLoop<ParametersType, DcaDataType, ClusterSolverType, dca::DistType::NONE>;

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    // Copy initial state from an aborted run.
    filesystem::copy_file(
        DCA_SOURCE_DIR "/test/system-level/dca/data.dca_sp_DCA+_mpi_test.hdf5.tmp",
        "./data.dca_sp_DCA+_mpi_test.hdf5.tmp", filesystem::copy_options::overwrite_existing);

    dca::util::GitVersion::print();
    dca::util::Modules::print();
    dca::config::CMakeOptions::print();

    dca::linalg::util::printInfoDevices();

    std::cout
        << "\n"
        << "********************************************************************************\n"
        << "**********                     DCA(+) Calculation                     **********\n"
        << "********************************************************************************\n"
        << "\n"
        << "Start time : " << dca::util::print_time() << "\n"
        << "\n"
        << "MPI-world set up: " << dca_test_env->concurrency.number_of_processors() << " processes."
        << "\n"
        << std::endl;
  }

  dca::linalg::util::initializeMagma();

  ParametersType parameters(dca::util::GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(dca_test_env->input_file_name);
  parameters.update_model();
  parameters.update_domains();

  DcaDataType dca_data(parameters);
  dca_data.initialize();

  DcaLoopType dca_loop(parameters, dca_data, dca_test_env->concurrency);
  dca_loop.initialize();
  dca_loop.execute();
  dca_loop.finalize();

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    std::cout << "\nProcessor " << dca_test_env->concurrency.id() << " is checking data "
              << std::endl;

    // Read self-energy from check_data file.
    DcaDataType::SpGreensFunction Sigma_check("Self_Energy");
    dca::io::HDF5Reader reader;
    reader.open_file(DCA_SOURCE_DIR "/test/system-level/dca/check_data.dca_sp_DCA+_mpi_test.hdf5");
    reader.open_group("functions");
    reader.execute(Sigma_check);
    reader.close_file();

    // Compare the computed self-energy with the expected result.
    auto diff = dca::func::util::difference(Sigma_check, dca_data.Sigma);

    EXPECT_NEAR(0, diff.l2, 1.e-12);

    std::cout << "\nProcessor " << dca_test_env->concurrency.id() << " is writing data." << std::endl;
    dca_loop.write();

    std::cout << "\nFinish time: " << dca::util::print_time() << "\n" << std::endl;
  }
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca::parallel::MPIConcurrency concurrency(argc, argv);
  dca_test_env = new dca::testing::DcaMpiTestEnvironment(
      concurrency, DCA_SOURCE_DIR "/test/system-level/dca/input.dca_sp_DCA+_mpi_test.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
