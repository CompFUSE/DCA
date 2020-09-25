// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// No-change test for a full DCA loop with a threaded ct-int solver.

#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/dca_loop/dca_loop.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

constexpr bool update_baseline = false;

dca::testing::DcaMpiTestEnvironment* dca_test_env;

TEST(dca_sp_DCAplus_mpi, Self_energy) {
  using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
  using DcaPointGroupType = dca::phys::domains::D4;
  using LatticeType = dca::phys::models::square_lattice<DcaPointGroupType>;
  using ModelType = dca::phys::models::TightBindingModel<LatticeType>;
  using Threading = dca::parallel::stdthread;
  using ParametersType =
      dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType, Threading,
                                    dca::profiling::NullProfiler, ModelType, RngType,
                                    dca::phys::solver::CT_INT>;
  using DcaDataType = dca::phys::DcaData<ParametersType>;
  using ClusterSolverType = dca::phys::solver::StdThreadQmciClusterSolver<
      dca::phys::solver::CtintClusterSolver<dca::linalg::CPU, ParametersType>>;
  using DcaLoopType = dca::phys::DcaLoop<ParametersType, DcaDataType, ClusterSolverType>;

  auto& concurrency = dca_test_env->concurrency;

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\nDCA main starting.\n"
              << "MPI-world set up: " << concurrency.number_of_processors() << " processes.\n"
              << std::endl;

    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  ParametersType parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(dca_test_env->input_file_name);
  parameters.update_model();
  parameters.update_domains();

  DcaDataType dca_data(parameters);
  dca_data.initialize();

  DcaLoopType dca_loop(parameters, dca_data, concurrency);
  dca_loop.initialize();
  dca_loop.execute();
  dca_loop.finalize();

  if (concurrency.id() == concurrency.first()) {
    const std::string baseline_name =
        DCA_SOURCE_DIR "/test/system-level/dca/baseline_ctint_sp_DCA_mpi_thread_test.hdf5";

    if (!update_baseline) {
      std::cout << "\nProcessor " << concurrency.id() << " is checking data " << std::endl;

      // Read self-energy from check_data file.
      decltype(dca_data.Sigma) Sigma_check(dca_data.Sigma.get_name());
      dca::io::HDF5Reader reader;
      reader.open_file(baseline_name);
      reader.open_group("functions");
      reader.execute(Sigma_check);
      reader.close_group();
      reader.close_file();

      // Compare the computed self-energy with the expected result.
      const auto err = dca::func::util::difference(Sigma_check, dca_data.Sigma);
      EXPECT_LE(err.l_inf, 1e-10);
    }
    else {
      //  Write results
      dca::io::HDF5Writer writer;
      writer.open_file(baseline_name);
      writer.open_group("functions");
      writer.execute(dca_data.Sigma);
      writer.close_group();
      writer.close_file();
    }
  }

  if (concurrency.id() == concurrency.last()) {
    std::cout << "\nProcessor " << concurrency.id() << " is writing data " << std::endl;
    dca_loop.write();

    std::cout << "\nDCA main ending.\n" << std::endl;
  }
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca::parallel::MPIConcurrency concurrency(argc, argv);
  dca_test_env = new dca::testing::DcaMpiTestEnvironment(
      concurrency, DCA_SOURCE_DIR "/test/system-level/dca/input.ctint_sp_DCA_mpi_thread_test.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
