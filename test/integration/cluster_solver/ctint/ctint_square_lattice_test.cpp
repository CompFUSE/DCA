// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// No-change test for CT-INT.
// Square lattice with single band and double occupancy repulsion U.

#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

dca::testing::DcaMpiTestEnvironment* dca_test_env;
const std::string input_dir = DCA_SOURCE_DIR "/test/integration/cluster_solver/ctint/";

constexpr bool update_baseline = false;

TEST(squareLattice_Nc4_nn, Self_Energy) {
  using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
  using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
  using Model = dca::phys::models::TightBindingModel<Lattice>;
  using Threading = dca::parallel::NoThreading;
  using Parameters =
      dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType, Threading,
                                    dca::profiling::NullProfiler, Model, RngType, dca::phys::solver::CT_INT>;
  using Data = dca::phys::DcaData<Parameters>;
  using QmcSolverType = dca::phys::solver::CtintClusterSolver<dca::linalg::CPU, Parameters, true>;

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  Parameters parameters(dca::util::GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(dca_test_env->input_file_name);
  parameters.update_model();
  parameters.update_domains();

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

  // Do one integration step.
  QmcSolverType qmc_solver(parameters, data);
  qmc_solver.initialize(0);
  qmc_solver.integrate();
  qmc_solver.finalize();

  EXPECT_NEAR(1.0, qmc_solver.computeDensity(), 1e-5);

  if (not update_baseline) {
    // Read and confront with previous run
    if (dca_test_env->concurrency.id() == 0) {
      Data::SpGreensFunction G_k_w_check(data.G_k_w.get_name());
      dca::io::HDF5Reader reader;
      reader.open_file(input_dir + "square_lattice_baseline.hdf5");
      reader.open_group("functions");
      reader.execute(G_k_w_check);
      reader.close_group(), reader.close_file();

      auto diff = dca::func::util::difference(G_k_w_check, data.G_k_w);
      EXPECT_GE(1e-6, diff.l2);
    }
  }
  else {
    //  Write results
    if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
      dca::io::HDF5Writer writer;
      writer.open_file(input_dir + "square_lattice_baseline.hdf5");
      writer.open_group("functions");
      writer.execute(data.G_k_w);
      writer.close_group(), writer.close_file();
    }
  }
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env =
      new dca::testing::DcaMpiTestEnvironment(argc, argv, input_dir + "square_lattice_input.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();
  return result;
}
