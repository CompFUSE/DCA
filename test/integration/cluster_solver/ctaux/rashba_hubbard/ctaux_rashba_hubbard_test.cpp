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
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/models/analytic_hamiltonians/rashba_hubbard.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

dca::testing::DcaMpiTestEnvironment* dca_test_env;
const std::string input_dir =
    DCA_SOURCE_DIR "/test/integration/cluster_solver/ctaux/rashba_hubbard/";

constexpr bool update_baseline = false;

TEST(CtauxSolverTest, RashaHubbardModel) {
  using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
  using Lattice = dca::phys::models::RashbaHubbard<dca::phys::domains::no_symmetry<2>>;
  using Model = dca::phys::models::TightBindingModel<Lattice>;
  using Threading = dca::parallel::stdthread;
  using Parameters =
      dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType, Threading,
                                    dca::profiling::NullProfiler, Model, RngType, dca::phys::solver::CT_AUX>;
  using Data = dca::phys::DcaData<Parameters>;
  using QmcBaseSolver = dca::phys::solver::CtauxClusterSolver<dca::linalg::GPU, Parameters, Data>;
  using QmcSolver = dca::phys::solver::StdThreadQmciClusterSolver<QmcBaseSolver>;

  dca::linalg::util::initializeMagma();

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  Parameters parameters(dca::util::GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir + "rashba_hubbard_input.json");
  parameters.update_model();
  parameters.update_domains();

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

  // Do one integration step.
  QmcSolver qmc_solver(parameters, data);
  qmc_solver.initialize(0);
  qmc_solver.integrate();

  dca::phys::DcaLoopData<Parameters> dca_loop_data;
  qmc_solver.finalize(dca_loop_data);

  const std::string baseline = input_dir + "rashba_hubbard_baseline.hdf5";
  if (!update_baseline) {
    // Read and confront with previous run
    if (dca_test_env->concurrency.id() == 0) {
      Data::SpGreensFunction G_k_w_check(data.G_k_w.get_name());
      Data::TpGreensFunction G4_check(data.get_G4()[0].get_name());
      dca::io::HDF5Reader reader;
      reader.open_file(baseline);
      reader.open_group("functions");
      reader.execute(G_k_w_check);
      reader.execute(G4_check);
      reader.close_group();
      reader.close_file();

      auto diff = dca::func::util::difference(G_k_w_check, data.G_k_w);
      auto diff_g4 = dca::func::util::difference(G4_check, data.get_G4()[0]);
      EXPECT_GE(1e-6, diff.l2);
      EXPECT_GE(1e-6, diff_g4.l2);
    }
  }
  else {
    //  Write results
    if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
      dca::io::HDF5Writer writer;
      writer.open_file(baseline);
      writer.open_group("functions");
      writer.execute(data.G_k_w);
      writer.execute(data.get_G4()[0]);
      writer.close_group(), writer.close_file();
    }
  }
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  dca::parallel::MPIConcurrency concurrency(argc, argv);
  dca_test_env =
      new dca::testing::DcaMpiTestEnvironment(concurrency, "");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();
  return result;
}
