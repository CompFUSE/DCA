// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// No-change test for CT-INT posix wrapper.

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
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

constexpr bool update_results = false;

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/cluster_solver/stdthread_qmci/";

using Concurrency = dca::parallel::NoConcurrency;
using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Threading = dca::parallel::stdthread;
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                                 Model, RngType, dca::phys::solver::CT_INT>;
using Data = dca::phys::DcaData<Parameters>;
using BaseSolver = dca::phys::solver::CtintClusterSolver<dca::linalg::CPU, Parameters>;
using QmcSolver = dca::phys::solver::StdThreadQmciClusterSolver<BaseSolver>;

void performTest(const std::string& input, const std::string& baseline) {
  static bool update_model = true;

  Concurrency concurrency(0, nullptr);

  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir + input);
  if (update_model) {
    parameters.update_model();
    parameters.update_domains();
  }
  update_model = false;

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

  // Do one integration step.
  QmcSolver qmc_solver(parameters, data);
  qmc_solver.initialize(0);
  qmc_solver.integrate();

  dca::phys::DcaLoopData<Parameters> loop_data;
  qmc_solver.finalize(loop_data);

  EXPECT_NEAR(1., qmc_solver.computeDensity(), 1e-2);

  if (not update_results) {
    // Read and confront with previous run.
    if (concurrency.id() == 0) {
      Data::SpGreensFunction G_k_w_check(data.G_k_w.get_name());
      dca::io::HDF5Reader reader;
      reader.open_file(input_dir + baseline);
      reader.open_group("functions");
      reader.execute(G_k_w_check);
      reader.close_group(), reader.close_file();

      const auto err_g = dca::func::util::difference(G_k_w_check, data.G_k_w);

      EXPECT_GE(5e-7, err_g.l_inf);
    }
  }
  else {
    //  Write results
    if (concurrency.id() == concurrency.first()) {
      dca::io::HDF5Writer writer;
      writer.open_file(input_dir + baseline);
      writer.open_group("functions");
      writer.execute(data.G_k_w);
      writer.close_group(), writer.close_file();
    }
  }
}

TEST(PosixCtintClusterSolverTest, NonShared) {
  performTest("stdthread_ctint_test_nonshared_input.json",
              "stdthread_ctint_test_nonshared_baseline.hdf5");
}

TEST(PosixCtintClusterSolverTest, Shared) {
  performTest("stdthread_ctint_test_shared_input.json",
              "stdthread_ctint_test_shared_baseline.hdf5");
}
