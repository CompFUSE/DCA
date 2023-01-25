// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// No-change test for CT-INT posix wrapper.

#include "dca/io/writer.hpp"

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
#include "dca/config/profiler.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

[[maybe_unused]] constexpr bool update_results = false;

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/cluster_solver/stdthread_qmci/";

using TestConcurrency = dca::parallel::MPIConcurrency;
using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using StdThreading = dca::parallel::stdthread;
using Parameters =
    dca::phys::params::Parameters<TestConcurrency, StdThreading, dca::profiling::NullProfiler,
                                  Model, RngType, dca::ClusterSolverId::CT_INT>;
using Data = dca::phys::DcaData<Parameters>;
using BaseSolver = dca::phys::solver::CtintClusterSolver<dca::linalg::CPU, Parameters>;
using QmcSolver = dca::phys::solver::StdThreadQmciClusterSolver<BaseSolver>;

// See below, life cycle issue with MPI
// std::unique_ptr<dca::parallel::MPIConcurrency> concurrency;
dca::parallel::MPIConcurrency* concurrency_ptr;

TEST(PosixCtintClusterSolverTest, PerMeasurementIO) {
  static bool update_model = true;

  TestConcurrency& concurrency = *concurrency_ptr;

  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(
      input_dir + "stdthread_ctint_test_per_sample_output.json");
  if (update_model) {
    parameters.update_model();
    parameters.update_domains();
  }
  update_model = false;

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

  auto writer = std::make_shared<dca::io::Writer<TestConcurrency>>(
      concurrency.get_adios(), std::ref(concurrency), parameters.get_output_format(), false);
  writer->open_file(parameters.get_filename_dca(), true);
  // Do one integration step.
  QmcSolver qmc_solver(parameters, data, writer);
  qmc_solver.initialize(0);
  qmc_solver.integrate();
  dca::phys::DcaLoopData<Parameters> loop_data;
  qmc_solver.finalize(loop_data);
  writer->close_file();

  /// \todo add checks of single measurement output here.
}

int main(int argc, char** argv) {
  // This results in a copy constructor beging called at somepoint,  resulting in an MPI_INIT after
  // the finalize. concurrency = std::make_unique<dca::parallel::MPIConcurrency>(argc, argv);
  concurrency_ptr = new dca::parallel::MPIConcurrency(argc, argv);

#ifdef DCA_HAVE_GPU
  dca::linalg::util::printInfoDevices();

  dca::linalg::util::initializeMagma();
#endif
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  // dca::linalg::util::printInfoDevices();

  if (concurrency_ptr->id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  int result = RUN_ALL_TESTS();

  delete concurrency_ptr;

  return result;
}
