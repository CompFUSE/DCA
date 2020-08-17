// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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

constexpr bool UPDATE_RESULTS = false;

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/cluster_solver/stdthread_qmci/gpu/";

using Concurrency = dca::parallel::NoConcurrency;
using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Threading = dca::parallel::stdthread;
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                                 Model, RngType, dca::phys::solver::CT_INT>;
using Data = dca::phys::DcaData<Parameters>;
template <dca::linalg::DeviceType device>
using BaseSolver = dca::phys::solver::CtintClusterSolver<device, Parameters, true>;
template <dca::linalg::DeviceType device>
using QmcSolver = dca::phys::solver::StdThreadQmciClusterSolver<BaseSolver<device>>;

TEST(StdthreadCtintGpuTest, GpuVsCpu) {
  Concurrency concurrency(0, nullptr);
  dca::linalg::util::initializeMagma();

  Parameters parameters("", concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir + "threaded_input.json");
  parameters.update_model();
  parameters.update_domains();

  // Initialize data with G0 computation.
  Data data_cpu(parameters), data_gpu(parameters);
  data_cpu.initialize();
  data_gpu.initialize();

  // Do one integration step.
  auto perform_integration = [&](auto& solver) {
    solver.initialize(0);
    solver.integrate();
    dca::phys::DcaLoopData<Parameters> loop_data;
    solver.finalize(loop_data);
  };

  QmcSolver<dca::linalg::GPU> qmc_solver_gpu(parameters, data_gpu);
  perform_integration(qmc_solver_gpu);

  RngType::resetCounter();  // Use the same seed for both solvers.
  QmcSolver<dca::linalg::CPU> qmc_solver_cpu(parameters, data_cpu);
  perform_integration(qmc_solver_cpu);

  const auto err_g = dca::func::util::difference(data_cpu.G_k_w, data_gpu.G_k_w);
  const auto err_g4 = dca::func::util::difference(data_cpu.get_G4()[0], data_gpu.get_G4()[0]);

  EXPECT_GE(5e-7, err_g.l_inf);
  EXPECT_GE(5e-7, err_g4.l_inf);
}
