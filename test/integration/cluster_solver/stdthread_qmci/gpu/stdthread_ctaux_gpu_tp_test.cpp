// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Confront the MC integration performed on the CPU and GPU over a square lattice with
// nearest-neighbour hopping and on site interaction. The results are expected to be the
// same up to numerical error.

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
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/cluster_solver/stdthread_qmci/gpu/";

using Concurrency = dca::parallel::NoConcurrency;
using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Threading = dca::parallel::stdthread;
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                                 Model, RngType, dca::phys::solver::CT_AUX>;
using Data = dca::phys::DcaData<Parameters>;

using BaseSolverGpu = dca::phys::solver::CtauxClusterSolver<dca::linalg::GPU, Parameters, Data>;
using QmcSolverGpu = dca::phys::solver::StdThreadQmciClusterSolver<BaseSolverGpu>;

using BaseSolverCpu = dca::phys::solver::CtauxClusterSolver<dca::linalg::CPU, Parameters, Data>;
using QmcSolverCpu = dca::phys::solver::StdThreadQmciClusterSolver<BaseSolverCpu>;

TEST(PosixCtauxClusterSolverTest, G_k_w) {
  dca::linalg::util::initializeMagma();
  Concurrency concurrency(0, nullptr);
  if (concurrency.id() == concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(
      input_dir + "stdthread_ctaux_gpu_tp_test_input.json");
  parameters.update_model();
  parameters.update_domains();

  // Initialize data with G0 computation.
  Data data_cpu(parameters), data_gpu(parameters);
  data_cpu.initialize();
  data_gpu.initialize();

  QmcSolverCpu qmc_solver_cpu(parameters, data_cpu);
  RngType::resetCounter();  // Use the same seed for both solvers.
  QmcSolverGpu qmc_solver_gpu(parameters, data_gpu);

  // Do one integration step.
  auto perform_integration = [&](auto& solver) {
    solver.initialize(0);
    solver.integrate();
    dca::phys::DcaLoopData<Parameters> loop_data;
    solver.finalize(loop_data);
  };
  perform_integration(qmc_solver_cpu);
  perform_integration(qmc_solver_gpu);

  const auto err_g = dca::func::util::difference(data_cpu.G_k_w, data_gpu.G_k_w);
  const auto err_g4 =
      dca::func::util::difference(data_cpu.get_G4(), data_gpu.get_G4());

  EXPECT_GE(5e-7, err_g.l_inf);
  EXPECT_GE(5e-7, err_g4.l_inf);
}
