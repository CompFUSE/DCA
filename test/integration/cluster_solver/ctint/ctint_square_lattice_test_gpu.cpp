// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Confront GPU and CPU runs with CT-INT.
// Model: square lattice with single band and double occupancy repulsion U.

#include <cuda_profiler_api.h>
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
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/cluster_solver/ctint/";

using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Threading = dca::parallel::NoThreading;
using Concurrency = dca::parallel::NoConcurrency;
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                                 Model, RngType, dca::phys::solver::CT_INT>;
using Data = dca::phys::DcaData<Parameters>;

TEST(SquareLatticeTest, GpuSolver) {
  dca::util::GitVersion::print();
  dca::util::Modules::print();

  Concurrency concurrency(1, NULL);
  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir +
                                                           "square_lattice_gpu_input.json");
  parameters.update_model();
  parameters.update_domains();

  Data data_gpu(parameters);
  data_gpu.initialize();
  dca::phys::solver::CtintClusterSolver<dca::linalg::GPU, Parameters, true> qmc_solver_gpu(
      parameters, data_gpu);
  qmc_solver_gpu.initialize(0);
  cudaProfilerStart();
  qmc_solver_gpu.integrate();
  cudaProfilerStop();
  qmc_solver_gpu.finalize();
  EXPECT_NEAR(1.0, qmc_solver_gpu.computeDensity(), 1e-3);

  // Confront with CPU run.
  Data data_cpu(parameters);
  data_cpu.initialize();
  RngType::resetCounter();  // Use the same random numbers.
  using GpuSolver = dca::phys::solver::CtintClusterSolver<dca::linalg::CPU, Parameters, true>;
  GpuSolver qmc_solver_cpu(parameters, data_cpu);
  qmc_solver_cpu.initialize();
  qmc_solver_cpu.integrate();
  qmc_solver_cpu.finalize();

  auto diff = dca::func::util::difference(data_cpu.G_k_w, data_gpu.G_k_w);
  auto tolerance = std::max(1e-10, 100 * std::numeric_limits<GpuSolver::Real>::epsilon());
  EXPECT_GE(tolerance, diff.l_inf);
}
