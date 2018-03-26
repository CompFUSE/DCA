// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// No-change test for CT-INT.
// Square lattice with single band and double occupancy repulsion U.

#ifndef DCA_HAVE_CUDA
#pragma error "CUDA must be enabled"
#endif

constexpr bool UPDATE_RESULTS = false;

#include <cuda_profiler_api.h>
#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/function/function.hpp"
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

const std::string input_dir =
    DCA_SOURCE_DIR "/test/integration/ctint/";

using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Threading = dca::parallel::NoThreading;
using Concurrency = dca::parallel::NoConcurrency;
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                                 Model, RngType, dca::phys::solver::CT_INT>;
using Data = dca::phys::DcaData<Parameters>;
using QmcSolverType = dca::phys::solver::CtintClusterSolver<dca::linalg::GPU, Parameters>;

TEST(SquareLatticeTest, GpuSolver) {
  dca::util::GitVersion::print();
  dca::util::Modules::print();

  Concurrency concurrency(1, NULL);
  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir + "square_lattice_input.json");
  parameters.update_model();
  parameters.update_domains();

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

  // Do one integration step.
  QmcSolverType qmc_solver_gpu(parameters, data);
  qmc_solver_gpu.initialize(0);
  cudaProfilerStart();
  qmc_solver_gpu.integrate();
  cudaProfilerStop();
  qmc_solver_gpu.finalize();

  EXPECT_NEAR(1.0, qmc_solver_gpu.computeDensity(), 1e-5);

  if (not UPDATE_RESULTS) {
    // Read and confront with previous run
    typeof(data.G_k_w) G_k_w_check(data.G_k_w.get_name());
    dca::io::HDF5Reader reader;
    reader.open_file(input_dir + "square_lattice_result_gpu.hdf5");
    reader.open_group("functions");
    reader.execute(G_k_w_check);
    reader.close_group(), reader.close_file();

    for (int i = 0; i < G_k_w_check.size(); i++) {
      EXPECT_NEAR(G_k_w_check(i).real(), data.G_k_w(i).real(), 1e-7);
      EXPECT_NEAR(G_k_w_check(i).imag(), data.G_k_w(i).imag(), 1e-7);
    }
  }
  else {
    //  Write results
    dca::io::HDF5Writer writer;
    writer.open_file("square_lattice_result_gpu.hdf5");
    writer.open_group("functions");
    writer.execute(data.G_k_w);
    writer.close_group(), writer.close_file();
  }
}
