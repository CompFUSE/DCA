// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Performance test for CT-INT.
// Bilayer lattice with two bands and 36 sites.

#include <iostream>
#include <string>
#ifdef DCA_HAVE_CUDA
#include <cuda_profiler_api.h>
#endif

#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/profiling/counting_profiler.hpp"
#include "dca/profiling/events/time_event.hpp"
#include "dca/util/ignore.hpp"

const std::string input_dir = DCA_SOURCE_DIR "/test/performance/phys/ctint/";

using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
using Lattice = dca::phys::models::bilayer_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Threading = dca::parallel::NoThreading;
using Concurrency = dca::parallel::NoConcurrency;
using Profiler = dca::profiling::CountingProfiler<dca::profiling::time_event<std::size_t>>;
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, Profiler, Model, RngType,
                                                 dca::phys::solver::CT_INT>;
using Data = dca::phys::DcaData<Parameters>;
template <dca::linalg::DeviceType device_t>
using QmcSolverType = dca::phys::solver::CtintClusterSolver<device_t, Parameters, true>;

int main(int argc, char** argv) {
  bool test_cpu(true), test_gpu(true);
  int submatrix_size = -1;
  dca::util::ignoreUnused(test_gpu);
  for (int i = 0; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--skip_cpu")
      test_cpu = false;
    else if (arg == "--skip_gpu")
      test_gpu = false;
    else if (arg == "--submatrix_size")
      submatrix_size = std::atoi(argv[i + 1]);
  }

  Concurrency concurrency(1, NULL);
  Parameters parameters("", concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir +
                                                           "bilayer_lattice_input.json");
  if (submatrix_size != -1)
    parameters.setMaxSubmatrixSize(submatrix_size);

  parameters.update_model();
  parameters.update_domains();

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

  auto printTime = [](const std::string& str, const auto& start, const auto& end) {
    dca::profiling::Duration time(end, start);
    std::cout << str << ": time taken: " << time.sec + 1e-6 * time.usec << std::endl;
  };

  std::cout << "Integrating with max-submatrix-size: " << parameters.getMaxSubmatrixSize()
            << std::endl;

  if (test_cpu) {
    std::cout << "\n\n  *********** CPU integration  ***************\n\n";
    std::cout.flush();

    // TODO: always start if the profiler supports the writing of multiple files.
    if (!test_gpu)
      Profiler::start();

    // Do one integration step.
    QmcSolverType<dca::linalg::CPU> qmc_solver(parameters, data);
    qmc_solver.initialize(0);
    // Timed section.
    dca::profiling::WallTime start_t;
    qmc_solver.integrate();
    dca::profiling::WallTime integration_t;

    //    qmc_solver.finalize();
    //    dca::profiling::WallTime finalize_t;

    std::cout << std::endl;
    printTime("Integration CPU", start_t, integration_t);
    //    printTime("Finalization", integration_t, finalize_t);

    if (!test_gpu)
      Profiler::stop(concurrency, "profile_cpu.txt");
  }

#ifdef DCA_HAVE_CUDA
  if (test_gpu) {
    std::cout << "\n\n  *********** GPU integration  ***************\n\n";
    std::cout.flush();

    Profiler::start();

    RngType::resetCounter();
    QmcSolverType<dca::linalg::GPU> solver_gpu(parameters, data);
    solver_gpu.initialize();

    // Timed section.
    cudaProfilerStart();
    dca::profiling::WallTime start_t;
    solver_gpu.integrate();
    dca::profiling::WallTime integration_t;
    cudaProfilerStop();

    Profiler::stop(concurrency, "profile_gpu.txt");

    std::cout << std::endl;
    printTime("Integration GPU", start_t, integration_t);
  }
#endif  // DCA_HAVE_CUDA
}
