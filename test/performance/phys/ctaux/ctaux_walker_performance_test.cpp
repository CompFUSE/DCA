// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Performance test for CT-AUX walker.
// Bilayer lattice with two bands and 36 sites.

#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_walker.hpp"

#include <iostream>
#include <string>

#include "dca/config/mc_options.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/profiling/null_profiler.hpp"

#ifdef DCA_HAVE_CUDA
#include <cuda_profiler_api.h>
#include "dca/profiling/cuda_profiler.hpp"
#endif

const std::string input_dir = DCA_SOURCE_DIR "/test/performance/phys/ctaux/";

#ifdef DCA_HAVE_CUDA
constexpr dca::linalg::DeviceType device = dca::linalg::GPU;
const std::string device_name = "GPU";
using Profiler = dca::profiling::CudaProfiler;
#else
constexpr dca::linalg::DeviceType device = dca::linalg::CPU;
const std::string device_name = "CPU";
using Profiler = dca::profiling::NullProfiler;
#endif  // DCA_HAVE_CUDA

using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using NoThreading = dca::parallel::NoThreading;
using Concurrency = dca::parallel::NoConcurrency;
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, Profiler, Model, RngType,
                                                 dca::phys::solver::CT_AUX>;
using Data = dca::phys::DcaData<Parameters>;
using Real = dca::config::McOptions::MCScalar;
using Walker = dca::phys::solver::ctaux::CtauxWalker<device, Parameters, Data, Real>;

int main(int argc, char** argv) {
  int submatrix_size = -1;
  int n_warmup = 30;
  int n_sweeps = 10;
  int n_walkers = -1;

  for (int i = 0; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--submatrix_size")
      submatrix_size = std::atoi(argv[i + 1]);
    else if (arg == "--n_sweeps")
      n_sweeps = std::atoi(argv[i + 1]);
    else if (arg == "--n_walkers")
      n_walkers = std::atoi(argv[i + 1]);
  }

  Concurrency concurrency(argc, argv);
  Parameters parameters("", concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir + "square_lattice_input.json");

  parameters.update_model();
  parameters.update_domains();

  if (submatrix_size != -1)
    parameters.set_max_submatrix_size(submatrix_size);
  if (n_walkers == -1)
    n_walkers = parameters.get_walkers();

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

  auto printTime = [](const std::string& str, const auto& start, const auto& end) {
    dca::profiling::Duration time(end, start);
    std::cout << str << ": time taken: " << time.sec + 1e-6 * time.usec << std::endl;
  };

  auto do_sweeps = [](auto& walker, int n, bool verbose) {
    for (int i = 0; i < n; ++i) {
      walker.doSweep();
      if (verbose)
        walker.updateShell(i, n);
    }
  };
  std::cout << "\n\n  *********** " + device_name + " integration  ***************\n";
  std::cout << "Nr walkers: " << n_walkers << "\n\n";

#ifdef DCA_HAVE_CUDA
  dca::linalg::util::initializeMagma();
  dca::linalg::util::resizeHandleContainer(n_walkers);
#endif  // DCA_HAVE_CUDA

  RngType::resetCounter();
  std::vector<RngType> rngs;
  std::vector<Walker> walkers;
  rngs.reserve(n_walkers);
  walkers.reserve(n_walkers);
  for (int i = 0; i < n_walkers; ++i) {
    rngs.emplace_back(0, 1, 0);
    walkers.emplace_back(parameters, data, rngs.back(), i);
  }

  std::vector<std::future<void>> fs;
  dca::parallel::ThreadPool pool(n_walkers);
  for (int i = 0; i < n_walkers; ++i) {
    fs.push_back(pool.enqueue([&do_sweeps, &walkers, i, n_warmup]() {
      walkers[i].initialize(0);
      do_sweeps(walkers[i], n_warmup, i == 0);
      walkers[i].markThermalized();
    }));
  }

  for (auto& f : fs)
    f.get();
  fs.clear();
  std::cout << "\n Warmed up.\n" << std::endl;

  // Timed section.
#ifdef DCA_HAVE_CUDA
  cudaProfilerStart();
#endif  // DCA_HAVE_CUDA
  Profiler::start();
  dca::profiling::WallTime start_t;

  for (int i = 0; i < n_walkers; ++i) {
    fs.push_back(pool.enqueue(
        [&do_sweeps, &walkers, i, n_sweeps]() { do_sweeps(walkers[i], n_sweeps, i == 0); }));
  }
  for (auto& f : fs)
    f.get();

  dca::profiling::WallTime integration_t;
  Profiler::stop("profile.txt");
#ifdef DCA_HAVE_CUDA
  cudaProfilerStop();
#endif  // DCA_HAVE_CUDA

  std::cout << std::endl;
  printTime("Integration " + device_name, start_t, integration_t);
}
