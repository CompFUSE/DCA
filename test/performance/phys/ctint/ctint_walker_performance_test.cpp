// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
#include "dca/phys/dca_step/cluster_solver/ctint/details/shrink_G0.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_choice.hpp"
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

using Real = float;

using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
using Lattice = dca::phys::models::bilayer_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using NoThreading = dca::parallel::NoThreading;
using Concurrency = dca::parallel::NoConcurrency;
using Profiler = dca::profiling::CountingProfiler<dca::profiling::time_event<std::size_t>>;
using Parameters = dca::phys::params::Parameters<Concurrency, NoThreading, Profiler, Model, RngType,
                                                 dca::phys::solver::CT_INT>;
using Data = dca::phys::DcaData<Parameters>;
template <dca::linalg::DeviceType device_t>

using Walker = dca::phys::solver::ctint::CtintWalkerChoice<device_t, Parameters, true, Real>;

using BDmn = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;
using RDmn = Parameters::RClusterDmn;
using BBRDmn = dca::func::dmn_variadic<BDmn, BDmn, RDmn>;

int main(int argc, char** argv) {
  bool test_cpu(true), test_gpu(true);
  int submatrix_size = -1;
  int n_sweeps = -1;
  dca::util::ignoreUnused(test_gpu);
  for (int i = 0; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--skip_cpu")
      test_cpu = false;
    else if (arg == "--skip_gpu")
      test_gpu = false;
    else if (arg == "--submatrix_size")
      submatrix_size = std::atoi(argv[i + 1]);
    else if (arg == "--n_sweeps")
      n_sweeps = std::atoi(argv[i + 1]);
  }

  Concurrency concurrency(1, NULL);
  Parameters parameters("", concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir +
                                                           "bilayer_lattice_input.json");
  if (submatrix_size != -1)
    parameters.setMaxSubmatrixSize(submatrix_size);
  if (n_sweeps == -1)
    n_sweeps = parameters.get_warm_up_sweeps();

  parameters.update_model();
  parameters.update_domains();

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

#ifdef DCA_HAVE_CUDA
  constexpr dca::linalg::DeviceType device = dca::linalg::GPU;
#else
  constexpr dca::linalg::DeviceType device = dca::linalg::CPU;
#endif  // DCA_HAVE_CUDA

  dca::phys::solver::ctint::G0Interpolation<device, Real> g0(
      dca::phys::solver::ctint::details::shrinkG0(data.G0_r_t));

  BBRDmn bbr_dmn;
  Walker<device>::setDMatrixBuilder(g0);
  Walker<device>::setDMatrixAlpha(parameters.getAlphas(), false);
  Walker<device>::setInteractionVertices(data, parameters);

  auto printTime = [](const std::string& str, const auto& start, const auto& end) {
    dca::profiling::Duration time(end, start);
    std::cout << str << ": time taken: " << time.sec + 1e-6 * time.usec << std::endl;
  };

  auto do_sweeps = [n_sweeps, &parameters](auto& walker) {
    walker.fixStepsPerSweep(parameters.getInitialConfigurationSize());
    for (int i = 0; i < n_sweeps; ++i) {
      walker.doSweep();
      walker.updateShell(i, n_sweeps);
    }
  };

  std::cout << "Integrating with max-submatrix-size: " << parameters.getMaxSubmatrixSize()
            << std::endl;

  if (test_cpu) {
    std::cout << "\n\n  *********** CPU integration  ***************\n" << std::endl;

    // TODO: always start if the profiler supports the writing of multiple files.
    if (!test_gpu)
      Profiler::start();

    // Do one integration step.
    RngType rng(0, 1, 0);
    Walker<dca::linalg::CPU> walker(parameters, data, rng);
    walker.initialize();

    // Timed section.
    dca::profiling::WallTime start_t;
    do_sweeps(walker);
    dca::profiling::WallTime integration_t;

    std::cout << std::endl;
    printTime("Integration CPU", start_t, integration_t);

    if (!test_gpu)
      Profiler::stop(concurrency, "profile_cpu.txt");
  }

#ifdef DCA_HAVE_CUDA
  if (test_gpu) {
    std::cout << "\n\n  *********** GPU integration  ***************\n\n";
    std::cout.flush();

    Profiler::start();

    RngType::resetCounter();
    RngType rng(0, 1, 0);
    Walker<dca::linalg::GPU> walker_gpu(parameters, data, rng, 0);
    walker_gpu.initialize();

    // Timed section.
    cudaProfilerStart();
    dca::profiling::WallTime start_t;
    do_sweeps(walker_gpu);
    dca::profiling::WallTime integration_t;
    cudaProfilerStop();

    Profiler::stop(concurrency, "profile_gpu.txt");

    std::cout << std::endl;
    printTime("Integration GPU", start_t, integration_t);
  }
#endif  // DCA_HAVE_CUDA
}
