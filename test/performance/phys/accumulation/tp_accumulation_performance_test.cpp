// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// Measure the runtime of the two particle accumulator.

#include "dca/config/haves_defines.hpp"

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator.hpp"

#include <array>
#include <vector>
#include <iostream>
#ifdef DCA_HAVE_CUDA
#include <cuda_profiler_api.h>
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"
#endif  // DCA_HAVE_CUDA

#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/counting_profiler.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/profiling/events/time_event.hpp"

struct ConfigElement {
  double get_tau() const {
    return tau_;
  }
  double get_left_band() const {
    return band_;
  }
  double get_right_band() const {
    return band_;
  }
  double get_left_site() const {
    return r_;
  }
  double get_right_site() const {
    return r_;
  }

  int band_;
  int r_;
  double tau_;
};

using dca::linalg::CPU;
using dca::linalg::GPU;

using Model =
    dca::phys::models::TightBindingModel<dca::phys::models::square_lattice<dca::phys::domains::D4>>;
using Concurrency = dca::parallel::NoConcurrency;
using Profiler = dca::profiling::CountingProfiler<dca::profiling::time_event<std::size_t>>;
using Parameters = dca::phys::params::Parameters<Concurrency, dca::parallel::NoThreading, Profiler,
                                                 Model, void, dca::phys::solver::CT_AUX>;
using Data = dca::phys::DcaData<Parameters>;

using Real = Parameters::TP_measurement_scalar_type;
template <dca::linalg::DeviceType device>
using MatrixPair = std::array<dca::linalg::Matrix<double, device>, 2>;
using Configuration = std::array<std::vector<ConfigElement>, 2>;

void prepareRandomConfig(Configuration& config, MatrixPair<CPU>& M, int n);

using BDmn = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;
using RDmn = typename Parameters::RClusterDmn;

int main(int argc, char** argv) {
  const int n = (argc > 1) ? std::atoi(argv[1]) : 3000;
  const bool skip_cpu = (argc > 2) ? std::atoi(argv[2]) : false;

  const std::string inputs_directory = DCA_SOURCE_DIR "/test/performance/phys/accumulation/";

  Concurrency concurrency(argc, argv);
  Parameters parameters("", concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(inputs_directory + "input.json");

  parameters.update_model();
  parameters.update_domains();

  Data data(parameters);
  data.initialize();

  MatrixPair<CPU> M;
  Configuration config;
  prepareRandomConfig(config, M, n);
  const int sign = 1;

  auto duration = [](dca::profiling::WallTime end, dca::profiling::WallTime start) {
    dca::profiling::Duration elapsed(end, start);
    return elapsed.sec + 1e-6 * elapsed.usec;
  };

  if (!skip_cpu) {
    dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::linalg::CPU> accumulator(
        data.G0_k_w_cluster_excluded, parameters);
    accumulator.resetAccumulation(0);

    // Allows memory to be assigned.
    accumulator.accumulate(M, config, sign);
    accumulator.resetAccumulation(0);

    Profiler::start();
    dca::profiling::WallTime start_time;
    accumulator.accumulate(M, config, sign);
    dca::profiling::WallTime end_time;
    Profiler::stop("tp_accumulation_profile.txt");

    const double time = duration(end_time, start_time);

    std::string precision("double");
    if (std::is_same<float, dca::phys::solver::accumulator::TpAccumulator<Parameters>::Real>::value)
      precision = "single";

    std::cout << "\nExpansion order:\t" << n;
    std::cout << "\nPrecision:\t" << precision;
    std::cout << "\nN positive frequencies:\t" << parameters.get_four_point_fermionic_frequencies();
    std::cout << "\nN bands:\t" << BDmn::dmn_size();
    std::cout << "\nN cluster sites:\t" << RDmn::dmn_size();
    std::cout << "\nType:\t" << dca::phys::toString(parameters.get_four_point_channels().at(0));
    std::cout << "\n\nTpAccumulation CPU time [sec]:\t " << time << "\n";
  }

#ifdef DCA_HAVE_CUDA
  dca::linalg::util::initializeMagma();
  dca::linalg::util::CudaEvent start_event;
  dca::linalg::util::CudaEvent stop_event;

  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::linalg::GPU> gpu_accumulator(
      data.G0_k_w_cluster_excluded, parameters);
  gpu_accumulator.resetAccumulation(0);
  MatrixPair<GPU> M_dev{M[0], M[1]};

  // Allow memory to be assigned.
  gpu_accumulator.accumulate(M_dev, config, sign);
  cudaDeviceSynchronize();
  gpu_accumulator.resetAccumulation(1);

  Profiler::start();
  cudaProfilerStart();

  // Time a single execution.
  start_event.record(*gpu_accumulator.get_stream());
  dca::profiling::WallTime host_start_time;
  gpu_accumulator.accumulate(M_dev, config, sign);
  dca::profiling::WallTime host_end_time;
  stop_event.record(*gpu_accumulator.get_stream());

  const double host_time = duration(host_end_time, host_start_time);
  const double dev_time = dca::linalg::util::elapsedTime(stop_event, start_event);

  std::cout << "\nTpAccumulation GPU: Host time [sec]:\t " << host_time;
  std::cout << "\nTpAccumulation GPU: Device time [sec]:\t " << dev_time << "\n\n";

  // Time a loop.
  cudaDeviceSynchronize();
  const dca::profiling::WallTime loop_start;
  constexpr int n_iters = 10;
  for (int i = 0; i < n_iters; ++i) {
    gpu_accumulator.accumulate(M, config, sign);
  }
  const dca::profiling::WallTime loop_end;
  std::cout << "\nTpAccumulation GPU loop: time per iteration [sec]:\t "
            << duration(loop_end, loop_start) / n_iters << "\n\n";

  cudaProfilerStop();
  Profiler::stop("tp_gpu_accumulation_profile.txt");
#endif  // DCA_HAVE_CUDA
}

void prepareRandomConfig(Configuration& config, MatrixPair<CPU>& M, const int n) {
  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 0);

  for (int s = 0; s < 2; ++s) {
    config[s].resize(n);
    M[s].resize(n);
    for (int i = 0; i < n; ++i) {
      const double tau = rng() - 0.5;
      const int r = rng() * RDmn::dmn_size();
      const int b = rng() * BDmn::dmn_size();
      config[s][i] = ConfigElement{b, r, tau};
    }

    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
        M[s](i, j) = 2 * rng() - 1.;
  }
}
