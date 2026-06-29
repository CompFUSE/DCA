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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "dca/config/cmake_options.hpp"
#include "dca/config/threading.hpp"

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/config/profiler.hpp"

using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/cluster_solver/stdthread_qmci/gpu/";

using TestConcurrency = dca::parallel::NoConcurrency;
using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
// This test will end up testing whatever lattice you've selected as DCA_LATTICE
//using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Parameters =
    dca::phys::params::Parameters<TestConcurrency, Threading, dca::profiling::NullProfiler, Model,
                                  RngType, dca::ClusterSolverId::CT_AUX,
                                  dca::NumericalTraits<dca::util::RealAlias<Scalar>, Scalar>>;
using Data = dca::phys::DcaData<Parameters>;

using BaseSolverGpu = dca::phys::solver::CtauxClusterSolver<dca::linalg::GPU, Parameters, Data>;
using QmcSolverGpu = dca::phys::solver::StdThreadQmciClusterSolver<BaseSolverGpu>;

using BaseSolverCpu = dca::phys::solver::CtauxClusterSolver<dca::linalg::CPU, Parameters, Data>;
using QmcSolverCpu = dca::phys::solver::StdThreadQmciClusterSolver<BaseSolverCpu>;

template <class Function>
dca::func::util::Difference differenceForTransfer(const Function& cpu, const Function& gpu,
                                                  const int q_index, const int w_index) {
  double l1 = 0.;
  double l2 = 0.;
  double linf = 0.;

  double l1_error = 0.;
  double l2_error = 0.;
  double linf_error = 0.;

  for (int i = 0; i < cpu.size(); ++i) {
    const auto subind = cpu.linind_2_subind(i);
    if (subind[8] != std::size_t(q_index) || subind[9] != std::size_t(w_index))
      continue;

    const double ref = std::abs(cpu(i));
    l1 += ref;
    l2 += ref * ref;
    linf = std::max(linf, ref);

    const double err = std::abs(cpu(i) - gpu(i));
    l1_error += err;
    l2_error += err * err;
    linf_error = std::max(linf_error, err);
  }

  const auto relative = [](const double numerator, const double denominator) {
    if (denominator > 0)
      return numerator / denominator;
    return numerator == 0 ? 0. : std::numeric_limits<double>::infinity();
  };

  return {relative(l1_error, l1), relative(std::sqrt(l2_error), std::sqrt(l2)),
          relative(linf_error, linf)};
}

TEST(PosixCtauxClusterSolverTest, G_k_w) {
  dca::linalg::util::initializeMagma();
  TestConcurrency concurrency(0, nullptr);
  if (concurrency.id() == concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  const char* input_override = std::getenv("DCA_CPU_GPU_TEST_INPUT");
  const std::string input_file =
      input_override ? input_override : input_dir + "threaded_input.json";
  std::cout << "Reading input " << input_file << "\n";
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
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
  std::cout << "Creating CPU solver\n";
  QmcSolverCpu qmc_solver_cpu(parameters, data_cpu, nullptr);
  perform_integration(qmc_solver_cpu);

  RngType::resetCounter();  // Use the same seed for both solvers.
  // i.e. assume that the consumption of random numbers is exactly the same in sequence for gpu and
  // cpu. This will not be true unless walker and accumulator share a thread.
  std::cout << "Creating GPU solver\n";
  QmcSolverGpu qmc_solver_gpu(parameters, data_gpu, nullptr);
  perform_integration(qmc_solver_gpu);

  // dca::func::util::difference reports relative CPU/GPU norm errors:
  // err_i = abs(cpu_i - gpu_i), ref_i = abs(cpu_i)
  // l1 = sum_i err_i / sum_i ref_i
  // l2 = sqrt(sum_i err_i^2 / sum_i ref_i^2)
  // l_inf = max_i err_i / max_i ref_i
  const auto err_g = dca::func::util::difference(data_cpu.G_k_w, data_gpu.G_k_w);
  std::cout << "CPU/GPU G_k_w relative differences:"
            << " l1=" << err_g.l1 << " l2=" << err_g.l2 << " l_inf=" << err_g.l_inf
            << "\n";

  EXPECT_GE(5e-7, err_g.l_inf);

  const auto& g4_cpu = data_cpu.get_G4();
  const auto& g4_gpu = data_gpu.get_G4();
  constexpr double g4_tolerance = 5e-5;
  ASSERT_EQ(g4_cpu.size(), g4_gpu.size());
  for (std::size_t channel = 0; channel < g4_cpu.size(); ++channel) {
    const auto err_g4 = dca::func::util::difference(g4_cpu[channel], g4_gpu[channel]);
    std::cout << "CPU/GPU G4 channel " << channel << " all-transfer relative differences:"
              << " l1=" << err_g4.l1 << " l2=" << err_g4.l2 << " l_inf=" << err_g4.l_inf
              << "\n";
    EXPECT_GE(g4_tolerance, err_g4.l_inf) << "G4 channel: " << channel;

    const auto& sizes = g4_cpu[channel].getDomainSizes();
    ASSERT_EQ(10, sizes.size());
    ASSERT_EQ(sizes, g4_gpu[channel].getDomainSizes());

    for (int w_index = 0; w_index < sizes[9]; ++w_index) {
      for (int q_index = 0; q_index < sizes[8]; ++q_index) {
        const auto err_g4_transfer =
            differenceForTransfer(g4_cpu[channel], g4_gpu[channel], q_index, w_index);
        if (err_g4_transfer.l_inf > g4_tolerance) {
          ADD_FAILURE() << "G4 channel: " << channel << ", q_index: " << q_index
                        << ", w_index: " << w_index
                        << " relative differences: l1=" << err_g4_transfer.l1
                        << " l2=" << err_g4_transfer.l2
                        << " l_inf=" << err_g4_transfer.l_inf;
        }
      }
    }
  }
}
