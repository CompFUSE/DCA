// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Anirudha Mirmira (anirudha.mirmira@gmail.com)
//
// Disorder stat test, Step 2 (finite disorder robustness): runs the full CT-AUX disorder ensemble
// average at a small finite potential (V = 0.1) for both the binary and box distributions, and
// checks that the disorder-averaged G_k_w stays finite and sane (no NaN/Inf, nonzero, bounded).
// The aim is NOT accuracy -- it is to confirm the disorder path does not blow up catastrophically
// once V != 0. It replicates the production disorder loop (DcaLoop::workTheClusters + averageGkw):
// generate the ensemble -> per config makeDisorderedG0 -> integrate -> accumulateGkw(weight) ->
// normalize -> translational average -> FT r->k.
//
// Only built when DCA_WITH_DISORDER is ON (see CMakeLists.txt), which is CPU-only.

#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "dca/config/profiler.hpp"
#include "dca/io/json/json_reader.hpp"

// The shared setup header uses a global Scalar and includes ctint_cluster_solver.hpp, which needs
// dca::config::McOptions. Provide both before pulling in the setup header (same pattern as the
// rashba stat test). The square lattice Green's function is real, so Scalar is double.
using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "test/integration/statistical_tests/square_lattice/square_lattice_setup.hpp"

// These pull in dca_shared_types.hpp / domain_aliases.hpp, so they must follow the setup header
// (which defines the time/frequency/cluster domains they alias).
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/dca_loop/disorder/disorder_average.hpp"
#include "dca/phys/dca_loop/disorder/makeDisorderConfigurations.hpp"

dca::testing::DcaMpiTestEnvironment* dca_test_env;

namespace {

using dca::ClusterSolverId;
using SolverParameters = dca::testing::ParametersType<ClusterSolverId::CT_AUX>;
using SolverData = dca::testing::DcaData<ClusterSolverId::CT_AUX>;

// Drive the full disorder ensemble average for the given input file and return the disorder-averaged
// G_k_w. Mirrors DcaLoop::workTheClusters (the disorder loop) followed by DcaLoop::averageGkw.
SolverData::SpGreensFunction runDisorderEnsemble(const std::string& input_file) {
  SolverParameters parameters(dca::util::GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
  parameters.update_model();
  parameters.update_domains();

  SolverData data(parameters);
  data.initialize();

  // Generate the disorder ensemble (binary or box, per the input's distribution).
  dca::phys::MakeDisorderConfigurations<SolverParameters> make_configs;
  std::vector<typename SolverData::DisorderConfiguration> configs;
  std::vector<double> weights;
  make_configs(parameters, configs, weights);
  EXPECT_GT(configs.size(), 0u) << "disorder generator produced no configurations";

  dca::testing::QuantumClusterSolver<ClusterSolverId::CT_AUX> solver(parameters, data, nullptr);

  // --- DcaLoop::workTheClusters disorder loop ---
  data.disorder_G_r_r_w = 0;
  for (std::size_t i = 0; i < configs.size(); ++i) {
    data.makeDisorderedG0(configs[i]);
    solver.initialize(0);
    solver.integrate();
    solver.accumulateGkw(weights[i]);
  }

  // Catastrophic-failure guard on the two-site accumulator, before the transforms, to localize a NaN.
  for (int i = 0; i < data.disorder_G_r_r_w.size(); ++i)
    EXPECT_TRUE(std::isfinite(data.disorder_G_r_r_w(i).real()) &&
                std::isfinite(data.disorder_G_r_r_w(i).imag()))
        << "non-finite disorder_G_r_r_w at linear index " << i;

  // --- DcaLoop::averageGkw ---
  const double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
  data.disorder_G_r_r_w /= total_weight;

  using NuDmn = SolverData::NuDmn;
  using RClusterDmn = SolverData::RClusterDmn;
  using KClusterDmn = SolverData::KClusterDmn;
  using WDmn = SolverData::WDmn;

  dca::phys::disorder::translationalAverage<NuDmn, RClusterDmn, WDmn>(data.disorder_G_r_r_w,
                                                                      data.G_r_w);
  dca::math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(data.G_r_w, data.G_k_w);

  return data.G_k_w;
}

// Hard robustness assertions: every element finite, the result is not identically zero, and its
// magnitude is bounded (a blow-up would show up as a huge or non-finite value).
void expectFiniteAndSane(const SolverData::SpGreensFunction& g_k_w) {
  double max_abs = 0.;
  for (int i = 0; i < g_k_w.size(); ++i) {
    ASSERT_TRUE(std::isfinite(g_k_w(i).real()) && std::isfinite(g_k_w(i).imag()))
        << "non-finite G_k_w at linear index " << i;
    max_abs = std::max(max_abs, std::abs(g_k_w(i)));
  }
  EXPECT_GT(max_abs, 0.) << "disorder-averaged G_k_w is identically zero";
  // |G(k,iw)| <= beta for a fermionic propagator; 1e3 is a very loose blow-up ceiling.
  EXPECT_LT(max_abs, 1e3) << "disorder-averaged G_k_w magnitude is implausibly large (blow-up)";
}

}  // namespace

TEST(DisorderFiniteVTest, BinaryDisorderIsFinite) {
  expectFiniteAndSane(
      runDisorderEnsemble(dca::testing::test_directory + "input_disorder_finite_v_binary.json"));
}

TEST(DisorderFiniteVTest, BoxDisorderIsFinite) {
  expectFiniteAndSane(
      runDisorderEnsemble(dca::testing::test_directory + "input_disorder_finite_v_box.json"));
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca::parallel::MPIConcurrency concurrency(argc, argv);
  dca_test_env = new dca::testing::DcaMpiTestEnvironment(
      concurrency, dca::testing::test_directory + "input_disorder_finite_v_binary.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
