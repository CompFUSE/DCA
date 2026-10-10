// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Anirudha Mirmira (anirudha.mirmira@gmail.com)
//
// Unit test for the DISORDERED_G0 disorder-averaging pipeline as wired in DcaLoop::workTheClusters
// / averageGkw, exercised on the real DcaData (2x4 cluster, Nc=8) WITHOUT running the QMC.
//
// With U = 0 the interaction vertices vanish, so M = 0 and the per-configuration Dyson collapses to
// G_c = G0_dis,c. We therefore drive the real pipeline kernels directly per configuration:
//   makeDisorderedG0(config) -> accumulateDisorderDyson(G0_dis, M=0, ...) -> /total_weight
//   -> translationalAverage -> FunctionTransform r->k,
// exactly as the loop does, and check:
//
//   1. Clean limit (zero disorder potential): the whole machinery is the identity on G0 - the
//      disorder-averaged G_k_w reproduces the clean cluster-excluded G0_k_w (the unfold, projection
//      and r<->k transforms compose to identity).
//   2. The weighted disorder average equals the explicit weighted mean of the per-configuration
//      disordered G0 (validates the accumulate + normalize wiring), and the result is finite and
//      actually shifted away from the clean G0 by the disorder.

#include "dca/phys/dca_loop/disorder/disorder_average.hpp"

#include <cmath>
#include <complex>
#include <memory>
#include <numeric>
#include <vector>

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/io/json/json_reader.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"

#include "dca/phys/parameters/parameters.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "test/unit/phys/dca_step/cluster_solver/stub_rng.hpp"

namespace {

using Scalar = double;
using Concurrency = dca::parallel::NoConcurrency;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Parameters =
    dca::phys::params::Parameters<Concurrency, dca::parallel::NoThreading, dca::profiling::NullProfiler,
                                  Model, dca::testing::StubRng, dca::ClusterSolverId::CT_AUX,
                                  dca::NumericalTraits<dca::util::RealAlias<Scalar>, Scalar>>;
using Data = dca::phys::DcaData<Parameters>;

using NuDmn = Data::NuDmn;
using RDmn = Data::RClusterDmn;  // 2x4 cluster, Nc=8
using KDmn = Data::KClusterDmn;
using WDmn = Data::WDmn;
using Complex = std::complex<double>;

constexpr char input_file[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_loop/disorder/input_disorder_loop.json";

std::unique_ptr<Data> makeData() {
  static Concurrency concurrency(0, nullptr);
  static Parameters parameters("", concurrency);
  static bool initialized = false;
  if (!initialized) {
    parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
    parameters.update_model();
    parameters.update_domains();
    initialized = true;
  }
  return std::make_unique<Data>(parameters);
}

// Clean, cluster-excluded G0 in frequency: distinct, diagonal-dominant (invertible) nu-blocks per k.
void fillCleanG0(Data& data) {
  auto& g0 = data.G0_k_w_cluster_excluded;
  for (int w = 0; w < WDmn::dmn_size(); ++w)
    for (int k = 0; k < KDmn::dmn_size(); ++k)
      for (int n2 = 0; n2 < NuDmn::dmn_size(); ++n2)
        for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1)
          g0(n1, n2, k, w) = (n1 == n2) ? Complex(3.0 + 0.1 * k + 0.01 * n1, 0.2 + 0.02 * w)
                                        : Complex(0.3, n1 < n2 ? 0.1 : -0.1);
}

// makeDisorderedG0's legacy time-domain loop inverts nu-blocks of G0_r_t_cluster_excluded; give it
// invertible (identity) blocks. Its output is unused by the DISORDERED_G0 path.
void fillG0rtIdentity(Data& data) {
  auto& g = data.G0_r_t_cluster_excluded;
  g = Scalar(0);
  for (int t = 0; t < Data::TDmn::dmn_size(); ++t)
    for (int r = 0; r < RDmn::dmn_size(); ++r)
      for (int n = 0; n < NuDmn::dmn_size(); ++n)
        g(n, n, r, t) = 1.0;
}

// Runs the real disorder-averaging pipeline (loop body of workTheClusters + averageGkw) with M = 0.
// beta is irrelevant here because the G0*M*G0 term vanishes for M = 0.
void runDisorderPipeline(Data& data, const std::vector<Data::DisorderConfiguration>& configs,
                         const std::vector<double>& weights) {
  Data::Disordered_G0_r_r_w m_zero("m_zero");
  m_zero = 0;
  data.disorder_G_r_r_w = 0;
  for (std::size_t c = 0; c < configs.size(); ++c) {
    data.makeDisorderedG0(configs[c]);
    dca::phys::disorder::accumulateDisorderDyson<NuDmn, RDmn, WDmn, Scalar>(
        data.disordered_G0_r_r_w_cl_exl, m_zero, /*beta=*/1.0, weights[c], data.disorder_G_r_r_w);
  }
  const double total = std::accumulate(weights.begin(), weights.end(), 0.0);
  data.disorder_G_r_r_w /= total;
  dca::phys::disorder::translationalAverage<NuDmn, RDmn, WDmn>(data.disorder_G_r_r_w, data.G_r_w);
  dca::math::transform::FunctionTransform<RDmn, KDmn>::execute(data.G_r_w, data.G_k_w);
}

}  // namespace

// Clean limit: zero disorder potential => G0_dis = clean G0, and the full disorder machinery
// (unfold -> dense identity Dyson -> disorder average -> translational average -> FT r->k) must
// reproduce the clean cluster-excluded G0_k_w.
TEST(DisorderLoopTest, CleanLimitReproducesClusterExcludedG0) {
  auto data = makeData();
  fillCleanG0(*data);
  fillG0rtIdentity(*data);

  // Three identical zero-potential configurations, equal weights.
  std::vector<Data::DisorderConfiguration> configs(3);  // default-constructed = zero potential
  std::vector<double> weights(3, 1.0 / 3.0);

  runDisorderPipeline(*data, configs, weights);

  constexpr double tol = 1e-9;
  for (int w = 0; w < WDmn::dmn_size(); ++w)
    for (int k = 0; k < KDmn::dmn_size(); ++k)
      for (int n2 = 0; n2 < NuDmn::dmn_size(); ++n2)
        for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1) {
          const Complex got = data->G_k_w(n1, n2, k, w);
          const Complex ref = data->G0_k_w_cluster_excluded(n1, n2, k, w);
          EXPECT_NEAR(got.real(), ref.real(), tol)
              << "(n1,n2,k,w)=(" << n1 << "," << n2 << "," << k << "," << w << ")";
          EXPECT_NEAR(got.imag(), ref.imag(), tol);
        }
}

// With distinct nonzero configurations and weights, the accumulated/normalized two-site G must equal
// the explicit weighted mean of the per-config disordered G0 (M = 0 => G_c = G0_dis,c). Also checks
// the averaged G_k_w is finite and actually shifted away from the clean G0 by the disorder.
TEST(DisorderLoopTest, DisorderAverageIsWeightedMeanOfConfigs) {
  auto data = makeData();
  fillCleanG0(*data);
  fillG0rtIdentity(*data);

  const int nu = NuDmn::dmn_size();
  const int nr = RDmn::dmn_size();

  // Three distinct, nonzero, site-local disorder configurations.
  std::vector<Data::DisorderConfiguration> configs(3);
  for (int c = 0; c < 3; ++c)
    for (int r = 0; r < nr; ++r)
      for (int n = 0; n < nu; ++n)
        configs[c](n, r) = 0.4 * (c + 1) * (1.0 + 0.1 * r) - 0.2 * n;
  const std::vector<double> weights = {0.2, 0.5, 0.3};
  const double total = std::accumulate(weights.begin(), weights.end(), 0.0);

  // Independent expected weighted mean of the per-config disordered G0 (M = 0 => G_c = G0_dis,c).
  Data::Disordered_G0_r_r_w expected("expected");
  expected = 0;
  for (int c = 0; c < 3; ++c) {
    data->makeDisorderedG0(configs[c]);
    for (int idx = 0; idx < expected.size(); ++idx)
      expected(idx) += weights[c] * data->disordered_G0_r_r_w_cl_exl(idx);
  }
  for (int idx = 0; idx < expected.size(); ++idx)
    expected(idx) /= total;

  runDisorderPipeline(*data, configs, weights);

  // The pipeline never modifies G0_k_w_cluster_excluded, so it still holds the clean G0.
  constexpr double tol = 1e-9;
  for (int idx = 0; idx < expected.size(); ++idx) {
    EXPECT_NEAR(data->disorder_G_r_r_w(idx).real(), expected(idx).real(), tol) << "idx=" << idx;
    EXPECT_NEAR(data->disorder_G_r_r_w(idx).imag(), expected(idx).imag(), tol) << "idx=" << idx;
  }

  // The averaged G_k_w is finite, and the disorder shifts it away from the clean G0.
  bool differs = false;
  for (int idx = 0; idx < data->G_k_w.size(); ++idx) {
    ASSERT_TRUE(std::isfinite(data->G_k_w(idx).real()) && std::isfinite(data->G_k_w(idx).imag()));
    if (std::abs(data->G_k_w(idx) - data->G0_k_w_cluster_excluded(idx)) > 1e-6)
      differs = true;
  }
  EXPECT_TRUE(differs) << "nonzero disorder must shift the averaged G_k_w away from the clean G0";
}
