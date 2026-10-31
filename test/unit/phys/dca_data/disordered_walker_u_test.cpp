// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Anirudha Mirmira (anirudha.mirmira@gmail.com)
//
// Step F (the U != 0 teeth) for the DISORDERED_G0 walker path.
//
// All other disorder tests run at U = 0, where M = 0 and the reconstructed G = G0_dis regardless of
// what the walker did -- so they cannot see whether the QMC walker actually samples in the disordered
// medium. This test drives the real CT-AUX walker + accumulator at U != 0 and inspects the sampled
// two-site M(R1,R2,w) directly:
//
//   * Determinism: same disorder config + same seed  ->  bitwise-identical M.
//   * Teeth:       different disorder potential V     ->  M MUST change.
//
// The teeth isolates the walker. Comparing the reconstructed local_G_k_w instead would NOT work:
// it differs with V through the Dyson's G0_dis even if the walker ignored V. M, by contrast, is
// produced solely by the walker's sampling, so a V-dependent M proves the walker's bare line is the
// disordered G0. Under the pre-fix bug (walker read the clean G0_r_t_cluster_excluded; V entered only
// the post-QMC Dyson) M would be V-independent and the sensitivity assertion would fail.

#include <algorithm>
#include <complex>
#include <memory>

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/io/json/json_reader.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"

using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"

namespace {

using Concurrency = dca::parallel::NoConcurrency;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
using Parameters =
    dca::phys::params::Parameters<Concurrency, dca::parallel::NoThreading, dca::profiling::NullProfiler,
                                  Model, RngType, dca::ClusterSolverId::CT_AUX,
                                  dca::NumericalTraits<dca::util::RealAlias<Scalar>, Scalar>>;
using Data = dca::phys::DcaData<Parameters>;
using Solver = dca::phys::solver::CtauxClusterSolver<dca::linalg::CPU, Parameters, Data>;
using Walker = Solver::Walker;
using Accumulator = Solver::Accumulator;

using NuDmn = Data::NuDmn;
using RDmn = Data::RClusterDmn;
using WDmn = Data::WDmn;

constexpr char input_file[] = DCA_SOURCE_DIR "/test/unit/phys/dca_data/input_disorder_walker_u.json";

// Domain of the accumulated two-site M(nu0,R0,nu1,R1,w) (mirrors the solver's NuRNuRClusterWDmn).
using MDmn = dca::func::dmn_variadic<NuDmn, RDmn, NuDmn, RDmn, WDmn>;
using MFunction = dca::func::function<std::complex<double>, MDmn>;

double maxAbs(const MFunction& f) {
  double m = 0.;
  for (int i = 0; i < f.size(); ++i)
    m = std::max(m, std::abs(f(i)));
  return m;
}

double maxAbsDiff(const MFunction& a, const MFunction& b) {
  double m = 0.;
  for (int i = 0; i < a.size(); ++i)
    m = std::max(m, std::abs(a(i) - b(i)));
  return m;
}

}  // namespace

TEST(DisorderedWalkerUTest, WalkerMSamplesDisorderedMedium) {
  Concurrency concurrency(0, nullptr);
  Parameters parameters("", concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
  parameters.update_model();
  parameters.update_domains();

  Data data(parameters);
  data.initialize();  // populates a valid bare cluster-excluded G0 from the model H_DCA.

  // Drive the real CT-AUX walker + accumulator for one disorder config and return the accumulated,
  // sign-weighted two-site M the walker sampled. Fresh walker/accumulator/rng each call (same seed),
  // so two calls with the same config are bit-for-bit identical.
  auto sampleM = [&](const typename Data::DisorderConfiguration& config) {
    data.makeDisorderedG0(config);  // sets disordered_G0_r_r_{w,t}_cl_exl; the walker reads the tau one.

    // StdRandomWrapper derives its engine seed from a static per-construction counter
    // (engine_seed = hash(counter_ + seedArg)), so a fixed seed argument alone does NOT reproduce a
    // stream. Reset the counter so every run draws the IDENTICAL random sequence; then the only thing
    // that can move M between runs is the disorder G0 the walker samples in.
    RngType::resetCounter();
    RngType rng(0, 1, parameters.get_seed());
    Walker walker(parameters, data, rng, 0);
    Accumulator accumulator(parameters, data, 0);

    walker.initialize(0);
    accumulator.initialize(0);

    for (int i = 0; i < parameters.get_warm_up_sweeps(); ++i)
      walker.doSweep();
    walker.markThermalized();

    const int n_meas = parameters.get_measurements()[0];
    for (int i = 0; i < n_meas; ++i) {
      walker.doSweep();
      accumulator.updateFrom(walker);
      accumulator.measure();
    }
    accumulator.finalize();

    return MFunction(accumulator.get_sign_times_M_r_w());
  };

  typename Data::DisorderConfiguration zero_config;  // V = 0 everywhere.

  typename Data::DisorderConfiguration nonzero_config;  // staggered site potential, breaks symmetry.
  for (int r = 0; r < RDmn::dmn_size(); ++r)
    for (int nu = 0; nu < NuDmn::dmn_size(); ++nu)
      nonzero_config(nu, r) = (r % 2 == 0) ? 0.6 : -0.4;

  const MFunction M_zero_a = sampleM(zero_config);
  const MFunction M_zero_b = sampleM(zero_config);
  const MFunction M_nonzero = sampleM(nonzero_config);

  // Sanity: U != 0, so the walker actually sampled a nontrivial M.
  EXPECT_GT(maxAbs(M_zero_a), 1e-6) << "expected a nontrivial interacting M at U != 0";

  // Determinism anchor: identical config + identical seed must reproduce M exactly.
  EXPECT_LT(maxAbsDiff(M_zero_a, M_zero_b), 1e-12)
      << "same disorder config and seed must give identical M";

  // The teeth: the walker's sampled M must respond to the disorder potential, i.e. the walker
  // propagates in the disordered G0. Under the old bug (walker in the clean G0) this difference is 0.
  EXPECT_GT(maxAbsDiff(M_zero_a, M_nonzero), 1e-6)
      << "walker M did not change with V: the walker is not sampling the disordered G0";
}
