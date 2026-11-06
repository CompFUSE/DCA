// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Anirudha Mirmira (anirudha.mirmira@gmail.com)
//
// Unit tests for the DISORDERED_G0 walker Green's-function path, covering three chunks:
//   A. DcaData::makeDisorderedG0 builds the dense disordered G0 in imaginary time via the w->t
//      subtract-the-clean-reference Fourier transform (replaces the old per-slice site-diagonal
//      loop). At zero potential the result must reproduce the clean tau G0 unfolded to two indices.
//   B. The CT-AUX walker's G0Interpolation stores the now two-site (R0,R1) disordered G0 in its
//      akima splines: the spline value at each tau node must equal the input G0.
//   C. build_G0_matrix keys the lookups by the two ABSOLUTE sites, not a single displacement: two
//      vertex pairs sharing a displacement but at different absolute sites must give different G0.

#include <complex>
#include <memory>
#include <vector>

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/io/json/json_reader.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
// vertex_singleton must precede the g0 interpolation header, which uses it without including it.
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_cpu.hpp"
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
using RDmn = Data::RClusterDmn;
using KDmn = Data::KClusterDmn;
using WDmn = Data::WDmn;
using TDmn = Data::TDmn;
using Complex = std::complex<double>;

using G0Interp = dca::phys::solver::ctaux::G0Interpolation<dca::linalg::CPU, Parameters>;
using Vertex = dca::phys::solver::ctaux::vertex_singleton;

constexpr char input_file[] = DCA_SOURCE_DIR "/test/unit/phys/dca_data/input_make_disorder.json";

// One set of globally initialized domains / parameters shared by all tests.
Parameters& sharedParameters() {
  static Concurrency concurrency(0, nullptr);
  static Parameters parameters("", concurrency);
  static bool initialized = false;
  if (!initialized) {
    parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
    parameters.update_model();
    parameters.update_domains();
    initialized = true;
  }
  return parameters;
}

std::unique_ptr<Data> makeData() {
  return std::make_unique<Data>(sharedParameters());
}

// Clean, k-diagonal, cluster-excluded G0 in frequency: invertible (diagonal-dominant) nu-blocks,
// distinct per k so the displacement structure is nontrivial. Drives makeDisorderedG0's freq path.
void fillCleanG0kw(Data& data) {
  auto& g0 = data.G0_k_w_cluster_excluded;
  for (int w = 0; w < WDmn::dmn_size(); ++w)
    for (int k = 0; k < KDmn::dmn_size(); ++k)
      for (int n2 = 0; n2 < NuDmn::dmn_size(); ++n2)
        for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1)
          g0(n1, n2, k, w) = (n1 == n2)
                                 ? Complex(3.0 + 0.1 * k + 0.01 * n1, 0.2 + 0.02 * w)
                                 : Complex(0.3, n1 < n2 ? 0.1 : -0.1);
}

// Clean cluster-excluded G0 in (single-displacement) imaginary time: the tau reference that the
// w->t step adds back. Structured in (nu, r, t) so the unfold/add-back has teeth.
void fillCleanG0rt(Data& data) {
  auto& g = data.G0_r_t_cluster_excluded;
  for (int t = 0; t < TDmn::dmn_size(); ++t)
    for (int r = 0; r < RDmn::dmn_size(); ++r)
      for (int n2 = 0; n2 < NuDmn::dmn_size(); ++n2)
        for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1)
          g(n1, n2, r, t) = 0.5 + 0.3 * r + 0.05 * t + 0.1 * (n1 == n2 ? 1.0 : -0.5);
}

// A synthetic two-site disordered G0 in tau that depends on BOTH absolute sites (R0 and R1), not
// just their displacement -- the property the walker rewrite must respect.
void fillSyntheticDisorderedG0rt(Data& data) {
  auto& g = data.disordered_G0_r_r_t_cl_exl;
  for (int t = 0; t < TDmn::dmn_size(); ++t)
    for (int r1 = 0; r1 < RDmn::dmn_size(); ++r1)
      for (int r0 = 0; r0 < RDmn::dmn_size(); ++r0)
        for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1)
          for (int n0 = 0; n0 < NuDmn::dmn_size(); ++n0)
            g(n0, r0, n1, r1, t) =
                (1.0 + 0.5 * r0 + 0.25 * r1 + 0.1 * (n0 + n1)) * (1.0 + 0.1 * t);
}

}  // namespace

// Chunk A: at zero potential the dense disordered tau G0 must equal the clean tau G0 unfolded with
// the same subtract(R2,R1) displacement (delta == 0, so only the added-back reference survives).
TEST(DisorderedWalkerG0Test, ImaginaryTimeUnfoldAtZeroPotential) {
  auto data = makeData();
  fillCleanG0kw(*data);
  fillCleanG0rt(*data);

  typename Data::DisorderConfiguration config;  // zero potential
  data->makeDisorderedG0(config);

  const int n_r = RDmn::dmn_size();
  constexpr double tol = 1e-8;
  for (int t = 0; t < TDmn::dmn_size(); ++t)
    for (int R1 = 0; R1 < n_r; ++R1)
      for (int R2 = 0; R2 < n_r; ++R2) {
        const int d = RDmn::parameter_type::subtract(R2, R1);
        for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1)
          for (int n2 = 0; n2 < NuDmn::dmn_size(); ++n2) {
            const double got = data->disordered_G0_r_r_t_cl_exl(n1, R1, n2, R2, t);
            const double ref = data->G0_r_t_cluster_excluded(n1, n2, d, t);
            EXPECT_NEAR(got, ref, tol)
                << "t=" << t << " R1=" << R1 << " R2=" << R2 << " n1=" << n1 << " n2=" << n2;
          }
      }
}

// Chunk B: the walker akima storage is now keyed by the two sites (R0,R1). The cubic's constant
// term (value at the left tau node) must equal the input disordered G0 at that node, per site pair.
// The two-site akima storage only exists under DISORDERED_G0; in the clean build the storage is
// keyed by a single displacement, so this property is not defined.
TEST(DisorderedWalkerG0Test, WalkerAkimaStoresTwoSiteG0) {
  auto data = makeData();
  fillSyntheticDisorderedG0rt(*data);

  G0Interp g0(0, sharedParameters());
  g0.initialize(*data);

  const auto& ak = g0.getAkimaCoefficients();
  const int n_r = RDmn::dmn_size();
  const int half = TDmn::dmn_size() / 2;
  constexpr double tol = 1e-9;

  // Negative-tau half: shifted node s in [0, half-2] holds the spline for [node s, node s+1], and
  // its constant coefficient a[0] equals the data value at full-tau index s.
  for (int s = 0; s < half - 1; ++s)
    for (int r1 = 0; r1 < n_r; ++r1)
      for (int r0 = 0; r0 < n_r; ++r0)
        for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1)
          for (int n0 = 0; n0 < NuDmn::dmn_size(); ++n0) {
            const double a0 = ak(0, n0, r0, n1, r1, s);
            const double ref = data->disordered_G0_r_r_t_cl_exl(n0, r0, n1, r1, s);
            EXPECT_NEAR(a0, ref, tol) << "s=" << s << " r0=" << r0 << " r1=" << r1;
          }
}

// Chunk C: build_G0_matrix must key by absolute sites. Two vertex pairs with the SAME displacement
// but at different absolute sites must yield different G0 (off-diagonal and site-resolved diagonal).
// Site-resolved G0 only exists under DISORDERED_G0; the clean build collapses both pairs to their
// shared displacement (giving identical G0), so this property is not defined there.
TEST(DisorderedWalkerG0Test, BuildG0MatrixUsesAbsoluteSites) {
  auto data = makeData();
  fillSyntheticDisorderedG0rt(*data);

  G0Interp g0(0, sharedParameters());
  g0.initialize(*data);

  const double beta = sharedParameters().get_beta();
  const double tau_a = 0.3 * beta;
  const double tau_b = 0.1 * beta;

  // Two pairs sharing a displacement but at different absolute sites (4-site chain cluster).
  const int rA0 = 0, rA1 = 1, rB0 = 1, rB1 = 2;
  ASSERT_EQ(RDmn::parameter_type::subtract(rA1, rA0), RDmn::parameter_type::subtract(rB1, rB0))
      << "test setup expects the two pairs to share a displacement";

  auto vtx = [](int nu, int r, double tau) {
    return Vertex(0, dca::phys::e_UP, nu, 0, r, 0, tau, dca::phys::solver::ctaux::HS_UP,
                  dca::phys::solver::ctaux::HS_FIELD_UP, 0);
  };
  std::vector<Vertex> cfg1 = {vtx(0, rA0, tau_a), vtx(0, rA1, tau_b)};
  std::vector<Vertex> cfg2 = {vtx(0, rB0, tau_a), vtx(0, rB1, tau_b)};

  dca::linalg::Matrix<Scalar, dca::linalg::CPU> G0_1, G0_2;
  g0.build_G0_matrix(cfg1, G0_1);
  g0.build_G0_matrix(cfg2, G0_2);

  // Off-diagonal: same displacement, different absolute sites -> different G0.
  EXPECT_GT(std::abs(G0_1(0, 1) - G0_2(0, 1)), 1e-6)
      << "build_G0_matrix collapsed to a displacement; disorder requires absolute sites";
  // Diagonal is site-resolved (on-site G0_dis(R,R) varies site to site).
  EXPECT_GT(std::abs(G0_1(0, 0) - G0_2(0, 0)), 1e-6);

  // Sanity: rebuilding the same configuration is deterministic.
  dca::linalg::Matrix<Scalar, dca::linalg::CPU> G0_1b;
  g0.build_G0_matrix(cfg1, G0_1b);
  EXPECT_DOUBLE_EQ(G0_1(0, 1), G0_1b(0, 1));
}
