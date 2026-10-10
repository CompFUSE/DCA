// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Anirudha Mirmira (anirudha.mirmira@gmail.com)
//
// Unit tests for the real-space (R1,R2) disordered-G0 construction in DcaData::makeDisorderedG0
// used by the DISORDERED_G0 path:
//   1. the single k->r FunctionTransform of the cluster-excluded clean G0,
//   2. unfolding the translationally invariant clean G0 into the full two space-index quantity
//      disordered_G0_r_r_w_cl_exl(nu1,R1,nu2,R2,w) = g(R1 - R2) via the cluster subtract table,
//   3. the per-frequency disorder Dyson in the dense (nu*N_R) basis (invert -> add a block ->
//      invert).
//
// The first two tests drive the real makeDisorderedG0 (the disorder block it injects is currently
// a zero dummy, so its output equals the clean unfolded G0 -- a clean inversion round trip). The
// third test exercises the invert -> add nonzero block -> invert pipeline that the real disorder
// potential will use, on a genuine clean block extracted from makeDisorderedG0's output.

#include <complex>
#include <memory>

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/io/json/json_reader.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/math/util/vector_operations.hpp"
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
using RDmn = Data::RClusterDmn;
using KDmn = Data::KClusterDmn;
using WDmn = Data::WDmn;
using TDmn = Data::TDmn;
using Complex = std::complex<double>;

constexpr char input_file[] = DCA_SOURCE_DIR "/test/unit/phys/dca_data/input_make_disorder.json";

// Builds a fresh DcaData sharing one set of (globally initialized) domains.
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

// Clean, k-diagonal, cluster-excluded G0 in frequency: diagonal-dominant (invertible) nu-blocks,
// distinct per k so the displacement structure is nontrivial.
void fillCleanG0(Data& data) {
  auto& g0 = data.G0_k_w_cluster_excluded;
  for (int w = 0; w < WDmn::dmn_size(); ++w)
    for (int k = 0; k < KDmn::dmn_size(); ++k)
      for (int n2 = 0; n2 < NuDmn::dmn_size(); ++n2)
        for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1)
          g0(n1, n2, k, w) = (n1 == n2)
                                 ? Complex(3.0 + 0.1 * k + 0.01 * n1, 0.2 + 0.02 * w)
                                 : Complex(0.3, n1 < n2 ? 0.1 : -0.1);
}

// The active (legacy) first loop in makeDisorderedG0 inverts nu-blocks of G0_r_t_cluster_excluded;
// give it invertible (identity) blocks so it is well defined. Its output is unused here.
void fillG0rtIdentity(Data& data) {
  auto& g = data.G0_r_t_cluster_excluded;
  g = Scalar(0);
  for (int t = 0; t < TDmn::dmn_size(); ++t)
    for (int r = 0; r < RDmn::dmn_size(); ++r)
      for (int n = 0; n < NuDmn::dmn_size(); ++n)
        g(n, n, r, t) = 1.0;
}

}  // namespace

// makeDisorderedG0 should FunctionTransform the clean k-space G0 to (r,w) and unfold it to the full
// two space-index quantity. With the (currently zero) disorder block, the invert/invert round trip
// returns the clean values, so disordered_G0_r_r_w_cl_exl must equal the manual inverse-DFT of the
// clean G0 unfolded via subtract(R2,R1).
TEST(MakeDisorderedG0Test, KtoRFtAndUnfold) {
  auto data = makeData();
  fillCleanG0(*data);
  fillG0rtIdentity(*data);

  typename Data::DisorderConfiguration config;  // zero potential
  data->makeDisorderedG0(config);

  const int n_k = KDmn::dmn_size();
  const int n_r = RDmn::dmn_size();
  constexpr double tol = 1e-9;

  // Manual clean g(d) = (1/N_k) sum_k G0_k exp(-i k.r_d) (FunctionTransform K->R convention).
  for (int w = 0; w < WDmn::dmn_size(); ++w)
    for (int R1 = 0; R1 < n_r; ++R1)
      for (int R2 = 0; R2 < n_r; ++R2) {
        const int d = RDmn::parameter_type::subtract(R2, R1);  // index of r_{R1} - r_{R2}
        const auto& r_d = RDmn::get_elements()[d];
        for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1)
          for (int n2 = 0; n2 < NuDmn::dmn_size(); ++n2) {
            Complex g_r(0);
            for (int k = 0; k < n_k; ++k) {
              const auto& k_val = KDmn::get_elements()[k];
              g_r += data->G0_k_w_cluster_excluded(n1, n2, k, w) *
                     std::exp(Complex(0, -dca::math::util::innerProduct(k_val, r_d)));
            }
            g_r /= static_cast<double>(n_k);

            const Complex got = data->disordered_G0_r_r_w_cl_exl(n1, R1, n2, R2, w);
            EXPECT_NEAR(got.real(), g_r.real(), tol)
                << "mismatch at w=" << w << " (n1,R1,n2,R2)=(" << n1 << "," << R1 << "," << n2 << ","
                << R2 << ")";
            EXPECT_NEAR(got.imag(), g_r.imag(), tol);
          }
      }
}

// The unfolded clean G0 must be translationally invariant: every (R1,R2) pair sharing a
// displacement subtract(R2,R1) carries the same nu-block.
TEST(MakeDisorderedG0Test, UnfoldingIsTranslationallyInvariant) {
  auto data = makeData();
  fillCleanG0(*data);
  fillG0rtIdentity(*data);

  typename Data::DisorderConfiguration config;
  data->makeDisorderedG0(config);

  const int n_r = RDmn::dmn_size();
  constexpr double tol = 1e-10;
  const int w = 0;

  for (int R1 = 0; R1 < n_r; ++R1)
    for (int R2 = 0; R2 < n_r; ++R2) {
      const int d = RDmn::parameter_type::subtract(R2, R1);
      // Reference pair for this displacement: (R1'=d, R2'=0) also has subtract(0,d) == d.
      for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1)
        for (int n2 = 0; n2 < NuDmn::dmn_size(); ++n2) {
          const Complex a = data->disordered_G0_r_r_w_cl_exl(n1, R1, n2, R2, w);
          const Complex b = data->disordered_G0_r_r_w_cl_exl(n1, d, n2, 0, w);
          EXPECT_NEAR(a.real(), b.real(), tol);
          EXPECT_NEAR(a.imag(), b.imag(), tol);
        }
    }
}

// Disorder Dyson on the dense (nu*N_R) basis: invert the clean block, add a nonzero disorder block
// V, invert again. Verifies the algebraic identity inverse(G_dis) - inverse(G_clean) == V (so the
// two inversions compose correctly through the contiguous block layout) and that V actually changes
// the Green's function.
TEST(MakeDisorderedG0Test, DisorderBlockInversionIdentity) {
  auto data = makeData();
  fillCleanG0(*data);
  fillG0rtIdentity(*data);

  typename Data::DisorderConfiguration config;
  data->makeDisorderedG0(config);  // disordered_G0_r_r_w_cl_exl now holds the clean unfolded G0.

  const int nu_r = NuDmn::dmn_size() * RDmn::dmn_size();
  const int w = 0;

  // Clean block for frequency w (same pointer idiom as makeDisorderedG0).
  dca::linalg::Matrix<Complex, dca::linalg::CPU> clean(nu_r);
  dca::linalg::matrixop::copyArrayToMatrix(
      nu_r, nu_r, &data->disordered_G0_r_r_w_cl_exl(0, 0, 0, 0, w), nu_r, clean);

  // Nonzero, modest disorder block (general, not site-diagonal).
  dca::linalg::Matrix<Complex, dca::linalg::CPU> V(nu_r);
  for (int j = 0; j < nu_r; ++j)
    for (int i = 0; i < nu_r; ++i)
      V(i, j) = Complex(0.05 * (i + 1) - 0.02 * j, 0.01 * (i - j));

  dca::linalg::Matrix<Complex, dca::linalg::CPU> inv_clean(clean);
  dca::linalg::matrixop::inverse(inv_clean);

  dca::linalg::Matrix<Complex, dca::linalg::CPU> disordered(inv_clean);
  for (int j = 0; j < nu_r; ++j)
    for (int i = 0; i < nu_r; ++i)
      disordered(i, j) += V(i, j);
  dca::linalg::matrixop::inverse(disordered);

  dca::linalg::Matrix<Complex, dca::linalg::CPU> inv_dis(disordered);
  dca::linalg::matrixop::inverse(inv_dis);

  constexpr double tol = 1e-8;
  bool differs = false;
  for (int j = 0; j < nu_r; ++j)
    for (int i = 0; i < nu_r; ++i) {
      const Complex recovered = inv_dis(i, j) - inv_clean(i, j);
      EXPECT_NEAR(recovered.real(), V(i, j).real(), tol) << "(" << i << "," << j << ")";
      EXPECT_NEAR(recovered.imag(), V(i, j).imag(), tol) << "(" << i << "," << j << ")";
      if (std::abs(disordered(i, j) - clean(i, j)) > 1e-6)
        differs = true;
    }
  EXPECT_TRUE(differs) << "the disorder block must change the Green's function";
}
