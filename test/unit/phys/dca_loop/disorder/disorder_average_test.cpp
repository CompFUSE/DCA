// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Anirudha Mirmira (anirudha.mirmira@gmail.com)
//
// Unit tests for the DISORDERED_G0 numeric kernels in
// include/dca/phys/dca_loop/disorder/disorder_average.hpp:
//
//   1. translationalAverage  - orientation / round-trip: unfolding a single-displacement g0 into
//      the two space-index quantity (using subtract(R2,R1)) and then projecting it back (using
//      subtract(R,R0)) must recover g0 exactly. Tested on three clusters: single site (Nc=1), a
//      2x2 cluster (Nc=4), and an irregular 3x4 cluster (Nc=12). The first two are self-inverse
//      (-r == r for every site) and CANNOT distinguish subtract(R,R0) from its negation; the 3x4
//      cluster has sites with -r != r, so it is the one that actually catches an add-vs-subtract
//      (orientation) flip in the projection.
//
//   3. accumulateDisorderDyson - the per-frequency dense (nu*N_R) Dyson
//      G = G0_dis - G0_dis*M*G0_dis/beta, accumulated weighted. Checked against an independent
//      naive matrix computation, and that repeated calls accumulate (+=) rather than overwrite.
//
// The tests use lightweight mock domains (a band/spin/frequency placeholder and a rectangular
// cluster mock that carries its own subtract table) so all three cluster shapes coexist in one
// executable - the real cluster_domain is a global singleton initialized once per binary.

#include "dca/phys/dca_loop/disorder/disorder_average.hpp"

#include <complex>
#include <string>
#include <vector>

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"

namespace {

using dca::func::dmn_0;
using dca::func::dmn_variadic;
using Complex = std::complex<double>;

// Minimal fixed-size leaf domain (band/spin, frequency placeholders). Satisfies the dmn_0
// parameter_type interface: get_size / element_type + get_elements / get_name.
template <int N>
struct FixedDmn {
  using element_type = int;
  static int get_size() {
    return N;
  }
  static int dmn_size() {
    return N;
  }
  static std::vector<int>& get_elements() {
    static std::vector<int> e = [] {
      std::vector<int> v(N);
      for (int i = 0; i < N; ++i)
        v[i] = i;
      return v;
    }();
    return e;
  }
  static std::string get_name() {
    return "FixedDmn";
  }
};

// Rectangular LX x LY cluster mock with periodic wrap. Site index = x + LX*y. Carries the cluster
// subtract table with the production convention subtract(i,j) = index of r_j - r_i. For LX or LY
// odd-or-not, sites with -r != r exist iff the cluster is not self-inverse (e.g. 3x4).
template <int LX, int LY>
struct RectCluster {
  using element_type = int;
  static int get_size() {
    return LX * LY;
  }
  static int dmn_size() {
    return LX * LY;
  }
  static std::vector<int>& get_elements() {
    static std::vector<int> e = [] {
      std::vector<int> v(LX * LY);
      for (int i = 0; i < LX * LY; ++i)
        v[i] = i;
      return v;
    }();
    return e;
  }
  static std::string get_name() {
    return "RectCluster";
  }
  // index of r_j - r_i, periodic in the LX x LY cell.
  static int subtract(int i, int j) {
    const int xi = i % LX, yi = i / LX;
    const int xj = j % LX, yj = j / LX;
    const int dx = ((xj - xi) % LX + LX) % LX;
    const int dy = ((yj - yi) % LY + LY) % LY;
    return dx + LX * dy;
  }
  // index of -r_i (used only by the test to reason about self-inverseness).
  static int negate(int i) {
    return subtract(i, 0);  // r_0 - r_i = -r_i
  }
};

using NuMock = dmn_0<FixedDmn<2>>;  // 1 band x 2 spin placeholder
using WMock = dmn_0<FixedDmn<3>>;   // 3 frequencies

// A distinct, displacement-sensitive value for the reference single-displacement g0.
Complex g0Value(int n1, int n2, int d, int w) {
  return Complex(1.0 + d + 0.1 * n1 + 0.01 * n2 + 0.001 * w, 0.5 * d - 0.2 * w + 0.03 * n1);
}

// Builds g0, unfolds it to the two space-index quantity via subtract(R2,R1) (as makeDisorderedG0
// does), runs translationalAverage, and asserts the projection recovers g0 exactly.
template <int LX, int LY>
void runTranslationalAverageRoundTrip() {
  using RDmn = dmn_0<RectCluster<LX, LY>>;
  using TwoR = dca::func::function<Complex, dmn_variadic<NuMock, RDmn, NuMock, RDmn, WMock>>;
  using OneR = dca::func::function<Complex, dmn_variadic<NuMock, NuMock, RDmn, WMock>>;

  const int N = RDmn::dmn_size();
  const int nu = NuMock::dmn_size();
  const int nw = WMock::dmn_size();

  OneR g0("g0");
  for (int w = 0; w < nw; ++w)
    for (int d = 0; d < N; ++d)
      for (int n2 = 0; n2 < nu; ++n2)
        for (int n1 = 0; n1 < nu; ++n1)
          g0(n1, n2, d, w) = g0Value(n1, n2, d, w);

  // Unfold: GG(n1,R1,n2,R2,w) = g0(n1,n2, subtract(R2,R1), w).
  TwoR gg("gg");
  for (int w = 0; w < nw; ++w)
    for (int R1 = 0; R1 < N; ++R1)
      for (int R2 = 0; R2 < N; ++R2) {
        const int d = RDmn::parameter_type::subtract(R2, R1);
        for (int n2 = 0; n2 < nu; ++n2)
          for (int n1 = 0; n1 < nu; ++n1)
            gg(n1, R1, n2, R2, w) = g0(n1, n2, d, w);
      }

  OneR projected("projected");
  dca::phys::disorder::translationalAverage<NuMock, RDmn, WMock>(gg, projected);

  constexpr double tol = 1e-12;
  for (int w = 0; w < nw; ++w)
    for (int R = 0; R < N; ++R)
      for (int n2 = 0; n2 < nu; ++n2)
        for (int n1 = 0; n1 < nu; ++n1) {
          const Complex got = projected(n1, n2, R, w);
          const Complex ref = g0(n1, n2, R, w);
          EXPECT_NEAR(got.real(), ref.real(), tol)
              << "Nc=" << N << " (n1,n2,R,w)=(" << n1 << "," << n2 << "," << R << "," << w << ")";
          EXPECT_NEAR(got.imag(), ref.imag(), tol);
        }
}

}  // namespace

TEST(DisorderAverageTest, TranslationalAverageRoundTripSingleSite) {
  runTranslationalAverageRoundTrip<1, 1>();  // Nc=1
}

TEST(DisorderAverageTest, TranslationalAverageRoundTripFourSite) {
  runTranslationalAverageRoundTrip<2, 2>();  // Nc=4, self-inverse
}

TEST(DisorderAverageTest, TranslationalAverageRoundTripIrregular3x4) {
  // Nc=12: this cluster has sites with -r != r, so it (unlike the 1- and 4-site clusters) fails if
  // the projection uses the wrong (negated) displacement orientation.
  using RDmn = RectCluster<3, 4>;
  bool has_non_self_inverse = false;
  for (int i = 0; i < RDmn::get_size(); ++i)
    if (RDmn::negate(i) != i)
      has_non_self_inverse = true;
  ASSERT_TRUE(has_non_self_inverse) << "the 3x4 cluster must not be self-inverse to be a teeth test";

  runTranslationalAverageRoundTrip<3, 4>();
}

// accumulateDisorderDyson on a small (Nc=2) cluster: compare to an independent naive computation of
// weight * (G0 - G0*M*G0/beta) on the dense (nu*N_R) block, and verify accumulation.
TEST(DisorderAverageTest, DisorderDysonMatchesNaiveAndAccumulates) {
  using RDmn = dmn_0<RectCluster<2, 1>>;  // Nc=2
  using TwoR = dca::func::function<Complex, dmn_variadic<NuMock, RDmn, NuMock, RDmn, WMock>>;

  const int N = RDmn::dmn_size();
  const int nu = NuMock::dmn_size();
  const int nw = WMock::dmn_size();
  const int dim = nu * N;  // matrix dimension
  const double beta = 2.5;
  const double weight = 0.75;

  // Controlled, invertible-ish G0_dis and a modest M.
  TwoR g0("g0_dis"), m("m"), g_acc("g_acc");
  auto block_index = [nu](int n, int R) { return n + nu * R; };  // row/col index in the dense block
  for (int w = 0; w < nw; ++w)
    for (int R2 = 0; R2 < N; ++R2)
      for (int n2 = 0; n2 < nu; ++n2)
        for (int R1 = 0; R1 < N; ++R1)
          for (int n1 = 0; n1 < nu; ++n1) {
            const int i = block_index(n1, R1), j = block_index(n2, R2);
            g0(n1, R1, n2, R2, w) =
                (i == j) ? Complex(2.0 + 0.1 * i + 0.05 * w, 0.3) : Complex(0.2 - 0.01 * (i + j), 0.1);
            m(n1, R1, n2, R2, w) = Complex(0.4 - 0.02 * i + 0.03 * j, 0.05 * (i - j) + 0.01 * w);
          }

  g_acc = 0;
  dca::phys::disorder::accumulateDisorderDyson<NuMock, RDmn, WMock, double>(g0, m, beta, weight, g_acc);

  // Independent naive computation per frequency: result = weight*(G0 - G0*M*G0/beta).
  constexpr double tol = 1e-11;
  for (int w = 0; w < nw; ++w) {
    std::vector<Complex> G0(dim * dim), M(dim * dim);
    auto at = [dim](int i, int j) { return i + dim * j; };
    for (int R2 = 0; R2 < N; ++R2)
      for (int n2 = 0; n2 < nu; ++n2)
        for (int R1 = 0; R1 < N; ++R1)
          for (int n1 = 0; n1 < nu; ++n1) {
            const int i = block_index(n1, R1), j = block_index(n2, R2);
            G0[at(i, j)] = g0(n1, R1, n2, R2, w);
            M[at(i, j)] = m(n1, R1, n2, R2, w);
          }
    std::vector<Complex> G0M(dim * dim, 0), G0MG0(dim * dim, 0);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        Complex s = 0;
        for (int k = 0; k < dim; ++k)
          s += G0[at(i, k)] * M[at(k, j)];
        G0M[at(i, j)] = s;
      }
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        Complex s = 0;
        for (int k = 0; k < dim; ++k)
          s += G0M[at(i, k)] * G0[at(k, j)];
        G0MG0[at(i, j)] = s;
      }
    for (int R2 = 0; R2 < N; ++R2)
      for (int n2 = 0; n2 < nu; ++n2)
        for (int R1 = 0; R1 < N; ++R1)
          for (int n1 = 0; n1 < nu; ++n1) {
            const int i = block_index(n1, R1), j = block_index(n2, R2);
            const Complex ref = weight * (G0[at(i, j)] - G0MG0[at(i, j)] / beta);
            const Complex got = g_acc(n1, R1, n2, R2, w);
            EXPECT_NEAR(got.real(), ref.real(), tol) << "w=" << w << " (i,j)=(" << i << "," << j << ")";
            EXPECT_NEAR(got.imag(), ref.imag(), tol);
          }
  }

  // Accumulation: a second identical call must double the accumulator (+= not =).
  TwoR after_one("after_one");
  after_one = g_acc;
  dca::phys::disorder::accumulateDisorderDyson<NuMock, RDmn, WMock, double>(g0, m, beta, weight, g_acc);
  for (int idx = 0; idx < g_acc.size(); ++idx) {
    EXPECT_NEAR(g_acc(idx).real(), 2.0 * after_one(idx).real(), tol);
    EXPECT_NEAR(g_acc(idx).imag(), 2.0 * after_one(idx).imag(), tol);
  }
}
