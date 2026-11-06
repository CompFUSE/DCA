// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Anirudha Mirmira (anirudha.mirmira@gmail.com)
//
// CPU test for the disordered (DISORDERED_G0) single-particle accumulator. It verifies that the
// two-site M function M(b1,s,r1,b2,s,r2,w) is filled with genuine two-site structure:
//   * the value M_ij from a configuration pair (i,j) lands in bin (r_i, r_j) -- ALL entries,
//     including the off-diagonal r1 != r2 ones, not only the diagonal;
//   * bins with no contributing pair stay exactly zero (correct, leak-free indexing);
//   * the two spin sectors are written into their own s slot.
//
// The test only exercises real code under DISORDERED_G0 (the macro that switches the accumulator's
// PDmn to <B,R,B,R> and MFunction to <Nu,R,Nu,R,W>). When the macro is off it degenerates to a
// trivial passing test so the target still builds in a non-disorder configuration.

#include <array>
#include <cmath>
#include <complex>
#include <vector>

#include "dca/testing/gtest_h_w_warning_blocking.h"

using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"

constexpr int n_bands = 1;
constexpr int n_sites = 4;
constexpr int n_frqs = 32;

using SpAccumulatorDisorderTest =
    dca::testing::AccumulationTest<Scalar, n_bands, n_sites, n_frqs, false>;

using Parameters = SpAccumulatorDisorderTest::Parameters;
using Configuration = SpAccumulatorDisorderTest::Configuration;  // std::array<vector<Vertex>, 2>
using MatrixPair = SpAccumulatorDisorderTest::Sample;            // std::array<Matrix, 2>
using Vertex = dca::testing::Vertex;

using HostAccumulator = dca::phys::solver::accumulator::SpAccumulator<Parameters, dca::linalg::CPU>;
using MFunction = typename HostAccumulator::MFunction;
using RDmn = typename HostAccumulator::RDmn;
using BDmn = typename HostAccumulator::BDmn;
using WDmn = typename HostAccumulator::WDmn;

namespace {

// Largest |M(b=0, s, r1, b=0, s, r2, w)| over all Matsubara frequencies for a given bin/spin.
double maxAbsOverW(const MFunction& f, const int s, const int r1, const int r2) {
  double m = 0.;
  for (int w = 0; w < WDmn::dmn_size(); ++w)
    m = std::max(m, std::abs(f(0, s, r1, 0, s, r2, w)));
  return m;
}

}  // namespace

// All sites carry band 0 (n_bands == 1); each spin's vertices sit on distinct sites with distinct
// imaginary times, so every config pair (i,j) maps to its own (r_i, r_j) bin with no cancellation.
// We then assert exactly which bins are populated.
TEST_F(SpAccumulatorDisorderTest, OffDiagonalIndexPlacement) {
  // Spin up lives on sites {0,1,2}; spin down on {2,3}. Site 3 is never touched by spin up, site 0
  // and 1 never by spin down -- this exposes both off-diagonal placement and leak-free indexing.
  const std::array<std::vector<int>, 2> sites{std::vector<int>{0, 1, 2}, std::vector<int>{2, 3}};

  Configuration config;
  MatrixPair M;
  for (int s = 0; s < 2; ++s) {
    const int n = sites[s].size();
    config[s].resize(n);
    M[s].resize(n);
    for (int k = 0; k < n; ++k)
      config[s][k] = Vertex{0, sites[s][k], 0.1 + 0.3 * k + 0.05 * s};  // distinct, in [0, beta]
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
        M[s](i, j) = 0.5 + i + 2. * j + 0.1 * s;  // never zero
  }

  HostAccumulator acc(parameters_);
  acc.resetAccumulation();
  acc.accumulate(M, config, 1.);  // factor / sign = 1
  acc.finalize();
  const MFunction& Mrw = acc.get_sign_times_M_r_w();

  constexpr double zero_tol = 1e-14;
  int off_diag_nonzero = 0;

  for (int s = 0; s < 2; ++s) {
    const auto& used = sites[s];
    auto is_used = [&](int r) { return std::find(used.begin(), used.end(), r) != used.end(); };

    for (int r1 = 0; r1 < RDmn::dmn_size(); ++r1)
      for (int r2 = 0; r2 < RDmn::dmn_size(); ++r2) {
        const double mag = maxAbsOverW(Mrw, s, r1, r2);
        if (is_used(r1) && is_used(r2)) {
          EXPECT_GT(mag, zero_tol) << "expected nonzero bin s=" << s << " (" << r1 << "," << r2 << ")";
          if (r1 != r2)
            ++off_diag_nonzero;
        }
        else {
          EXPECT_EQ(0.0, mag) << "expected exactly-zero bin s=" << s << " (" << r1 << "," << r2 << ")";
        }
      }
  }

  // Spin up contributes 3x3 - 3 = 6 off-diagonal bins, spin down 2x2 - 2 = 2 -> 8 total.
  EXPECT_EQ(8, off_diag_nonzero) << "off-diagonal (r1 != r2) entries must be populated, not skipped";
}

// Exact, NFFT-independent check that the VALUE M(i,j) is placed in bin (r_i, r_j) for the full
// matrix (every off-diagonal included). With all imaginary times equal, every bin shares one
// frequency profile g(w), so M_rw(bin_ij, w) * M(0,0) == M_rw(bin_00, w) * M(i,j) for all w; g(w)
// (whatever the NFFT produces) cancels, leaving an exact proportionality of the placed values.
TEST_F(SpAccumulatorDisorderTest, ValuePlacementFullMatrix) {
  const std::vector<int> sites{0, 1, 2};  // distinct -> one pair per bin
  const int n = sites.size();

  Configuration config;
  MatrixPair M;
  for (int s = 0; s < 2; ++s) {
    config[s].resize(n);
    M[s].resize(n);
    for (int k = 0; k < n; ++k)
      config[s][k] = Vertex{0, sites[k], 0.0};  // identical tau -> shared frequency profile
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
        M[s](i, j) = 1.0 + i + 3. * j;  // real, nonzero, distinct per (i,j)
  }

  HostAccumulator acc(parameters_);
  acc.resetAccumulation();
  acc.accumulate(M, config, 1.);
  acc.finalize();
  const MFunction& Mrw = acc.get_sign_times_M_r_w();

  const int s = 0;
  const int r0 = sites[0];
  const double w00 = M[s](0, 0);  // reference weight, nonzero

  // bin_00 must be nonzero so the proportionality reference is meaningful.
  ASSERT_GT(maxAbsOverW(Mrw, s, r0, r0), 1e-14);

  constexpr double tol = 1e-10;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
      const int ri = sites[i];  // right vertex i -> first R index
      const int rj = sites[j];  // left  vertex j -> second R index
      const double wij = M[s](i, j);
      for (int w = 0; w < WDmn::dmn_size(); ++w) {
        const std::complex<double> lhs = Mrw(0, s, ri, 0, s, rj, w) * w00;
        const std::complex<double> rhs = Mrw(0, s, r0, 0, s, r0, w) * wij;
        EXPECT_NEAR(lhs.real(), rhs.real(), tol)
            << "value misplaced for pair (" << i << "," << j << ") bin (" << ri << "," << rj
            << ") at w=" << w;
        EXPECT_NEAR(lhs.imag(), rhs.imag(), tol);
      }
    }
}
