// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// Integration tests for the cached_ndft class.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/function_utils.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/profiling/events/time.hpp"

using dca::func::function;
using dca::func::dmn_variadic;
class PositiveFrq {
public:
  using element_type = double;
  static std::size_t get_size() {
    return size_;
  }
  static void initialize(const int n_frq) {
    size_ = n_frq;
  }

private:
  static int size_;
};
int PositiveFrq::size_ = -1;

using RDmn = dca::func::dmn_0<dca::func::dmn<2, int>>;

using FreqDmn = dca::func::dmn_0<dca::phys::domains::frequency_domain>;
using FreqDmnPos = dca::func::dmn_0<PositiveFrq>;
using BDmn = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;

struct Vertex {
  double get_tau() const {
    return tau_;
  }
  double get_left_band() const {
    return b_;
  }
  double get_right_band() const {
    return b_;
  }
  double get_left_site() const {
    return r_;
  }
  double get_right_site() const {
    return r_;
  }

  int b_;
  int r_;
  double tau_;
};
using Configuration = std::vector<Vertex>;
using Matrix = dca::linalg::Matrix<double, dca::linalg::CPU>;

template <class W2Dmn>
double computeWithFastDNFT(
    const Configuration& config, const Matrix& M,
    function<std::complex<double>, dmn_variadic<BDmn, BDmn, RDmn, RDmn, FreqDmn, FreqDmn>>& f_w);
void computeWithDft(
    const Configuration& config, const Matrix& M,
    function<std::complex<double>, dmn_variadic<BDmn, BDmn, RDmn, RDmn, FreqDmn, FreqDmn>>& f_w);
void prepareConfiguration(Configuration& config, Matrix& M, int n);

namespace global {
// Global data, to be initialized once.
Configuration configuration;
Matrix M;
// Stores f(w) computed with the direct definition of the discrete Fourier transform (DFT).
std::unique_ptr<function<std::complex<double>, dmn_variadic<BDmn, BDmn, RDmn, RDmn, FreqDmn, FreqDmn>>> f_w_dft;
}  // global

TEST(CachedNdftTest, prepareConfigAndDomains) {
  const int positive_frequencies = 8;
  const int config_size = 40;

  const int n_bands = 2;

  // Initialize time and frequency domains.
  dca::phys::domains::frequency_domain::initialize(0.5, positive_frequencies);
  PositiveFrq::initialize(positive_frequencies);
  int mock_parameter = 0;

  BDmn::parameter_type::initialize(
      mock_parameter, n_bands, std::vector<int>(),
      std::vector<std::vector<double>>(n_bands, std::vector<double>(n_bands, 0)));

  prepareConfiguration(global::configuration, global::M, config_size);
  global::f_w_dft.reset(
      new function<std::complex<double>, dmn_variadic<BDmn, BDmn, RDmn, RDmn, FreqDmn, FreqDmn>>(
          "f_w_dft"));
  computeWithDft(global::configuration, global::M, *global::f_w_dft);
}

TEST(CachedNdftTest, execute) {
  function<std::complex<double>, dmn_variadic<BDmn, BDmn, RDmn, RDmn, FreqDmn, FreqDmn>> f_w_fast(
      "f_w_fast");

  const double time = computeWithFastDNFT<FreqDmnPos>(global::configuration, global::M, f_w_fast);

  // Check errors.
  const auto err = dca::func::utils::difference(*global::f_w_dft, f_w_fast);
  EXPECT_LT(err.l_inf, 1e-14);

  std::cout << "\nTrimmed cached nft time [sec]:\t " << time << "\n";
}

template <class W2Dmn>
double computeWithFastDNFT(
    const Configuration& config, const Matrix& M,
    function<std::complex<double>, dmn_variadic<BDmn, BDmn, RDmn, RDmn, FreqDmn, FreqDmn>>& f_w) {
  function<std::complex<double>, dmn_variadic<BDmn, BDmn, RDmn, RDmn, W2Dmn, FreqDmn>> f_b_b_r_r_w_w;
  dca::phys::solver::accumulator::CachedNdft<double, RDmn, FreqDmn, FreqDmnPos, dca::linalg::CPU> nft_obj;

  dca::profiling::WallTime start_time;
  nft_obj.execute(config, M, f_b_b_r_r_w_w);
  dca::profiling::WallTime end_time;

  const int n_w = PositiveFrq::get_size();
  auto invert_w = [=](const int w) { return 2 * n_w - 1 - w; };
  for (int b2 = 0; b2 < BDmn::dmn_size(); ++b2)
    for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1)
      for (int r2 = 0; r2 < RDmn::dmn_size(); ++r2)
        for (int r1 = 0; r1 < RDmn::dmn_size(); ++r1)
          for (int w2 = 0; w2 < FreqDmn::dmn_size(); ++w2)
            for (int w1 = 0; w1 < n_w; ++w1) {
              f_w(b1, b2, r1, r2, w1 + n_w, w2) = f_b_b_r_r_w_w(b1, b2, r1, r2, w1, w2);
              f_w(b1, b2, r1, r2, invert_w(w1 + n_w), invert_w(w2)) =
                  std::conj(f_b_b_r_r_w_w(b1, b2, r1, r2, w1, w2));
            }

  dca::profiling::Duration duration(end_time, start_time);
  return duration.sec + 1.e-6 * duration.usec;
}

void computeWithDft(
    const Configuration& config, const Matrix& M,
    function<std::complex<double>, dmn_variadic<BDmn, BDmn, RDmn, RDmn, FreqDmn, FreqDmn>>& f_w) {
  const std::complex<double> img(0, 1);

  for (int w_ind2 = 0; w_ind2 < FreqDmn::dmn_size(); ++w_ind2) {
    const double w_val2 = FreqDmn::get_elements()[w_ind2];
    for (int w_ind1 = 0; w_ind1 < FreqDmn::dmn_size(); ++w_ind1) {
      const double w_val1 = FreqDmn::get_elements()[w_ind1];
      for (int j = 0; j < config.size(); ++j) {
        const auto t_val2 = config[j].get_tau();
        const int b2 = config[j].b_;
        const int r2 = config[j].r_;
        for (int i = 0; i < config.size(); ++i) {
          const auto t_val1 = config[i].get_tau();
          const int b1 = config[i].b_;
          const int r1 = config[i].r_;

          const auto f_t = M(i, j);
          f_w(b1, b2, r1, r2, w_ind1, w_ind2) +=
              f_t * std::exp(img * (t_val1 * w_val1 - t_val2 * w_val2));
        }
      }
    }
  }
}

void prepareConfiguration(Configuration& config, Matrix& M, const int n) {
  config.resize(n);
  M.resize(n);
  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 0);

  for (int i = 0; i < n; ++i) {
    const double tau = rng() - 0.5;
    const int b = rng() * BDmn::dmn_size();
    const int r = rng() * RDmn::dmn_size();
    config[i] = Vertex{b, r, tau};
  }

  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
      const double t1 = config[i].get_tau();
      const double t2 = config[j].get_tau();
      M(i, j) = std::sin(2 * M_PI * t1) * std::sin(6 * M_PI * t2);
    }
}
