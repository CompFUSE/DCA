// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file provides a common setup to the CachedNdft tests.

#ifndef TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_TEST_HPP
#define TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_TEST_HPP

#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

namespace dca {
namespace testing {

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

template <int n_samples, int n_bands, int n_frqs>
class CachedNdftTest : public ::testing::Test {
public:
  using RDmn = dca::func::dmn_0<dca::func::dmn<2, int>>;
  using FreqDmn = dca::func::dmn_0<dca::phys::domains::frequency_domain>;
  using PosFreqDmn = dca::func::dmn_0<PositiveFrq>;
  using BDmn = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;

  using Configuration = std::vector<Vertex>;
  using Matrix = dca::linalg::Matrix<double, dca::linalg::CPU>;
  using F_w_w =
      dca::func::function<std::complex<double>,
                          dca::func::dmn_variadic<BDmn, BDmn, RDmn, RDmn, FreqDmn, FreqDmn>>;

protected:
  static void SetUpTestCase() {
    // Initialize time and frequency domains.
    dca::phys::domains::frequency_domain::initialize(0.5, n_frqs);
    PositiveFrq::initialize(n_frqs);
    // Initialize band domain.
    int mock_parameter = 0;
    BDmn::parameter_type::initialize(
        mock_parameter, n_bands, std::vector<int>(),
        std::vector<std::vector<double>>(n_bands, std::vector<double>(n_bands, 0)));
  }

  void SetUp() {
    prepareConfiguration(configuration_, M_, n_samples);
    computeWithDft(configuration_, M_, f_baseline_);
  }

protected:
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

  void computeWithDft(const Configuration& config, const Matrix& M, F_w_w& f_w) {
    const std::complex<double> imag(0, 1);

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
                f_t * std::exp(imag * (t_val1 * w_val1 - t_val2 * w_val2));
          }
        }
      }
    }
  }

protected:
  Configuration configuration_;
  Matrix M_;
  F_w_w f_baseline_;
};

}  // testing
}  // dca

#endif  // TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_TEST_HPP
