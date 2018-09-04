// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the Dnfft1DGpu class by comparing its accumulation results with the CPU version.

#include "dca/math/nfft/dnfft_1d_gpu.hpp"

#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

using dca::func::function;
using dca::func::dmn_variadic;

class MockRDmn {
public:
  using element_type = short;

  static void initialize(const int n) {
    n_sites_ = n;
    auto& sub_matrix = get_matrix();
    sub_matrix.resizeNoCopy(n);
    for (int j = 0; j < n_sites_; ++j)
      for (int i = 0; i < n_sites_; ++i)
        sub_matrix(i, j) = (j - i + n_sites_) % n_sites_;
  }
  static const auto& get_subtract_matrix() {
    return get_matrix();
  }
  static int get_size() {
    return n_sites_;
  }
  static int subtract(const int i, const int j) {
    return get_matrix()(i, j);
  }

private:
  static int n_sites_;

  static inline dca::linalg::Matrix<int, dca::linalg::CPU>& get_matrix() {
    static dca::linalg::Matrix<int, dca::linalg::CPU> sub_matrix;
    return sub_matrix;
  }
};
int MockRDmn::n_sites_ = -1;

struct ConfigElement {
  double get_tau() const {
    return tau_;
  }
  int get_left_band() const {
    return band_;
  }
  int get_right_band() const {
    return band_;
  }
  int get_left_site() const {
    return r_;
  }
  int get_right_site() const {
    return r_;
  }

  int band_;
  int r_;
  double tau_;
};

using FreqDmn = dca::func::dmn_0<dca::phys::domains::frequency_domain>;
using BDmn = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;
using RDmn = dca::func::dmn_0<MockRDmn>;
using LabelDmn = dmn_variadic<BDmn, BDmn, RDmn>;
using Configuration = std::vector<ConfigElement>;

template <typename DnfftType>
void computeWithCpuDnfft(dca::linalg::Matrix<double, dca::linalg::CPU>& M, Configuration& config,
                         DnfftType& dnfft_obj,
                         function<std::complex<double>, dmn_variadic<FreqDmn, LabelDmn>>& f_w);
void prepareConfiguration(dca::linalg::Matrix<double, dca::linalg::CPU>& M, Configuration& config,
                          int n);
constexpr double beta = 10.;

TEST(Dnfft1DGpuTest, Accumulate) {
  // Initialize time and frequency domains.
  const int positive_frequencies = 25;
  const int n_bands = 2;
  const int n_sites = 3;

  dca::phys::domains::frequency_domain::initialize(beta, positive_frequencies);
  int mock_par(0);
  dca::phys::domains::electron_band_domain::initialize(
      mock_par, n_bands, std::vector<int>(),
      std::vector<std::vector<double>>(n_bands, std::vector<double>(n_bands, 0)));
  MockRDmn::initialize(n_sites);

  // Prepare random samples.
  const int samples = 62;
  dca::linalg::Matrix<double, dca::linalg::CPU> M;
  Configuration config;
  prepareConfiguration(M, config, samples);

  // Compute f(w) using the delayed-NFFT algorithm on the CPU.
  constexpr int oversampling = 8;
  dca::math::nfft::Dnfft1D<double, FreqDmn, LabelDmn, oversampling, dca::math::nfft::CUBIC> cpu_dnfft_obj;
  function<std::complex<double>, dmn_variadic<FreqDmn, LabelDmn>> f_w_dnfft_cpu("f_w_dnfft_cpu");
  computeWithCpuDnfft(M, config, cpu_dnfft_obj, f_w_dnfft_cpu);

  // Compute f(w) using the delayed-NFFT algorithm on the GPU.
  cudaStream_t stream;
  cudaStreamCreate(&stream);
  dca::math::nfft::Dnfft1DGpu<double, FreqDmn, RDmn, oversampling, dca::math::nfft::CUBIC> gpu_dnfft_obj(
      beta, stream);
  function<std::complex<double>, dmn_variadic<FreqDmn, LabelDmn>> f_w_dnfft_gpu("f_w_dnfft_gpu");
  gpu_dnfft_obj.accumulate(M, config, 1);
  gpu_dnfft_obj.finalize(f_w_dnfft_gpu);
  cudaStreamDestroy(stream);

  // Check errors.
  const auto err = dca::func::util::difference(f_w_dnfft_cpu, f_w_dnfft_gpu);
  EXPECT_LT(err.l_inf, 1.e-9);
}

template <typename DnfftType>
void computeWithCpuDnfft(dca::linalg::Matrix<double, dca::linalg::CPU>& M, Configuration& config,
                         DnfftType& dnfft_obj,
                         function<std::complex<double>, dmn_variadic<FreqDmn, LabelDmn>>& f_w) {
  dnfft_obj.initialize();
  const static LabelDmn bbr_dmn;
  const int n = config.size();
  const double scale = 0.5 / beta;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
      const int delta_r =
          RDmn::parameter_type::subtract(config[j].get_left_site(), config[i].get_right_site());
      const int index = bbr_dmn(config[i].get_right_band(), config[j].get_left_band(), delta_r);
      const double delta_t = (config[i].get_tau() - config[j].get_tau()) * scale;
      dnfft_obj.accumulate(index, delta_t, M(i, j));
    }

  dnfft_obj.finalize(f_w);
}

void prepareConfiguration(dca::linalg::Matrix<double, dca::linalg::CPU>& M, Configuration& config,
                          const int n) {
  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 42);
  config.resize(n);
  M.resize(n);
  for (int i = 0; i < n; ++i) {
    const double tau = beta * rng();
    const int r = rng() * RDmn::dmn_size();
    const int b = rng() * BDmn::dmn_size();
    config[i] = ConfigElement{b, r, tau};
  }

  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      M(i, j) = 2 * rng() - 1.;
}
