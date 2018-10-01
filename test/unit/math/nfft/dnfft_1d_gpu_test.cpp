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
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/single_sector_accumulation_test.hpp"

using dca::func::function;
using dca::func::dmn_variadic;

constexpr int n_bands = 3;
constexpr int n_sites = 5;
constexpr int n_frequencies = 64;
using Dnfft1DGpuTest =
    dca::testing::SingleSectorAccumulationTest<double, n_bands, n_sites, n_frequencies>;

using FreqDmn = typename Dnfft1DGpuTest::FreqDmn;
using BDmn = typename Dnfft1DGpuTest::BDmn;
using RDmn = typename Dnfft1DGpuTest::RDmn;
using LabelDmn = dmn_variadic<BDmn, BDmn, RDmn>;
using Configuration = typename Dnfft1DGpuTest::Configuration;

template <typename DnfftType>
void computeWithCpuDnfft(dca::linalg::Matrix<double, dca::linalg::CPU>& M, Configuration& config,
                         DnfftType& dnfft_obj,
                         function<std::complex<double>, dmn_variadic<FreqDmn, LabelDmn>>& f_w);

TEST_F(Dnfft1DGpuTest, Accumulate) {
  prepareConfiguration(configuration_, M_, 128);
  dca::linalg::Matrix<double, dca::linalg::GPU> M_dev(M_);

  constexpr int oversampling = 8;
  // Compute f(w) using the delayed-NFFT algorithm on the CPU.
  dca::math::nfft::Dnfft1D<double, FreqDmn, LabelDmn, oversampling, dca::math::nfft::CUBIC> cpu_dnfft_obj;
  function<std::complex<double>, dmn_variadic<FreqDmn, LabelDmn>> f_w_dnfft_cpu("f_w_dnfft_cpu");
  computeWithCpuDnfft(M_, configuration_, cpu_dnfft_obj, f_w_dnfft_cpu);

  // Compute f(w) using the delayed-NFFT algorithm on the GPU.
  cudaStream_t stream;
  cudaStreamCreate(&stream);
  dca::math::nfft::Dnfft1DGpu<double, FreqDmn, RDmn, oversampling, dca::math::nfft::CUBIC> gpu_dnfft_obj(
      beta_, stream);
  function<std::complex<double>, dmn_variadic<FreqDmn, LabelDmn>> f_w_dnfft_gpu("f_w_dnfft_gpu");

  gpu_dnfft_obj.resetAccumulation();
  gpu_dnfft_obj.accumulate(M_dev, configuration_, 1);
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
  const double beta = Dnfft1DGpuTest::get_beta();
  dnfft_obj.resetAccumulation();
  const static LabelDmn bbr_dmn;
  const int n = config.size();
  const double scale = 1. / (2. * beta);
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
