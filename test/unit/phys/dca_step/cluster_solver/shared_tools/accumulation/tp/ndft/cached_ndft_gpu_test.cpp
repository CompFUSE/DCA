// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// Unit tests for the device version of the cached_ndft class.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft_gpu.hpp"

#include <complex>
#include <limits>

#include "gtest/gtest.h"

#include "dca/function/util/difference.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/reshapable_matrix.hpp"
#include "dca/linalg/util/util_cublas.hpp"
#include "dca/profiling/events/time.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/single_sector_accumulation_test.hpp"

constexpr int n_bands = 2;
constexpr int n_sites = 3;
constexpr int n_frqs = 7;

template <typename Real>
using CachedNdftGpuTest = dca::testing::SingleSectorAccumulationTest<Real, n_bands, n_sites, n_frqs>;

template <typename Real>
double computeWithFastNDFT(const typename CachedNdftGpuTest<Real>::Configuration& config,
                           const typename CachedNdftGpuTest<Real>::Matrix& M,
                           typename CachedNdftGpuTest<Real>::F_w_w& f_w);

using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_CASE(CachedNdftGpuTest, TestTypes);

// Compare the result provided by the GPU version of CachedNdft::execute with the definition of the
// NDFT f(w1, w2) = \sum_{t1, t2} f(t1, t2) exp(i * t1 * w1 - t2 w2) stored in f_baseline_.
TYPED_TEST(CachedNdftGpuTest, Execute) {
  constexpr int n_samples = 31;
  TestFixture::prepareConfiguration(TestFixture::configuration_, TestFixture::M_, n_samples);

  using Real = TypeParam;
  typename TestFixture::F_w_w f_w_fast("f_w_fast");

  // Compute the NDFT with the CachedNdft class and rearrange the result with the same order as
  // f_baseline_.
  const double time =
      computeWithFastNDFT<Real>(TestFixture::configuration_, TestFixture::M_, f_w_fast);

  auto f_baseline = TestFixture::compute2DFTBaseline();
  const auto err = dca::func::util::difference(f_baseline, f_w_fast);
  EXPECT_LT(err.l_inf, 100 * std::numeric_limits<Real>::epsilon());

  std::cout << "\nCached GPU ndft time [sec]:\t " << time << "\n";
}

template <typename Real>
double computeWithFastNDFT(const typename CachedNdftGpuTest<Real>::Configuration& config,
                           const typename CachedNdftGpuTest<Real>::Matrix& M,
                           typename CachedNdftGpuTest<Real>::F_w_w& f_w) {
  using Complex = std::complex<Real>;
  dca::linalg::util::initializeMagma();
  magma_queue_t queue;
  magma_queue_create(&queue);

  using BDmn = typename CachedNdftGpuTest<Real>::BDmn;
  using RDmn = typename CachedNdftGpuTest<Real>::RDmn;
  using FreqDmn = typename CachedNdftGpuTest<Real>::FreqDmn;
  using PosFreqDmn = typename CachedNdftGpuTest<Real>::PosFreqDmn;

  using NftType =
      dca::phys::solver::accumulator::CachedNdft<Real, RDmn, FreqDmn, PosFreqDmn, dca::linalg::GPU>;
  NftType nft_obj(queue);
  EXPECT_EQ(magma_queue_get_cuda_stream(queue), nft_obj.get_stream());

  dca::linalg::Matrix<double, dca::linalg::GPU> M_dev(M);
  typename NftType::RMatrix result_device(64);

  dca::profiling::WallTime start_time;
  nft_obj.execute(config, M_dev, result_device);
  cudaStreamSynchronize(nft_obj.get_stream());
  dca::profiling::WallTime end_time;

  dca::linalg::ReshapableMatrix<Complex, dca::linalg::CPU> result_host(result_device);

  // Rearrange the output from a function of (r1, b1, w1, r2, b2, w2) to a function of (b1, b2, r1,
  // r2, w1, w2).
  const int nb = BDmn::dmn_size();
  const int nr = RDmn::dmn_size();
  const int n_w = PosFreqDmn::dmn_size();
  auto invert_w = [=](const int w) { return 2 * n_w - 1 - w; };
  for (int b2 = 0; b2 < nb; ++b2)
    for (int b1 = 0; b1 < nb; ++b1)
      for (int r2 = 0; r2 < nr; ++r2)
        for (int r1 = 0; r1 < nr; ++r1)
          for (int w2 = 0; w2 < FreqDmn::dmn_size(); ++w2)
            for (int w1 = 0; w1 < n_w; ++w1) {
              const auto val = result_host(r1 + b1 * nr + w1 * nr * nb, r2 + b2 * nr + w2 * nr * nb);
              f_w(b1, b2, r1, r2, w1 + n_w, w2) = val;
              f_w(b1, b2, r1, r2, invert_w(w1 + n_w), invert_w(w2)) = std::conj(val);
            }

  magma_queue_destroy(queue);

  dca::profiling::Duration duration(end_time, start_time);
  return duration.sec + 1.e-6 * duration.usec;
}
