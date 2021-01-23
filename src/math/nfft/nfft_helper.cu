// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Helper class for computing interpolation indices during the single particle accumulation.

#include "dca/math/nfft/nfft_helper.cuh"

#include <mutex>

namespace dca {
namespace math {
namespace nfft {
namespace details {
// dca::math::nfft::details::

__device__ __constant__ NfftHelper nfft_helper;

void NfftHelper::set(int nb, const int nc, const int *add_r, int lda, const int *sub_r, int lds,
                     Scalar t0, Scalar delta_t, const Scalar t0_window, const Scalar delta_t_window,
                     const Scalar beta) {
  static std::once_flag flag;
  std::call_once(flag, [=]() {
    dca::phys::solver::details::ClusterHelper::set(nc, add_r, lda, sub_r, lds, false);

    NfftHelper host_helper;

    host_helper.nb_ = nb;
    host_helper.tau_scale_ = 0.5 / beta;
    host_helper.t0_ = t0;
    host_helper.delta_t_ = delta_t;
    host_helper.t0_window_ = t0_window;
    host_helper.delta_t_window_ = delta_t_window;
    host_helper.one_div_delta_ = 1. / delta_t;
    host_helper.one_div_delta_t_window_ = 1. / delta_t_window;

    cudaMemcpyToSymbol(nfft_helper, &host_helper, sizeof(NfftHelper));
  });
}

}  // namespace details
}  // namespace nfft
}  // namespace math
}  // namespace dca
