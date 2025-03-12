// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// Helper class for computing interpolation indices during the single particle accumulation.

#include "dca/math/nfft/nfft_helper.cuh"

#include <mutex>

namespace dca {
namespace math {
namespace nfft {
namespace details {
// dca::math::nfft::details::

__device__ __constant__ NfftHelper<double> nfft_helper_double;
__device__ __constant__ NfftHelper<float> nfft_helper_float;

template void NfftHelper<double>::set(int nb, const int nc, const int* add_r, int lda,
                                      const int* sub_r, int lds, double t0, double delta_t,
                                      const double t0_window, const double delta_t_window,
                                      const double beta);
template void NfftHelper<float>::set(int nb, const int nc, const int* add_r, int lda,
                                     const int* sub_r, int lds, float t0, float delta_t,
                                     const float t0_window, const float delta_t_window,
                                     const float beta);

template <typename REAL>
void NfftHelper<REAL>::set(int nb, const int nc, const int* add_r, int lda, const int* sub_r,
                           int lds, REAL t0, REAL delta_t, const REAL t0_window,
                           const REAL delta_t_window, const REAL beta) {
  static std::once_flag flag;
  std::call_once(flag, [=]() {
    // Initialize real space cluster \todo this should be done somewhere definite instead of via call_once from multiple places
    // other is SolverHelper.
    dca::phys::solver::details::ClusterHelper::set(nc, add_r, lda, sub_r, lds);

    NfftHelper host_helper;

    host_helper.nb_ = nb;
    host_helper.tau_scale_ = 0.5 / beta;
    host_helper.t0_ = t0;
    host_helper.delta_t_ = delta_t;
    host_helper.t0_window_ = t0_window;
    host_helper.delta_t_window_ = delta_t_window;
    host_helper.one_div_delta_ = 1. / delta_t;
    host_helper.one_div_delta_t_window_ = 1. / delta_t_window;

    if constexpr (std::is_same_v<Real, double>)
      checkRC(cudaMemcpyToSymbol(HIP_SYMBOL(nfft_helper_double), &host_helper, sizeof(NfftHelper<double>)));
    else if constexpr (std::is_same_v<Real, float>)
      checkRC(cudaMemcpyToSymbol(HIP_SYMBOL(nfft_helper_float), &host_helper, sizeof(NfftHelper<float>)));
  });
}

}  // namespace details
}  // namespace nfft
}  // namespace math
}  // namespace dca
