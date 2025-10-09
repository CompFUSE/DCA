// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file provides the interface to the kernels for the Dnfft1DGpu accumulator

#ifndef DCA_MATH_NFFT_KERNELS_INTERFACE_HPP
#define DCA_MATH_NFFT_KERNELS_INTERFACE_HPP

#include "dca/config/haves_defines.hpp"
#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu.h"
#else
#pragma error "This file requires GPU."
#endif


namespace dca {
namespace math {
namespace nfft {
namespace details {
// dca::math::nfft::details::

struct ConfigElem {
  int band;
  int site;
};

template <int oversampling, int window_sampling, typename Scalar, typename Real>
void accumulateOnDevice(const Scalar* M, const int ldm, const Scalar factor, Scalar* out,
                        Scalar* out_sqr, const int ldo, const ConfigElem* config_left,
                        const ConfigElem* config_right, const Real* tau,
                        const Real* cubic_coeff, const int size, cudaStream_t stream_);

template <typename ScalarType>
void sum(const ScalarType* in, int ldi, ScalarType* out, int ldo, int n, int m, cudaStream_t stream);

template<typename REAL>
void initializeNfftHelper(int nb, int nc, const int* add_r, int lda, const int* sub_r, int lds,
                          REAL t0, REAL delta_t, REAL t0_window, REAL delta_t_window,
                          REAL beta);

extern template void initializeNfftHelper<double>(int nb, int nc, const int* add_r, int lda, const int* sub_r, int lds,
                          double t0, double delta_t, double t0_window, double delta_t_window,
                          double beta);
// extern template void initializeNfftHelper<float>(int nb, int nc, const int* add_r, int lda, const int* sub_r, int lds,
//                           float t0, float delta_t, float t0_window, float delta_t_window,
//                           float beta);
  
}  // namespace details
}  // namespace nfft
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_NFFT_KERNELS_INTERFACE_HPP
