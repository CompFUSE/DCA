// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file provides the interface to the kernels for the Dnfft1DGpu accumulator

#ifndef DCA_MATH_NFFT_KERNELS_INTERFACE_HPP
#define DCA_MATH_NFFT_KERNELS_INTERFACE_HPP

#ifndef DCA_HAVE_CUDA
#pragma error "This file requires CUDA."
#endif

#include <cuda.h>

namespace dca {
namespace math {
namespace nfft {
namespace details {
// dca::math::nfft::details::

struct ConfigElem {
  int band;
  int site;
};

template <typename ScalarType>
void accumulateOnDevice(const ScalarType* M, int ldm, ScalarType sign, ScalarType* out,
                        ScalarType* out_sqr, const int ldo, const ConfigElem* config_left,
                        const ConfigElem* config_right, const ScalarType* tau,
                        const ScalarType* cubic_coeff, int size, cudaStream_t stream_);

template <typename ScalarType>
void sum(const ScalarType* in, int ldi, ScalarType* out, int ldo, int n, int m, cudaStream_t stream);

template <typename ScalarType>
void initializeNfftHelper(int nb, int nr, const int* sub_r, int lds, int oversampling,
                          int window_sampling, ScalarType t0, ScalarType delta_t,
                          ScalarType t0_window, ScalarType delta_t_window, ScalarType beta);

}  // details
}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_KERNELS_INTERFACE_HPP
