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

void initializeNfftHelper(int nb, int nc, const int *add_r, int lda, const int *sub_r, int lds,
                          double t0, double delta_t, double t0_window, double delta_t_window,
                          double beta);

}  // details
}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_KERNELS_INTERFACE_HPP
