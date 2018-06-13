// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the GPU implementation of the laset function.

#ifndef DCA_LINALG_LAPACK_LASET_GPU_HPP
#define DCA_LINALG_LAPACK_LASET_GPU_HPP

#include <complex>
#include <cuComplex.h>
#include "dca/linalg/util/cast_cuda.hpp"

namespace dca {
namespace linalg {
namespace lapack {
// dca::linalg::lapack::

// Sets the diagonal elements of the matrix to diag, and the off diagonal elements to offdiag.
// Preconditions: lda >= m.
// Type can be float, double, cuComplex, cuDoubleComplex, std::complex<float>, std::complex<double>.
template <typename Type>
void laset_gpu(int m, int n, Type offdiag, Type diag, Type* a, int lda, int thread_id, int stream_id);
template <typename Type>
inline void laset_gpu(int m, int n, std::complex<Type> offdiag, std::complex<Type> diag,
                      std::complex<Type>* a, int lda, int thread_id, int stream_id) {
  auto cu_offdiag = linalg::util::castCudaComplex(offdiag);
  auto cu_diag = linalg::util::castCudaComplex(&diag);
  auto cu_a = linalg::util::castCudaComplex(a);
  laset_gpu(m, n, *cu_offdiag, *cu_diag, cu_a, lda, thread_id, stream_id);
}

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_LASET_GPU_HPP
