// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a class to store and use the arguments of a batched gemm operation.

#ifndef DCA_LINALG_UTIL_MAGMA_BATCHED_GEMM_HPP
#define DCA_LINALG_UTIL_MAGMA_BATCHED_GEMM_HPP
#ifdef DCA_HAVE_CUDA

#include <cassert>
#include <vector>

#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/util/allocators/vectors_typedefs.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <typename ScalarType>
class MagmaBatchedGemm {
public:
  // Creates a plan for a batched gemm.
  MagmaBatchedGemm(magma_queue_t queue);
  // Creates a plan for a batched gemm and allocates the memory for the arguments of `size`
  // multiplications.
  MagmaBatchedGemm(int size, magma_queue_t queue);

  // Allocates the memory for the arguments of `size` multiplications.
  void reserve(int size);

  // This methods sets the argument for a single gemm operation.
  // See the documentation of the MAGMA library for a description of the arguments.
  void addGemm(const ScalarType* a, const ScalarType* b, ScalarType* c);

  // Transfers the arguments to the device and performs the batched gemm of matrices of equal size.
  // See the documentation of the MAGMA library for a description of the arguments.
  void execute(char transa, char transb, int m, int n, int k, ScalarType alpha, ScalarType beta,
               int lda, int ldb, int ldc);

  // Synchronizes the copy of the arguments to the device and clears the arguments of the batched
  // gemm before reusing the object for a new batched operation.
  void synchronizeCopy();

private:
  magma_queue_t queue_;
  const cudaStream_t stream_;
  CudaEvent copied_;

  linalg::util::HostVector<const ScalarType *> a_ptr_, b_ptr_;
  linalg::util::HostVector<ScalarType*> c_ptr_;

  linalg::Vector<const ScalarType *, linalg::GPU> a_ptr_dev_, b_ptr_dev_;
  linalg::Vector<ScalarType*, linalg::GPU> c_ptr_dev_;
};

template <typename ScalarType>
MagmaBatchedGemm<ScalarType>::MagmaBatchedGemm(magma_queue_t queue)
    : queue_(queue), stream_(magma_queue_get_cuda_stream(queue_)) {}

template <typename ScalarType>
MagmaBatchedGemm<ScalarType>::MagmaBatchedGemm(const int size, magma_queue_t queue)
    : MagmaBatchedGemm(queue) {
  reserve(size);
}

template <typename ScalarType>
void MagmaBatchedGemm<ScalarType>::reserve(int size) {
  a_ptr_.reserve(size), b_ptr_.reserve(size), c_ptr_.reserve(size);
}

template <typename ScalarType>
void MagmaBatchedGemm<ScalarType>::synchronizeCopy() {
  copied_.block();

  a_ptr_.resize(0);
  b_ptr_.resize(0);
  c_ptr_.resize(0);
}

template <typename ScalarType>
void MagmaBatchedGemm<ScalarType>::addGemm(const ScalarType* a, const ScalarType* b, ScalarType* c) {
  a_ptr_.push_back(a);
  b_ptr_.push_back(b);
  c_ptr_.push_back(c);
}

template <typename ScalarType>
void MagmaBatchedGemm<ScalarType>::execute(const char transa, const char transb, const int m,
                                           const int n, const int k, const ScalarType alpha,
                                           const ScalarType beta, const int lda, const int ldb,
                                           const int ldc) {
  a_ptr_dev_.setAsync(a_ptr_, stream_);
  b_ptr_dev_.setAsync(b_ptr_, stream_);
  c_ptr_dev_.setAsync(c_ptr_, stream_);
  copied_.record(stream_);

  const int n_batched = a_ptr_.size();
  magma::magmablas_gemm_batched(transa, transb, m, n, k, alpha, a_ptr_dev_.ptr(), lda,
                                b_ptr_dev_.ptr(), ldb, beta, c_ptr_dev_.ptr(), ldc, n_batched,
                                queue_);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

}  // util
}  // linalg
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_LINALG_UTIL_MAGMA_BATCHED_GEMM_HPP
