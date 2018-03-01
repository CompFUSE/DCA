// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a class to store and use the arguments of a batched gemm operation.

#ifndef DCA_LINALG_UTIL_MAGMA_VBATCHED_GEMM_HPP
#define DCA_LINALG_UTIL_MAGMA_VBATCHED_GEMM_HPP
#ifdef DCA_HAVE_CUDA

#include <cassert>
#include <vector>

#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/util/allocators.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <typename ScalarType>
class MagmaVBatchedGemm {
public:
  // Creates a plan for a batched gemm with variable size.
  MagmaVBatchedGemm(magma_queue_t queue);
  // Creates a plan and allocates the memory for the arguments of `size`
  // multiplications.
  MagmaVBatchedGemm(int size, magma_queue_t queue);

  // Allocates the memory for the arguments of `size` multiplications.
  void reserve(int size);

  // This methods sets the argument for a single gemm operation.
  // See the documentation of the MAGMA library for a description of the arguments.
  void addGemm(int n, int m, int k, const ScalarType* a, int lda, const ScalarType* b, int ldb,
               ScalarType* c, int ldc);

  // Transfers the arguments to the device and performs the batched gemm of matrices of unequal
  // size.
  void execute(char transa = 'N', char transb = 'N');

  // Synchronizes the copy of the arguments to the device and clears the arguments of the batched
  // gemm before reusing the object for a new batched operation.
  void synchronizeCopy();

private:
  magma_queue_t queue_;
  const cudaStream_t stream_;
  CudaEvent copied_;

  linalg::util::HostVector<const ScalarType *> a_ptr_, b_ptr_;
  linalg::util::HostVector<ScalarType*> c_ptr_;
  linalg::util::HostVector<int> n_, m_, k_, lda_, ldb_, ldc_;

  linalg::Vector<const ScalarType *, linalg::GPU> a_ptr_dev_, b_ptr_dev_;
  linalg::Vector<ScalarType*, linalg::GPU> c_ptr_dev_;
  linalg::Vector<int, linalg::GPU> n_dev_, m_dev_, k_dev_, lda_dev_, ldc_dev_, ldb_dev_;
};

template <typename ScalarType>
MagmaVBatchedGemm<ScalarType>::MagmaVBatchedGemm(magma_queue_t queue)
    : queue_(queue), stream_(magma_queue_get_cuda_stream(queue_)) {}

template <typename ScalarType>
MagmaVBatchedGemm<ScalarType>::MagmaVBatchedGemm(const int size, magma_queue_t queue)
    : MagmaVBatchedGemm(queue) {
  reserve(size);
}

template <typename ScalarType>
void MagmaVBatchedGemm<ScalarType>::reserve(int size) {
  a_ptr_.reserve(size), b_ptr_.reserve(size), c_ptr_.reserve(size);
  n_.reserve(size + 1), m_.reserve(size + 1), k_.reserve(size + 1);
  lda_.reserve(size + 1), ldb_.reserve(size + 1), ldc_.reserve(size + 1);
}

template <typename ScalarType>
void MagmaVBatchedGemm<ScalarType>::addGemm(const int n, const int m, const int k,
                                            const ScalarType* a, const int lda, const ScalarType* b,
                                            const int ldb, ScalarType* c, const int ldc) {
  n_.push_back(n);
  m_.push_back(m);
  k_.push_back(k);
  lda_.push_back(lda);
  ldb_.push_back(ldb);
  ldc_.push_back(ldc);
  a_ptr_.push_back(a);
  b_ptr_.push_back(b);
  c_ptr_.push_back(c);
}

template <typename ScalarType>
void MagmaVBatchedGemm<ScalarType>::execute(const char transa, const char transb) {
  n_.push_back(0);
  m_.push_back(0);
  k_.push_back(0);
  lda_.push_back(0);
  ldb_.push_back(0);
  ldc_.push_back(0);

  // TODO: store in a buffer if the performance gain is neccessair.
  a_ptr_dev_.setAsync(a_ptr_, stream_);
  b_ptr_dev_.setAsync(b_ptr_, stream_);
  c_ptr_dev_.setAsync(c_ptr_, stream_);
  lda_dev_.setAsync(lda_, stream_);
  ldb_dev_.setAsync(ldb_, stream_);
  ldc_dev_.setAsync(ldc_, stream_);
  n_dev_.setAsync(n_, stream_);
  k_dev_.setAsync(k_, stream_);
  m_dev_.setAsync(m_, stream_);

  copied_.record(stream_);

  const int n_batched = a_ptr_.size();
  magma::magmablas_gemm_vbatched(transa, transb, n_dev_.ptr(), m_dev_.ptr(), k_dev_.ptr(),
                                 ScalarType(1), a_ptr_dev_.ptr(), lda_dev_.ptr(), b_ptr_dev_.ptr(),
                                 ldb_dev_.ptr(), ScalarType(0), c_ptr_dev_.ptr(), ldc_dev_.ptr(),
                                 n_batched, queue_);

  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <typename ScalarType>
void MagmaVBatchedGemm<ScalarType>::synchronizeCopy() {
  copied_.block();

  a_ptr_.resize(0);
  b_ptr_.resize(0);
  c_ptr_.resize(0);
  lda_.resize(0);
  ldb_.resize(0);
  ldc_.resize(0);
  n_.resize(0);
  m_.resize(0);
  k_.resize(0);
}

}  // util
}  // linalg
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_LINALG_UTIL_MAGMA_VBATCHED_GEMM_HPP
