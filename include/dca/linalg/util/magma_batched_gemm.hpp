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

#ifndef DCA_LINALG_UTIL_MAGMA_BATCHED_GEMM_HPP
#define DCA_LINALG_UTIL_MAGMA_BATCHED_GEMM_HPP
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

template <typename Scalar>
class MagmaBatchedGemm {
public:
  // Creates a plan for a batched gemm.
  MagmaBatchedGemm(magma_queue_t queue);
  // Creates a plan for a batched gemm and allocates the memory for the arguments of `size`
  // multiplications.
  MagmaBatchedGemm(int size, magma_queue_t queue);

  // Allocates the memory for the arguments of `size` multiplications.
  void resize(int size);

  // Transfer the arguments to the device and performs the batched gemm of matrices of equal size.
  // See the documentation of the MAGMA library for a description of the arguments.
  void executeBatched(char transa, char transb, int n, int m, int k, Scalar alpha, Scalar beta,
                      int lda, int ldb, int ldc);
  // Transfer the arguments to the device and performs the batched gemm of matrices of unequal size.
  void executeVBatched(char transa = 'N', char transb = 'N');

  // Synchronizes the copy of the arguments to the device before reusing the object for a new
  // batched operation.
  void synchronizeCopy();

  // This methods set the argument for a single gemm operation relative to the index idx.
  // set_lda, set_ldb, set_ldc, set_n, set_k and set_m beeds to be called only before executeVGemm.
  // See the documentation of the MAGMA library for a description of the arguments.
  inline void set_a(int idx, const Scalar* a);
  inline void set_b(int idx, const Scalar* b);
  inline void set_c(int idx, Scalar* c);
  inline void set_lda(int idx, int lda);
  inline void set_ldb(int idx, int ldb);
  inline void set_ldc(int idx, int ldc);
  inline void set_n(int idx, int n);
  inline void set_m(int idx, int m);
  inline void set_k(int idx, int k);

private:
  int n_batched_;
  magma_queue_t queue_;
  const cudaStream_t stream_;
  CudaEvent copied_;

  linalg::util::HostVector<const Scalar *> a_ptr_, b_ptr_;
  linalg::util::HostVector<Scalar*> c_ptr_;
  linalg::util::HostVector<int> n_, m_, k_, lda_, ldb_, ldc_;

  linalg::Vector<const Scalar *, linalg::GPU> a_ptr_dev_, b_ptr_dev_;
  linalg::Vector<Scalar*, linalg::GPU> c_ptr_dev_;
  linalg::Vector<int, linalg::GPU> n_dev_, m_dev_, k_dev_, lda_dev_, ldc_dev_, ldb_dev_;
};

template <typename Scalar>
MagmaBatchedGemm<Scalar>::MagmaBatchedGemm(magma_queue_t queue)
    : n_batched_(0), queue_(queue), stream_(magma_queue_get_cuda_stream(queue_)) {}

template <typename Scalar>
MagmaBatchedGemm<Scalar>::MagmaBatchedGemm(const int size, magma_queue_t queue)
    : MagmaBatchedGemm(queue) {
  resize(size);
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::resize(int size) {
  a_ptr_.resize(size), b_ptr_.resize(size), c_ptr_.resize(size);
  n_.resize(size + 1), m_.resize(size + 1), k_.resize(size + 1);
  lda_.resize(size + 1), ldb_.resize(size + 1), ldc_.resize(size + 1);

  n_batched_ = size;
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::executeVBatched(const char transa, const char transb) {
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

  magma::magmablas_gemm_vbatched(transa, transb, n_dev_.ptr(), m_dev_.ptr(), k_dev_.ptr(), Scalar(1),
                                 a_ptr_dev_.ptr(), lda_dev_.ptr(), b_ptr_dev_.ptr(), ldb_dev_.ptr(),
                                 Scalar(0), c_ptr_dev_.ptr(), ldc_dev_.ptr(), n_batched_, queue_);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::executeBatched(const char transa, const char transb, const int n,
                                              const int m, const int k, const Scalar alpha,
                                              const Scalar beta, const int lda, const int ldb,
                                              const int ldc) {
  a_ptr_dev_.setAsync(a_ptr_, stream_);
  b_ptr_dev_.setAsync(b_ptr_, stream_);
  c_ptr_dev_.setAsync(c_ptr_, stream_);
  copied_.record(stream_);
  magma::magmablas_gemm_batched(transa, transb, n, m, k, alpha, a_ptr_dev_.ptr(), lda,
                                b_ptr_dev_.ptr(), ldb, beta, c_ptr_dev_.ptr(), ldc, n_batched_,
                                queue_);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::synchronizeCopy() {
  copied_.block();
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::set_a(const int idx, const Scalar* a) {
  assert(0 <= idx && idx < n_batched_);
  a_ptr_[idx] = a;
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::set_b(const int idx, const Scalar* b) {
  assert(0 <= idx && idx < n_batched_);
  b_ptr_[idx] = b;
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::set_c(const int idx, Scalar* c) {
  assert(0 <= idx && idx < n_batched_);
  c_ptr_[idx] = c;
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::set_lda(const int idx, const int lda) {
  assert(0 <= idx && idx < n_batched_);
  lda_[idx] = lda;
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::set_ldb(const int idx, const int ldb) {
  assert(0 <= idx && idx < n_batched_);
  ldb_[idx] = ldb;
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::set_ldc(const int idx, const int ldc) {
  assert(0 <= idx && idx < n_batched_);
  ldc_[idx] = ldc;
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::set_n(const int idx, const int n) {
  assert(0 <= idx && idx < n_batched_);
  n_[idx] = n;
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::set_m(const int idx, const int m) {
  assert(0 <= idx && idx < n_batched_);
  m_[idx] = m;
}

template <typename Scalar>
void MagmaBatchedGemm<Scalar>::set_k(const int idx, const int k) {
  assert(0 <= idx && idx < n_batched_);
  k_[idx] = k;
}

}  // util
}  // linalg
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_LINALG_UTIL_MAGMA_BATCHED_GEMM_HPP
