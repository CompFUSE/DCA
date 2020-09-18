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

#ifndef DCA_LINALG_UTIL_MAGMA_VBATCHED_GEMM_HPP
#define DCA_LINALG_UTIL_MAGMA_VBATCHED_GEMM_HPP
#ifdef DCA_HAVE_CUDA

#include <cassert>
#include <vector>

#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/util/allocators/vectors_typedefs.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/util/magma_queue.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <typename ScalarType>
class MagmaVBatchedGemm {
public:
  // Creates a plan for a batched gemm with variable size.
  MagmaVBatchedGemm(const linalg::util::MagmaQueue& queue);
  // Creates a plan and allocates the memory for the arguments of `size`
  // multiplications.
  MagmaVBatchedGemm(int size, const linalg::util::MagmaQueue& queue);

  // Allocates the memory for the arguments of `size` multiplications.
  void reserve(int size);

  // This methods sets the argument for a single gemm operation.
  // See the documentation of the MAGMA library for a description of the arguments.
  void addGemm(int m, int n, int k, const ScalarType* a, int lda, const ScalarType* b, int ldb,
               ScalarType* c, int ldc);

  // Transfers the arguments to the device and performs the batched gemm of matrices of unequal
  // size.
  void execute(char transa = 'N', char transb = 'N');

  // Synchronizes the copy of the arguments to the device and clears the arguments of the batched
  // gemm before reusing the object for a new batched operation.
  void synchronizeCopy();

private:
  const linalg::util::MagmaQueue& queue_;
  CudaEvent copied_;

  linalg::util::HostVector<const ScalarType*> a_ptr_, b_ptr_;
  linalg::util::HostVector<ScalarType*> c_ptr_;
  linalg::util::HostVector<int> m_, n_, k_, lda_, ldb_, ldc_;
  int m_max_, n_max_, k_max_;

  linalg::Vector<const ScalarType*, linalg::GPU> a_ptr_dev_, b_ptr_dev_;
  linalg::Vector<ScalarType*, linalg::GPU> c_ptr_dev_;
  linalg::Vector<int, linalg::GPU> m_dev_, n_dev_, k_dev_, lda_dev_, ldc_dev_, ldb_dev_;
};

template <typename ScalarType>
MagmaVBatchedGemm<ScalarType>::MagmaVBatchedGemm(const linalg::util::MagmaQueue& queue)
    : queue_(queue), m_max_(0), n_max_(0), k_max_(0) {}

template <typename ScalarType>
MagmaVBatchedGemm<ScalarType>::MagmaVBatchedGemm(const int size, const linalg::util::MagmaQueue& queue)
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
void MagmaVBatchedGemm<ScalarType>::addGemm(const int m, const int n, const int k,
                                            const ScalarType* a, const int lda, const ScalarType* b,
                                            const int ldb, ScalarType* c, const int ldc) {
  m_.push_back(m);
  n_.push_back(n);
  k_.push_back(k);

  m_max_ = std::max(m, m_max_);
  n_max_ = std::max(n, n_max_);
  k_max_ = std::max(k, k_max_);

  lda_.push_back(lda);
  ldb_.push_back(ldb);
  ldc_.push_back(ldc);
  a_ptr_.push_back(a);
  b_ptr_.push_back(b);
  c_ptr_.push_back(c);
}

template <typename ScalarType>
void MagmaVBatchedGemm<ScalarType>::execute(const char transa, const char transb) {
  m_.push_back(0);
  n_.push_back(0);
  k_.push_back(0);
  lda_.push_back(0);
  ldb_.push_back(0);
  ldc_.push_back(0);

  // TODO: store in a buffer if the performance gain is necessary.
  const auto& stream = queue_.getStream();
  a_ptr_dev_.setAsync(a_ptr_, stream);
  b_ptr_dev_.setAsync(b_ptr_, stream);
  c_ptr_dev_.setAsync(c_ptr_, stream);
  lda_dev_.setAsync(lda_, stream);
  ldb_dev_.setAsync(ldb_, stream);
  ldc_dev_.setAsync(ldc_, stream);
  m_dev_.setAsync(m_, stream);
  n_dev_.setAsync(n_, stream);
  k_dev_.setAsync(k_, stream);

  copied_.record(stream);

  const int n_batched = a_ptr_.size();
  magma::magmablas_gemm_vbatched_max_nocheck(
      transa, transb, m_dev_.ptr(), n_dev_.ptr(), k_dev_.ptr(), ScalarType(1), a_ptr_dev_.ptr(),
      lda_dev_.ptr(), b_ptr_dev_.ptr(), ldb_dev_.ptr(), ScalarType(0), c_ptr_dev_.ptr(),
      ldc_dev_.ptr(), n_batched, m_max_, n_max_, k_max_, queue_);

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
  m_.resize(0);
  n_.resize(0);
  k_.resize(0);

  m_max_ = n_max_ = k_max_ = 0;
}

}  // namespace util
}  // namespace linalg
}  // namespace dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_LINALG_UTIL_MAGMA_VBATCHED_GEMM_HPP
