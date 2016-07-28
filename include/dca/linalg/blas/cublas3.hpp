// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the declaration of the BLAS level 3 routines and implements C++ wrappers.

#ifndef DCA_LINALG_CUBLAS3_HPP
#define DCA_LINALG_CUBLAS3_HPP

#include <complex>
#include <cublas_v2.h>

#include "comp_library/linalg/basic_cuda_functions.h"
#include "comp_library/linalg/basic_cublas_functions.h"

// C++ wrappers
namespace dca {
namespace linalg {
namespace cublas {
// dca::linalg::cublas::

using LIN_ALG::cublas_side_type;
using LIN_ALG::cublas_triangle_type;
using LIN_ALG::cublas_operation_type;
using LIN_ALG::cublas_diagonal_type;

inline void gemm(cublasHandle_t handle, const char* transa, const char* transb, int m, int n, int k,
                 float alpha, const float* a, int lda, const float* b, int ldb, float beta,
                 float* c, int ldc) {
  cublasStatus_t status =
      cublasSgemm(handle, cublas_operation_type(*transa), cublas_operation_type(*transb), m, n, k,
                  &alpha, a, lda, b, ldb, &beta, c, ldc);
  cublasCheckReturnCode(status);
}
inline void gemm(cublasHandle_t handle, const char* transa, const char* transb, int m, int n, int k,
                 double alpha, const double* a, int lda, const double* b, int ldb, double beta,
                 double* c, int ldc) {
  cublasStatus_t status =
      cublasDgemm(handle, cublas_operation_type(*transa), cublas_operation_type(*transb), m, n, k,
                  &alpha, a, lda, b, ldb, &beta, c, ldc);
  cublasCheckReturnCode(status);
}
inline void gemm(cublasHandle_t handle, const char* transa, const char* transb, int m, int n, int k,
                 std::complex<float> alpha, const std::complex<float>* a, int lda,
                 const std::complex<float>* b, int ldb, std::complex<float> beta,
                 std::complex<float>* c, int ldc) {
  const cuComplex* cu_alpha = reinterpret_cast<const cuComplex*>(&alpha);
  const cuComplex* cu_a = reinterpret_cast<const cuComplex*>(a);
  const cuComplex* cu_b = reinterpret_cast<const cuComplex*>(b);
  const cuComplex* cu_beta = reinterpret_cast<const cuComplex*>(&beta);
  cuComplex* cu_c = reinterpret_cast<cuComplex*>(c);

  cublasStatus_t status =
      cublasCgemm(handle, cublas_operation_type(*transa), cublas_operation_type(*transb), m, n, k,
                  cu_alpha, cu_a, lda, cu_b, ldb, cu_beta, cu_c, ldc);
  cublasCheckReturnCode(status);
}
inline void gemm(cublasHandle_t handle, const char* transa, const char* transb, int m, int n, int k,
                 std::complex<double> alpha, const std::complex<double>* a, int lda,
                 const std::complex<double>* b, int ldb, std::complex<double> beta,
                 std::complex<double>* c, int ldc) {
  const cuDoubleComplex* cu_alpha = reinterpret_cast<const cuDoubleComplex*>(&alpha);
  const cuDoubleComplex* cu_a = reinterpret_cast<const cuDoubleComplex*>(a);
  const cuDoubleComplex* cu_b = reinterpret_cast<const cuDoubleComplex*>(b);
  const cuDoubleComplex* cu_beta = reinterpret_cast<const cuDoubleComplex*>(&beta);
  cuDoubleComplex* cu_c = reinterpret_cast<cuDoubleComplex*>(c);

  cublasStatus_t status =
      cublasZgemm(handle, cublas_operation_type(*transa), cublas_operation_type(*transb), m, n, k,
                  cu_alpha, cu_a, lda, cu_b, ldb, cu_beta, cu_c, ldc);
  cublasCheckReturnCode(status);
}

inline void trsm(cublasHandle_t handle, const char* side, const char* uplo, const char* transa,
                 const char* diag, int m, int n, float alpha, const float* a, int lda, float* b,
                 int ldb) {
  cublasStatus_t status = cublasStrsm(handle, cublas_side_type(*side), cublas_triangle_type(*uplo),
                                      cublas_operation_type(*transa), cublas_diagonal_type(*diag),
                                      m, n, &alpha, a, lda, b, ldb);
  cublasCheckReturnCode(status);
}
inline void trsm(cublasHandle_t handle, const char* side, const char* uplo, const char* transa,
                 const char* diag, int m, int n, double alpha, const double* a, int lda, double* b,
                 int ldb) {
  cublasStatus_t status = cublasDtrsm(handle, cublas_side_type(*side), cublas_triangle_type(*uplo),
                                      cublas_operation_type(*transa), cublas_diagonal_type(*diag),
                                      m, n, &alpha, a, lda, b, ldb);
  cublasCheckReturnCode(status);
}
inline void trsm(cublasHandle_t handle, const char* side, const char* uplo, const char* transa,
                 const char* diag, int m, int n, std::complex<float> alpha,
                 const std::complex<float>* a, int lda, std::complex<float>* b, int ldb) {
  const cuComplex* cu_alpha = reinterpret_cast<const cuComplex*>(&alpha);
  const cuComplex* cu_a = reinterpret_cast<const cuComplex*>(a);
  cuComplex* cu_b = reinterpret_cast<cuComplex*>(b);

  cublasStatus_t status = cublasCtrsm(handle, cublas_side_type(*side), cublas_triangle_type(*uplo),
                                      cublas_operation_type(*transa), cublas_diagonal_type(*diag),
                                      m, n, cu_alpha, cu_a, lda, cu_b, ldb);
  cublasCheckReturnCode(status);
}
inline void trsm(cublasHandle_t handle, const char* side, const char* uplo, const char* transa,
                 const char* diag, int m, int n, std::complex<double> alpha,
                 const std::complex<double>* a, int lda, std::complex<double>* b, int ldb) {
  const cuDoubleComplex* cu_alpha = reinterpret_cast<const cuDoubleComplex*>(&alpha);
  const cuDoubleComplex* cu_a = reinterpret_cast<const cuDoubleComplex*>(a);
  cuDoubleComplex* cu_b = reinterpret_cast<cuDoubleComplex*>(b);

  cublasStatus_t status = cublasZtrsm(handle, cublas_side_type(*side), cublas_triangle_type(*uplo),
                                      cublas_operation_type(*transa), cublas_diagonal_type(*diag),
                                      m, n, cu_alpha, cu_a, lda, cu_b, ldb);
  cublasCheckReturnCode(status);
}
}  // cublas
}  // linalg
}  // dca
#endif  // DCA_LINALG_BLAS3_HPP
