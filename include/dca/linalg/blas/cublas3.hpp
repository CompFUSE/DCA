// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the C++ wrappers for some of the level 3 CUBLAS functions.

#ifndef DCA_LINALG_BLAS_CUBLAS3_HPP
#define DCA_LINALG_BLAS_CUBLAS3_HPP

#include <complex>
#include <cublas_v2.h>

#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/linalg/blas/cublas_conversion_char_types.hpp"
#include "dca/linalg/util/error_cublas.hpp"
#include "dca/linalg/util/error_cuda.hpp"

// C++ wrappers
namespace dca {
namespace linalg {
namespace cublas {
// dca::linalg::cublas::

inline void gemm(cublasHandle_t handle, const char* transa, const char* transb, int m, int n, int k,
                 float alpha, const float* a, int lda, const float* b, int ldb, float beta,
                 float* c, int ldc) {
  checkErrorsCudaDebug();
  cublasStatus_t status =
      cublasSgemm(handle, getCublasTransValue(*transa), getCublasTransValue(*transb), m, n, k,
                  &alpha, a, lda, b, ldb, &beta, c, ldc);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void gemm(cublasHandle_t handle, const char* transa, const char* transb, int m, int n, int k,
                 double alpha, const double* a, int lda, const double* b, int ldb, double beta,
                 double* c, int ldc) {
  checkErrorsCudaDebug();
  cublasStatus_t status =
      cublasDgemm(handle, getCublasTransValue(*transa), getCublasTransValue(*transb), m, n, k,
                  &alpha, a, lda, b, ldb, &beta, c, ldc);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void gemm(cublasHandle_t handle, const char* transa, const char* transb, int m, int n, int k,
                 std::complex<float> alpha, const std::complex<float>* a, int lda,
                 const std::complex<float>* b, int ldb, std::complex<float> beta,
                 std::complex<float>* c, int ldc) {
  checkErrorsCudaDebug();
  const cuComplex* cu_alpha = util::castCudaComplex(alpha);
  const cuComplex* cu_a = util::castCudaComplex(a);
  const cuComplex* cu_b = util::castCudaComplex(b);
  const cuComplex* cu_beta = util::castCudaComplex(beta);
  cuComplex* cu_c = util::castCudaComplex(c);

  cublasStatus_t status =
      cublasCgemm(handle, getCublasTransValue(*transa), getCublasTransValue(*transb), m, n, k,
                  cu_alpha, cu_a, lda, cu_b, ldb, cu_beta, cu_c, ldc);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void gemm(cublasHandle_t handle, const char* transa, const char* transb, int m, int n, int k,
                 std::complex<double> alpha, const std::complex<double>* a, int lda,
                 const std::complex<double>* b, int ldb, std::complex<double> beta,
                 std::complex<double>* c, int ldc) {
  checkErrorsCudaDebug();
  const cuDoubleComplex* cu_alpha = util::castCudaComplex(alpha);
  const cuDoubleComplex* cu_a = util::castCudaComplex(a);
  const cuDoubleComplex* cu_b = util::castCudaComplex(b);
  const cuDoubleComplex* cu_beta = util::castCudaComplex(beta);
  cuDoubleComplex* cu_c = util::castCudaComplex(c);

  cublasStatus_t status =
      cublasZgemm(handle, getCublasTransValue(*transa), getCublasTransValue(*transb), m, n, k,
                  cu_alpha, cu_a, lda, cu_b, ldb, cu_beta, cu_c, ldc);
  checkRC(status);
  checkErrorsCudaDebug();
}

inline void trsm(cublasHandle_t handle, const char* side, const char* uplo, const char* transa,
                 const char* diag, int m, int n, float alpha, const float* a, int lda, float* b,
                 int ldb) {
  checkErrorsCudaDebug();
  cublasStatus_t status = cublasStrsm(handle, getCublasSideValue(*side), getCublasUploValue(*uplo),
                                      getCublasTransValue(*transa), getCublasDiagValue(*diag), m, n,
                                      &alpha, a, lda, b, ldb);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void trsm(cublasHandle_t handle, const char* side, const char* uplo, const char* transa,
                 const char* diag, int m, int n, double alpha, const double* a, int lda, double* b,
                 int ldb) {
  checkErrorsCudaDebug();
  cublasStatus_t status = cublasDtrsm(handle, getCublasSideValue(*side), getCublasUploValue(*uplo),
                                      getCublasTransValue(*transa), getCublasDiagValue(*diag), m, n,
                                      &alpha, a, lda, b, ldb);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void trsm(cublasHandle_t handle, const char* side, const char* uplo, const char* transa,
                 const char* diag, int m, int n, std::complex<float> alpha,
                 const std::complex<float>* a, int lda, std::complex<float>* b, int ldb) {
  checkErrorsCudaDebug();
  const cuComplex* cu_alpha = util::castCudaComplex(alpha);
  const cuComplex* cu_a = util::castCudaComplex(a);
  cuComplex* cu_b = util::castCudaComplex(b);

  cublasStatus_t status = cublasCtrsm(handle, getCublasSideValue(*side), getCublasUploValue(*uplo),
                                      getCublasTransValue(*transa), getCublasDiagValue(*diag), m, n,
                                      cu_alpha, cu_a, lda, cu_b, ldb);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void trsm(cublasHandle_t handle, const char* side, const char* uplo, const char* transa,
                 const char* diag, int m, int n, std::complex<double> alpha,
                 const std::complex<double>* a, int lda, std::complex<double>* b, int ldb) {
  checkErrorsCudaDebug();
  const cuDoubleComplex* cu_alpha = util::castCudaComplex(alpha);
  const cuDoubleComplex* cu_a = util::castCudaComplex(a);
  cuDoubleComplex* cu_b = util::castCudaComplex(b);

  cublasStatus_t status = cublasZtrsm(handle, getCublasSideValue(*side), getCublasUploValue(*uplo),
                                      getCublasTransValue(*transa), getCublasDiagValue(*diag), m, n,
                                      cu_alpha, cu_a, lda, cu_b, ldb);
  checkRC(status);
  checkErrorsCudaDebug();
}
}  // cublas
}  // linalg
}  // dca

#endif  // DCA_LINALG_BLAS_CUBLAS3_HPP
