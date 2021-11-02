// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the C++ wrappers for some of the level 1 CUBLAS functions.

#ifndef DCA_LINALG_BLAS_CUBLAS1_HPP
#define DCA_LINALG_BLAS_CUBLAS1_HPP

#include <complex>
#include "dca/platform/dca_gpu_blas.h"

#include "dca/linalg/util/cast_gpu.hpp"
#include "dca/linalg/util/error_gpuBLAS.hpp"

// C++ wrappers
namespace dca {
namespace linalg {
namespace cublas {
// dca::linalg::cublas::

inline void swap(cublasHandle_t handle, int n, float* x, int incx, float* y, int incy) {
  cublasStatus_t status = cublasSswap(handle, n, x, incx, y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void swap(cublasHandle_t handle, int n, double* x, int incx, double* y, int incy) {
  cublasStatus_t status = cublasDswap(handle, n, x, incx, y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void swap(cublasHandle_t handle, int n, std::complex<float>* x, int incx,
                 std::complex<float>* y, int incy) {
  cublasComplex* cu_x = util::castCUBLASComplex(x);
  cublasComplex* cu_y = util::castCUBLASComplex(y);
  cublasStatus_t status = cublasCswap(handle, n, cu_x, incx, cu_y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void swap(cublasHandle_t handle, int n, std::complex<double>* x, int incx,
                 std::complex<double>* y, int incy) {
  cublasDoubleComplex* cu_x = util::castCUBLASComplex(x);
  cublasDoubleComplex* cu_y = util::castCUBLASComplex(y);
  cublasStatus_t status = cublasZswap(handle, n, cu_x, incx, cu_y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}

inline void scal(cublasHandle_t handle, int n, float alpha, float* x, int incx) {
  cublasStatus_t status = cublasSscal(handle, n, &alpha, x, incx);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void scal(cublasHandle_t handle, int n, double alpha, double* x, int incx) {
  cublasStatus_t status = cublasDscal(handle, n, &alpha, x, incx);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void scal(cublasHandle_t handle, int n, std::complex<float> alpha, std::complex<float>* x,
                 int incx) {
  const cublasComplex* cu_alpha = util::castCUBLASComplex(alpha);
  cublasComplex* cu_x = util::castCUBLASComplex(x);
  cublasStatus_t status = cublasCscal(handle, n, cu_alpha, cu_x, incx);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void scal(cublasHandle_t handle, int n, std::complex<double> alpha, std::complex<double>* x,
                 int incx) {
  const cublasDoubleComplex* cu_alpha = util::castCUBLASComplex(alpha);
  cublasDoubleComplex* cu_x = util::castCUBLASComplex(x);
  cublasStatus_t status = cublasZscal(handle, n, cu_alpha, cu_x, incx);
  checkRC(status);
  checkErrorsCudaDebug();
}

inline void copy(cublasHandle_t handle, int n, const float* x, int incx, float* y, int incy) {
  cublasStatus_t status = cublasScopy(handle, n, x, incx, y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void copy(cublasHandle_t handle, int n, const double* x, int incx, double* y, int incy) {
  cublasStatus_t status = cublasDcopy(handle, n, x, incx, y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void copy(cublasHandle_t handle, int n, const std::complex<float>* x, int incx,
                 std::complex<float>* y, int incy) {
  const cublasComplex* cu_x = util::castCUBLASComplex(x);
  cublasComplex* cu_y = util::castCUBLASComplex(y);
  cublasStatus_t status = cublasCcopy(handle, n, cu_x, incx, cu_y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void copy(cublasHandle_t handle, int n, const std::complex<double>* x, int incx,
                 std::complex<double>* y, int incy) {
  const cublasDoubleComplex* cu_x = util::castCUBLASComplex(x);
  cublasDoubleComplex* cu_y = util::castCUBLASComplex(y);
  cublasStatus_t status = cublasZcopy(handle, n, cu_x, incx, cu_y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}

inline void axpy(cublasHandle_t handle, int n, float alpha, const float* x, int incx, float* y,
                 int incy) {
  cublasStatus_t status = cublasSaxpy(handle, n, &alpha, x, incx, y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void axpy(cublasHandle_t handle, int n, double alpha, const double* x, int incx, double* y,
                 int incy) {
  cublasStatus_t status = cublasDaxpy(handle, n, &alpha, x, incx, y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void axpy(cublasHandle_t handle, int n, std::complex<float> alpha,
                 const std::complex<float>* x, int incx, std::complex<float>* y, int incy) {
  const cublasComplex* cu_alpha = util::castCUBLASComplex(alpha);
  const cublasComplex* cu_x = util::castCUBLASComplex(x);
  cublasComplex* cu_y = util::castCUBLASComplex(y);
  cublasStatus_t status = cublasCaxpy(handle, n, cu_alpha, cu_x, incx, cu_y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}
inline void axpy(cublasHandle_t handle, int n, std::complex<double> alpha,
                 const std::complex<double>* x, int incx, std::complex<double>* y, int incy) {
  const cublasDoubleComplex* cu_alpha = util::castCUBLASComplex(alpha);
  const cublasDoubleComplex* cu_x = util::castCUBLASComplex(x);
  cublasDoubleComplex* cu_y = util::castCUBLASComplex(y);
  cublasStatus_t status = cublasZaxpy(handle, n, cu_alpha, cu_x, incx, cu_y, incy);
  checkRC(status);
  checkErrorsCudaDebug();
}
}  // cublas
}  // linalg
}  // dca

#endif  // DCA_LINALG_BLAS_CUBLAS1_HPP
