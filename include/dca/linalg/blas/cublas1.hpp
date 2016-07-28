// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the implementation of the C++ wrappers for the CUBLAS level 1 functions.

#ifndef DCA_LINALG_CUBLAS1_HPP
#define DCA_LINALG_CUBLAS1_HPP

#include <complex>
#include <cublas_v2.h>

#include "comp_library/linalg/basic_cuda_functions.h"
#include "comp_library/linalg/basic_cublas_functions.h"

// C++ wrappers
namespace dca {
namespace linalg {
namespace cublas {
// dca::linalg::cublas::

inline void axpy(cublasHandle_t handle, int n, float alpha, const float* x, int incx, float* y,
                 int incy) {
  cublasStatus_t status = cublasSaxpy(handle, n, &alpha, x, incx, y, incy);
  cublasCheckReturnCode(status);
}
inline void axpy(cublasHandle_t handle, int n, double alpha, const double* x, int incx, double* y,
                 int incy) {
  cublasStatus_t status = cublasDaxpy(handle, n, &alpha, x, incx, y, incy);
  cublasCheckReturnCode(status);
}
inline void axpy(cublasHandle_t handle, int n, std::complex<float> alpha,
                 const std::complex<float>* x, int incx, std::complex<float>* y, int incy) {
  const cuComplex* cu_alpha = reinterpret_cast<const cuComplex*>(&alpha);
  const cuComplex* cu_x = reinterpret_cast<const cuComplex*>(x);
  cuComplex* cu_y = reinterpret_cast<cuComplex*>(y);
  cublasStatus_t status = cublasCaxpy(handle, n, cu_alpha, cu_x, incx, cu_y, incy);
  cublasCheckReturnCode(status);
}
inline void axpy(cublasHandle_t handle, int n, std::complex<double> alpha,
                 const std::complex<double>* x, int incx, std::complex<double>* y, int incy) {
  const cuDoubleComplex* cu_alpha = reinterpret_cast<const cuDoubleComplex*>(&alpha);
  const cuDoubleComplex* cu_x = reinterpret_cast<const cuDoubleComplex*>(x);
  cuDoubleComplex* cu_y = reinterpret_cast<cuDoubleComplex*>(y);
  cublasStatus_t status = cublasZaxpy(handle, n, cu_alpha, cu_x, incx, cu_y, incy);
  cublasCheckReturnCode(status);
}
}  // cublas
}  // linalg
}  // dca

#endif  // DCA_LINALG_BLAS1_HPP
