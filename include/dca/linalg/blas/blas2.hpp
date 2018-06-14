// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the declaration of the BLAS level 2 routines and implements C++ wrappers.

#ifndef DCA_LINALG_BLAS_BLAS2_HPP
#define DCA_LINALG_BLAS_BLAS2_HPP

#include <complex>

// Declaration of the BLAS functions. Do not use them in the code but use the provided wrappers.
extern "C" {
void sgemv_(const char* trans, const int* m, const int* n, const float* alpha, const float* a,
            const int* lda, const float* x, const int* incx, const float* beta, float* y,
            const int* incy);
void dgemv_(const char* trans, const int* m, const int* n, const double* alpha, const double* a,
            const int* lda, const double* x, const int* incx, const double* beta, double* y,
            const int* incy);
void cgemv_(const char* trans, const int* m, const int* n, const std::complex<float>* alpha,
            const std::complex<float>* a, const int* lda, const std::complex<float>* x,
            const int* incx, const std::complex<float>* beta, std::complex<float>* y,
            const int* incy);
void zgemv_(const char* trans, const int* m, const int* n, const std::complex<double>* alpha,
            const std::complex<double>* a, const int* lda, const std::complex<double>* x,
            const int* incx, const std::complex<double>* beta, std::complex<double>* y,
            const int* incy);

void chemv_(const char* uplo, const int* n, const std::complex<float>* alpha,
            const std::complex<float>* a, const int* lda, const std::complex<float>* x,
            const int* incx, const std::complex<float>* beta, std::complex<float>* y,
            const int* incy);
void zhemv_(const char* uplo, const int* n, const std::complex<double>* alpha,
            const std::complex<double>* a, const int* lda, const std::complex<double>* x,
            const int* incx, const std::complex<double>* beta, std::complex<double>* y,
            const int* incy);

void ssymv_(const char* uplo, const int* n, const float* alpha, const float* a, const int* lda,
            const float* x, const int* incx, const float* beta, float* y, const int* incy);
void dsymv_(const char* uplo, const int* n, const double* alpha, const double* a, const int* lda,
            const double* x, const int* incx, const double* beta, double* y, const int* incy);

void strmv_(const char* uplo, const char* transa, const char* diag, const int* n, const float* a,
            const int* lda, float* x, const int* incx);
void dtrmv_(const char* uplo, const char* transa, const char* diag, const int* n, const double* a,
            const int* lda, double* x, const int* incx);
void ctrmv_(const char* uplo, const char* transa, const char* diag, const int* n,
            const std::complex<float>* a, const int* lda, std::complex<float>* x, const int* incx);
void ztrmv_(const char* uplo, const char* transa, const char* diag, const int* n,
            const std::complex<double>* a, const int* lda, std::complex<double>* x, const int* incx);

void strsv_(const char* uplo, const char* transa, const char* diag, const int* n, const float* a,
            const int* lda, float* x, const int* incx);
void dtrsv_(const char* uplo, const char* transa, const char* diag, const int* n, const double* a,
            const int* lda, double* x, const int* incx);
void ctrsv_(const char* uplo, const char* transa, const char* diag, const int* n,
            const std::complex<float>* a, const int* lda, std::complex<float>* x, const int* incx);
void ztrsv_(const char* uplo, const char* transa, const char* diag, const int* n,
            const std::complex<double>* a, const int* lda, std::complex<double>* x, const int* incx);

void sger_(const int* m, const int* n, const float* alpha, const float* x, const int* incx,
           const float* y, const int* incy, float* a, const int* lda);
void dger_(const int* m, const int* n, const double* alpha, const double* x, const int* incx,
           const double* y, const int* incy, double* a, const int* lda);

void cgeru_(const int* m, const int* n, const std::complex<float>* alpha,
            const std::complex<float>* x, const int* incx, const std::complex<float>* y,
            const int* incy, std::complex<float>* a, const int* lda);
void zgeru_(const int* m, const int* n, const std::complex<double>* alpha,
            const std::complex<double>* x, const int* incx, const std::complex<double>* y,
            const int* incy, std::complex<double>* a, const int* lda);

void cgerc_(const int* m, const int* n, const std::complex<float>* alpha,
            const std::complex<float>* x, const int* incx, const std::complex<float>* y,
            const int* incy, std::complex<float>* a, const int* lda);
void zgerc_(const int* m, const int* n, const std::complex<double>* alpha,
            const std::complex<double>* x, const int* incx, const std::complex<double>* y,
            const int* incy, std::complex<double>* a, const int* lda);

void ssyr_(const char* uplo, const int* n, const float* alpha, const float* x, const int* incx,
           float* a, const int* lda);
void dsyr_(const char* uplo, const int* n, const double* alpha, const double* x, const int* incx,
           double* a, const int* lda);

void cher_(const char* uplo, const int* n, const float* alpha, const std::complex<float>* x,
           const int* incx, std::complex<float>* a, const int* lda);
void zher_(const char* uplo, const int* n, const double* alpha, const std::complex<double>* x,
           const int* incx, std::complex<double>* a, const int* lda);

void ssyr2_(const char* uplo, const int* n, const float* alpha, const float* x, const int* incx,
            const float* y, const int* incy, float* a, const int* lda);
void dsyr2_(const char* uplo, const int* n, const double* alpha, const double* x, const int* incx,
            const double* y, const int* incy, double* a, const int* lda);

void cher2_(const char* uplo, const int* n, const std::complex<float>* alpha,
            const std::complex<float>* x, const int* incx, const std::complex<float>* y,
            const int* incy, std::complex<float>* a, const int* lda);
void zher2_(const char* uplo, const int* n, const std::complex<double>* alpha,
            const std::complex<double>* x, const int* incx, const std::complex<double>* y,
            const int* incy, std::complex<double>* a, const int* lda);
}

// C++ wrappers
namespace dca {
namespace linalg {
namespace blas {
// dca::linalg::blas::
inline void gemv(const char* trans, int m, int n, float alpha, const float* a, int lda,
                 const float* x, int incx, float beta, float* y, int incy) {
  sgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}
inline void gemv(const char* trans, int m, int n, double alpha, const double* a, int lda,
                 const double* x, int incx, double beta, double* y, int incy) {
  dgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}
inline void gemv(const char* trans, int m, int n, std::complex<float> alpha,
                 const std::complex<float>* a, int lda, const std::complex<float>* x, int incx,
                 std::complex<float> beta, std::complex<float>* y, int incy) {
  cgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}
inline void gemv(const char* trans, int m, int n, std::complex<double> alpha,
                 const std::complex<double>* a, int lda, const std::complex<double>* x, int incx,
                 std::complex<double> beta, std::complex<double>* y, int incy) {
  zgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

inline void hemv(const char* uplo, int n, std::complex<float> alpha, const std::complex<float>* a,
                 int lda, const std::complex<float>* x, int incx, std::complex<float> beta,
                 std::complex<float>* y, int incy) {
  chemv_(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}
inline void hemv(const char* uplo, int n, std::complex<double> alpha, const std::complex<double>* a,
                 int lda, const std::complex<double>* x, int incx, std::complex<double> beta,
                 std::complex<double>* y, int incy) {
  zhemv_(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

inline void symv(const char* uplo, int n, float alpha, const float* a, int lda, const float* x,
                 int incx, float beta, float* y, int incy) {
  ssymv_(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}
inline void symv(const char* uplo, int n, double alpha, const double* a, int lda, const double* x,
                 int incx, double beta, double* y, int incy) {
  dsymv_(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

inline void trmv(const char* uplo, const char* transa, const char* diag, int n, const float* a,
                 int lda, float* x, int incx) {
  strmv_(uplo, transa, diag, &n, a, &lda, x, &incx);
}
inline void trmv(const char* uplo, const char* transa, const char* diag, int n, const double* a,
                 int lda, double* x, int incx) {
  dtrmv_(uplo, transa, diag, &n, a, &lda, x, &incx);
}
inline void trmv(const char* uplo, const char* transa, const char* diag, int n,
                 const std::complex<float>* a, int lda, std::complex<float>* x, int incx) {
  ctrmv_(uplo, transa, diag, &n, a, &lda, x, &incx);
}
inline void trmv(const char* uplo, const char* transa, const char* diag, int n,
                 const std::complex<double>* a, int lda, std::complex<double>* x, int incx) {
  ztrmv_(uplo, transa, diag, &n, a, &lda, x, &incx);
}

inline void trsv(const char* uplo, const char* transa, const char* diag, int n, const float* a,
                 int lda, float* x, int incx) {
  strsv_(uplo, transa, diag, &n, a, &lda, x, &incx);
}
inline void trsv(const char* uplo, const char* transa, const char* diag, int n, const double* a,
                 int lda, double* x, int incx) {
  dtrsv_(uplo, transa, diag, &n, a, &lda, x, &incx);
}
inline void trsv(const char* uplo, const char* transa, const char* diag, int n,
                 const std::complex<float>* a, int lda, std::complex<float>* x, int incx) {
  ctrsv_(uplo, transa, diag, &n, a, &lda, x, &incx);
}
inline void trsv(const char* uplo, const char* transa, const char* diag, int n,
                 const std::complex<double>* a, int lda, std::complex<double>* x, int incx) {
  ztrsv_(uplo, transa, diag, &n, a, &lda, x, &incx);
}

inline void ger(int m, int n, float alpha, const float* x, int incx, const float* y, int incy,
                float* a, int lda) {
  sger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}
inline void ger(int m, int n, double alpha, const double* x, int incx, const double* y, int incy,
                double* a, int lda) {
  dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

inline void geru(int m, int n, std::complex<float> alpha, const std::complex<float>* x, int incx,
                 const std::complex<float>* y, int incy, std::complex<float>* a, int lda) {
  cgeru_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}
inline void geru(int m, int n, std::complex<double> alpha, const std::complex<double>* x, int incx,
                 const std::complex<double>* y, int incy, std::complex<double>* a, int lda) {
  zgeru_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

inline void gerc(int m, int n, std::complex<float> alpha, const std::complex<float>* x, int incx,
                 const std::complex<float>* y, int incy, std::complex<float>* a, int lda) {
  cgerc_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}
inline void gerc(int m, int n, std::complex<double> alpha, const std::complex<double>* x, int incx,
                 const std::complex<double>* y, int incy, std::complex<double>* a, int lda) {
  zgerc_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

inline void syr(const char* uplo, int n, float alpha, const float* x, int incx, float* a, int lda) {
  ssyr_(uplo, &n, &alpha, x, &incx, a, &lda);
}
inline void syr(const char* uplo, int n, double alpha, const double* x, int incx, double* a, int lda) {
  dsyr_(uplo, &n, &alpha, x, &incx, a, &lda);
}

inline void her(const char* uplo, int n, float alpha, const std::complex<float>* x, int incx,
                std::complex<float>* a, int lda) {
  cher_(uplo, &n, &alpha, x, &incx, a, &lda);
}
inline void her(const char* uplo, int n, double alpha, const std::complex<double>* x, int incx,
                std::complex<double>* a, int lda) {
  zher_(uplo, &n, &alpha, x, &incx, a, &lda);
}

inline void syr2(const char* uplo, int n, float alpha, const float* x, int incx, const float* y,
                 int incy, float* a, int lda) {
  ssyr2_(uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}
inline void syr2(const char* uplo, int n, double alpha, const double* x, int incx, const double* y,
                 int incy, double* a, int lda) {
  dsyr2_(uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

inline void her2(const char* uplo, int n, std::complex<float> alpha, const std::complex<float>* x,
                 int incx, const std::complex<float>* y, int incy, std::complex<float>* a, int lda) {
  cher2_(uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}
inline void her2(const char* uplo, int n, std::complex<double> alpha, const std::complex<double>* x,
                 int incx, const std::complex<double>* y, int incy, std::complex<double>* a, int lda) {
  zher2_(uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}
}  // blas
}  // linalg
}  // dca

#endif  // DCA_LINALG_BLAS_BLAS2_HPP
