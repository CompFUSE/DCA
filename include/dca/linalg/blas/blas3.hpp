// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the declaration of the BLAS level 3 routines and implements C++ wrappers.

#ifndef DCA_LINALG_BLAS_BLAS3_HPP
#define DCA_LINALG_BLAS_BLAS3_HPP

#include <complex>

// Declaration of the BLAS functions. Do not use them in the code but use the provided wrappers.
extern "C" {
void sgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
            const float* alpha, const float* a, const int* lda, const float* b, const int* ldb,
            const float* beta, float* c, const int* ldc);
void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
            const double* alpha, const double* a, const int* lda, const double* b, const int* ldb,
            const double* beta, double* c, const int* ldc);
void cgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
            const std::complex<float>* alpha, const std::complex<float>* a, const int* lda,
            const std::complex<float>* b, const int* ldb, const std::complex<float>* beta,
            std::complex<float>* c, const int* ldc);
void zgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
            const std::complex<double>* alpha, const std::complex<double>* a, const int* lda,
            const std::complex<double>* b, const int* ldb, const std::complex<double>* beta,
            std::complex<double>* c, const int* ldc);

void ssymm_(const char* side, const char* uplo, const int* m, const int* n, const float* alpha,
            const float* a, const int* lda, const float* b, const int* ldb, const float* beta,
            float* c, const int* ldc);
void dsymm_(const char* side, const char* uplo, const int* m, const int* n, const double* alpha,
            const double* a, const int* lda, const double* b, const int* ldb, const double* beta,
            double* c, const int* ldc);
void csymm_(const char* side, const char* uplo, const int* m, const int* n,
            const std::complex<float>* alpha, const std::complex<float>* a, const int* lda,
            const std::complex<float>* b, const int* ldb, const std::complex<float>* beta,
            std::complex<float>* c, const int* ldc);
void zsymm_(const char* side, const char* uplo, const int* m, const int* n,
            const std::complex<double>* alpha, const std::complex<double>* a, const int* lda,
            const std::complex<double>* b, const int* ldb, const std::complex<double>* beta,
            std::complex<double>* c, const int* ldc);

void chemm_(const char* side, const char* uplo, const int* m, const int* n,
            const std::complex<float>* alpha, const std::complex<float>* a, const int* lda,
            const std::complex<float>* b, const int* ldb, const std::complex<float>* beta,
            std::complex<float>* c, const int* ldc);
void zhemm_(const char* side, const char* uplo, const int* m, const int* n,
            const std::complex<double>* alpha, const std::complex<double>* a, const int* lda,
            const std::complex<double>* b, const int* ldb, const std::complex<double>* beta,
            std::complex<double>* c, const int* ldc);

void ssyrk_(const char* uplo, const char* trans, const int* n, const int* k, const float* alpha,
            const float* a, const int* lda, const float* beta, float* c, const int* ldc);
void dsyrk_(const char* uplo, const char* trans, const int* n, const int* k, const double* alpha,
            const double* a, const int* lda, const double* beta, double* c, const int* ldc);
void csyrk_(const char* uplo, const char* trans, const int* n, const int* k,
            const std::complex<float>* alpha, const std::complex<float>* a, const int* lda,
            const std::complex<float>* beta, std::complex<float>* c, const int* ldc);
void zsyrk_(const char* uplo, const char* trans, const int* n, const int* k,
            const std::complex<double>* alpha, const std::complex<double>* a, const int* lda,
            const std::complex<double>* beta, std::complex<double>* c, const int* ldc);

void cherk_(const char* uplo, const char* trans, const int* n, const int* k, const float* alpha,
            const std::complex<float>* a, const int* lda, const float* beta, std::complex<float>* c,
            const int* ldc);
void zherk_(const char* uplo, const char* trans, const int* n, const int* k, const double* alpha,
            const std::complex<double>* a, const int* lda, const double* beta,
            std::complex<double>* c, const int* ldc);

void ssyr2k_(const char* uplo, const char* trans, const int* n, const int* k, const float* alpha,
             const float* a, const int* lda, const float* b, const int* ldb, const float* beta,
             float* c, const int* ldc);
void dsyr2k_(const char* uplo, const char* trans, const int* n, const int* k, const double* alpha,
             const double* a, const int* lda, const double* b, const int* ldb, const double* beta,
             double* c, const int* ldc);
void csyr2k_(const char* uplo, const char* trans, const int* n, const int* k,
             const std::complex<float>* alpha, const std::complex<float>* a, const int* lda,
             const std::complex<float>* b, const int* ldb, const std::complex<float>* beta,
             std::complex<float>* c, const int* ldc);
void zsyr2k_(const char* uplo, const char* trans, const int* n, const int* k,
             const std::complex<double>* alpha, const std::complex<double>* a, const int* lda,
             const std::complex<double>* b, const int* ldb, const std::complex<double>* beta,
             std::complex<double>* c, const int* ldc);

void cher2k_(const char* uplo, const char* trans, const int* n, const int* k,
             const std::complex<float>* alpha, const std::complex<float>* a, const int* lda,
             const std::complex<float>* b, const int* ldb, const float* beta,
             std::complex<float>* c, const int* ldc);
void zher2k_(const char* uplo, const char* trans, const int* n, const int* k,
             const std::complex<double>* alpha, const std::complex<double>* a, const int* lda,
             const std::complex<double>* b, const int* ldb, const double* beta,
             std::complex<double>* c, const int* ldc);

void strmm_(const char* side, const char* uplo, const char* transa, const char* diag, const int* m,
            const int* n, const float* alpha, const float* a, const int* lda, float* b,
            const int* ldb);
void dtrmm_(const char* side, const char* uplo, const char* transa, const char* diag, const int* m,
            const int* n, const double* alpha, const double* a, const int* lda, double* b,
            const int* ldb);
void ctrmm_(const char* side, const char* uplo, const char* transa, const char* diag, const int* m,
            const int* n, const std::complex<float>* alpha, const std::complex<float>* a,
            const int* lda, std::complex<float>* b, const int* ldb);
void ztrmm_(const char* side, const char* uplo, const char* transa, const char* diag, const int* m,
            const int* n, const std::complex<double>* alpha, const std::complex<double>* a,
            const int* lda, std::complex<double>* b, const int* ldb);

void strsm_(const char* side, const char* uplo, const char* transa, const char* diag, const int* m,
            const int* n, const float* alpha, const float* a, const int* lda, float* b,
            const int* ldb);
void dtrsm_(const char* side, const char* uplo, const char* transa, const char* diag, const int* m,
            const int* n, const double* alpha, const double* a, const int* lda, double* b,
            const int* ldb);
void ctrsm_(const char* side, const char* uplo, const char* transa, const char* diag, const int* m,
            const int* n, const std::complex<float>* alpha, const std::complex<float>* a,
            const int* lda, std::complex<float>* b, const int* ldb);
void ztrsm_(const char* side, const char* uplo, const char* transa, const char* diag, const int* m,
            const int* n, const std::complex<double>* alpha, const std::complex<double>* a,
            const int* lda, std::complex<double>* b, const int* ldb);
}

// C++ wrappers
namespace dca {
namespace linalg {
namespace blas {
// dca::linalg::blas::
inline void gemm(const char* transa, const char* transb, int m, int n, int k, float alpha,
                 const float* a, int lda, const float* b, int ldb, float beta, float* c, int ldc) {
  sgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void gemm(const char* transa, const char* transb, int m, int n, int k, double alpha,
                 const double* a, int lda, const double* b, int ldb, double beta, double* c, int ldc) {
  dgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void gemm(const char* transa, const char* transb, int m, int n, int k,
                 std::complex<float> alpha, const std::complex<float>* a, int lda,
                 const std::complex<float>* b, int ldb, std::complex<float> beta,
                 std::complex<float>* c, int ldc) {
  cgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void gemm(const char* transa, const char* transb, int m, int n, int k,
                 std::complex<double> alpha, const std::complex<double>* a, int lda,
                 const std::complex<double>* b, int ldb, std::complex<double> beta,
                 std::complex<double>* c, int ldc) {
  zgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

inline void symm(const char* side, const char* uplo, int m, int n, float alpha, const float* a,
                 int lda, const float* b, int ldb, float beta, float* c, int ldc) {
  ssymm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void symm(const char* side, const char* uplo, int m, int n, double alpha, const double* a,
                 int lda, const double* b, int ldb, double beta, double* c, int ldc) {
  dsymm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void symm(const char* side, const char* uplo, int m, int n, std::complex<float> alpha,
                 const std::complex<float>* a, int lda, const std::complex<float>* b, int ldb,
                 std::complex<float> beta, std::complex<float>* c, int ldc) {
  csymm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void symm(const char* side, const char* uplo, int m, int n, std::complex<double> alpha,
                 const std::complex<double>* a, int lda, const std::complex<double>* b, int ldb,
                 std::complex<double> beta, std::complex<double>* c, int ldc) {
  zsymm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

inline void hemm(const char* side, const char* uplo, int m, int n, std::complex<float> alpha,
                 const std::complex<float>* a, int lda, const std::complex<float>* b, int ldb,
                 std::complex<float> beta, std::complex<float>* c, int ldc) {
  chemm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void hemm(const char* side, const char* uplo, int m, int n, std::complex<double> alpha,
                 const std::complex<double>* a, int lda, const std::complex<double>* b, int ldb,
                 std::complex<double> beta, std::complex<double>* c, int ldc) {
  zhemm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

inline void syrk(const char* uplo, const char* trans, int n, int k, float alpha, const float* a,
                 int lda, float beta, float* c, int ldc) {
  ssyrk_(uplo, trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}
inline void syrk(const char* uplo, const char* trans, int n, int k, double alpha, const double* a,
                 int lda, double beta, double* c, int ldc) {
  dsyrk_(uplo, trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}
inline void syrk(const char* uplo, const char* trans, int n, int k, std::complex<float> alpha,
                 const std::complex<float>* a, int lda, std::complex<float> beta,
                 std::complex<float>* c, int ldc) {
  csyrk_(uplo, trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}
inline void syrk(const char* uplo, const char* trans, int n, int k, std::complex<double> alpha,
                 const std::complex<double>* a, int lda, std::complex<double> beta,
                 std::complex<double>* c, int ldc) {
  zsyrk_(uplo, trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

inline void herk(const char* uplo, const char* trans, int n, int k, float alpha,
                 const std::complex<float>* a, int lda, float beta, std::complex<float>* c, int ldc) {
  cherk_(uplo, trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}
inline void herk(const char* uplo, const char* trans, int n, int k, double alpha,
                 const std::complex<double>* a, int lda, double beta, std::complex<double>* c,
                 int ldc) {
  zherk_(uplo, trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

inline void syr2k(const char* uplo, const char* trans, int n, int k, float alpha, const float* a,
                  int lda, const float* b, int ldb, float beta, float* c, int ldc) {
  ssyr2k_(uplo, trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void syr2k(const char* uplo, const char* trans, int n, int k, double alpha, const double* a,
                  int lda, const double* b, int ldb, double beta, double* c, int ldc) {
  dsyr2k_(uplo, trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void syr2k(const char* uplo, const char* trans, int n, int k, std::complex<float> alpha,
                  const std::complex<float>* a, int lda, const std::complex<float>* b, int ldb,
                  std::complex<float> beta, std::complex<float>* c, int ldc) {
  csyr2k_(uplo, trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void syr2k(const char* uplo, const char* trans, int n, int k, std::complex<double> alpha,
                  const std::complex<double>* a, int lda, const std::complex<double>* b, int ldb,
                  std::complex<double> beta, std::complex<double>* c, int ldc) {
  zsyr2k_(uplo, trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

inline void her2k(const char* uplo, const char* trans, int n, int k, std::complex<float> alpha,
                  const std::complex<float>* a, int lda, const std::complex<float>* b, int ldb,
                  float beta, std::complex<float>* c, int ldc) {
  cher2k_(uplo, trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void her2k(const char* uplo, const char* trans, int n, int k, std::complex<double> alpha,
                  const std::complex<double>* a, int lda, const std::complex<double>* b, int ldb,
                  double beta, std::complex<double>* c, int ldc) {
  zher2k_(uplo, trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

inline void trmm(const char* side, const char* uplo, const char* transa, const char* diag, int m,
                 int n, float alpha, const float* a, int lda, float* b, int ldb) {
  strmm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
}
inline void trmm(const char* side, const char* uplo, const char* transa, const char* diag, int m,
                 int n, double alpha, const double* a, int lda, double* b, int ldb) {
  dtrmm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
}
inline void trmm(const char* side, const char* uplo, const char* transa, const char* diag, int m,
                 int n, std::complex<float> alpha, const std::complex<float>* a, int lda,
                 std::complex<float>* b, int ldb) {
  ctrmm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
}
inline void trmm(const char* side, const char* uplo, const char* transa, const char* diag, int m,
                 int n, std::complex<double> alpha, const std::complex<double>* a, int lda,
                 std::complex<double>* b, int ldb) {
  ztrmm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

inline void trsm(const char* side, const char* uplo, const char* transa, const char* diag, int m,
                 int n, float alpha, const float* a, int lda, float* b, int ldb) {
  strsm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
}
inline void trsm(const char* side, const char* uplo, const char* transa, const char* diag, int m,
                 int n, double alpha, const double* a, int lda, double* b, int ldb) {
  dtrsm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
}
inline void trsm(const char* side, const char* uplo, const char* transa, const char* diag, int m,
                 int n, std::complex<float> alpha, const std::complex<float>* a, int lda,
                 std::complex<float>* b, int ldb) {
  ctrsm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
}
inline void trsm(const char* side, const char* uplo, const char* transa, const char* diag, int m,
                 int n, std::complex<double> alpha, const std::complex<double>* a, int lda,
                 std::complex<double>* b, int ldb) {
  ztrsm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
}
}  // blas
}  // linalg
}  // dca

#endif  // DCA_LINALG_BLAS_BLAS3_HPP
