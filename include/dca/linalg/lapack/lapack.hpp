// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the declaration of some of the LAPACK routines and implements C++ wrappers.

#ifndef DCA_LINALG_LAPACK_LAPACK_HPP
#define DCA_LINALG_LAPACK_LAPACK_HPP

#include <complex>

// Declaration of the LAPACK functions. Do not use them in the code but use the provided wrappers.
extern "C" {
void slaset_(const char* uplo, const int* m, const int* n, const float* alpha, const float* beta,
             float* a, const int* lda);
void dlaset_(const char* uplo, const int* m, const int* n, const double* alpha, const double* beta,
             double* a, const int* lda);
void claset_(const char* uplo, const int* m, const int* n, const std::complex<float>* alpha,
             const std::complex<float>* beta, std::complex<float>* a, const int* lda);
void zlaset_(const char* uplo, const int* m, const int* n, const std::complex<double>* alpha,
             const std::complex<double>* beta, std::complex<double>* a, const int* lda);

void sgetrf_(const int* m, const int* n, float* a, const int* lda, int* ipiv, int* info);
void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
void cgetrf_(const int* m, const int* n, std::complex<float>* a, const int* lda, int* ipiv,
             int* info);
void zgetrf_(const int* m, const int* n, std::complex<double>* a, const int* lda, int* ipiv,
             int* info);

void sgetri_(const int* n, float* a, const int* lda, int* ipiv, float* work, const int* lwork,
             int* info);
void dgetri_(const int* n, double* a, const int* lda, int* ipiv, double* work, const int* lwork,
             int* info);
void cgetri_(const int* n, std::complex<float>* a, const int* lda, int* ipiv,
             std::complex<float>* work, const int* lwork, int* info);
void zgetri_(const int* n, std::complex<double>* a, const int* lda, int* ipiv,
             std::complex<double>* work, const int* lwork, int* info);
}

// C++ wrappers
namespace dca {
namespace linalg {
namespace lapack {
// dca::linalg::lapack::

inline void laset(const char* uplo, int m, int n, float alpha, float beta, float* a, int lda) {
  slaset_(uplo, &m, &n, &alpha, &beta, a, &lda);
}
inline void laset(const char* uplo, int m, int n, double alpha, double beta, double* a, int lda) {
  dlaset_(uplo, &m, &n, &alpha, &beta, a, &lda);
}
inline void laset(const char* uplo, int m, int n, std::complex<float> alpha,
                  std::complex<float> beta, std::complex<float>* a, int lda) {
  claset_(uplo, &m, &n, &alpha, &beta, a, &lda);
}
inline void laset(const char* uplo, int m, int n, std::complex<double> alpha,
                  std::complex<double> beta, std::complex<double>* a, int lda) {
  zlaset_(uplo, &m, &n, &alpha, &beta, a, &lda);
}

inline void getrf(int m, int n, float* a, int lda, int* ipiv, int* info) {
  sgetrf_(&m, &n, a, &lda, ipiv, info);
}
inline void getrf(int m, int n, double* a, int lda, int* ipiv, int* info) {
  dgetrf_(&m, &n, a, &lda, ipiv, info);
}
inline void getrf(int m, int n, std::complex<float>* a, int lda, int* ipiv, int* info) {
  cgetrf_(&m, &n, a, &lda, ipiv, info);
}
inline void getrf(int m, int n, std::complex<double>* a, int lda, int* ipiv, int* info) {
  zgetrf_(&m, &n, a, &lda, ipiv, info);
}

inline void getri(int n, float* a, int lda, int* ipiv, float* work, int lwork, int* info) {
  sgetri_(&n, a, &lda, ipiv, work, &lwork, info);
}
inline void getri(int n, double* a, int lda, int* ipiv, double* work, int lwork, int* info) {
  dgetri_(&n, a, &lda, ipiv, work, &lwork, info);
}
inline void getri(int n, std::complex<float>* a, int lda, int* ipiv, std::complex<float>* work,
                  int lwork, int* info) {
  cgetri_(&n, a, &lda, ipiv, work, &lwork, info);
}
inline void getri(int n, std::complex<double>* a, int lda, int* ipiv, std::complex<double>* work,
                  int lwork, int* info) {
  zgetri_(&n, a, &lda, ipiv, work, &lwork, info);
}

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_LAPACK_HPP
