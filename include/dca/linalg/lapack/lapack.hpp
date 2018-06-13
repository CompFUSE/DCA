// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the declaration of some of the LAPACK routines and implements C++ wrappers.

#ifndef DCA_LINALG_LAPACK_LAPACK_HPP
#define DCA_LINALG_LAPACK_LAPACK_HPP

#include <complex>
#include <dca/linalg/util/util_lapack.hpp>

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

void slacpy_(const char* uplo, const int* m, const int* n, const float* a, const int* lda, float* b,
             const int* ldb);
void dlacpy_(const char* uplo, const int* m, const int* n, const double* a, const int* lda,
             double* b, const int* ldb);
void clacpy_(const char* uplo, const int* m, const int* n, const std::complex<float>* a,
             const int* lda, std::complex<float>* b, const int* ldb);
void zlacpy_(const char* uplo, const int* m, const int* n, const std::complex<double>* a,
             const int* lda, std::complex<double>* b, const int* ldb);

void sgesv_(const int* n, const int* nrhs, float* a, const int* lda, int* ipiv, float* b,
            const int* ldb, int* info);
void dgesv_(const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* b,
            const int* ldb, int* info);
void cgesv_(const int* n, const int* nrhs, std::complex<float>* a, const int* lda, int* ipiv,
            std::complex<float>* b, const int* ldb, int* info);
void zgesv_(const int* n, const int* nrhs, std::complex<double>* a, const int* lda, int* ipiv,
            std::complex<double>* b, const int* ldb, int* info);

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

void sgeev_(const char* job_vl, const char* job_vr, const int* n, float* a, const int* lda,
            float* wr, float* wi, float* vl, const int* ldvl, float* vr, const int* ldvr,
            float* work, const int* lwork, int* info);
void dgeev_(const char* job_vl, const char* job_vr, const int* n, double* a, const int* lda,
            double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr,
            double* work, const int* lwork, int* info);

void cgeev_(const char* job_vl, const char* job_vr, const int* n, std::complex<float>* a,
            const int* lda, std::complex<float>* w, std::complex<float>* vl, const int* ldvl,
            std::complex<float>* vr, const int* ldvr, std::complex<float>* work, const int* lwork,
            float* rwork, int* info);
void zgeev_(const char* job_vl, const char* job_vr, const int* n, std::complex<double>* a,
            const int* lda, std::complex<double>* w, std::complex<double>* vl, const int* ldvl,
            std::complex<double>* vr, const int* ldvr, std::complex<double>* work, const int* lwork,
            double* rwork, int* info);

int spotrf_(const char* uplo, const int* n, float* a, const int* lda, int* info);
int dpotrf_(const char* uplo, const int* n, double* a, const int* lda, int* info);
int cpotrf_(const char* uplo, const int* n, std::complex<float>* a, const int* lda, int* info);
int zpotrf_(const char* uplo, const int* n, std::complex<double>* a, const int* lda, int* info);

void spocon_(const char* uplo, const int* n, const float* a, const int* lda, const float* norm,
             float* rcond, float* work, int* iwork, int* info);
void dpocon_(const char* uplo, const int* n, const double* a, const int* lda, const double* norm,
             double* rcond, double* work, int* iwork, int* info);
void cpocon_(const char* uplo, const int* n, const std::complex<float>* a, const int* lda,
             const float* norm, float* rcond, std::complex<float>* work, float* rwork, int* info);
void zpocon_(const char* uplo, const int* n, const std::complex<double>* a, const int* lda,
             const double* norm, double* rcond, std::complex<double>* work, double* rwork, int* info);

void ssyevd_(const char* job_v, const char* uplo, const int* n, float* a, const int* lda, float* w,
             float* work, const int* lwork, int* iwork, const int* liwork, int* info);
void dsyevd_(const char* job_v, const char* uplo, const int* n, double* a, const int* lda,
             double* w, double* work, const int* lwork, int* iwork, const int* liwork, int* info);

void cheevd_(const char* job_v, const char* uplo, const int* n, std::complex<float>* a,
             const int* lda, float* w, std::complex<float>* work, const int* lwork, float* rwork,
             const int* lrwork, int* iwork, const int* liwork, int* info);
void zheevd_(const char* job_v, const char* uplo, const int* n, std::complex<double>* a,
             const int* lda, double* w, std::complex<double>* work, const int* lwork, double* rwork,
             const int* lrwork, int* iwork, const int* liwork, int* info);
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

inline void lacpy(const char* uplo, int m, int n, const float* a, int lda, float* b, int ldb) {
  slacpy_(uplo, &m, &n, a, &lda, b, &ldb);
}
inline void lacpy(const char* uplo, int m, int n, const double* a, int lda, double* b, int ldb) {
  dlacpy_(uplo, &m, &n, a, &lda, b, &ldb);
}
inline void lacpy(const char* uplo, int m, int n, const std::complex<float>* a, int lda,
                  std::complex<float>* b, int ldb) {
  clacpy_(uplo, &m, &n, a, &lda, b, &ldb);
}
inline void lacpy(const char* uplo, int m, int n, const std::complex<double>* a, int lda,
                  std::complex<double>* b, int ldb) {
  zlacpy_(uplo, &m, &n, a, &lda, b, &ldb);
}

inline void gesv(int n, int nrhs, float* a, int lda, int* ipiv, float* b, int ldb) {
  int info = 0;
  sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  checkLapackInfo(info);
}
inline void gesv(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb) {
  int info = 0;
  dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  checkLapackInfo(info);
}
inline void gesv(int n, int nrhs, std::complex<float>* a, int lda, int* ipiv,
                 std::complex<float>* b, int ldb) {
  int info = 0;
  cgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  checkLapackInfo(info);
}
inline void gesv(int n, int nrhs, std::complex<double>* a, int lda, int* ipiv,
                 std::complex<double>* b, int ldb) {
  int info = 0;
  zgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  checkLapackInfo(info);
}

inline void getrf(int m, int n, float* a, int lda, int* ipiv) {
  int info = 0;
  sgetrf_(&m, &n, a, &lda, ipiv, &info);
  checkLapackInfo(info);
}
inline void getrf(int m, int n, double* a, int lda, int* ipiv) {
  int info = 0;
  dgetrf_(&m, &n, a, &lda, ipiv, &info);
  checkLapackInfo(info);
}
inline void getrf(int m, int n, std::complex<float>* a, int lda, int* ipiv) {
  int info = 0;
  cgetrf_(&m, &n, a, &lda, ipiv, &info);
  checkLapackInfo(info);
}
inline void getrf(int m, int n, std::complex<double>* a, int lda, int* ipiv) {
  int info = 0;
  zgetrf_(&m, &n, a, &lda, ipiv, &info);
  checkLapackInfo(info);
}

inline void getri(int n, float* a, int lda, int* ipiv, float* work, int lwork) {
  int info = 0;
  sgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
  checkLapackInfo(info);
}
inline void getri(int n, double* a, int lda, int* ipiv, double* work, int lwork) {
  int info = 0;
  dgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
  checkLapackInfo(info);
}
inline void getri(int n, std::complex<float>* a, int lda, int* ipiv, std::complex<float>* work,
                  int lwork) {
  int info = 0;
  cgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
  checkLapackInfo(info);
}
inline void getri(int n, std::complex<double>* a, int lda, int* ipiv, std::complex<double>* work,
                  int lwork) {
  int info = 0;
  zgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
  checkLapackInfo(info);
}

inline void geev(const char* job_vl, const char* job_vr, int n, float* a, int lda, float* wr,
                 float* wi, float* vl, int ldvl, float* vr, int ldvr, float* work, int lwork) {
  int info = 0;
  sgeev_(job_vl, job_vr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
  checkLapackInfo(info);
}
inline void geev(const char* job_vl, const char* job_vr, int n, double* a, int lda, double* wr,
                 double* wi, double* vl, int ldvl, double* vr, int ldvr, double* work, int lwork) {
  int info = 0;
  dgeev_(job_vl, job_vr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
  checkLapackInfo(info);
}

inline void geev(const char* job_vl, const char* job_vr, int n, std::complex<float>* a, int lda,
                 std::complex<float>* w, std::complex<float>* vl, int ldvl, std::complex<float>* vr,
                 int ldvr, std::complex<float>* work, int lwork, float* rwork) {
  int info = 0;
  cgeev_(job_vl, job_vr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
  checkLapackInfo(info);
}
inline void geev(const char* job_vl, const char* job_vr, int n, std::complex<double>* a, int lda,
                 std::complex<double>* w, std::complex<double>* vl, int ldvl,
                 std::complex<double>* vr, int ldvr, std::complex<double>* work, int lwork,
                 double* rwork) {
  int info = 0;
  zgeev_(job_vl, job_vr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
  checkLapackInfo(info);
}

inline int potrf(const char* uplo, int n, float* a, int lda) {
  int info = 0;
  int res = spotrf_(uplo, &n, a, &lda, &info);
  checkLapackInfo(info);
  return res;
}
inline int potrf(const char* uplo, int n, double* a, int lda) {
  int info = 0;
  int res = dpotrf_(uplo, &n, a, &lda, &info);
  checkLapackInfo(info);
  return res;
}
inline int potrf(const char* uplo, int n, std::complex<float>* a, int lda) {
  int info = 0;
  int res = cpotrf_(uplo, &n, a, &lda, &info);
  checkLapackInfo(info);
  return res;
}
inline int potrf(const char* uplo, int n, std::complex<double>* a, int lda) {
  int info = 0;
  int res = zpotrf_(uplo, &n, a, &lda, &info);
  checkLapackInfo(info);
  return res;
}

inline float pocon(const char* uplo, int n, const float* a, int lda, float norm, float* work,
                   int* iwork) {
  float cond;
  int info = 0;
  spocon_(uplo, &n, a, &lda, &norm, &cond, work, iwork, &info);
  checkLapackInfo(info);
  return cond;
}
inline double pocon(const char* uplo, int n, const double* a, int lda, double norm, double* work,
                    int* iwork) {
  double cond;
  int info = 0;
  dpocon_(uplo, &n, a, &lda, &norm, &cond, work, iwork, &info);
  checkLapackInfo(info);
  return cond;
}
inline float pocon(const char* uplo, int n, const std::complex<float>* a, int lda, float norm,
                   std::complex<float>* work, float* rwork) {
  float cond;
  int info = 0;
  cpocon_(uplo, &n, a, &lda, &norm, &cond, work, rwork, &info);
  checkLapackInfo(info);
  return cond;
}
inline double pocon(const char* uplo, int n, const std::complex<double>* a, int lda, double norm,
                    std::complex<double>* work, double* rwork) {
  double cond;
  int info = 0;
  zpocon_(uplo, &n, a, &lda, &norm, &cond, work, rwork, &info);
  checkLapackInfo(info);
  return cond;
}

inline void syevd(const char* job_v, const char* uplo, int n, float* a, int lda, float* w,
                  float* work, int lwork, int* iwork, int liwork) {
  int info = 0;
  ssyevd_(job_v, uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
  checkLapackInfo(info);
}
inline void syevd(const char* job_v, const char* uplo, int n, double* a, int lda, double* w,
                  double* work, int lwork, int* iwork, int liwork) {
  int info = 0;
  dsyevd_(job_v, uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
  checkLapackInfo(info);
}

inline void heevd(const char* job_v, const char* uplo, int n, std::complex<float>* a, int lda,
                  float* w, std::complex<float>* work, int lwork, float* rwork, int lrwork,
                  int* iwork, int liwork) {
  int info = 0;
  cheevd_(job_v, uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  checkLapackInfo(info);
}
inline void heevd(const char* job_v, const char* uplo, int n, std::complex<double>* a, int lda,
                  double* w, std::complex<double>* work, int lwork, double* rwork, int lrwork,
                  int* iwork, int liwork) {
  int info = 0;
  zheevd_(job_v, uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  checkLapackInfo(info);
}

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_LAPACK_HPP
