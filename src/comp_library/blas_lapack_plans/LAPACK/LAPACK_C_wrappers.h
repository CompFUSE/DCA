
//-*-C++-*-
// ****************************************************************************
// * C++ wrapper for LAPACK                                                   *
// *                                                                          *
// * Thomas Schulthess, ORNL, October 1999                                    *
// * Peter Staar, ETHZ, December 2011                                         *
// ****************************************************************************

#ifndef LAPACK_C_WRAPPERS
#define LAPACK_C_WRAPPERS

/** \file LAPACK
 *  \author Thomas C. Schulthess, Michael S. Summers, Peter Staar
 */

#include <complex>

namespace LAPACK {

// ============================================================================


extern "C" float slamch_(char* JOBZ);
extern "C" double dlamch_(char* JOBZ);
// solve
extern "C" void sgesv_(int*, int*, float*, int*, int*, float*, int*, int*);
extern "C" void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
extern "C" void cgesv_(int*, int*, std::complex<float>*, int*, int*, std::complex<float>*, int*,
                       int*);
extern "C" void zgesv_(int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*,
                       int*);

// eigenvalue decomposition

extern "C" void sgeev_(char* JOBVL, char* JOBVR, int* N, float* A, int* LDA, float* WR, float* WI,
                       float* VL, int* LDVL, float* VR, int* LDVR, float* WORK, int* LWORK,
                       int* INFO);
extern "C" void dgeev_(char* JOBVL, char* JOBVR, int* N, double* A, int* LDA, double* WR,
                       double* WI, double* VL, int* LDVL, double* VR, int* LDVR, double* WORK,
                       int* LWORK, int* INFO);

extern "C" void cgeev_(char* JOBVL, char* JOBVR, int* N, std::complex<float>* A, int* LDA,
                       std::complex<float>* W, std::complex<float>* VL, int* LDVL,
                       std::complex<float>* VR, int* LDVR, std::complex<float>* WORK, int* LWORK,
                       float* RWORK, int* INFO);
extern "C" void zgeev_(char* JOBVL, char* JOBVR, int* N, std::complex<double>* A, int* LDA,
                       std::complex<double>* W, std::complex<double>* VL, int* LDVL,
                       std::complex<double>* VR, int* LDVR, std::complex<double>* WORK, int* LWORK,
                       double* RWORK, int* INFO);

extern "C" void ssyev_(char* JOBZ, char* UPLO, int* N, float* Matrix, int* LDA, float* eigenvalues,
                       float* WORK, int* LWORK, int* INFO);
extern "C" void dsyev_(char* JOBZ, char* UPLO, int* N, double* Matrix, int* LDA,
                       double* eigenvalues, double* WORK, int* LWORK, int* INFO);
extern "C" void cheev_(char* jobz, char* uplo, int* n, std::complex<float>* a, int* lda, float* w,
                       std::complex<float>* work, int* lwork, float* rwork, int* info);
extern "C" void zheev_(char* jobz, char* uplo, int* n, std::complex<double>* a, int* lda, double* w,
                       std::complex<double>* work, int* lwork, double* rwork, int* info);
extern "C" void ssyevd_(char* JOBZ, char* UPLO, int* N, float* Matrix, int* LDA, float* eigenvalues,
                        float* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
extern "C" void dsyevd_(char* JOBZ, char* UPLO, int* N, double* Matrix, int* LDA, double* eigenvalues,
                        double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

extern "C" void cheevd_(char* jobz, char* uplo, int* n, std::complex<float>* a, int* lda, float* w,
                        std::complex<float>* work, int* lwork, float* rwork, int* LRWORK,
                        int* IWORK, int* LIWORK, int* info);
extern "C" void zheevd_(char* jobz, char* uplo, int* n, std::complex<double>* a, int* lda,
                        double* w, std::complex<double>* work, int* lwork, double* rwork,
                        int* LRWORK, int* IWORK, int* LIWORK, int* info);
extern "C" void ssyevx_(char* JOBZ, char* RANGE, char* UPLO, int* N, float* A, int* LDA, float* VL,
                        float* VU, int* IL, int* UL, float* ABSTOL, int* M, float* W, float* Z,
                        int* LDZ, float* WORK, int* LWORK, int* IWORK, int* IFAIL, int* INFO);

extern "C" void dsyevx_(char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA, double* VL,
                        double* VU, int* IL, int* UL, double* ABSTOL, int* M, double* W, double* Z,
                        int* LDZ, double* WORK, int* LWORK, int* IWORK, int* IFAIL, int* INFO);

extern "C" void cheevx_(char* JOBZ, char* RANGE, char* UPLO, int* N, std::complex<float>* A,
                        int* LDA, float* VL, float* VU, int* IL, int* UL, float* ABSTOL, int* M,
                        float* W, std::complex<float>* Z, int* LDZ, std::complex<float>* WORK,
                        int* LWORK, float* RWORK, int* LRWORK, int* IWORK, int* IFAIL, int* INFO);

extern "C" void zheevx_(char* JOBZ, char* RANGE, char* UPLO, int* N, std::complex<double>* A,
                        int* LDA, double* VL, double* VU, int* IL, int* UL, double* ABSTOL, int* M,
                        double* W, std::complex<double>* Z, int* LDZ, std::complex<double>* WORK,
                        int* LWORK, double* RWORK, int* LRWORK, int* IWORK, int* IFAIL, int* INFO);
extern "C" void ssyevr_(char* JOBZ, char* RANGE, char* UPLO, int* N, float* A, int* LDA, float* VL,
                        float* VU, int* IL, int* UL, float* ABSTOL, int* M, float* W, float* Z,
                        int* LDZ, int* ISUPPZ, float* WORK, int* LWORK, int* IWORK, int* LIWORK,
                        int* INFO);

extern "C" void dsyevr_(char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA,
                        double* VL, double* VU, int* IL, int* UL, double* ABSTOL, int* M, double* W,
                        double* Z, int* LDZ, int* ISUPPZ, double* WORK, int* LWORK, int* IWORK,
                        int* LIWORK, int* INFO);

extern "C" void cheevr_(char* JOBZ, char* RANGE, char* UPLO, int* N, std::complex<float>* A, int* LDA,
                        float* VL, float* VU, int* IL, int* UL, float* ABSTOL, int* M, float* W,
                        std::complex<float>* Z, int* LDZ, int* ISUPPZ, std::complex<float>* WORK,
                        int* LWORK, float* RWORK, int* LRWORK, int* IWORK, int* LIWORK, int* INFO);

extern "C" void zheevr_(char* JOBZ, char* RANGE, char* UPLO, int* N, std::complex<double>* A,
                        int* LDA, double* VL, double* VU, int* IL, int* UL, double* ABSTOL, int* M,
                        double* W, std::complex<double>* Z, int* LDZ, int* ISUPPZ,
                        std::complex<double>* WORK, int* LWORK, double* RWORK, int* LRWORK,
                        int* IWORK, int* LIWORK, int* INFO);
} /* namespace LAPACK */

#endif  // PSIMAG_LAPACK_H
