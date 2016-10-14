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

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_LAPACK_HPP
