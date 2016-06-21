//-*-C++-*-

/*
 * invert_plan.h
 *
 *      Author: peter staar
 */

#ifndef INVERT_PLAN_H_
#define INVERT_PLAN_H_

#include <cassert>
#include <complex>

#include "comp_library/blas_lapack_plans/LAPACK/LAPACK_C_wrappers.h"

/*
*******************************************
***       REAL MATRICES                 ***
*******************************************
*/

template <typename real_scalartype>
class invert_plan {
public:
  invert_plan(int n);
  invert_plan(int n, int lda);
  ~invert_plan();

  void execute_plan();

private:
  void reset_inverted_matrix();

public:
  int N;
  int LDA;
  int LWORK;
  int INFO;

  real_scalartype* Matrix;
  real_scalartype* inverted_matrix;

private:
  int* IPIV;
};

template <typename real_scalartype>
invert_plan<real_scalartype>::invert_plan(int n) : N(n), LDA(N) {
  Matrix = new real_scalartype[N * N];
  memset(Matrix, 0, sizeof(real_scalartype) * N * N);

  inverted_matrix = new real_scalartype[N * N];
  memset(inverted_matrix, 0, sizeof(real_scalartype) * N * N);

  for (int i = 0; i < N; i++)
    inverted_matrix[i + N * i] = real_scalartype(1);

  IPIV = new int[N];
}

template <typename real_scalartype>
invert_plan<real_scalartype>::invert_plan(int n, int lda) : N(n), LDA(lda) {
  Matrix = new real_scalartype[LDA * LDA];
  memset(Matrix, 0, sizeof(real_scalartype) * LDA * LDA);

  inverted_matrix = new real_scalartype[LDA * LDA];
  memset(inverted_matrix, 0, sizeof(real_scalartype) * LDA * LDA);

  for (int i = 0; i < LDA; i++)
    inverted_matrix[i + LDA * i] = real_scalartype(1);

  IPIV = new int[N];
}

template <typename real_scalartype>
invert_plan<real_scalartype>::~invert_plan() {
  delete[] Matrix;
  delete[] inverted_matrix;

  delete[] IPIV;
}

template <typename real_scalartype>
void invert_plan<real_scalartype>::reset_inverted_matrix() {
  int lda = std::max(N, LDA);

  memset(inverted_matrix, 0, sizeof(real_scalartype) * lda * lda);

  for (int i = 0; i < lda; i++)
    inverted_matrix[i + lda * i] = real_scalartype(1.);
}

template <typename real_scalartype>
void invert_plan<real_scalartype>::execute_plan() {
  throw std::logic_error(__PRETTY_FUNCTION__);
  assert(false);
}

template <>
void invert_plan<float>::execute_plan() {
  reset_inverted_matrix();
  LAPACK::sgesv_(&N, &N, Matrix, &LDA, IPIV, inverted_matrix, &LDA, &INFO);

  if (INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

template <>
void invert_plan<double>::execute_plan() {
  reset_inverted_matrix();
  LAPACK::dgesv_(&N, &N, Matrix, &LDA, IPIV, inverted_matrix, &LDA, &INFO);

  if (INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

/*
**********************************************
***       COMPLEX  MATRICES                ***
**********************************************
*/

template <typename real_scalartype>
class invert_plan<std::complex<real_scalartype>> {
public:
  invert_plan(int n);
  invert_plan(int n, int LDA);

  ~invert_plan();

  void execute_plan();

private:
  void reset_inverted_matrix();

public:
  int N;
  int LDA;
  int LWORK;
  int INFO;

  std::complex<real_scalartype>* Matrix;
  std::complex<real_scalartype>* inverted_matrix;

private:
  int* IPIV;

  //   const static bool check = true;
};

template <typename real_scalartype>
invert_plan<std::complex<real_scalartype>>::invert_plan(int n) : N(n), LDA(N) {
  Matrix = new std::complex<real_scalartype>[N * N];
  memset(Matrix, 0, sizeof(std::complex<real_scalartype>) * N * N);

  inverted_matrix = new std::complex<real_scalartype>[N * N];
  memset(inverted_matrix, 0, sizeof(std::complex<real_scalartype>) * N * N);

  for (int i = 0; i < N; i++)
    inverted_matrix[i + N * i] = std::complex<real_scalartype>(1.);

  IPIV = new int[N];
}

template <typename real_scalartype>
invert_plan<std::complex<real_scalartype>>::invert_plan(int n, int lda) : N(n), LDA(lda) {
  Matrix = new std::complex<real_scalartype>[LDA * LDA];
  memset(Matrix, 0, sizeof(std::complex<real_scalartype>) * LDA * LDA);

  inverted_matrix = new std::complex<real_scalartype>[LDA * LDA];
  memset(inverted_matrix, 0, sizeof(std::complex<real_scalartype>) * LDA * LDA);

  for (int i = 0; i < LDA; i++)
    inverted_matrix[i + LDA * i] = std::complex<real_scalartype>(1.);

  IPIV = new int[LDA];
}

template <typename real_scalartype>
invert_plan<std::complex<real_scalartype>>::~invert_plan() {
  delete[] Matrix;
  delete[] inverted_matrix;
  delete[] IPIV;
}

template <typename real_scalartype>
void invert_plan<std::complex<real_scalartype>>::reset_inverted_matrix() {
  int lda = std::max(N, LDA);

  memset(inverted_matrix, 0, sizeof(std::complex<real_scalartype>) * lda * lda);

  for (int i = 0; i < lda; i++)
    inverted_matrix[i + lda * i] = std::complex<real_scalartype>(1.);
}

template <typename real_scalartype>
void invert_plan<std::complex<real_scalartype>>::execute_plan() {
  throw std::logic_error(__PRETTY_FUNCTION__);
}

template <>
void invert_plan<std::complex<float>>::execute_plan() {
  reset_inverted_matrix();

  LAPACK::cgesv_(&N, &N, Matrix, &LDA, IPIV, inverted_matrix, &LDA, &INFO);

  assert(INFO == 0);
}

template <>
void invert_plan<std::complex<double>>::execute_plan() {
  reset_inverted_matrix();

  LAPACK::zgesv_(&N, &N, Matrix, &LDA, IPIV, inverted_matrix, &LDA, &INFO);

  assert(INFO == 0);
}

#endif
