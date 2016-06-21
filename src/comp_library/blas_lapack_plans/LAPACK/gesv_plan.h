//-*-C++-*-

/*
 * invert_plan.h
 *
 *      Author: peter staar
 */

#ifndef SOLVE_PLAN_H_
#define SOLVE_PLAN_H_

#include <complex>
#include <stdexcept>

#include "comp_library/blas_lapack_plans/LAPACK/LAPACK_C_wrappers.h"

template <typename real_scalartype>
class solve_plan {
public:
  solve_plan(int n, int nrhs);
  solve_plan(int n, int lda, int nrhs);
  ~solve_plan();

  int execute_plan();

public:
  int N;
  int NRHS;
  int LDA;
  int LWORK;
  int INFO;

  real_scalartype* matrix;
  real_scalartype* solved_matrix;

private:
  int* IPIV;
};

template <typename real_scalartype>
solve_plan<real_scalartype>::solve_plan(int n, int nrhs) : N(n), NRHS(nrhs), LDA(N) {
  matrix = new real_scalartype[N * N];
  memset(matrix, 0, sizeof(real_scalartype) * N * N);

  solved_matrix = new real_scalartype[N * NRHS];
  memset(solved_matrix, 0, sizeof(real_scalartype) * N * NRHS);

  IPIV = new int[N];
}

template <typename real_scalartype>
solve_plan<real_scalartype>::solve_plan(int n, int lda, int /*nrhs*/) : N(n), LDA(lda) {
  matrix = new real_scalartype[LDA * LDA];
  memset(matrix, 0, sizeof(real_scalartype) * LDA * LDA);

  solved_matrix = new real_scalartype[LDA * NRHS];
  memset(solved_matrix, 0, sizeof(real_scalartype) * LDA * NRHS);

  IPIV = new int[N];
}

template <typename real_scalartype>
solve_plan<real_scalartype>::~solve_plan() {
  delete[] matrix;
  delete[] solved_matrix;

  delete[] IPIV;
}

template <typename real_scalartype>
int solve_plan<real_scalartype>::execute_plan() {
  throw std::logic_error(__PRETTY_FUNCTION__);
}

template <>
int solve_plan<float>::execute_plan() {
  LAPACK::sgesv_(&N, &NRHS, matrix, &LDA, IPIV, solved_matrix, &LDA, &INFO);
  return INFO;
}

template <>
int solve_plan<double>::execute_plan() {
  LAPACK::dgesv_(&N, &NRHS, matrix, &LDA, IPIV, solved_matrix, &LDA, &INFO);
  return INFO;
}

template <>
int solve_plan<std::complex<float>>::execute_plan() {
  LAPACK::cgesv_(&N, &NRHS, matrix, &LDA, IPIV, solved_matrix, &LDA, &INFO);
  return INFO;
}

template <>
int solve_plan<std::complex<double>>::execute_plan() {
  LAPACK::zgesv_(&N, &NRHS, matrix, &LDA, IPIV, solved_matrix, &LDA, &INFO);
  return INFO;
}

#endif
