//-*-C++-*-

/*
 * matrix_vector_multiplication_plan.h
 *
 *      Author: peter staar
 */

#ifndef BANDED_MATRIX_VECTOR_PLAN_H_
#define BANDED_MATRIX_VECTOR_PLAN_H_

template<typename scalartype>
class gbmv_plan
{

 public:

  gbmv_plan(int n);
  gbmv_plan(int m, int n);
  ~gbmv_plan();

  void execute_plan();

public:

  char TRANS;
  int M; // --> matrix-rows 
  int N; // --> matrix-columns 
  int KL; // number of sub-diagonals --> 0 for diagonal matrix
  int KU; // number of super-diagonals --> 0 for diagonal matrix
 
  scalartype   alpha;
  scalartype*  A; // -> matrix
  int LDA;

  scalartype* X; // -> source-vector
  int INCX;

  scalartype   beta;
  scalartype* Y; // -> target-vector
  int INCY;
};

template<typename scalartype>
gbmv_plan<scalartype>::gbmv_plan(int n):
  TRANS('N'),
  M(n),
  N(n),
  KL(0),
  KU(0),
  alpha(1.),
  LDA(1),
  INCX(1),
  beta(0.),
  INCY(1)
{}

template<typename scalartype>
gbmv_plan<scalartype>::gbmv_plan(int m, int n):
  TRANS('N'),
  M(m),
  N(n),
  KL(0),
  KU(0),
  alpha(1.),
  LDA(1),
  INCX(1),
  beta(0.),
  INCY(1)
{}

template<typename scalartype>
gbmv_plan<scalartype>::~gbmv_plan()
{}

template<typename scalartype>
void gbmv_plan<scalartype>::execute_plan()
{
  cout << __PRETTY_FUNCTION__ << endl;
  assert(false);
}

template<>
void gbmv_plan<float>::execute_plan()
{
  BLAS::sgbmv_(&TRANS, &M, &N, &KL, &KU, &alpha, A, &LDA, X, &INCX, &beta, Y, &INCY);
}

template<>
void gbmv_plan<double>::execute_plan()
{
  BLAS::dgbmv_(&TRANS, &M, &N, &KL, &KU, &alpha, A, &LDA, X, &INCX, &beta, Y, &INCY);
}

template<>
void gbmv_plan<std::complex<float> >::execute_plan()
{
  BLAS::cgbmv_(&TRANS, &M, &N, &KL, &KU, &alpha, A, &LDA, X, &INCX, &beta, Y, &INCY);
}

template<>
void gbmv_plan<std::complex<double> >::execute_plan()
{
  BLAS::zgbmv_(&TRANS, &M, &N, &KL, &KU, &alpha, A, &LDA, X, &INCX, &beta, Y, &INCY);
}

#endif
