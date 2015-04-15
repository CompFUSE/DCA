//-*-C++-*-

#ifndef MATRIX_VECTOR_PLAN_H_
#define MATRIX_VECTOR_PLAN_H_

/*!
 *  \author: peter staar
 */
template<typename scalartype>
class gemv_plan
{

 public:

  gemv_plan(int m, int n);
  ~gemv_plan();

  void execute_plan();

 public:

  char TRANS;
  int M; // --> matrix-rows 
  int N; // --> matrix-columns 
 
  int LDA;
  int INCX;
  int INCY;

  scalartype   alpha;
  scalartype   beta;

  scalartype*  matrix;
  scalartype*  vector_source;
  scalartype*  vector_target;
};

template<typename scalartype>
gemv_plan<scalartype>::gemv_plan(int m, int n):
  TRANS('n'),
  M(m),
  N(n),
  LDA(m),
  INCX(1),
  INCY(1),
  alpha(1.),
  beta(0.)
{}

template<typename scalartype>
gemv_plan<scalartype>::~gemv_plan()
{}

template<typename scalartype>
void gemv_plan<scalartype>::execute_plan()
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  assert(false);
}

template<>
void gemv_plan<float>::execute_plan()
{
  BLAS::sgemv_(&TRANS, &M, &N, &alpha, matrix, &LDA, vector_source, &INCX, &beta, vector_target , &INCY);
}

template<>
void gemv_plan<double>::execute_plan()
{
  BLAS::dgemv_(&TRANS, &M, &N, &alpha, matrix, &LDA, vector_source, &INCX, &beta, vector_target , &INCY);
}

template<>
void gemv_plan<std::complex<float> >::execute_plan()
{
  BLAS::cgemv_(&TRANS, &M, &N, &alpha, matrix, &LDA, vector_source, &INCX, &beta, vector_target , &INCY);
}

template<>
void gemv_plan<std::complex<double> >::execute_plan()
{
  BLAS::zgemv_(&TRANS, &M, &N, &alpha, matrix, &LDA, vector_source, &INCX, &beta, vector_target , &INCY);
}

#endif
