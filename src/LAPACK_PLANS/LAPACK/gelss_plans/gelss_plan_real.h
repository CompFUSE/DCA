//-*-C++-*-

#ifndef LEAST_SQUARE_PLAN_REAL_H
#define LEAST_SQUARE_PLAN_REAL_H

/*!
 *   \author: peter staar
 */
template<typename scalartype>
class least_square_plan
{

 public:

  least_square_plan(int m, int n, int nrhs, scalartype* a, scalartype* b);
  ~least_square_plan();

  int execute_plan();

 public:

  int M;
  int N;

  int NRHS;

  scalartype* A;
  int LDA;

  scalartype* B;
  int LDB;

  scalartype* S;

  scalartype RCOND;
  int RANK;

  scalartype* WORK;
  int LWORK;

  int INFO;

};


template<typename scalartype>
least_square_plan<scalartype>::least_square_plan(int m, int n, int nrhs, scalartype* a, scalartype* b):
  M(m),
  N(n),
  NRHS(nrhs),

  A(a),
  LDA(m),

  B(b),
  LDB(m),

  RCOND(-1),

  LWORK(10*(3*std::min(M,N) + std::max(std::max(2*std::min(M,N), std::max(M,N)), NRHS)))
{
  //cout << __FUNCTION__ << endl;

  S     = new scalartype[min(N,M)];
  WORK  = new scalartype[LWORK];
}

template<typename scalartype>
least_square_plan<scalartype>::~least_square_plan()
{
  delete [] S;
  delete [] WORK;
}

template<typename scalartype>
int least_square_plan<scalartype>::execute_plan()
{
  throw std::logic_error(__FUNCTION__);
}

template<>
int least_square_plan<float>::execute_plan()
{
  LAPACK::sgelss_( &M, &N, &NRHS, A, &LDA, B, &LDB, S, &RCOND, &RANK, WORK, &LWORK, &INFO);
  return INFO;
}

template<>
int least_square_plan<double>::execute_plan()
{
  LAPACK::dgelss_( &M, &N, &NRHS, A, &LDA, B, &LDB, S, &RCOND, &RANK, WORK, &LWORK, &INFO);
  return INFO;
}

#endif
