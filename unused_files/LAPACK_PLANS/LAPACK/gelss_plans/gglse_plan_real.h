//-*-C++-*-

#ifndef GENERALIZED_LEAST_SQUARE_PLAN_REAL_H
#define GENERALIZED_LEAST_SQUARE_PLAN_REAL_H

/*!
 *   \author: peter staar
 */
template<typename scalartype>
class generalized_least_square_plan
{

 public:

  generalized_least_square_plan(int m, int n, int p);
  ~generalized_least_square_plan();

  int execute_plan();

 public:

  int M;
  int N;
  int P;

  scalartype* A;
  int LDA;

  scalartype* B;
  int LDB;

  scalartype* C;
  scalartype* D;

  scalartype* X;

  scalartype* WORK;
  int LWORK;

  int INFO;
};

template<typename scalartype>
generalized_least_square_plan<scalartype>::generalized_least_square_plan(int m, int n, int p):
  M(m),
  N(n),
  P(p),

  LDA(m),
  LDB(p),

  LWORK(std::max(1,M+N+P))
{
  WORK  = new scalartype[LWORK];
  LWORK = -1;

  execute_plan();

  LWORK = WORK[0];
  delete [] WORK;
  WORK  = new scalartype[LWORK];
}

template<typename scalartype>
generalized_least_square_plan<scalartype>::~generalized_least_square_plan()
{
  delete [] WORK;
}

template<typename scalartype>
int generalized_least_square_plan<scalartype>::execute_plan()
{
  throw std::logic_error(__FUNCTION__);
}

template<>
int generalized_least_square_plan<float>::execute_plan()
{
  LAPACK::sgglse_( &M, &N, &P, A, &LDA, B, &LDB, C, D, X, WORK, &LWORK, &INFO);
  return INFO;
}

template<>
int generalized_least_square_plan<double>::execute_plan()
{
  LAPACK::dgglse_( &M, &N, &P, A, &LDA, B, &LDB, C, D, X, WORK, &LWORK, &INFO);
  return INFO;
}



#endif
