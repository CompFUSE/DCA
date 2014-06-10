//-*-C++-*-

#ifndef MATRIX_MATRIX_PLAN_H
#define MATRIX_MATRIX_PLAN_H

/*!
 *  \author P. Staar
 */
template<typename scalartype>
class gemm_plan
{
public:
    
  gemm_plan();
    
  gemm_plan(int n);
    
  gemm_plan(int m, int k, int n);
    
  ~gemm_plan();
    
  inline void execute_plan(LINEAR_ALGEBRA_LIBRARY_TYPE LAL_t = BLAS_LIBRARY);
    
public:
    
  char TRANSA;
  char TRANSB;
    
  int M;
  int N;
  int K;
  
  int LDA;
  int LDB;
  int LDC;
    
  scalartype alpha;
  scalartype beta;
    
  scalartype* A;
  scalartype* B;
  scalartype* C;
};

template<typename scalartype>
gemm_plan<scalartype>::gemm_plan():
  TRANSA('N'),
  TRANSB('N'),

  alpha(1.),
  beta(0.)
{}

template<typename scalartype>
gemm_plan<scalartype>::gemm_plan(int n):
  TRANSA('N'),
  TRANSB('N'),

  M(n),
  N(n),
  K(n),

  LDA(n),
  LDB(n),
  LDC(n),

  alpha(1.),
  beta(0.)
{}

template<typename scalartype>
gemm_plan<scalartype>::gemm_plan(int m, int k, int n):
  TRANSA('N'),
  TRANSB('N'),

  M(m),
  N(n),
  K(k),

  LDA(m),
  LDB(k),
  LDC(m),

  alpha(1.),
  beta(0.)
{}

template<typename scalartype>
gemm_plan<scalartype>::~gemm_plan()
{}

template<typename scalartype>
void gemm_plan<scalartype>::execute_plan(LINEAR_ALGEBRA_LIBRARY_TYPE LAL_t /*= BLAS_LIBRARY*/)
{
  switch(LAL_t)
    {
    case BLAS_LIBRARY:
      {
	blas_gemm_t<scalartype>(TRANSA, TRANSB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
      }
      break;
      
    default:
      throw std::logic_error(__FUNCTION__);
    }
}


#endif
