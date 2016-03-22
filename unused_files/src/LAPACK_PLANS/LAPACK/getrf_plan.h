//-*-C++-*-

#ifndef LU_DECOMPOSITION_PLAN_H
#define LU_DECOMPOSITION_PLAN_H

/*!
 *  \class   lu_decomposition_plan
 *  \author  Peter Staar
 *  \version 1.0
 */
template<typename scalartype>
class LU_decomposition_plan
{
public:

  LU_decomposition_plan(int m, int n):
    M(m),
    N(n),

    A(NULL),
    LDA(m),

    IPIV(NULL)
  {
    A    = new scalartype[LDA*N];
    IPIV = new int[std::min(M,N)]; 
  }

  ~LU_decomposition_plan()
  {
    delete [] A;
    delete [] IPIV;
  }
  
  int execute_plan();

  void get_row_swap_indices(int* T)
  {
    for(int i=0; i<M; ++i)
      if(i<std::min(M,N))
	T[i] = IPIV[i]-1; 
      else
	T[i] = i; 

    for(int i=0; i<M; ++i)
      if(i==T[i])
	for(int j=0; j<i; ++j)
	  if(i==T[j])
	    T[i]=j;    
  }

public:

  int M;
  int N;

  scalartype* A;
  int         LDA;

  int* IPIV;
};

template<typename scalartype>
int LU_decomposition_plan<scalartype>::execute_plan()
{
  throw std::logic_error(__FUNCTION__);
}

template<>
int LU_decomposition_plan<double>::execute_plan()
{
  int INFO;
  LAPACK::dgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
  return INFO;
}

template<>
int LU_decomposition_plan<std::complex<double> >::execute_plan()
{
  int INFO;
  LAPACK::zgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
  return INFO;
}

#endif
