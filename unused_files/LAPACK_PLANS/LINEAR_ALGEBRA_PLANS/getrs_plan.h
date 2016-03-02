//-*-C++-*-

#ifndef ANY_GETRS_PLAN_H
#define ANY_GETRS_PLAN_H

/*!
 *   \author Peter Staar
 */
template<typename scalartype>
class getrs_plan
{
public:
    
  getrs_plan(int n, int nrhs, int lda, int ldb);
  ~getrs_plan();
    
  inline void execute_plan(LINEAR_ALGEBRA_LIBRARY_TYPE LAL_t = LAPACK_LIBRARY);
    
private:

  template<typename whatever_t>
  inline void set_values(whatever_t& whatever_ref);
    
public:
    
  char TRANS;
  int  N;
  int  NRHS;
  int  LDA;
  int  LDB;
  int  INFO;

  scalartype*  Matrix_A;
  scalartype*  Matrix_B;
  int*             IPIV;

};

template<typename scalartype>
getrs_plan<scalartype>::getrs_plan(int n, int nrhs, int lda, int ldb):
  TRANS('N'),
  N(n),
  NRHS(nrhs),
  LDA(lda),
  LDB(ldb)
{
  IPIV = new int[N];
  
  for(int i=0; i<N; i++)
    IPIV[i] = i+1;
}

template<typename scalartype>
getrs_plan<scalartype>::~getrs_plan()
{
  delete [] IPIV;
}

template<typename scalartype>
void getrs_plan<scalartype>::execute_plan(LINEAR_ALGEBRA_LIBRARY_TYPE LAL_t /*= LAPACK_LIBRARY*/)
{
  switch(LAL_t)
    {
    case LAPACK_LIBRARY:
      {
	lapack_getrs(TRANS, N, NRHS, Matrix_A, LDA, IPIV, Matrix_B, LDB, INFO);
      }
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
}

#endif
