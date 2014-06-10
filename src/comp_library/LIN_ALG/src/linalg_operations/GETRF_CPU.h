//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GETRF_CPU_H
#define LINALG_GETRF_CPU_H

namespace LIN_ALG {

  template<>
  class GETRF<CPU>
  {
    public:

      template<typename scalartype>
      static void execute(matrix<scalartype, CPU>& A, int* IPIV){
	  
	  int M = A.get_current_size().first;
	  int N = A.get_current_size().second;
	  
	  int LDA = A.get_global_size().first;
	  
	  int INFO=0;
	  execute(M, N, A.get_ptr(), LDA, IPIV, INFO);

	  if(INFO!=0)
	    throw std::logic_error(__FUNCTION__);
	}

      static void execute( int M, int N, float* A, int LDA, int* IPIV, int INFO){
	  LAPACK::sgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
      }

      static void execute( int M, int N, double* A, int LDA, int* IPIV, int INFO){
	  LAPACK::dgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
      }

      static void execute( int M, int N, std::complex<float>* A, int LDA, int* IPIV, int INFO){
	  LAPACK::cgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
      }

      static void execute( int M, int N, std::complex<double>* A, int LDA, int* IPIV, int INFO){
	  LAPACK::zgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
      }
  };

}

#endif
