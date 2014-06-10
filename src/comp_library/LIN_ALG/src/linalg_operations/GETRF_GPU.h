//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GETRF_GPU_H
#define LINALG_GETRF_GPU_H

namespace LIN_ALG {

  namespace GPU_KERNELS_GETRF {

    void sgetrf(int M, int N, float* A, int LDA, int* IPIV, int* INFO);
    void dgetrf(int M, int N, double* A, int LDA, int* IPIV, int* INFO);
  }

  template<>
  class GETRF<GPU>
  {
    public:

      	template<typename scalartype>
	static void execute(matrix<scalartype, GPU>& A, int* IPIV){
	  
	  int M = A.get_current_size().first;
	  int N = A.get_current_size().second;

	  int LDA = A.get_global_size().first;

	  int INFO=0;
	  execute(M, N, A.get_ptr(), LDA, IPIV, INFO);
	}

      static void execute( int M, int N, float* A, int LDA, int* IPIV, int& INFO){
	  GPU_KERNELS_GETRF::sgetrf(M, N, A, LDA, IPIV, &INFO);

	  if(INFO!=0)
	    throw std::logic_error(__PRETTY_FUNCTION__);
      }

      static void execute( int M, int N, double* A, int LDA, int* IPIV, int& INFO){
	  GPU_KERNELS_GETRF::dgetrf(M, N, A, LDA, IPIV, &INFO);
	  
	  if(INFO!=0)
	      throw std::logic_error(__PRETTY_FUNCTION__);
      }
  };
    
}

#endif
