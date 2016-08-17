//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GETRI_GPU_H
#define LINALG_GETRI_GPU_H

namespace LIN_ALG {

  namespace GPU_KERNELS_GETRI {

      int get_work_space(int n);

      void sgetri(int N, float* A, int LDA, int* IPIV, float* WORK, int LWORK, int* INFO);
      void dgetri(int N, double* A, int LDA, int* IPIV, double* WORK, int LWORK, int* INFO);
  }

  template<>
  class GETRI<GPU>
  {
    public:

      	template<typename scalartype>
	static void execute(matrix<scalartype, GPU>& A, int* IPIV){
	  
	  int M = A.size().first;
	  int N = A.size().second;

	  if(N != M)
	    throw std::logic_error(__FUNCTION__);

	  int LDA = A.leadingDimension();

	  int LWORK        = GPU_KERNELS_GETRI::get_work_space(N)*N;
	  dca::linalg::Vector<scalartype, GPU> WORK(LWORK);

	  int INFO=0;
	  execute(N, A.ptr(), LDA, IPIV, WORK.ptr(), LWORK, INFO);
	}

      static void execute(int N, float* A, int LDA, int* IPIV, float* WORK, int LWORK, int& INFO){
	  GPU_KERNELS_GETRI::sgetri(N, A, LDA, IPIV, WORK, LWORK, &INFO);

	  if(INFO!=0)
	    throw std::logic_error(__FUNCTION__);
      }

      static void execute(int N, double* A, int LDA, int* IPIV, double* WORK, int LWORK, int& INFO){
	  GPU_KERNELS_GETRI::dgetri(N, A, LDA, IPIV, WORK, LWORK, &INFO);

	  if(INFO!=0)
	    throw std::logic_error(__FUNCTION__);
      }
  };

}

#endif
