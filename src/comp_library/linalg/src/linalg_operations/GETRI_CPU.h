//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GETRI_CPU_H
#define LINALG_GETRI_CPU_H

namespace LIN_ALG {

  template<>
  class GETRI<CPU>
  {
  public:

    template<typename scalartype>
    static void execute(matrix<scalartype, CPU>& A, int* IPIV){
	  
      int M = A.get_current_size().first;
      int N = A.get_current_size().second;

      if(N != M)
	throw std::logic_error(__FUNCTION__);

      int LDA = A.get_global_size().first;

      int LWORK        = 16*N;
      scalartype* WORK = new scalartype[LWORK];

      int INFO=0;
      execute(N, A.get_ptr(), LDA, IPIV, WORK, LWORK, INFO);

      delete [] WORK;

      if(INFO!=0)
	throw std::logic_error(__FUNCTION__);
    }

  private:

    static void execute(int N, float* A, int LDA, int* IPIV, float* WORK, int LWORK, int INFO){
      LAPACK::sgetri_(&N, A, &LDA, IPIV, WORK, &LWORK, &INFO);
    }

    static void execute(int N, double* A, int LDA, int* IPIV, double* WORK, int LWORK, int INFO){
      LAPACK::dgetri_(&N, A, &LDA, IPIV, WORK, &LWORK, &INFO);
    }

    static void execute(int N, std::complex<float>* A, int LDA, int* IPIV, std::complex<float>* WORK, int LWORK, int INFO){
      LAPACK::cgetri_(&N, A, &LDA, IPIV, WORK, &LWORK, &INFO);
    }

    static void execute(int N, std::complex<double>* A, int LDA, int* IPIV, std::complex<double>* WORK, int LWORK, int INFO){
      LAPACK::zgetri_(&N, A, &LDA, IPIV, WORK, &LWORK, &INFO);
    }
  };

}

#endif
