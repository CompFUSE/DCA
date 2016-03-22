//-*-C++-*-

#ifndef LINALG_GETRS_GPU_H
#define LINALG_GETRS_GPU_H

#include <vector>

namespace LIN_ALG {

  namespace GPU_KERNELS_GETRS {

    void sgetrs(char TRANS, int N, int NRHS, float* Matrix_A, int LDA, int* IPIV, float* Matrix_B, int LDB, int* INFO);
    void dgetrs(char TRANS, int N, int NRHS, double* Matrix_A, int LDA, int* IPIV, double* Matrix_B, int LDB, int* INFO);

  }

  template<>
  class GETRS<GPU>
  {
  public:

    template<typename scalartype>
    static void execute(matrix<scalartype, GPU>& A, matrix<scalartype, GPU>& X){

      int m_A = A.get_current_size().first;
      int n_A = A.get_current_size().second;
      int LDA = A.get_global_size().first;

      if(m_A != n_A)
        throw std::logic_error(__FUNCTION__);

      int m_X  = X.get_current_size().first;
      int NRHS = X.get_current_size().second;
      int LDX  = X.get_global_size().first;

      if(m_A != m_X)
        throw std::logic_error(__FUNCTION__);

      char TRANS('N');

      int INFO=0;

      std::vector<int> IPIV(m_A); 
      
      for(int i=0; i<m_A; i++)
        IPIV[i] = i+1;

      execute(TRANS, m_A, NRHS, A.get_ptr(), LDA, &IPIV[0], X.get_ptr(), LDX, INFO);
    }

    static void execute(char& TRANS, int& N, int& NRHS, float* Matrix_A, int& LDA, int* IPIV, float* Matrix_B, int& LDB, int& INFO)
    {
      GPU_KERNELS_GETRS::sgetrs(TRANS, N, NRHS, Matrix_A, LDA, IPIV, Matrix_B, LDB, &INFO);

      if(INFO!=0)
        throw std::logic_error(__PRETTY_FUNCTION__);
    }

    static void execute(char& TRANS, int& N, int& NRHS, double* Matrix_A, int& LDA, int* IPIV, double* Matrix_B, int& LDB, int& INFO)
    {
      GPU_KERNELS_GETRS::dgetrs(TRANS, N, NRHS, Matrix_A, LDA, IPIV, Matrix_B, LDB, &INFO);

      if(INFO!=0)
        throw std::logic_error(__PRETTY_FUNCTION__);
    }

  };

}

#endif
