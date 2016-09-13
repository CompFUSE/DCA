//-*-C++-*-

#ifndef LINALG_GETRS_CPU_H
#define LINALG_GETRS_CPU_H

namespace LIN_ALG 
{

  template<>
  class GETRS<CPU>
  {
  public:

    template<typename scalartype>
    static void execute(matrix<scalartype, CPU>& A, matrix<scalartype, CPU>& X){

      int m_A = A.size().first;
      int n_A = A.size().second;
      int LDA = A.leadingDimension();

      if(m_A != n_A)
        throw std::logic_error(__FUNCTION__);

      int m_X  = X.size().first;
      int NRHS = X.size().second;
      int LDX  = X.leadingDimension();

      if(m_A != m_X)
        throw std::logic_error(__FUNCTION__);

      char TRANS('N');

      int INFO=0;

      //int IPIV[m_A];
      int* IPIV = new int[m_A];

      for(int i=0; i<m_A; i++)
        IPIV[i] = i+1;

      execute(TRANS, m_A, NRHS, A.ptr(), LDA, IPIV, X.ptr(), LDX, INFO);

      delete [] IPIV;
    }

    template<typename scalartype>
    static void execute(matrix<scalartype, CPU>& A, dca::linalg::Vector<scalartype, CPU>& X){

      int m_A = A.size().first;
      int n_A = A.size().second;
      int LDA = A.leadingDimension();

      if(m_A != n_A)
        throw std::logic_error(__FUNCTION__);

      int m_X  = X.size();
      int NRHS = 1;//X.size().second;
      int LDX  = X.capacity();

      if(m_A != m_X)
        throw std::logic_error(__FUNCTION__);

      char TRANS('N');

      int INFO=0;

      //int IPIV[m_A];
      int* IPIV = new int[m_A];

      for(int i=0; i<m_A; i++)
        IPIV[i] = i+1;

      execute(TRANS, m_A, NRHS, A.ptr(), LDA, IPIV, X.ptr(), LDX, INFO);

      delete [] IPIV;
    }

  private:

    static void execute(char& TRANS, int& N, int& NRHS, float* Matrix_A, int& LDA, int* IPIV, float* Matrix_B, int& LDB, int& INFO)
    {
      LAPACK::sgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);

      if(INFO!=0)
        throw std::logic_error(__PRETTY_FUNCTION__);
    }

    static void execute(char& TRANS, int& N, int& NRHS, double* Matrix_A, int& LDA, int* IPIV, double* Matrix_B, int& LDB, int& INFO)
    {
      LAPACK::dgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);

      if(INFO!=0)
        throw std::logic_error(__PRETTY_FUNCTION__);
    }

    static void execute(char& TRANS, int& N, int& NRHS, std::complex<float>* Matrix_A, int& LDA, int* IPIV, std::complex<float>* Matrix_B, int& LDB, int& INFO)
    {
      LAPACK::cgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);

      if(INFO!=0)
        throw std::logic_error(__PRETTY_FUNCTION__);
    }

    static void execute(char& TRANS, int& N, int& NRHS, std::complex<double>* Matrix_A, int& LDA, int* IPIV, std::complex<double>* Matrix_B, int& LDB, int& INFO)
    {
      LAPACK::zgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);

      if(INFO!=0)
        throw std::logic_error(__PRETTY_FUNCTION__);
    }

  };

}

#endif
