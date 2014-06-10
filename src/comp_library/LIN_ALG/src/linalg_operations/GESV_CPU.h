//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GESV_CPU_H
#define LINALG_GESV_CPU_H

namespace LIN_ALG {

  template<>
  class GESV<LIN_ALG::CPU>
  {
  public:

    template<typename scalar_type>
    static void execute(int N, scalar_type* A, scalar_type* X);

    template<typename scalar_type>
    static void execute(int N, int RHS, scalar_type* A, scalar_type* X);

  private:

    inline static void execute(int N, int RHS, float* A, int LDA, int* IPIV, float*X, int LDX, int INFO);

    inline static void execute(int N, int RHS, double* A, int LDA, int* IPIV, double* X, int LDX, int INFO);

    inline static void execute(int N, int RHS, std::complex<float>* A, int LDA, int* IPIV, std::complex<float>* X, int LDX, int INFO);

    inline static void execute(int N, int RHS, std::complex<double>* A, int LDA, int* IPIV, std::complex<double>* X, int LDX, int INFO);
  };

  template<typename scalar_type>
  void GESV<LIN_ALG::CPU>::execute(int N, scalar_type* A, scalar_type* X)
  {
    int IPIV[N];        
    execute(N, 1, A, N, IPIV, X, N, 0);
  }

  template<typename scalar_type>
  void GESV<LIN_ALG::CPU>::execute(int N, int RHS, scalar_type* A, scalar_type* X)
  {
    int IPIV[N];    
    execute(N, RHS, A, N, IPIV, X, N, 0);
  }

  void GESV<LIN_ALG::CPU>::execute(int N, int RHS, float* A, int LDA, int* IPIV, float*X, int LDX, int INFO)
  {
    LAPACK::sgesv_(&N, &RHS, A, &LDA, IPIV, X, &LDX, &INFO);
    assert(INFO==0);
  }
  
  void GESV<LIN_ALG::CPU>::execute(int N, int RHS, double* A, int LDA, int* IPIV, double* X, int LDX, int INFO)
  {
    LAPACK::dgesv_(&N, &RHS, A, &LDA, IPIV, X, &LDX, &INFO);  
    assert(INFO==0);
  }
    
  void GESV<LIN_ALG::CPU>::execute(int N, int RHS, std::complex<float>* A, int LDA, int* IPIV, std::complex<float>* X, int LDX, int INFO)
  {
    LAPACK::cgesv_(&N, &RHS, A, &LDA, IPIV, X, &LDX, &INFO);
    assert(INFO==0);
  }
  
  void GESV<LIN_ALG::CPU>::execute(int N, int RHS, std::complex<double>* A, int LDA, int* IPIV, std::complex<double>* X, int LDX, int INFO)
  {
    LAPACK::zgesv_(&N, &RHS, A, &LDA, IPIV, X, &LDX, &INFO);
    assert(INFO==0);
  }
}

#endif
