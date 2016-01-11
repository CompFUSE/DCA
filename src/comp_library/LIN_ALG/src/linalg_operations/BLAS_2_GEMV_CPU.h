//-*-C++-*-

#ifndef LINALG_GEMV_CPU_H
#define LINALG_GEMV_CPU_H

namespace LIN_ALG {

  template<>
  class GEMV<CPU>
  {
  public:

    template<typename scalartype>
    static void execute(matrix<scalartype, CPU>& A,
                        vector<scalartype, CPU>& B,
                        vector<scalartype, CPU>& C,
                        int thread_id=0, int stream_id=0);

    template<typename scalartype>
    static void execute(int M, int N,
                        scalartype* A,
                        scalartype* B,
                        scalartype* C);

    template<typename scalartype>
    static void execute(char TRANS,
                        matrix<scalartype, CPU>& A, vector<scalartype, CPU>& B, vector<scalartype, CPU>& C,
                        int thread_id=0, int stream_id=0);

    template<typename scalartype>
    static void execute(scalartype a, matrix<scalartype, CPU>& A, vector<scalartype, CPU>& B,
                        scalartype b, vector<scalartype, CPU>& C,
                        int thread_id=0, int stream_id=0);

    template<typename scalartype>
    static void execute(char TRANS,
                        scalartype a, matrix<scalartype, CPU>& A, vector<scalartype, CPU>& B,
                        scalartype b, vector<scalartype, CPU>& C,
                        int thread_id=0, int stream_id=0);


    template<typename scalartype>
    static void execute(char TRANS, int M, int N,
                        scalartype alpha,
                        scalartype* A, int LDA,
                        scalartype* B, int LDB,
                        scalartype  beta,
                        scalartype* C, int LDC,
                        int thread_id=0, int stream_id=0);

  private:

    inline static void execute_gemv(char TRANS, int M, int N,
                                    float alpha,
                                    float* A, int LDA,
                                    float* B, int LDB,
                                    float  beta,
                                    float* C, int LDC);

    inline static void execute_gemv(char TRANS, int M, int N,
                                    double alpha,
                                    double* A, int LDA,
                                    double* B, int LDB,
                                    double  beta,
                                    double* C, int LDC);

    inline static void execute_gemv(char TRANS, int M, int N,
                                    std::complex<float> alpha,
                                    std::complex<float>* A, int LDA,
                                    std::complex<float>* B, int LDB,
                                    std::complex<float>  beta,
                                    std::complex<float>* C, int LDC);

    inline static void execute_gemv(char TRANS, int M, int N,
                                    std::complex<double> alpha,
                                    std::complex<double>* A, int LDA,
                                    std::complex<double>* B, int LDB,
                                    std::complex<double>  beta,
                                    std::complex<double>* C, int LDC);
  };

  template<typename scalartype>
  void GEMV<CPU>::execute(int M, int N,
                          scalartype* A,
                          scalartype* B,
                          scalartype* C)
  {
    char TRANS = 'N';

    scalartype alpha = 1;
    scalartype beta  = 0;

    int LDA = M;
    int LDB = 1;//N;
    int LDC = 1;//N;

    execute_gemv(TRANS, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
  }

  template<typename scalartype>
  void GEMV<CPU>::execute(char TRANS,
                          matrix<scalartype, CPU>& A,
                          vector<scalartype, CPU>& B,
                          vector<scalartype, CPU>& C,
                          int /*thread_id*/, int /*stream_id*/)
  {
    scalartype ONE (1);
    scalartype ZERO(0);

    int M = A.get_current_size().first;
    int N = A.get_current_size().second;

    if(TRANS=='N'){
      assert(M == C.get_current_size());
      assert(N == B.get_current_size());
    }

    if(TRANS=='T' or TRANS=='H'){
      assert(M == B.get_current_size());
      assert(N == C.get_current_size());
    }

    int LDA = A.get_global_size().first;
    int LDB = 1;//B.get_global_size();
    int LDC = 1;//C.get_global_size();

    execute_gemv(TRANS, M, N, ONE, &A(0,0), LDA, &B[0], LDB, ZERO, &C[0], LDC);
  }

  template<typename scalartype>
  void GEMV<CPU>::execute(scalartype alpha, matrix<scalartype, CPU>& A, vector<scalartype, CPU>& B,
                          scalartype beta , vector<scalartype, CPU>& C,
                          int /*thread_id*/, int /*stream_id*/)
  {
    assert(A.get_current_size().second  == B.get_current_size());
    assert(B.get_current_size()         == C.get_current_size());

    char TRANS = 'N';

    int M = A.get_current_size().first;
    int N = A.get_current_size().second;

    int LDA = A.get_global_size().first;
    int LDB = 1;//B.get_global_size();
    int LDC = 1;//C.get_global_size();

    execute_gemv(TRANS, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
  }

  template<typename scalartype>
  void GEMV<CPU>::execute(char TRANS,
                          scalartype alpha, matrix<scalartype, CPU>& A, vector<scalartype, CPU>& B,
                          scalartype beta , vector<scalartype, CPU>& C,
                          int /*thread_id*/, int /*stream_id*/)
  {
    int M = A.get_current_size().first;
    int N = A.get_current_size().second;

    if(TRANS=='N'){
      assert(M == C.get_current_size());
      assert(N == B.get_current_size());
    }

    if(TRANS=='T' or TRANS=='H'){
      assert(M == B.get_current_size());
      assert(N == C.get_current_size());
    }

    int LDA = A.get_global_size().first;
    int LDB = 1;//B.get_global_size();
    int LDC = 1;//C.get_global_size();

    execute_gemv(TRANS, M, N, alpha, &A(0,0), LDA, &B[0], LDB, beta, &C[0], LDC);
  }


  template<typename scalartype>
  void GEMV<CPU>::execute(char TRANS, int M, int N,
                          scalartype alpha,
                          scalartype* A, int LDA,
                          scalartype* B, int LDB,
                          scalartype  beta,
                          scalartype* C, int LDC,
                          int /*thread_id*/, int /*stream_id*/)
  {
    execute_gemv(TRANS, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
  }

  void GEMV<CPU>::execute_gemv(char TRANS, int M, int N,
                               float alpha,
                               float* A, int LDA,
                               float* B, int LDB,
                               float  beta,
                               float* C, int LDC)
  {
    BLAS::sgemv_(&TRANS, &M, &N, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
  }

  void GEMV<CPU>::execute_gemv(char TRANS, int M, int N,
                               double alpha,
                               double* A, int LDA,
                               double* B, int LDB,
                               double  beta,
                               double* C, int LDC)
  {
    BLAS::dgemv_(&TRANS, &M, &N, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
  }

  void GEMV<CPU>::execute_gemv(char TRANS, int M, int N,
                               std::complex<float> alpha,
                               std::complex<float>* A, int LDA,
                               std::complex<float>* B, int LDB,
                               std::complex<float>  beta,
                               std::complex<float>* C, int LDC)
  {
    BLAS::cgemv_(&TRANS, &M, &N, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
  }

  void GEMV<CPU>::execute_gemv(char TRANS, int M, int N,
                               std::complex<double> alpha,
                               std::complex<double>* A, int LDA,
                               std::complex<double>* B, int LDB,
                               std::complex<double>  beta,
                               std::complex<double>* C, int LDC)
  {
    BLAS::zgemv_(&TRANS, &M, &N, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
  }

}

#endif
