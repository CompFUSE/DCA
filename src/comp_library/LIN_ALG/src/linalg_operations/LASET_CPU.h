//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_LASET_CPU_H
#define LINALG_LASET_CPU_H

namespace LIN_ALG 
{
  template<>
  class LASET<CPU>
  {
  public:
    
    static void set_zero(int M, int N, double* A, int LDA, int thread_id, int stream_id)
    {
      char   UPLO = 'A';
      double ALPHA = 0;
      double BETA  = 0;

      LAPACK::dlaset_(&UPLO, &M, &N, &ALPHA, &BETA, A, &LDA);
    }

    static void set_unity(int M, int N, double* A, int LDA, int thread_id, int stream_id)
    {
      char   UPLO  = 'A';
      double ALPHA = 0;
      double BETA  = 1;

      LAPACK::dlaset_(&UPLO, &M, &N, &ALPHA, &BETA, A, &LDA);
    }
  };
}

#endif
