//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_LASET_GPU_H
#define LINALG_LASET_GPU_H

namespace LIN_ALG 
{
  namespace GPU_KERNEL_LASET 
  {
    void set_zero (char UPLO, int M, int N, double* A, int LDA);
    void set_unity(           int M, int N, double* A, int LDA);

    void set_zero (int M, int N, double* A, int LDA, int thread_id, int stream_id);
    void set_unity(int M, int N, double* A, int LDA, int thread_id, int stream_id);
  }

  template<>
  class LASET<GPU>
  {
  public:

    static void set_zero(int M, int N, double* A, int LDA, int thread_id, int stream_id)
    {
      if(false)
	{
	  const char UPLO = 'A';
	  GPU_KERNEL_LASET::set_zero(UPLO, M, N, A, LDA);
	}
      else
	{
	  GPU_KERNEL_LASET::set_zero(M, N, A, LDA, thread_id, stream_id);
	}
    }

    static void set_unity(int M, int N, double* A, int LDA, int thread_id, int stream_id)
    {
      if(false)
	GPU_KERNEL_LASET::set_unity(M, N, A, LDA);
      else
	GPU_KERNEL_LASET::set_unity(M, N, A, LDA, thread_id, stream_id);
    }

  };
}

#endif
