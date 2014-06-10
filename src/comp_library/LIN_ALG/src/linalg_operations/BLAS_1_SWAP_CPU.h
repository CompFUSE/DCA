//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_SWAP_CPU_H
#define LINALG_SWAP_CPU_H

namespace LIN_ALG {

  template<>
  class SWAP<CPU>
  {
  public:

    template<typename scalartype>
    static void row(matrix<scalartype, CPU>& M, int i, int j,
			   int thread_id, int stream_id)
    {
      
      
      assert(i>-1 && i<M.get_current_size().first);
      assert(j>-1 && j<M.get_current_size().first);
     
      execute(M.get_current_size().second, &M(i,0), M.get_global_size().first, &M(j,0), M.get_global_size().first, thread_id, stream_id);
    }

    template<typename scalartype>
    static void col(matrix<scalartype, CPU>& M, int i, int j,
		    int thread_id, int stream_id)
    {

      assert(i>-1 && i<M.get_current_size().second);
      assert(j>-1 && j<M.get_current_size().second);

      execute(M.get_current_size().first, &M(0,i), 1, &M(0,j), 1, thread_id, stream_id);
    }

    template<typename scalartype>
    static void row_and_column(matrix<scalartype, CPU>& M, int i, int j,
			   int thread_id, int stream_id)
    {
      SWAP<CPU>::row( M, i, j, thread_id, stream_id);
      SWAP<CPU>::col( M, i, j, thread_id, stream_id);
    }

  private:

    template<typename scalar_type>
    static void execute(int length, scalar_type* a, int inc_a, scalar_type* b, int inc_b,
			   int thread_id, int stream_id)
    {
      throw std::logic_error(__FUNCTION__);
    }

    inline static void execute(int length, float* a, int inc_a, float* b, int inc_b,
			       int thread_id, int stream_id)
    {
      BLAS::sswap_(&length, a, &inc_a, b, &inc_b);
    }
    
    inline static void execute(int length, double* a, int inc_a, double* b, int inc_b,
			   int thread_id, int stream_id)
    {
      BLAS::dswap_(&length, a, &inc_a, b, &inc_b);
    }
    
    inline static void execute(int length, std::complex<float>* a, int inc_a, std::complex<float>* b, int inc_b,
			       int thread_id, int stream_id)
    {
      BLAS::cswap_(&length, a, &inc_a, b, &inc_b);
    }
    
    inline static void execute(int length, std::complex<double>* a, int inc_a, std::complex<double>* b, int inc_b,
			       int thread_id, int stream_id)
    {
      BLAS::zswap_(&length, a, &inc_a, b, &inc_b);
    }
  };

}

#endif
