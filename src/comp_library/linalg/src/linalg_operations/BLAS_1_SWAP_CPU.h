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
      
      
      assert(i>-1 && i<M.size().first);
      assert(j>-1 && j<M.size().first);
     
      execute(M.size().second, &M(i,0), M.leadingDimension(), &M(j,0), M.leadingDimension(), thread_id, stream_id);
    }

    template<typename scalartype>
    static void col(matrix<scalartype, CPU>& M, int i, int j,
		    int thread_id, int stream_id)
    {

      assert(i>-1 && i<M.size().second);
      assert(j>-1 && j<M.size().second);

      execute(M.size().first, &M(0,i), 1, &M(0,j), 1, thread_id, stream_id);
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
                        int /*thread_id*/, int /*stream_id*/)
    {
      dca::linalg::blas::swap(length, a, inc_a, b, inc_b);
    }
  };

}

#endif
