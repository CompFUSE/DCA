//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_SCALE_CPU_H
#define LINALG_SCALE_CPU_H

namespace LIN_ALG {

  template<>
  class SCALE<CPU>
  {
  public:

    template<typename scalartype>
    inline static void row(matrix<scalartype, CPU>& M, scalartype val, int i,
                           int /*thread_id*/, int /*stream_id*/)
    {
      execute(M.size().second, val, &M(i,0), M.leadingDimension());
    }

    template<typename scalartype>
    inline static void col(matrix<scalartype, CPU>& M, scalartype val, int i,
                           int /*thread_id*/, int /*stream_id*/)
    {
      execute(M.size().first, val, &M(0,i), M.capacity().second);
    }

    template<typename scalartype>
    inline static void row_and_col(matrix<scalartype, CPU>& M, scalartype val, int i,
                                   int /*thread_id*/, int /*stream_id*/)
    {
      SCALE<CPU>::row(M, val, i);
      SCALE<CPU>::col(M, val, i);
    }

    template<typename scalartype>
    inline static void execute(int length, scalartype f, scalartype* a, int inc_a,
                               int /*thread_id*/, int /*stream_id*/)
    {
      dca::linalg::blas::scal(&length, &f, a, &inc_a);
    }

    inline static void many_rows(int Nc, int Ni, int* r_i, double* alpha, double* A, int LD,
                                 int /*thread_id*/, int /*stream_id*/)
    {
      for(int j=0; j<Nc; ++j)
        for(int l=0; l<Ni; ++l)
          A[r_i[l]+j*LD] *= alpha[l]; 
    }
  };

}

#endif
