//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_SCALE_GPU_H
#define LINALG_SCALE_GPU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_SCALE {

    void dscale(int length, double f, double* a, int inc_a);
    void dscale(int length, double f, double* a, int inc_a, int id);

    void many_rows(int Nc, int Ni, int* r_i, double* alpha, double* A, int LD,
		   int thread_id, int stream_id);
  }

  template<>
  class SCALE<GPU>
  {
  public:

    template<typename scalartype>
    inline static void row(matrix<scalartype, GPU>& M, scalartype val, int i,
                           int /*thread_id*/, int stream_id)
    {
      assert(stream_id==0);
      execute(M.get_current_size().second, val, M.get_ptr(i,0), M.get_global_size().first);
    }

    template<typename scalartype>
    inline static void col(matrix<scalartype, GPU>& M, scalartype val, int i,
                           int /*thread_id*/, int stream_id)
    {
      assert(stream_id==0);
      execute(M.get_current_size().first, val, M.get_ptr(0,i), val, M.get_global_size().second);
    }

    template<typename scalartype>
    inline static void row_and_col(matrix<scalartype, GPU>& M, scalartype val, int i,
                                   int /*thread_id*/, int stream_id)
    {
      assert(stream_id==0);
      SCALE<GPU>::row(M, val, i);
      SCALE<GPU>::col(M, val, i);
    }

    /*
    inline static void execute(int length, float f, float* a, int inc_a){
      GPU_KERNEL_SCALE::sscale(length, f, a, inc_a);
    }
    */

    inline static void execute(int length, double f, double* a, int inc_a,
			       int thread_id, int stream_id)
    {
      assert(stream_id==0);
      GPU_KERNEL_SCALE::dscale(length, f, a, inc_a, thread_id);
    }

    /*
    inline static void execute(int length, float f, cuComplex* a, int inc_a){
      cublasCsscal(length, f, a, inc_a);
    }

    inline static void execute(int length, double f, cuDoubleComplex* a, int inc_a){
      cublasZdscal(length, f, a, inc_a);
    }
    */

    inline static void many_rows(int Nc, int Ni, int* r_i, double* alpha, double* A, int LD,
				 int thread_id, int stream_id)
    {
      assert(stream_id==0);
      GPU_KERNEL_SCALE::many_rows(Nc, Ni, r_i, alpha, A, LD, thread_id, stream_id);
    }
  };

}

#endif
