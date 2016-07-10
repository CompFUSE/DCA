//-*-C++-*-

#ifndef LINALG_COPY_GPU_H
#define LINALG_COPY_GPU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_COPY {

    void dcopy(int length, double* x, int inc_x, double* y, int inc_y);
    void dcopy(int length, double* x, int inc_x, double* y, int inc_y, int thread_id);

    void many_row_copies   (int N_x, int N_i, int* i_x, double* x, int inc_x, int* i_y, double* y, int inc_y, int thread_id, int stream_id);
    void many_column_copies(int N_x, int N_i, int* i_x, double* x, int inc_x, int* i_y, double* y, int inc_y, int thread_id, int stream_id);
  }

  template<>
  class COPY<GPU>
  {
  public:
      
    template<typename scalartype>
    static void row(matrix<scalartype, GPU>& X, int xi, matrix<scalartype, GPU>& Y, int yi,
                    int /*thread_id*/, int /*stream_id*/)
    {

      assert(xi>-1 && xi<X.get_current_size().first);
      assert(yi>-1 && yi<Y.get_current_size().first);

      assert(X.get_current_size().second == Y.get_current_size().second);
	  
      execute(X.get_current_size().second, X.get_ptr(xi,0), X.get_leading_dimension(), Y.get_ptr(yi,0), Y.get_leading_dimension());
    }
      
    template<typename scalartype>
    static void col(matrix<scalartype, GPU>& X, int xi, matrix<scalartype, GPU>& Y, int yi,
                    int /*thread_id*/, int /*stream_id*/)
    {
      assert(xi>-1 && xi<X.get_current_size().second);
      assert(yi>-1 && yi<Y.get_current_size().second);
      assert(X.get_current_size().first == Y.get_current_size().first);
	  
      execute(X.get_number_of_rows(), X.get_ptr(0,xi), 1, Y.get_ptr(0,yi), 1);
    }
     
    inline static void execute(int length, double* x, int inc_x, double* y, int inc_y,
			       int thread_id, int /*stream_id*/)
    {
      // assert(stream_id==0);
      GPU_KERNEL_COPY::dcopy(length, x, inc_x, y, inc_y, thread_id);
    }

    template<typename scalartype>
    inline static void many_rows(int N_x, int N_i, int* i_x, scalartype* x, int inc_x, int* i_y, scalartype* y, int inc_y,
				 int thread_id, int stream_id)
    {
      GPU_KERNEL_COPY::many_row_copies(N_x, N_i, i_x, x, inc_x, i_y, y, inc_y, thread_id, stream_id);
    }

    template<typename scalartype>
    inline static void many_columns(int N_x, int N_i, int* i_x, scalartype* x, int inc_x, int* i_y, scalartype* y, int inc_y,
				    int thread_id, int stream_id)
    {
      GPU_KERNEL_COPY::many_column_copies(N_x, N_i, i_x, x, inc_x, i_y, y, inc_y, thread_id, stream_id);
    }
  };

}

#endif
