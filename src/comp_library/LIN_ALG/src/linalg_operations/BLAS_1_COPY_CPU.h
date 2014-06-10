//-*-C++-*-                                                                                                                                                                             

#ifndef LINALG_COPY_CPU_H
#define LINALG_COPY_CPU_H

namespace LIN_ALG {

  template<>
  class COPY<CPU>
  {
  public:

    template<typename scalartype>
    static void row(matrix<scalartype, CPU>& X, int xi, matrix<scalartype, CPU>& Y, int yi,
		    int thread_id, int stream_id)
    {

      assert(xi>-1 && xi<X.get_current_size().first);
      assert(yi>-1 && yi<Y.get_current_size().first);

      assert(X.get_current_size().second == Y.get_current_size().second);

      execute(X.get_current_size().second, X.get_ptr(xi,0), X.get_global_size().first, Y.get_ptr(yi,0), Y.get_global_size().first);
    }

    template<typename scalartype>
    static void col(matrix<scalartype, CPU>& X, int xi, matrix<scalartype, CPU>& Y, int yi,
		    int thread_id, int stream_id)
    {
      assert(xi>-1 && xi<X.get_current_size().second);
      assert(yi>-1 && yi<Y.get_current_size().second);
      assert(X.get_current_size().first == Y.get_current_size().first);

      execute(X.get_current_size().first, X.get_ptr(0,xi), 1, Y.get_ptr(0,yi), 1);
    }

    inline static void execute(int length, float* x, int inc_x, float* y, int inc_y,
			       int thread_id, int stream_id)
    {
      BLAS::scopy_(&length, x, &inc_x, y, &inc_y);
    }

    inline static void execute(int length, double* x, int inc_x, double* y, int inc_y,
			       int thread_id, int stream_id)
    {
      BLAS::dcopy_(&length, x, &inc_x, y, &inc_y);
    }

    inline static void execute(int length, std::complex<float>* x, int inc_x, std::complex<float>* y, int inc_y,
			       int thread_id, int stream_id)
    {
      BLAS::ccopy_(&length, x, &inc_x, y, &inc_y);
    }

    inline static void execute(int length, std::complex<double>* x, int inc_x, std::complex<double>* y, int inc_y,
			       int thread_id, int stream_id)
    {
      BLAS::zcopy_(&length, x, &inc_x, y, &inc_y);
    }


    template<typename scalartype>
    inline static void many_rows(int N_x, int N_i, int* i_x, scalartype* x, int inc_x, int* i_y, scalartype* y, int inc_y,
				 int thread_id, int stream_id)
    {
      for(int i=0; i<N_x; ++i)
	for(int l=0; l<N_i; ++l)
	  y[i_y[l]+i*inc_y] = x[i_x[l]+i*inc_x];
    }

    template<typename scalartype>
    inline static void many_columns(int N_x, int N_i, int* i_x, scalartype* x, int inc_x, int* i_y, scalartype* y, int inc_y,
				    int thread_id, int stream_id)
    {
      for(int l=0; l<N_i; ++l)
	for(int i=0; i<N_x; ++i)
	  y[i+i_y[l]*inc_y] = x[i+i_x[l]*inc_x];
    }
  };

}

#endif

