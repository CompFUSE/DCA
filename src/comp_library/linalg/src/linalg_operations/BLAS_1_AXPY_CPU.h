//-*-C++-*-

#ifndef LINALG_AXPY_CPU_H
#define LINALG_AXPY_CPU_H

namespace LIN_ALG {

  /*!
   *  y <- \alpha x + y
   */
  template<>
  class AXPY<CPU>
  {
  public:

    template<typename scalartype>
    static void execute(int N, scalartype& alpha, scalartype* x, scalartype* y, int /*thread_id*/, int /*stream_id*/){
      execute(N, alpha, x, 1, y, 1);
    }

    inline static void execute(int length, float alpha, float* a, int inc_a, float* b, int inc_b, int /*thread_id*/, int /*stream_id*/)
    {
      BLAS::saxpy_(&length, &alpha, a, &inc_a, b, &inc_b);
    }

    inline static void execute(int length, double alpha, double* a, int inc_a, double* b, int inc_b, int /*thread_id*/, int /*stream_id*/)
    {
      BLAS::daxpy_(&length, &alpha, a, &inc_a, b, &inc_b);
    }

    inline static void execute(int length, std::complex<float> alpha, std::complex<float>* a, int inc_a, std::complex<float>* b, int inc_b, int /*thread_id*/, int /*stream_id*/)
    {
      BLAS::caxpy_(&length, &alpha, a, &inc_a, b, &inc_b);
    }

    inline static void execute(int length, std::complex<double> alpha, std::complex<double>* a, int inc_a, std::complex<double>* b, int inc_b, int /*thread_id*/, int /*stream_id*/)
    {
      BLAS::zaxpy_(&length, &alpha, a, &inc_a, b, &inc_b);
    }
  };

}

#endif
