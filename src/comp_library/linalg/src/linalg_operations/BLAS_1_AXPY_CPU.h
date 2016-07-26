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
    inline static void execute(int n, scalartype alpha, scalartype* x, int incx, scalartype* y, int incy, int /*thread_id*/, int /*stream_id*/){
      dca::linalg::axpy(n, alpha, x, incx, y, incy);
    }
  };
}

#endif
