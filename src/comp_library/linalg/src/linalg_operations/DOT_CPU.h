//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_DOT_CPU_H
#define LINALG_DOT_CPU_H

namespace LIN_ALG {

  /*!
   *  y <-  x*y
   */
  template<>
  class DOT<CPU>
  {
  public:

    template<typename scalartype>
    static void execute(int N, scalartype& alpha, scalartype* x, scalartype* y){     
      execute(N, alpha, x, 1, y, 1);
    }

    inline static float execute(int length, float* a, int inc_a, float* b, int inc_b)
    {
      return BLAS::sdot_(&length, a, &inc_a, b, &inc_b);
    }
    
    inline static double execute(int length, double* a, int inc_a, double* b, int inc_b)
    {
      return BLAS::ddot_(&length, a, &inc_a, b, &inc_b);
    }
    
//     inline static std::complex<float> execute(int length, std::complex<float>* a, int inc_a, std::complex<float>* b, int inc_b)
//     {
//       return BLAS::cdot_(&length, a, &inc_a, b, &inc_b);
//     }
    
//     inline static std::complex<double> execute(int length, std::complex<double>* a, int inc_a, std::complex<double>* b, int inc_b)
//     {
//       return BLAS::zdot_(&length, a, &inc_a, b, &inc_b);
//     }
  };

}

#endif
