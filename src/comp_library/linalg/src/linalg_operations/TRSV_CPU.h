//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_TRSV_CPU_H
#define LINALG_TRSV_CPU_H

namespace LIN_ALG {

  template<>
  class TRSV<CPU>
  {
  public:
    
    template<typename scalartype>
    inline static void execute(char uplo, char trans, char diag, int n, scalartype* A, int LDA, scalartype* X, int incx){
      dca::linalg::trsv(&uplo, &trans, &diag, n, A, LDA , X, incx);
    }
  };
}

#endif
