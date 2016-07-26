//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_TRSM_CPU_H
#define LINALG_TRSM_CPU_H

namespace LIN_ALG {

  template<>
  class TRSM<CPU>
  {
  public:
    
    template<typename scalartype>
    inline static void execute(char uplo, char diag, matrix<scalartype, CPU>& A, matrix<scalartype, CPU>& X,  int /*thread_id*/, int /*stream_id*/)
    {      
      assert(uplo=='U' or uplo=='L');
      assert(diag=='U' or diag=='N');
      
      dca::linalg::trsm("L", &uplo, "N", &diag, X.get_number_of_rows(), X.get_number_of_cols(), scalartype(1), A.get_ptr(), A.get_leading_dimension(), X.get_ptr(), X.get_leading_dimension());
    }
  };
}

#endif
