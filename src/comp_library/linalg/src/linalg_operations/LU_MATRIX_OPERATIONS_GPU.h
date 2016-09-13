//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_LU_MATRIX_OPERATIONS_GPU_H
#define LINALG_LU_MATRIX_OPERATIONS_GPU_H

namespace LIN_ALG {

  namespace LU_MATRIX_OPERATIONS_GPU {
    
    template<class scalartype>
    scalartype determinant_tridiagonal(int n, scalartype* A, int lda);
    
    template<class scalartype>
    std::pair<scalartype,scalartype> minmax_diagonal(int n, scalartype* A, int lda);
  }

  template<>
  class LU_MATRIX_OPERATIONS<GPU>
  {
  public:

    template<typename scalartype>
    static scalartype det(matrix<scalartype, GPU>& M){

      assert(M.is_square());

      scalartype determinant = LU_MATRIX_OPERATIONS_GPU::determinant_tridiagonal(M.nrRows(), M.ptr(0,0), M.leadingDimension());

      return determinant;
    }

    template<typename scalartype>
    static scalartype ratio(matrix<scalartype, GPU>& M){

      assert(M.is_square());

      std::pair<scalartype, scalartype> p = LU_MATRIX_OPERATIONS_GPU::minmax_diagonal(M.nrRows(), M.ptr(0,0), M.leadingDimension());

      return p.first/p.second;
    }

  };

}

#endif
