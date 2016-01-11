//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_LU_MATRIX_OPERATIONS_CPU_H
#define LINALG_LU_MATRIX_OPERATIONS_CPU_H

namespace LIN_ALG {

  template<>
  class LU_MATRIX_OPERATIONS<CPU>
  {
  public:

    template<typename scalartype>
    static double det(matrix<scalartype, CPU>& M){

      assert(M.is_square());

      scalartype determinant = 1.;

      for(int i=0; i<M.get_current_size().first; i++)
        determinant *= M(i,i);

      return determinant;
    }

    template<typename scalartype>
    static scalartype ratio(matrix<scalartype, CPU>& M){

      assert(M.is_square());

      scalartype MIN = std::abs(M(0,0));
      scalartype MAX = std::abs(M(0,0));

      for(int i=0; i<M.get_current_size().first; i++)
	{
	  if(std::fabs(M(i,i)) < MIN ) MIN = std::abs(M(i,i));
	  if(std::fabs(M(i,i)) > MAX ) MAX = std::abs(M(i,i));
	}

      return MAX/MIN;
    }

  };

}

#endif
