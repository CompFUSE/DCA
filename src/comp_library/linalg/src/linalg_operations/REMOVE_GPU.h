//-*-C++-*- 

#ifndef LINALG_REMOVE_GPU_H
#define LINALG_REMOVE_GPU_H

namespace LIN_ALG {

    template<>
    struct REMOVE<GPU>
    {
      template<typename scalartype>
      static void execute(dca::linalg::Vector<scalartype, GPU>& /*v*/, int /*index*/)
      {
        throw std::logic_error(__FUNCTION__);
      }

      template<typename scalartype>
      static void row(matrix<scalartype, GPU>& M, int i){

	int Nr = M.size().first;
	int Nc = M.size().second;
	
	assert(i>-1 && i<Nr);
	
	if(Nc>0 && (Nr-i-1)>0){
	  int LD = M.leadingDimension();
	  MEMORY_MANAGEMENT<GPU>::remove_first_row(Nr-i, Nc, M.ptr(i,0), LD);
	}

        auto size = M.size();
        size.first--;
        M.resize(size);
      }
	
      template<typename scalartype>
      static void col(matrix<scalartype, GPU>& M, int i){
	
	int Nr = M.size().first;
	int Nc = M.size().second;
	
	assert(i>-1 && i<Nc);
	
	if(Nr>0 && (Nc-i-1)>0){
	  int LD = M.leadingDimension();
	  MEMORY_MANAGEMENT<GPU>::remove_first_col(Nr, Nc-i, M.ptr(0,i), LD);
	}
	
        auto size = M.size();
        size.second--;
        M.resize(size);
      }
	
      template<typename scalartype>
      static void row_and_column(matrix<scalartype, GPU>& M, int i){
	REMOVE<GPU>::col(M,i);
	REMOVE<GPU>::row(M,i);
      }
    };
}

#endif
