//-*-C++-*- 

#ifndef LINALG_REMOVE_GPU_H
#define LINALG_REMOVE_GPU_H

namespace LIN_ALG {

    template<>
    struct REMOVE<GPU>
    {
      template<typename scalartype>
      static void execute(vector<scalartype, GPU>& /*v*/, int /*index*/)
      {
        throw std::logic_error(__FUNCTION__);
      }

      template<typename scalartype>
      static void row(matrix<scalartype, GPU>& M, int i){

	int Nr = M.get_current_size().first;
	int Nc = M.get_current_size().second;
	
	assert(i>-1 && i<Nr);
	
	if(Nc>0 && (Nr-i-1)>0){
	  int LD = M.get_global_size().first;
	  MEMORY_MANAGEMENT<GPU>::remove_first_row(Nr-i, Nc, M.get_ptr(i,0), LD);
	}

	M.get_current_size().first--;
      }
	
      template<typename scalartype>
      static void col(matrix<scalartype, GPU>& M, int i){
	
	int Nr = M.get_current_size().first;
	int Nc = M.get_current_size().second;
	
	assert(i>-1 && i<Nc);
	
	if(Nr>0 && (Nc-i-1)>0){
	  int LD = M.get_global_size().first;
	  MEMORY_MANAGEMENT<GPU>::remove_first_col(Nr, Nc-i, M.get_ptr(0,i), LD);
	}
	
	M.get_current_size().second--;
      }
	
      template<typename scalartype>
      static void row_and_column(matrix<scalartype, GPU>& M, int i){
	REMOVE<GPU>::col(M,i);
	REMOVE<GPU>::row(M,i);
      }
    };
}

#endif
