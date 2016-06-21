//-*-C++-*- 

#ifndef LINALG_REMOVE_CPU_H
#define LINALG_REMOVE_CPU_H

namespace LIN_ALG {

  template<>
  struct REMOVE<CPU>
  {
    template<typename scalartype>
    static void execute(vector<scalartype, CPU>& v, int index)
    {
      memmove(v.get_ptr(index), v.get_ptr(index+1), sizeof(scalartype)*(v.size()-index-1));
      v.size() = v.size()-1;
    }

    template<typename scalartype>
    static void row(matrix<scalartype, CPU>& M, int index)
    {
      assert(index>-1 && index<M.get_current_size().first);

      int L = M.get_number_of_rows()-index-1;

      if(M.get_number_of_cols()>0 && L>0){

// 	scalartype* data = M.get_ptr(0,0);
	
// 	int N = M.get_global_size().first;
// 	int L = M.get_current_size().first-index-1;
	
// 	for(int i=0; i<M.get_current_size().second; i++)
// 	  memmove(&data[index + i*N], &data[(index+1) + i*N], sizeof(scalartype)*L);

// 	int L = M.get_number_of_rows()-index-1;

	for(int i=0; i<M.get_number_of_cols(); i++)
	  memmove(M.get_ptr(index,i), M.get_ptr((index+1),i), sizeof(scalartype)*L);
      }

      M.get_current_size().first--;
    }

    template<typename scalartype>
    static void col(matrix<scalartype, CPU>& M, int index)
    {
      assert(index>-1 && index<M.get_current_size().second);

      int L = M.get_number_of_cols()-index-1;

      if(M.get_number_of_rows()>0 && L>0){

// 	scalartype* data = &M(0,0);

// 	int N = M.get_global_size().first;
// 	int L = M.get_current_size().second-index-1;
	
// 	memmove(&data[index*N], &data[(index+1)*N], sizeof(scalartype)*N*L);

	int N = M.get_leading_dimension();
	memmove(M.get_ptr(0,index), M.get_ptr(0,index+1), sizeof(scalartype)*N*L);
      }

      M.get_current_size().second--;
    }

    template<typename scalartype>
    static void row_and_column(matrix<scalartype, CPU>& M, int i)
    {
      REMOVE<CPU>::row(M, i);
      REMOVE<CPU>::col(M, i);
    }
  };

}

#endif
