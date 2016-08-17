//-*-C++-*- 

#ifndef LINALG_REMOVE_CPU_H
#define LINALG_REMOVE_CPU_H

namespace LIN_ALG {

  template<>
  struct REMOVE<CPU>
  {
    template<typename scalartype>
    static void execute(dca::linalg::Vector<scalartype, CPU>& v, int index)
    {
      memmove(v.ptr(index), v.ptr(index+1), sizeof(scalartype)*(v.size()-index-1));
      v.size() = v.size()-1;
    }

    template<typename scalartype>
    static void row(matrix<scalartype, CPU>& M, int index)
    {
      assert(index>-1 && index<M.size().first);

      int L = M.nrRows()-index-1;

      if(M.nrCols()>0 && L>0){

// 	scalartype* data = M.ptr(0,0);
	
// 	int N = M.leadingDimension();
// 	int L = M.size().first-index-1;
	
// 	for(int i=0; i<M.size().second; i++)
// 	  memmove(&data[index + i*N], &data[(index+1) + i*N], sizeof(scalartype)*L);

// 	int L = M.nrRows()-index-1;

	for(int i=0; i<M.nrCols(); i++)
	  memmove(M.ptr(index,i), M.ptr((index+1),i), sizeof(scalartype)*L);
      }

      auto size = M.size();
      size.first--;
      M.resize(size);
    }

    template<typename scalartype>
    static void col(matrix<scalartype, CPU>& M, int index)
    {
      assert(index>-1 && index<M.size().second);

      int L = M.nrCols()-index-1;

      if(M.nrRows()>0 && L>0){

// 	scalartype* data = &M(0,0);

// 	int N = M.leadingDimension();
// 	int L = M.size().second-index-1;
	
// 	memmove(&data[index*N], &data[(index+1)*N], sizeof(scalartype)*N*L);

	int N = M.leadingDimension();
	memmove(M.ptr(0,index), M.ptr(0,index+1), sizeof(scalartype)*N*L);
      }

      auto size = M.size();
      size.second--;
      M.resize(size);
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
