//-*-C++-*-

#ifndef LINALG_GEMD_CPU_H
#define LINALG_GEMD_CPU_H

namespace LIN_ALG {

  template<>
  class GEMD<CPU>
  {
  public:

    template<typename scalartype>
    static void execute(matrix<scalartype, CPU>& M, scalartype* D, matrix<scalartype, CPU>& A,
                        int /*thread_id*/=0, int /*stream_id*/=0)
    {      
      int m = M.size().first;
      int n = M.size().second;

      int LDM = M.leadingDimension();
      int LDA = A.leadingDimension();

      scalartype* A_ptr = A.ptr();
      scalartype* M_ptr = M.ptr();

      for(int j=0; j<n; ++j){

        scalartype  D_j     = D[j];
        scalartype* A_ptr_j = &A_ptr[j*LDA];
        scalartype* M_ptr_j = &M_ptr[j*LDM];

        for(int i=0; i<m; ++i)
          A_ptr_j[i] = M_ptr_j[i]*D_j;
      }
    }

    template<typename scalartype>
    static void execute(scalartype* D, matrix<scalartype, CPU>& M, matrix<scalartype, CPU>& A,
                        int /*thread_id*/=0, int /*stream_id*/=0)
    {
      int m = M.size().first;
      int n = M.size().second;

      int LDM = M.leadingDimension();
      int LDA = A.leadingDimension();

      scalartype* A_ptr = A.ptr();
      scalartype* M_ptr = M.ptr();

      for(int j=0; j<n; ++j){

        scalartype* A_ptr_j = &A_ptr[j*LDA];
        scalartype* M_ptr_j = &M_ptr[j*LDM];

        for(int i=0; i<m; ++i)
          A_ptr_j[i] = M_ptr_j[i]*D[i];
      }
    }


    template<typename scalartype>
    static void execute(std::pair<int, int> current_size,
                        scalartype* M_ptr, int LDM,
                        scalartype* D,
                        scalartype* A_ptr, int LDA,
                        int /*thread_id*/=0, int /*stream_id*/=0)
    {
      int m = current_size.first;
      int n = current_size.second;

      for(int j=0; j<n; ++j){

        scalartype  D_j     = D[j];
        scalartype* A_ptr_j = &A_ptr[j*LDA];
        scalartype* M_ptr_j = &M_ptr[j*LDM];

        for(int i=0; i<m; ++i)
          A_ptr_j[i] = M_ptr_j[i]*D_j;
      }
    }

  };

}

#endif
