//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GEMD_GPU_H
#define LINALG_GEMD_GPU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_GEMD {

    void execute(int m, int n, double* M, int LDM, double* D, double* A, int LDA,
		 int thread_id, int stream_id);
  }

  template<>
  class GEMD<GPU>
  {
  public:

    template<typename scalartype>
    static void execute(matrix<scalartype, GPU>& M, scalartype* D, matrix<scalartype, GPU>& A,
			int thread_id, int stream_id)
    {
      int m = M.get_current_size().first;
      int n = M.get_current_size().second;

      int LDM = M.get_global_size().first;
      int LDA = A.get_global_size().first;

      scalartype* M_ptr = M.get_ptr();
      scalartype* A_ptr = A.get_ptr();

      GPU_KERNEL_GEMD::execute(m, n, M_ptr, LDM, D, A_ptr, LDA, thread_id, stream_id);
    }

    template<typename scalartype>
    static void execute(std::pair<int, int> current_size, 
			scalartype* M_ptr, int LDM, 
			scalartype* D, 
			scalartype* A_ptr, int LDA,
			int thread_id, int stream_id)
    {
      GPU_KERNEL_GEMD::execute(current_size.first, current_size.second, M_ptr, LDM, D, A_ptr, LDA, thread_id, stream_id);
    }
  };
}

#endif
