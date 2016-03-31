//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_BENNET_GPU_H
#define LINALG_BENNET_GPU_H

namespace LIN_ALG {

   namespace GPU_KERNEL_BENNETT {

     
     void standard_Bennet(int N, int LD, double* M, double* c, double* r);
   }

    template<>
    class BENNET<GPU>
    {
    public:

	template<typename scalartype>
	static void execute(matrix<scalartype, GPU>& M, scalartype* c, scalartype* r){
	    
	    assert(M.get_current_size().first == M.get_current_size().second);

	    int N  = M.get_current_size().first;
	    int LD = M.get_global_size() .first;

	    GPU_KERNEL_BENNETT::standard_Bennet(N, LD, M.get_ptr(), &c[0], &r[0]);
	}

      template<typename scalartype>
      static void standard_Bennet(int N, int LD, scalartype* M, scalartype* c, scalartype* r){
	GPU_KERNEL_BENNETT::standard_Bennet(N, LD, M, c, r);
      }
    };

}

#endif
