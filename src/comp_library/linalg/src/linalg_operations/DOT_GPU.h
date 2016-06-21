//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_DOT_GPU_H
#define LINALG_DOT_GPU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_DOT {

    float  sdot(int length, float* a, int inc_a, float* b, int inc_b);
    double ddot(int length, double* a, int inc_a, double* b, int inc_b);
  }

  template<>
  class DOT<GPU>
    {
    public:

      inline static float execute(int length, float* a, int inc_a, float* b, int inc_b){
	return GPU_KERNEL_DOT::sdot(length, a, inc_a, b, inc_b);
      }
      
      inline static double execute(int length, double* a, int inc_a, double* b, int inc_b){
	return GPU_KERNEL_DOT::ddot(length, a, inc_a, b, inc_b);
      }
  
      /*    
      inline static void execute(int length, cuComplex* a, int inc_a, cuComplex* b, int inc_b){
	    cublasCswap(length, a, inc_a, b, inc_b);
      }
	
      inline static void execute(int length, cuDoubleComplex* a, int inc_a, cuDoubleComplex* b, int inc_b){
	cublasZswap(length, a, inc_a, b, inc_b);
      }
      */
  };
  
}

#endif
