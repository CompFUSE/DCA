//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_AXPY_GPU_H
#define LINALG_AXPY_GPU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_AXPY {

    void daxpy(int length, double alpha, double* a, int inc_a, double* b, int inc_b);
    void daxpy(int length, double alpha, double* a, int inc_a, double* b, int inc_b, int id);

  }

  template<>
  class AXPY<GPU>
    {
    public:
	
      /*
      inline static void execute(int length, float alpha, float* a, int inc_a, float* b, int inc_b){
	  GPU_KERNEL_AXPY::saxpy(length, alpha, a, inc_a, b, inc_b);
	}
      */

      inline static void execute(int length, double alpha, double* a, int inc_a, double* b, int inc_b, int thread_id, int stream_id)
      {
	assert(stream_id==0);
	GPU_KERNEL_AXPY::daxpy(length, alpha, a, inc_a, b, inc_b, thread_id);
      }

      inline static void execute(int length, double alpha, double* a, int inc_a, double* b, int inc_b, int id, int thread_id, int stream_id)
      {
	assert(stream_id==0);
	GPU_KERNEL_AXPY::daxpy(length, alpha, a, inc_a, b, inc_b, thread_id);
      }

      /*	
	inline static void execute(int length, cuComplex* a, int inc_a, cuComplex* b, int inc_b){
	  GPU_KERNEL_AXPY::execute(length, a, inc_a, b, inc_b);
	}
	
	inline static void execute(int length, cuDoubleComplex* a, int inc_a, cuDoubleComplex* b, int inc_b){
	  GPU_KERNEL_AXPY::execute(length, a, inc_a, b, inc_b);
	}
      */
    };

}

#endif
