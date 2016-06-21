//-*-C++-*-

#ifndef CUBLAS_MALLOC_H
#define CUBLAS_MALLOC_H

namespace CUDA_GPU
{
  /*!
   *   \author R. Solca, P. Staar
   */
 
  template<typename cuda_scalartype>
  void allocate_gpu(cuda_scalartype** A, int LD, int D)
  {
    cudaError_t ret = cudaMalloc( (void**)A, D*LD*sizeof(cuda_scalartype) );
    
    if( ret != cudaSuccess){
      std::cout << "NOT ENOUGH GPU MEMORY " << (D*LD)*sizeof(cuda_scalartype) 
		<< " ret code: " << (ret != cudaSuccess) << std::endl;
      abort();
    }
  }

  template<typename cuda_scalartype>
  void deallocate_gpu(cuda_scalartype* A)
  {
    cudaFree(A);
  }
}

#endif
