//-*-C++-*-

#ifndef QMC_MATRIX_KERNELS_H
#define QMC_MATRIX_KERNELS_H

template<typename scalartype>
__global__ void set_to_zero_kernel(scalartype* A, int size) 
{ 
    int i = threadIdx.x; 

    A[i] = scalartype(0); 
} 

#endif
