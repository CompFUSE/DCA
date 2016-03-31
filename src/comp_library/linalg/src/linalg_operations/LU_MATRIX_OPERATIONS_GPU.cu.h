//-*-C++-*-

#ifndef LINALG_LU_MATRIX_OPERATIONS_GPU_CU_H
#define LINALG_LU_MATRIX_OPERATIONS_GPU_CU_H

namespace LIN_ALG {

  namespace LU_MATRIX_OPERATIONS_GPU {

    template<class T>
    __global__ void det_tri_small_size(int n, T* A, int lda, T* res)
    {
      extern __shared__ double2 shared_array[];
      T* diag = (T*)(shared_array);

      int id = threadIdx.x;
      int id2 = blockDim.x+id;

      diag[id] = A[id*(lda+1)];
      if(id2 < n)
	diag[id2] = A[id2*(lda+1)];
      __syncthreads();

      int i = n;
      while(n > 1){
	i = n / 2;
	n -= i;
	if (id < i)
	  diag[id] *= diag[id+n];
	__syncthreads();
      }

      if (id==0)
	*res = diag[id];
    }

    template<class T>
    T determinant_tridiagonal(int n, T* A, int lda)
    {
      T* res;
      cudaMalloc(&res, sizeof(T));
      int nt = (n+1)/2;

      if(nt <= get_device_properties().maxThreadsPerBlock)
	det_tri_small_size <<< 1, nt, n * sizeof(T)>>>(n, A, lda, res);
      else
	throw; 

      T ret;
      cudaMemcpy(&ret, res, sizeof(T), cudaMemcpyDeviceToHost);
      cudaFree(res);
      return ret;
    }

    template float  determinant_tridiagonal(int n, float * A, int lda);
    template double determinant_tridiagonal(int n, double* A, int lda);

    template<class T>
    __global__ void minmax_diag_small_size(int n, T* A, int lda, T* res)
    {
      extern __shared__ double2 shared_array[];

      T* diag = (T*)(shared_array);
      T* diag2 = diag + n;

      int id = threadIdx.x;
      int id2 = blockDim.x+id;

      diag2[id] = abs(A[id*(lda+1)]);
      diag[id] = diag2[id];

      if(id2 < n){
	diag2[id2] = abs(A[id2*(lda+1)]);
	diag[id2] = diag2[id2];
      }
      __syncthreads();

      int i = n;
      while(n>1){
	i = n/2;
	n -= i;
	if (id < i){
	  if (diag[id+n] > diag[id])
	    diag[id] = diag[id+n];
	  if (diag2[id+n] < diag2[id])
	    diag2[id] = diag2[id+n];
	}
	__syncthreads();
      }

      if (id==0){
	res[0] = diag[id];
	res[1] = diag2[id];
      }
    }

    template<class T>
    std::pair<T,T> minmax_diagonal(int n, T* A, int lda)
    {
      T* res;

      cudaMalloc(&res, 2*sizeof(T));
      int nt = (n+1)/2;

      if(nt <= get_device_properties().maxThreadsPerBlock)
	minmax_diag_small_size <<< 1, nt, 2*n*sizeof(T)>>>(n, A, lda, res);
      else
	throw;  

      std::pair<T,T> ret;
      cudaMemcpy(&ret.first, res, 2*sizeof(T), cudaMemcpyDeviceToHost);
      cudaFree(res);

      return ret;
    }

    template std::pair<float ,float > minmax_diagonal(int n, float*  A, int lda);
    template std::pair<double,double> minmax_diagonal(int n, double* A, int lda);

  }

}

#endif
