//-*-C++-*-

#ifndef GAMMA_MATRIX_TOOLS_GPU_CU_H
#define GAMMA_MATRIX_TOOLS_GPU_CU_H

namespace QMC {
  
  namespace GAMMA_MATRIX_TOOLS_GPU_KERNELS {

    __global__ void compute_beta_kernel(int n, double* Gamma_LU_ptr, int LD, double corr_term){

      int l = threadIdx.x;

      __shared__ double cache[256];

      cache[l] = 0.;

      if(l<n-1)
	cache[l] = Gamma_LU_ptr[(n-1)+l*LD]*Gamma_LU_ptr[l+(n-1)*LD];
      __syncthreads();

      int i = blockDim.x/2;

      while(i!=0) {

	if(l<i)
	  cache[l] += cache[l+i];
	
	__syncthreads();
	
	i /= 2; 
      }

      if(l==0)
	Gamma_LU_ptr[(n-1)+(n-1)*LD] -= (cache[0]+corr_term);
    }

    void compute_beta(int n, double* Gamma_LU_ptr, int LD, double corr_term){

      assert(n<256);
      assert(n<get_number_of_threads());
      
      compute_beta_kernel<<<1,256>>>(n, Gamma_LU_ptr, LD, corr_term);
    }
  }

  /*
  __device__ double solve_Gamma_test(int n, double* Gamma_LU, int LD, double exp_delta_V)
  {
    int id = threadIdx.x;
    
    {
      {
	double* y = Gamma_LU.get_ptr(0,n-1);

	for(int i=0; i<n-1; i++)
	  for(int j=0; j<i; j++)
	    y[i] -= Gamma_LU(i,j)*y[j];
      }

      {
	double* x = Gamma_LU.get_ptr(n-1,0);

	for(int j=0; j<n-1; j++){	

	  for(int i=0; i<j; i++)
	    x[j*LD] -= x[i*LD]*Gamma_LU(i,j);

	  x[j*LD] /= Gamma_LU(j,j);
	}
      }

      for(int i=0; i<n-1; i++)
	Gamma_LU(n-1, n-1) -= Gamma_LU(n-1,i)*Gamma_LU(i,n-1);
    }

    if(id==0)
      {
	{
	  double max = fabs(Gamma_LU(0,0));
	  double min = fabs(Gamma_LU(0,0));
	  
	  for(int i=1; i<n; i++){
	    max = fabs(Gamma_LU(i,i))>max? fabs(Gamma_LU(i,i)) : max;
	    min = fabs(Gamma_LU(i,i))<min? fabs(Gamma_LU(i,i)) : min;
	  }
	  
	  if( (max/min) > 1.e6)
	    return 0;
	}
	
	double phani_gamma       = exp_delta_V-1.;
	double determinant_ratio = -phani_gamma*Gamma_LU(i,i);
	
	return determinant_ratio;
      }
  }
  */
}

#endif
