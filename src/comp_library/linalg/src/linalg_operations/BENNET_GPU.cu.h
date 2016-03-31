//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_BENNET_GPU_CU_H
#define LINALG_BENNET_GPU_CU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_BENNETT {

    __global__ void Bennet_kernel(int N, int LD, double* M, double* c, double* r){

	int b_id_i = blockIdx.x;
	int b_id_j = blockIdx.y;

	int t_id_i = threadIdx.x;
	int t_id_j = threadIdx.y;

	int I = t_id_i+b_id_i*blockDim.x;
	int J = t_id_j+b_id_j*blockDim.y;

	for(int i=0; i<N; ++i){

	    if(I==i && J==i){
		M[I+I*LD] += c[I]*r[I];
		r[I]      /= M[I+I*LD];
	    }
	    syncthreads();
    
	    if(I>i && J==i && I<N && J<N){
		c[I]      -= M[I+J*LD]*c[J];
		M[I+J*LD] += c[I]*r[J]; 
	    }
	    //syncthreads();

	    if(I==i && J>i && I<N && J<N){
		M[I+J*LD] += c[I]*r[J]; 
		r[J]      -= r[I]*M[I+J*LD];
	    }
	    //syncthreads();
	}
    }

    __global__ void update_the_diagonal(int I, int LD, double* M, double* c, double* r){

	M[I+I*LD] += c[I]*r[I];
	r[I]      /= M[I+I*LD];
       
    }

    __global__ void update_the_row(int i, int N, int LD, double* M, double* c, double* r){
      
      int I = threadIdx.x + blockIdx.x*blockDim.x;
      int J = i;
      
      if(I>i && I<N){
	c[I]      -= M[I+J*LD]*c[J];
	M[I+J*LD] += c[I]*r[J]; 
      }
    }

    __global__ void update_the_col(int i, int N, int LD, double* M, double* c, double* r){

	int I=i;
	int J = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(J>i && J<N){
	    M[I+J*LD] += c[I]*r[J]; 
	    r[J]      -= r[I]*M[I+J*LD];
	}
    }

    void standard_Bennet(int N, int LD, double* M, double* c, double* r)
    {
	if(N<32){
	
	    int Nr_t = 32;
	    int Nr_b = 1;
	    
	    dim3 blocks (Nr_b,Nr_b);
	    dim3 threads(Nr_t,Nr_t);
	
	    Bennet_kernel<<<blocks,threads>>>(N, LD, M, c, r);
	}
	else{

	    int Nr_t = 32;
	    int Nr_b = N/Nr_t;
	    
	    if(N>Nr_b*Nr_t)
		Nr_b += 1;

	    dim3 b_d(1);
	    dim3 t_d(1);

	    dim3 b_r(Nr_b);
	    dim3 t_r(Nr_t);

	    dim3 b_c(Nr_b);
	    dim3 t_c(Nr_t);

	    for(int i=0; i<N; ++i){

	      update_the_diagonal<<<b_d,t_d>>>(i, LD, M, c, r);

	      update_the_row<<<b_r,t_r>>>(i, N, LD, M, c, r);

	      update_the_col<<<b_c,t_c>>>(i, N, LD, M, c, r);
	    }
	}
    }//end: standard_Bennet
    
  }//namespace LIN_ALG 
}//namespace GPU_KERNEL_BENNETT


#endif
