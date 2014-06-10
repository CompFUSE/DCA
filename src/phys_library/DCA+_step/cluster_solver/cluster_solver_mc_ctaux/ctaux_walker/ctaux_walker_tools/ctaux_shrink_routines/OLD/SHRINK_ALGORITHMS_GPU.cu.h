//-*-C++-*-

#ifndef SHRINK_CU_H
#define SHRINK_CU_H

namespace QMC {
  
  namespace SHRINK_TOOLS_ALGORITHMS_GPU {

    /*
    const static int NUMBER_OF_THREADS = 128;

    __global__ void swap_rows_kernel(int N_i, int* i_s_ptr, int* i_t_ptr, double* A_ptr, int A_r, int A_LD)
    {
      int I = threadIdx.x + blockIdx.x*blockDim.x;

      if(I<A_r)
	{
	  for(int l=0; l<N_i; ++l){
	    int i0_j = i_s_ptr[l]+I*A_LD;
	    int i1_j = i_t_ptr[l]+I*A_LD;
	    
	    double t    = A_ptr[i0_j];
	    A_ptr[i0_j] = A_ptr[i1_j];
	    A_ptr[i1_j] = t;	
	  }
	}
    }

    void swap_rows(int Ns, int* i_s_ptr, int* i_t_ptr, double* A_ptr, int A_r, int A_LD)
    {
      int Nr_t = NUMBER_OF_THREADS; 
      int Nr_b = get_number_of_blocks(A_r, Nr_t);

      dim3 block(Nr_t);
      dim3 grid (Nr_b);

      swap_rows_kernel<<<grid, block>>>(Ns, i_s_ptr, i_t_ptr, A_ptr, A_r, A_LD);
    }

    __global__ void swap_cols_kernel(int N_i, int* i_s_ptr, int* i_t_ptr, double* A_ptr, int A_r, int A_LD)
    {
      int I = threadIdx.x + blockIdx.x*blockDim.x;

      if(I<A_r)
	{
	  for(int l=0; l<N_i; ++l){
	    int i_j0 = I+i_s_ptr[l]*A_LD;
	    int i_j1 = I+i_t_ptr[l]*A_LD;
	    
	    double t    = A_ptr[i_j0];
	    A_ptr[i_j0] = A_ptr[i_j1];
	    A_ptr[i_j1] = t;	
	  }
	}
    }

    void swap_cols(int Ns, int* i_s_ptr, int* i_t_ptr, double* A_ptr, int A_r, int A_LD)
    {
      int Nr_t = NUMBER_OF_THREADS; 
      int Nr_b = get_number_of_blocks(A_r, Nr_t);

      dim3 block(Nr_t);
      dim3 grid (Nr_b);

      swap_cols_kernel<<<grid, block>>>(Ns, i_s_ptr, i_t_ptr, A_ptr, A_r, A_LD);
    }
    */


    /*
    __global__ void swap_rows_kernel_2(int Ns, int* i_s_ptr, int* i_t_ptr, int M, int N, double* A_ptr, int A_LD)
    {
      int BLOCK_SIZE = blockDim.x;

      int l = threadIdx.x;
      int I = threadIdx.x + (blockIdx.x+0)*BLOCK_SIZE;

      int J_min = (blockIdx.y+0)*BLOCK_SIZE;
      int J_max = (blockIdx.y+1)*BLOCK_SIZE;

      I_max = min(I_max, Ns);
      J_max = min(J_max, N );

      if(I<Ns)
	{
	  int source_index[BLOCK_SIZE];
	  int target_index[BLOCK_SIZE];

	  __shared__ double T_ptr[BLOCK_SIZE*BLOCK_SIZE];

	  source_index[l] = i_s_ptr[I];
	  target_index[l] = i_t_ptr[I];

	  for(int j=J_min; j<J_max; ++j)
	    T_ptr[l+(j-J_min)*BLOCK_SIZE] = A_ptr[source_index[l]+j*A_LD];

	  for(int j=J_min; j<J_max; ++j)
	    A_ptr[target_index[l]+j*A_LD] = T_ptr[l+(j-J_min)*BLOCK_SIZE];
	  }
	}
    }

    void swap_rows_2(int Ns, int* i_s_ptr, int* i_t_ptr, int M, int N, double* A_ptr, int A_LD, int stream_id)
    {
      int Nr_t = 32; 

      int Nr_b_x = get_number_of_blocks(Ns, Nr_t);
      int Nr_b_y = get_number_of_blocks(N , Nr_t);

      dim3 block(Nr_t);
      dim3 grid (Nr_b_x, Nr_b_y);

      swap_rows_kernel_2<<<grid, block>>>(Ns, i_s_ptr, i_t_ptr, M, N, A_ptr, A_LD);
    }

    __global__ void swap_cols_kernel_2(int Ns, int* i_s_ptr, int* i_t_ptr, int M, int N, double* A_ptr, int A_LD)
    {
      int BLOCK_SIZE = blockDim.x;

      int I = threadIdx.x + (blockIdx.x+0)*BLOCK_SIZE;

      int L_min = (blockIdx.y+0)*BLOCK_SIZE;
      int L_max = (blockIdx.y+1)*BLOCK_SIZE;

      L_max = min(L_max, Ns);

      if(I<Ns)
	{
	  __shared__ int source_index[BLOCK_SIZE];
	  __shared__ int target_index[BLOCK_SIZE];

	  __shared__ double T_ptr[BLOCK_SIZE*BLOCK_SIZE];

	  source_index[l] = i_s_ptr[J];
	  target_index[l] = i_t_ptr[J];

	  for(int l=L_min; l<L_max; ++l)
	    T_ptr[l+(j-J_min)*BLOCK_SIZE] = A_ptr[i+source_index[l]*A_LD];

	  for(int j=J_min; j<J_max; ++j)
	    A_ptr[target_index[l]+j*A_LD] = T_ptr[l+(j-J_min)*BLOCK_SIZE];
	}
    }

    void swap_cols_2(int Ns, int* i_s_ptr, int* i_t_ptr, int M, int N, double* A_ptr, int A_LD, int stream_id)
    {
      int Nr_t = 32; 

      int Nr_b_x = get_number_of_blocks(M , Nr_t);
      int Nr_b_y = get_number_of_blocks(Ns, Nr_t);

      dim3 block(Nr_t);
      dim3 grid (Nr_b_x, Nr_b_y);

      swap_rows_kernel_2<<<grid, block>>>(Ns, i_s_ptr, i_t_ptr, M, N, A_ptr, A_LD);
    }
    */

    /*
    const static int BLOCK_SIZE_x = 32;
    const static int BLOCK_SIZE_y = 8;

    __global__ void swap_rows_kernel_1(int N_i, int* i_s_ptr, int* i_t_ptr, 
				       double* N_ptr , int N_r, int N_LD,
				       double* G0_ptr, int G0_r, int G0_LD)
    {
      int I = threadIdx.x + blockIdx.x*blockDim.x;
      
      if(I<N_r)
	{
	  {// rows
	    for(int l=0; l<N_i; ++l){
	      int i0_j = i_s_ptr[l]+I*N_LD;
	      int i1_j = i_t_ptr[l]+I*N_LD;

	      double t    = N_ptr[i0_j];
	      N_ptr[i0_j] = N_ptr[i1_j];
	      N_ptr[i1_j] = t;	
	    }
	  }

	  {// rows
	    for(int l=0; l<N_i; ++l){
	      int i0_j = i_s_ptr[l]+I*G0_LD;
	      int i1_j = i_t_ptr[l]+I*G0_LD;

	      double t     = G0_ptr[i0_j];
	      G0_ptr[i0_j] = G0_ptr[i1_j];
	      G0_ptr[i1_j] = t;	
	    }
	  }
	}
    }

    __global__ void swap_rows_kernel_2(int N_i, int* i_s_ptr, int* i_t_ptr, 
				       double* N_ptr , int N_r, int N_LD,
				       double* G0_ptr, int G0_r, int G0_LD)
    {
      assert(blockDim.x == BLOCK_SIZE_x);

      int J = threadIdx.x + blockIdx.x*BLOCK_SIZE_x;

      int l_MIN = BLOCK_SIZE_y*(blockIdx.y+0);
      int l_MAX = BLOCK_SIZE_y*(blockIdx.y+1);

      l_MIN = max(l_MIN, 0);
      l_MAX = min(l_MAX, N_i);

//       if(threadIdx.x==0 and blockIdx.x==1)
// 	printf("\t block %d, %d, %d \n", blockIdx.y, l_MIN, l_MAX);

      if(J<N_r)
	{
	  {// rows
	    for(int l=l_MIN; l<l_MAX; ++l)
	      {
		int i0_j = i_s_ptr[l]+J*N_LD;
		int i1_j = i_t_ptr[l]+J*N_LD;

		double t    = N_ptr[i0_j];
		N_ptr[i0_j] = N_ptr[i1_j];
		N_ptr[i1_j] = t;	
	      }
	  }

	  {// rows
	    //for(int l=0; l<N_i; ++l){
	    for(int l=l_MIN; l<l_MAX; ++l)
	      {
		int i0_j = i_s_ptr[l]+J*G0_LD;
		int i1_j = i_t_ptr[l]+J*G0_LD;
		
		double t     = G0_ptr[i0_j];
		G0_ptr[i0_j] = G0_ptr[i1_j];
		G0_ptr[i1_j] = t;	
	      }
	  }
	}
    }

    __global__ void swap_cols_kernel_1(int N_i, int* i_s_ptr, int* i_t_ptr, 
				       double* N_ptr , int N_r, int N_LD,
				       double* G0_ptr, int G0_r, int G0_LD)
    {
      int I = threadIdx.x + blockIdx.x*blockDim.x;

      if(I<N_r)
	{
	  {// columns
	    for(int l=0; l<N_i; ++l){
	      int i_j0 = I+i_s_ptr[l]*N_LD;
	      int i_j1 = I+i_t_ptr[l]*N_LD;
	      
	      double t    = N_ptr[i_j0];
	      N_ptr[i_j0] = N_ptr[i_j1];
	      N_ptr[i_j1] = t;	
	    }
	  }

	  {// columns
	    for(int l=0; l<N_i; ++l){
	      int i_j0 = I+i_s_ptr[l]*G0_LD;
	      int i_j1 = I+i_t_ptr[l]*G0_LD;
	      
	      double t     = G0_ptr[i_j0];
	      G0_ptr[i_j0] = G0_ptr[i_j1];
	      G0_ptr[i_j1] = t;	
	    }
	  }
	}
    }

    __global__ void swap_cols_kernel_2(int N_i, int* i_s_ptr, int* i_t_ptr, 
				       double* N_ptr , int N_r, int N_LD,
				       double* G0_ptr, int G0_r, int G0_LD)
    {
      assert(blockDim.x == BLOCK_SIZE_x);

      int I = threadIdx.x + blockIdx.x*BLOCK_SIZE_x;

      int l_MIN = BLOCK_SIZE_y*(blockIdx.y+0);
      int l_MAX = BLOCK_SIZE_y*(blockIdx.y+1);

      l_MIN = max(l_MIN, 0);
      l_MAX = min(l_MAX, N_i);

      if(I<N_r)
	{
	  {// columns
	    //for(int l=0; l<N_i; ++l){
	    for(int l=l_MIN; l<l_MAX; ++l)
	      {
		int i_j0 = I+i_s_ptr[l]*N_LD;
		int i_j1 = I+i_t_ptr[l]*N_LD;
		
		double t    = N_ptr[i_j0];
		N_ptr[i_j0] = N_ptr[i_j1];
		N_ptr[i_j1] = t;	
	      }
	  }

	  {// columns
	    //for(int l=0; l<N_i; ++l){
	    for(int l=l_MIN; l<l_MAX; ++l)  
	      {
		int i_j0 = I+i_s_ptr[l]*G0_LD;
		int i_j1 = I+i_t_ptr[l]*G0_LD;
		
		double t     = G0_ptr[i_j0];
		G0_ptr[i_j0] = G0_ptr[i_j1];
		G0_ptr[i_j1] = t;	
	      }
	  }
	}
    }

    void execute(int Ns, int* i_s_ptr, int* i_t_ptr, 
		 double* N_ptr , int N_r , int N_LD,
		 double* G0_ptr, int G0_r, int G0_LD)
    {
      assert(N_r==G0_r);

//       if(false)
// 	{
// 	  int Nr_t = BLOCK_SIZE_x;//NUMBER_OF_THREADS; 
// 	  int Nr_b = get_number_of_blocks(N_r, Nr_t);
	  
// 	  dim3 threads(Nr_t);
// 	  dim3 blocks (Nr_b);

// 	  swap_rows_kernel_1<<<blocks,threads>>>(Ns, i_s_ptr, i_t_ptr, 
// 						 N_ptr , N_r , N_LD,
// 						 G0_ptr, G0_r, G0_LD);
// 	}
//       else
	{
	  if(N_r>0 and Ns>0)
	    {
	      int bl_x = get_number_of_blocks(N_r, BLOCK_SIZE_x);
	      int bl_y = get_number_of_blocks(Ns , BLOCK_SIZE_y);
	      
	      dim3 threads(BLOCK_SIZE_x);
	      dim3 blocks (bl_x, bl_y);
	      
	      swap_rows_kernel_2<<<blocks,threads>>>(Ns, i_s_ptr, i_t_ptr, 
						     N_ptr , N_r , N_LD,
						     G0_ptr, G0_r, G0_LD);
	    }
	}
      
      //assert(cuda_check_for_errors(__FUNCTION__));

//       if(false)
// 	{
// 	  int Nr_t = BLOCK_SIZE_x;//NUMBER_OF_THREADS; 
// 	  int Nr_b = get_number_of_blocks(N_r, Nr_t);
	  
// 	  dim3 threads(Nr_t);
// 	  dim3 blocks (Nr_b);
	  
// 	  swap_cols_kernel_1<<<blocks,threads>>>(Ns, i_s_ptr, i_t_ptr, 
// 					       N_ptr , N_r , N_LD,
// 					       G0_ptr, G0_r, G0_LD);
// 	}
//       else
	{
	  if(N_r>0 and Ns>0)
	    {
	      int bl_x = get_number_of_blocks(N_r, BLOCK_SIZE_x);
	      int bl_y = get_number_of_blocks(Ns , BLOCK_SIZE_y);
	      
	      dim3 threads(BLOCK_SIZE_x);
	      dim3 blocks (bl_x, bl_y);
	      
	      swap_cols_kernel_2<<<blocks,threads>>>(Ns, i_s_ptr, i_t_ptr, 
						     N_ptr , N_r , N_LD,
						     G0_ptr, G0_r, G0_LD);
	      
	      //assert(cuda_check_for_errors(__FUNCTION__));
	    }
	}
    }
    */


    /*****************************************************
     *****************************************************
     ***
     ***       EFFICIENT SWAP IMPLEMENTATION           ***     
     ***
     *****************************************************
     *****************************************************/

    /*
    const static int BLOCK_SIZE_i = 32; // rows
    const static int BLOCK_SIZE_j = 8;  // cols

    __global__ void swap_rows_kernel_3(int N_i, int* i_s_ptr, int* i_t_ptr, double* A_ptr, int A_r, int A_LD)
    {
      assert(blockDim.x == BLOCK_SIZE_i);

      __shared__ double T[BLOCK_SIZE_i*BLOCK_SIZE_j];

      int t_id    = threadIdx.x;
      int i_array = threadIdx.x + blockIdx.x*BLOCK_SIZE_i; // index of the index-array

      int i_MIN, i_MAX, J_MIN, J_MAX;

      {// min/max of array
	i_MIN = BLOCK_SIZE_i*(blockIdx.x+0);
	i_MAX = BLOCK_SIZE_i*(blockIdx.x+1);
	
	i_MIN = max(i_MIN, 0  );
	i_MAX = min(i_MAX, N_i);
      }

      {// min/max of column-index
	J_MIN = BLOCK_SIZE_j*(blockIdx.y+0);
	J_MAX = BLOCK_SIZE_j*(blockIdx.y+1);
	
	J_MIN = max(J_MIN, 0  );
	J_MAX = min(J_MAX, A_r);
      }
      
      if(i_array<i_MAX)//i_MIN+t_id<i_MAX)
	{
	  int i_s = i_s_ptr[i_array];
	  int i_t = i_t_ptr[i_array];

	  for(int j=J_MIN; j<J_MAX; ++j)
	    T[t_id+(j-J_MIN)*BLOCK_SIZE_j] = A_ptr[i_s+j*A_LD];

	  for(int j=J_MIN; j<J_MAX; ++j)
	    A_ptr[i_s+j*A_LD] = A_ptr[i_t+j*A_LD];

	  for(int j=J_MIN; j<J_MAX; ++j)
	    A_ptr[i_t+j*A_LD] = T[t_id+(j-J_MIN)*BLOCK_SIZE_j];
	}
    }

    __global__ void swap_cols_kernel_3(int N_i, int* i_s_ptr, int* i_t_ptr, double* A_ptr, int A_r, int A_LD)
    {
      int I = threadIdx.x + blockIdx.x*BLOCK_SIZE_i;

      int l_MIN = BLOCK_SIZE_j*(blockIdx.y+0);
      int l_MAX = BLOCK_SIZE_j*(blockIdx.y+1);

      l_MIN = max(l_MIN, 0);
      l_MAX = min(l_MAX, N_i);

      if(I<A_r)
	{
	  for(int l=l_MIN; l<l_MAX; ++l)
	    {
	      int i_j0 = I+i_s_ptr[l]*A_LD;
	      int i_j1 = I+i_t_ptr[l]*A_LD;
	      
	      double t    = A_ptr[i_j0];
	      A_ptr[i_j0] = A_ptr[i_j1];
	      A_ptr[i_j1] = t;	
	    }
	}
    }

    void swap_rows_and_cols(int N_s, int* i_s_ptr, int* i_t_ptr, double* A_ptr , int A_r , int A_LD, int thread_id, int stream_id)
    {
      if(A_r>0 and N_s>0)
	{
	  int bl_x = get_number_of_blocks(N_s, BLOCK_SIZE_i);
	  int bl_y = get_number_of_blocks(A_r, BLOCK_SIZE_j);
	  
	  dim3 threads(BLOCK_SIZE_i);
	  dim3 blocks (bl_x, bl_y);
	  
	  swap_rows_kernel_3<<<blocks,threads>>>(N_s, i_s_ptr, i_t_ptr, 
						 A_ptr, A_r , A_LD);
	  
	  //assert(cuda_check_for_errors(__FUNCTION__));
	}
      
      if(A_r>0 and N_s>0)
	{
	  int bl_x = get_number_of_blocks(A_r, BLOCK_SIZE_i);
	  int bl_y = get_number_of_blocks(N_s, BLOCK_SIZE_j);
	  
	  dim3 threads(BLOCK_SIZE_i);
	  dim3 blocks (bl_x, bl_y);
	  
	  swap_cols_kernel_3<<<blocks,threads>>>(N_s, i_s_ptr, i_t_ptr, 
						 A_ptr , A_r , A_LD);
	  
	  //assert(cuda_check_for_errors(__FUNCTION__));
	}
    }
    */
  }
}

#endif
