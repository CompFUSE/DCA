//-*-C++-*-

#ifndef BASIC_CUBLAS_FUNCTIONS_H
#define BASIC_CUBLAS_FUNCTIONS_H

void initialize_cublas(int N_threads)
{
  for(int l=0; l<N_threads; ++l)
    create_cublas_handle(l);
}

void finalize_cublas(int N_threads)
{
  for(int l=0; l<N_threads; ++l)
    destroy_cublas_handle(l);
}

void create_cublas_handle(int thread_id)
{
  if(CUDA_SUCCESS != cublasCreate(&get_thread_handle(thread_id)))
    {                        
      fprintf(stderr, "CUBLAS: Not initialized\n"); 
      exit(-1);               
    }    
  
  int* version=NULL;
  if(CUDA_SUCCESS != cublasGetVersion(&get_thread_handle(thread_id), version))
    {                        
      fprintf(stderr, "CUBLAS-verssion: Not initialized\n"); 
      exit(-1);               
    }    
  else
    cout << "\n\n\t\t CUBLAS-VERSION : " << *version << endl;
}

void destroy_cublas_handle(int thread_id)
{
  if(CUDA_SUCCESS != cublasDestroy(&get_thread_handle(thread_id)))
    {                        
      fprintf(stderr, "CUBLAS: Not destroyed\n"); exit(-1);               
    }
}

cublasHandle_t& get_thread_handle(int thread_id)
{
  switch(thread_id)
    {
    case 0 :
      {
	static cublasHandle_t handle;
	return handle;
      }
      break;

    case 1 :
      {
	static cublasHandle_t handle;
	return handle;
      }
      break;

    case 2 :
      {
	static cublasHandle_t handle;
	return handle;
      }
      break;

    case 3 :
      {
	static cublasHandle_t handle;
	return handle;
      }
      break;

    case 4 :
      {
	static cublasHandle_t handle;
	return handle;
      }
      break;

    case 5 :
      {
	static cublasHandle_t handle;
	return handle;
      }
      break;

    case 6 :
      {
	static cublasHandle_t handle;
	return handle;
      }
      break;

    case 7 :
      {
	static cublasHandle_t handle;
	return handle;
      }
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
  
  static cublasHandle_t handle;
  return handle;
};


cublasOperation_t cublas_operation_type(char ch)
{
  cublasOperation_t op;

  switch(ch)
    {
    case 'N' :
      op = CUBLAS_OP_N;
      break;

    case 'T' :
      op = CUBLAS_OP_T;
      break;

    case 'C' :
      op = CUBLAS_OP_C;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
  
  return op;
}

cublasFillMode_t cublas_triangle_type(char ch)
{
  cublasFillMode_t fill;

  switch(ch)
    {
    case 'L' :
      fill = CUBLAS_FILL_MODE_LOWER;
      break;

    case 'U' :
      fill = CUBLAS_FILL_MODE_UPPER;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
  
  return fill;
}

cublasDiagType_t cublas_diagonal_type(char ch)
{
  cublasDiagType_t diag;

  switch(ch)
    {
    case 'U' :
      diag = CUBLAS_DIAG_UNIT;
      break;

    case 'N' :
      diag = CUBLAS_DIAG_NON_UNIT;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
  
  return diag;  
}

cublasSideMode_t cublas_side_type(char ch)
{
  cublasSideMode_t side;

  switch(ch)
    {
    case 'L' :
      side = CUBLAS_SIDE_LEFT;
      break;

    case 'R' :
      side = CUBLAS_SIDE_RIGHT;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
  
  return diag;  
}

#endif
