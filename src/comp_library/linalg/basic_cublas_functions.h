//-*-C++-*-

#ifndef BASIC_CUBLAS_FUNCTIONS_H
#define BASIC_CUBLAS_FUNCTIONS_H

namespace LIN_ALG {

  // cuBLAS API errors
  const char* cublas_error_str(cublasStatus_t error);

  void cublas_error_msg(cublasStatus_t error, std::string function_name, std::string file_name, int line);

  int get_version(int thread_id);  

  /******************************************
   ***   IMPLEMENTATION CUBLAS-HANDLES    ***
   ******************************************/

  cublasHandle_t& get_thread_handle(int thread_id);

  void create_cublas_handle(int thread_id);

  void destroy_cublas_handle(int thread_id);

  /******************************************
   ***   IMPLEMENTATION STREAM-HANDLES    ***
   ******************************************/

  cudaStream_t&   get_stream_handle(int thread_id, int stream_id);

  void create_stream_handle(int thread_id);

  void destroy_stream_handle(int thread_id);
    
  void synchronize_stream_handle(int thread_id, int stream_id);

  /******************************************
   ***   IMPLEMENTATION MATRIX-TYPES      ***
   ******************************************/

  cublasOperation_t cublas_operation_type(char ch);
  
  cublasFillMode_t cublas_triangle_type(char ch);

  cublasDiagType_t cublas_diagonal_type(char ch);
  
  cublasSideMode_t cublas_side_type(char ch);

  /****************************
   ***   IMPLEMENTATION     ***
   *****************************/

  const char* cublas_error_str(cublasStatus_t error)
  {
    switch (error)
      {
      case CUBLAS_STATUS_SUCCESS:                    return "CUBLAS_STATUS_SUCCESS";
      case CUBLAS_STATUS_NOT_INITIALIZED:            return "CUBLAS_STATUS_NOT_INITIALIZED";
      case CUBLAS_STATUS_ALLOC_FAILED:               return "CUBLAS_STATUS_ALLOC_FAILED";
      case CUBLAS_STATUS_INVALID_VALUE:              return "CUBLAS_STATUS_INVALID_VALUE";
      case CUBLAS_STATUS_ARCH_MISMATCH:              return "CUBLAS_STATUS_ARCH_MISMATCH";
      case CUBLAS_STATUS_MAPPING_ERROR:              return "CUBLAS_STATUS_MAPPING_ERROR";
      case CUBLAS_STATUS_EXECUTION_FAILED:           return "CUBLAS_STATUS_EXECUTION_FAILED";
      case CUBLAS_STATUS_INTERNAL_ERROR:             return "CUBLAS_STATUS_INTERNAL_ERROR";
    }
    return "<unknown>";
  }

  void cublas_error_msg(cublasStatus_t error, std::string function_name, std::string file_name, int line)
  {
    std::stringstream ss;

    ss << "\n\n\t error in function : " << function_name 
       << "\t file-name : "             << file_name 
       << "\t at line : "               << line 
       << "\t error : " << cublas_error_str(error) << "\n\n";

    std::cout << ss.str();
  }

  void create_cublas_handle(int thread_id)
  {
    if(CUBLAS_STATUS_SUCCESS != cublasCreate(&get_thread_handle(thread_id)))
      {
	fprintf(stderr, "CUBLAS: Not initialized\n");
	exit(-1);
      }
  }

  void destroy_cublas_handle(int thread_id)
  {
    if(CUBLAS_STATUS_SUCCESS != cublasDestroy(get_thread_handle(thread_id)))
      {
	fprintf(stderr, "CUBLAS: Not destroyed\n"); exit(-1);
      }
  }

  int get_version(int thread_id)
  {
    int version;
    if(CUBLAS_STATUS_SUCCESS != cublasGetVersion(get_thread_handle(thread_id), &version))
      {
	fprintf(stderr, "CUBLAS-verssion: Not initialized\n");
	exit(-1);
      }
    
    return version;
  }

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
    
    return side;
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
	
      case 1 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}
	
      case 2 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}

      case 3 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}
	
      case 4 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}
	
      case 5 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}

      case 6 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}
	
      case 7 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}

      case 8 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}
	
      case 9 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}
	
      case 10 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}

      case 11 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}
	
      case 12 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}
	
      case 13 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}

      case 14 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}
	
      case 15 :
	{
	  static cublasHandle_t handle;
	  return handle;
	}
	
      default:
	throw std::logic_error(__FUNCTION__);
      }
    
  }

  cudaStream_t& get_stream_handle(int thread_id, int stream_id)
  {
    assert(thread_id>-1 and thread_id<16);
    assert(stream_id>-1 and stream_id<8);

    switch(thread_id)
      {
      case 0 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      case 1 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      case 2 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}

      case 3 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      case 4 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      case 5 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}

      case 6 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      case 7 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}

      case 8 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      case 9 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      case 10 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}

      case 11 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      case 12 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      case 13 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}

      case 14 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      case 15 :
	{
	  static cudaStream_t stream[8];
	  return stream[stream_id];
	}
	
      default:
	throw std::logic_error(__FUNCTION__);
      }
  }

  void create_stream_handle(int thread_id)
  {
    for(int stream_id=0; stream_id<8; ++stream_id)
      cudaStreamCreate(&get_stream_handle(thread_id, stream_id));
  }

  void destroy_stream_handle(int thread_id)
  {
    for(int stream_id=0; stream_id<8; ++stream_id)
      cudaStreamDestroy(get_stream_handle(thread_id, stream_id));
  }

  void synchronize_stream_handle(int thread_id, int stream_id)
  {
    cudaStreamSynchronize(get_stream_handle(thread_id, stream_id));
  }

  void link_thread_stream_to_cublas_handle(int thread_id, int stream_id)
  {
    cublasStatus_t error = cublasSetStream(get_thread_handle(thread_id), 
					   get_stream_handle(thread_id, stream_id));

    if(error != CUBLAS_STATUS_SUCCESS)
      {
	std::cout << cublas_error_str(error) << std::endl;
	throw std::logic_error(__FUNCTION__);
      }      
  }
}

#endif











