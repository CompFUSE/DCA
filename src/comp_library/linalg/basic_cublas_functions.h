//-*-C++-*-

#ifndef BASIC_CUBLAS_FUNCTIONS_H
#define BASIC_CUBLAS_FUNCTIONS_H

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cublas_v2.h>
#include "comp_library/linalg/basic_cuda_functions.h"

namespace LIN_ALG {

  void cublasCheckReturnCodeInternal(cublasStatus_t return_code, std::string function_name, std::string file_name, int line);
#define cublasCheckReturnCode(return_code) LIN_ALG::cublasCheckReturnCodeInternal(return_code, __FUNCTION__, __FILE__, __LINE__); 

  const char* cublas_error_str(cublasStatus_t error);

  void cublas_error_msg(cublasStatus_t error, std::string function_name, std::string file_name, int line);

  void create_cublas_handle(int thread_id);

  void destroy_cublas_handle(int thread_id);

  int get_version(int thread_id);

  cublasOperation_t cublas_operation_type(char ch);

  cublasFillMode_t cublas_triangle_type(char ch);

  cublasDiagType_t cublas_diagonal_type(char ch);
  
  cublasSideMode_t cublas_side_type(char ch);

  cublasHandle_t& get_thread_handle(int thread_id);

  cudaStream_t& get_stream_handle(int thread_id, int stream_id);

  void create_stream_handle(int thread_id);

  void destroy_stream_handle(int thread_id);

  void synchronize_stream_handle(int thread_id, int stream_id);

  void link_thread_stream_to_cublas_handle(int thread_id, int stream_id);
}

#endif
