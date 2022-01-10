// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file implements cublas related utilities.

#include "dca/linalg/util/error_gpuBLAS.hpp"
#include <stdexcept>
#include <string>

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

std::string errorStringCublas(cublasStatus_t error) {
  switch (error) {
    case CUBLAS_STATUS_SUCCESS:
      return "CUBLAS_STATUS_SUCCESS";
    case CUBLAS_STATUS_NOT_INITIALIZED:
      return "CUBLAS_STATUS_NOT_INITIALIZED";
    case CUBLAS_STATUS_ALLOC_FAILED:
      return "CUBLAS_STATUS_ALLOC_FAILED";
    case CUBLAS_STATUS_INVALID_VALUE:
      return "CUBLAS_STATUS_INVALID_VALUE";
    case CUBLAS_STATUS_ARCH_MISMATCH:
      return "CUBLAS_STATUS_ARCH_MISMATCH";
    case CUBLAS_STATUS_MAPPING_ERROR:
      return "CUBLAS_STATUS_MAPPING_ERROR";
    case CUBLAS_STATUS_EXECUTION_FAILED:
      return "CUBLAS_STATUS_EXECUTION_FAILED";
    case CUBLAS_STATUS_INTERNAL_ERROR:
      return "CUBLAS_STATUS_INTERNAL_ERROR";
    case CUBLAS_STATUS_NOT_SUPPORTED:
      return "CUBLAS_STATUS_NOT_SUPPORTED";
    default:
      return "UNKNOWN_CUBLAS_ERROR";
  }
}

void printErrorMessage(cublasStatus_t error, std::string function_name, std::string file_name,
                       int line, std::string extra_error_string) {
  auto cuda_error = cudaGetLastError();
  std::string cuda_error_str(cudaGetErrorString(cuda_error));
  cuda_error_str += extra_error_string;
  printErrorMessage(errorStringCublas(error) + " (" + std::to_string(error) + ")", function_name,
                    file_name, line, cuda_error_str);
}

}  // util
}  // linalg
}  // dca
