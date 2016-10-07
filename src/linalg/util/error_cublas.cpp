// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file implements cublas related utilities.

#include "dca/linalg/util/error_cublas.hpp"
#include <cublas_v2.h>
#include <stdexcept>
#include <string>
#include "dca/linalg/util/error_cuda.hpp"

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
  printErrorMessage(errorStringCublas(error) + " (" + std::to_string(error) + ")", function_name,
                    file_name, line, extra_error_string);
}

}  // util
}  // linalg
}  // dca
