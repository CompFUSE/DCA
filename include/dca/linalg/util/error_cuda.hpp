// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides cuda related utilities to
// - return code checking,
// - error message printing.

#ifndef DCA_LINALG_UTIL_ERROR_CUDA_HPP
#define DCA_LINALG_UTIL_ERROR_CUDA_HPP

#include <cuda_runtime.h>
#include <string>
#include <stdexcept>

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

// Performs extra cuda code checking when debugging.
// Remark: Use the macro checkErrorsCudaDebug instead of the function in the code to avoid overhead
// when DEBUG_CUDA is not defined.
void checkErrorsCudaDebugInternal(std::string function_name, std::string file_name, int line);
#ifdef DEBUG_CUDA
#define checkErrorsCudaDebug() \
  dca::linalg::util::checkErrorsCudaDebugInternal(__FUNCTION__, __FILE__, __LINE__)
#else
#define checkErrorsCudaDebug()
#endif

// Prints an error message containing error, function_name, file_name, line and extra_error_string.
void printErrorMessage(std::string error, std::string function_name, std::string file_name,
                       int line, std::string extra_error_string = "");

// Prints an error message containing function_name, file_name, line, extra_error_string, and the
// error message and error code related to the cudaError_t error.
void printErrorMessage(cudaError_t error, std::string function_name, std::string file_name,
                       int line, std::string extra_error_string = "");

// Prints an error message and throws a std::logic_error if the return code of a cuda function is
// not cudaSuccess.
// The macros provide the interfaces that automatically pass the function name, the filename, and
// the line to the function call.
#define checkRC(return_code) \
  dca::linalg::util::checkRCInternal(return_code, __FUNCTION__, __FILE__, __LINE__)
#define checkRCMsg(return_code, extra_error_string)                                 \
  dca::linalg::util::checkRCInternal(return_code, __FUNCTION__, __FILE__, __LINE__, \
                                     extra_error_string)
inline void checkRCInternal(cudaError_t return_code, std::string function_name,
                            std::string file_name, int line, std::string extra_error_string = "") {
  if (return_code != cudaSuccess) {
    printErrorMessage(return_code, function_name, file_name, line, extra_error_string);
    throw std::logic_error(function_name);
  }
}

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_ERROR_CUDA_HPP
