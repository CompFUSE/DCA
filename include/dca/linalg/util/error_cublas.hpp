// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides cublas related utilities to
// - return code checking,
// - error message printing.

#ifndef DCA_LINALG_UTIL_ERROR_CUBLAS_HPP
#define DCA_LINALG_UTIL_ERROR_CUBLAS_HPP

#include <cublas_v2.h>
#include <stdexcept>
#include <string>
#include "dca/linalg/util/error_cuda.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

// Returns the error string related to error.
std::string errorStringCublas(cublasStatus_t error);

// Prints an error message containing function_name, file_name, line, and the error message and
// error code related to error.
void printErrorMessage(cublasStatus_t error, std::string function_name, std::string file_name,
                       int line, std::string extra_error_string = "");

// Prints an error message and throws a std::logic_error if the return code of a cuda function is
// not CUBLAS_STATUS_SUCCESS.
// This function can be invoked with the macros checkRC and checkRCMsg (defined in error_cuda.hpp)
// that automatically include the function name, the filename, and the line to the function call.
inline void checkRCInternal(cublasStatus_t return_code, std::string function_name,
                            std::string file_name, int line, std::string extra_error_string = "") {
  if (return_code != CUBLAS_STATUS_SUCCESS) {
    printErrorMessage(return_code, function_name, file_name, line, extra_error_string);
    throw std::logic_error(function_name);
  }
}

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_ERROR_CUBLAS_HPP
