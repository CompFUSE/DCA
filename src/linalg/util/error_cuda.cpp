// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file implements error_cuda functions.

#include "dca/linalg/util/error_cuda.hpp"
#include <cuda_runtime.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

void checkErrorsCudaDebugInternal(std::string function_name, std::string file_name, int line) {
  // cudaDeviceSynchronize();

  cudaError_t ret = cudaGetLastError();

  if (ret != cudaSuccess) {
    std::stringstream s;

    s << std::endl;
    s << "DEBUG_CUDA: error in function: " << function_name;
    s << " (" << file_name << ":" << line << ")" << std::endl;
    s << "cudaGetLastError returned: " << cudaGetErrorString(ret) << " (" << ret << ")" << std::endl;

    std::cout << s.str() << std::endl;

    throw std::logic_error(s.str());
  }
}

void printErrorMessage(cudaError_t error, std::string function_name, std::string file_name,
                       int line, std::string extra_error_string) {
  printErrorMessage(std::string(cudaGetErrorString(error)) + " (" + std::to_string(error) + ")",
                    function_name, file_name, line, extra_error_string);
}

void printErrorMessage(std::string error, std::string function_name, std::string file_name,
                       int line, std::string extra_error_string) {
  std::stringstream s;

  s << "Error in function: " << function_name;
  s << " (" << file_name << ":" << line << ")" << std::endl;
  s << "The function returned: " << error << std::endl;
  if (extra_error_string != "")
    s << extra_error_string << std::endl;

  std::cout << s.str() << std::endl;
}
}  // util
}  // linalg
}  // dca
