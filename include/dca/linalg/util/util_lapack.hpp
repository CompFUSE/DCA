// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides utils for LAPACK.

#ifndef DCA_LINALG_UTIL_UTIL_LAPACK_HPP
#define DCA_LINALG_UTIL_UTIL_LAPACK_HPP

#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sstream>
#include "dca/linalg/util/lapack_exception.hpp"

namespace dca {
namespace linalg {
namespace lapack {
namespace util {
// dca::linalg::util::

// Performs a safe cast to int for the size of workspaces returned by lapack routines.
inline int getWorkSize(float tmp) {
  return static_cast<int>(tmp * (1 + std::numeric_limits<float>::epsilon()));
}
inline int getWorkSize(double tmp) {
  return static_cast<int>(tmp * (1 + std::numeric_limits<double>::epsilon()));
}
template <typename Type>
inline int getWorkSize(std::complex<Type> tmp) {
  return getWorkSize(tmp.real());
}

// Prints an error message and throws a std::logic_error if info is not zero
// The macro provides the interface that automatically passes the function name, the filename, and
// the line to the function call.
#define checkLapackInfo(info) \
  dca::linalg::lapack::util::checkLapackInfoInternal(info, __FUNCTION__, __FILE__, __LINE__)
inline void checkLapackInfoInternal(int info, std::string function_name, std::string file_name,
                                    int line) {
  if (info != 0) {
    std::stringstream s;
    s << "Error in function: " << function_name << " (" << file_name << ":" << line << ")"
      << "\n";
    s << "The Lapack function returned info = " << info << std::endl;

    throw LapackException(s.str(), info);
  }
}

}  // util
}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_UTIL_LAPACK_HPP
