// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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

// Prints an error message and throws a std::logic_error if info in not zero
// The macro provides the interface that automatically pass the function name, the filename, and
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
