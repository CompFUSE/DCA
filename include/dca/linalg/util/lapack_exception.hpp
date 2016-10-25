// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides LapackException.

#ifndef DCA_LINALG_UTIL_LAPACK_EXCEPTION_HPP
#define DCA_LINALG_UTIL_LAPACK_EXCEPTION_HPP

#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sstream>

namespace dca {
namespace linalg {
namespace lapack {
namespace util {
// dca::linalg::util::

class LapackException : public std::runtime_error {
public:
  LapackException(std::string msg, int info) : runtime_error(msg), info_(info) {}

  int info() const {
    return info_;
  }

private:
  int info_;
};

}  // util
}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_LAPACK_EXCEPTION_HPP
