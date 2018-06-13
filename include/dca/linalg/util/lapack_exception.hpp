// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
// dca::linalg::lapack::util::

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
