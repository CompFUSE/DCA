// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// TODO: This file has to be modified when the LAPACK functions are cleaned up.

#ifndef DCA_LINALG_UTIL_UTIL_CUBLAS_HPP
#define DCA_LINALG_UTIL_UTIL_CUBLAS_HPP

#include <cublas_v2.h>
#include <stdexcept>
#include <string>
#include "dca/linalg/util/error_cuda.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

int getCublasVersion();

void initializeMagma();

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_UTIL_CUBLAS_HPP
