// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides cuda related utilities to device specs printing.

#ifndef DCA_LINALG_UTIL_INFO_CUDA_HPP
#define DCA_LINALG_UTIL_INFO_CUDA_HPP

#if defined(DCA_HAVE_CUDA)
#include <cuda_runtime.h>
#elif defined(DCA_HAVE_HIP)
#include <hip/hip_runtime.h>
#include "dca/util/cuda2hip.h"
#endif
#include <string>
#include <stdexcept>

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

// Prints the specs of the cuda devices found by the application.
void printInfoDevices();

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_INFO_CUDA_HPP
