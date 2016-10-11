// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides cuda related utilities to device specs printing.

#ifndef DCA_LINALG_UTIL_INFO_CUDA_HPP
#define DCA_LINALG_UTIL_INFO_CUDA_HPP

#include <cuda_runtime.h>
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
