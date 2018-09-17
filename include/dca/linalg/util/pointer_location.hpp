// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides an utility to find if a pointer belongs to the host or the device.

#ifndef DCA_LINALG_UTIL_POINTER_LOCATION_HPP
#define DCA_LINALG_UTIL_POINTER_LOCATION_HPP

#ifndef DCA_HAVE_CUDA
#error "This file requires CUDA support."
#endif  // DCA_HAVE_CUDA

#include <cuda_runtime.h>

#include "dca/linalg/device_type.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <typename T>
DeviceType pointerLocation(const T* ptr) {
  cudaPointerAttributes attributes;
  cudaPointerGetAttributes(&attributes, ptr);

  return attributes.memoryType == cudaMemoryTypeHost ? CPU : GPU;
}

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_POINTER_LOCATION_HPP
