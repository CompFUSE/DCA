// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides some utilities to test simple GPU operations.

#ifndef DCA_TEST_UNIT_LINALG_GPU_TEST_UTIL_HPP
#define DCA_TEST_UNIT_LINALG_GPU_TEST_UTIL_HPP

#include <cuda_runtime.h>

namespace testing {
template <typename ScalarType>
cudaMemoryType PointerType(const ScalarType* ptr) {
  cudaPointerAttributes attributes;
  // TODO: check return code.
  cudaError_t ret = cudaPointerGetAttributes(&attributes, reinterpret_cast<const void*>(ptr));
  // If the pointer is not managed by CUDA, assume it is a Host memory pointer.
  if (ret == cudaErrorInvalidValue)
    return cudaMemoryTypeHost;
  return attributes.memoryType;
}

template <typename ScalarType>
bool isDevicePointer(const ScalarType* ptr) {
  return cudaMemoryTypeDevice == PointerType(ptr);
}
template <typename ScalarType>
bool isHostPointer(const ScalarType* ptr) {
  return cudaMemoryTypeHost == PointerType(ptr);
}

template <typename ScalarType>
void setOnDevice(ScalarType* ptr, ScalarType value) {
  // TODO: check return code.
  cudaMemcpy(ptr, &value, sizeof(ScalarType), cudaMemcpyDefault);
}

template <typename ScalarType>
ScalarType getFromDevice(ScalarType* ptr) {
  ScalarType value;
  // TODO: check return code.
  cudaMemcpy(&value, ptr, sizeof(ScalarType), cudaMemcpyDefault);
  return value;
}

}  // testing

#endif  // DCA_TEST_UNIT_LINALG_GPU_TEST_UTIL_HPP
