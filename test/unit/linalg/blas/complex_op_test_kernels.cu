// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests whether the platform has working assignment operators *= etc.

#include <cassert>
#include "complex_op_test_kernels.hpp"
#include <magma_operators.h>
#include "dca/util/cuda2hip.h"

template<typename T>
__global__ void gpu_operator_opmult_kernel(T* a, const T* b) {
  if (threadIdx.x + blockIdx.x == 0)
    a[0] *= b[0];
}

template <typename Type>
void gpu_operator_opmult(Type* a, const Type* b, int thread_id, int stream_id) {
  const int block_size = 16;
  int threads = 1;
  int blocks = 1;
  cudaStream_t stream = dca::linalg::util::getStream(thread_id, stream_id);
  gpu_operator_opmult_kernel<<<dim3(blocks), dim3(threads), 0, 0>>> (a, b);
}

template void gpu_operator_opmult(magmaFloatComplex* a, const magmaFloatComplex* b, int thread_id , int stream_id);
template void gpu_operator_opmult(magmaDoubleComplex* a, const magmaDoubleComplex* b, int thread_id , int stream_id);
