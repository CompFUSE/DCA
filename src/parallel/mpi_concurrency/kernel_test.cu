// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Test the execution of a GPU kernel.

#include "dca/parallel/mpi_concurrency/kernel_test.hpp"

#include <array>
#include <cuda.h>
#include <cuda_runtime.h>

__global__ void kernel(int* out) {
  out[threadIdx.x] = threadIdx.x;
}

namespace dca {
namespace parallel {
// dca::parallel::

bool kernelTest() {
  int* dev;
  std::array<int, 32> host;

  cudaMalloc(&dev, sizeof(int) * 32);
  kernel<<<1, 32>>>(dev);
  cudaMemcpy(host.data(), dev, sizeof(int) * 32, cudaMemcpyDeviceToHost);
  cudaFree(dev);

  for (int i = 0; i < 32; ++i)
    if (host[i] != i)
      return false;

  return cudaSuccess == cudaPeekAtLastError();
}

}  // parallel
}  // dca
