// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#include "dca/phys/dca_step/cluster_solver/ctint/walker/kernels_interface.hpp"

#include <array>
#include <cassert>
#include <cuda.h>
#include <cuda_runtime.h>

#include "dca/util/integer_division.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {

std::array<dim3, 2> getBlockSize(const int i, const int j) {
  assert(i > 0 && j > 0);
  const int n_threads_i = std::min(32, i);
  const int n_threads_j = std::min(32, j);
  const int n_blocks_i = util::ceilDiv(i, n_threads_i);
  const int n_blocks_j = util::ceilDiv(j, n_threads_j);

  return std::array<dim3, 2>{dim3(n_blocks_i, n_blocks_j), dim3(n_threads_i, n_threads_j)};
}

__global__ void setRightSectorToIdKernel(double* m, const int ldm, const int n0, const int n_max) {
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y + n0;

  if (i >= n_max || j >= n_max)
    return;

  m[i + ldm * j] = (i == j) ? 1. : 0.;
}

void setRightSectorToId(double* m, const int ldm, const int n0, const int n_max, cudaStream_t stream) {
  auto blocks = getBlockSize(n_max, n_max - n0);

  setRightSectorToIdKernel<<<blocks[0], blocks[1], 0, stream>>>(m, ldm, n0, n_max);
}

}  // details
}  // ctint
}  // solver
}  // phys
}  // dca