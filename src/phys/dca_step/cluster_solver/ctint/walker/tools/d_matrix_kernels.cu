// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the device methods of DMatrixBuilder.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"

#include <cuda.h>

#include "dca/linalg/util/error_cuda.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/device_helper/ctint_helper.cuh"
#include "dca/util/cuda_blocks.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

__global__ void buildG0MatrixKernel(linalg::MatrixView<double, linalg::GPU> G0, const int n_init,
                                    const bool right_section, DeviceConfiguration config,
                                    DeviceInterpolationData g0_interp) {
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;
  if (id_i >= G0.nrRows() || id_j >= G0.nrCols())
    return;
  int i(id_i);
  int j(id_j);

  if (right_section)
    j += n_init;
  else
    i += n_init;

  const int b_i = config.getLeftB(i);
  const double tau_i = config.getTau(i);

  const int b_j = config.getRightB(j);
  const double tau_j = config.getTau(j);

  const int label = ctint_helper.index(b_i, b_j, config.getLeftR(i), config.getRightR(j));

  G0(id_i, id_j) = g0_interp(tau_i - tau_j, label);
}

void buildG0Matrix(linalg::MatrixView<double, linalg::GPU> G0, const int n_init,
                   const bool right_section, DeviceConfiguration config,
                   DeviceInterpolationData g0_interp, cudaStream_t stream) {
  // assert(CtintHelper::is_initialized());
  const auto blocks = dca::util::getBlockSize(G0.nrRows(), G0.nrCols());

  buildG0MatrixKernel<<<blocks[0], blocks[1], 0, stream>>>(G0, n_init, right_section, config,
                                                           g0_interp);
  checkErrorsCudaDebug();
}

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
