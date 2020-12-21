// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the device methods of G0Interpolation<GPU>.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"

#include <cuda_runtime.h>

#include "dca/linalg/util/error_cuda.hpp"
#include "dca/util/cuda_blocks.hpp"
#include "../../../../../../../include/dca/phys/dca_step/cluster_solver/ctint/walker/tools/device_interpolation_data.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details

template <typename Real>
__global__ void interpolateSlowKernel(Real tau, const int lindex, DeviceInterpolationData<Real> g0,
                                      Real* result) {
  *result = g0(tau, lindex);
}

template <typename Real>
Real interpolateSlow(Real tau, int lindex, const DeviceInterpolationData<Real>& g0) {
  Real* d_result;
  Real result;
  cudaMalloc((void**)&d_result, sizeof(Real));

  interpolateSlowKernel<<<1, 1>>>(tau, lindex, g0, d_result);

  assert(cudaSuccess == cudaPeekAtLastError());
  cudaMemcpy(&result, d_result, sizeof(Real), cudaMemcpyDeviceToHost);
  cudaFree(d_result);
  return result;
}

template float interpolateSlow(float, int, const DeviceInterpolationData<float>&);
template double interpolateSlow(double, int, const DeviceInterpolationData<double>&);

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
