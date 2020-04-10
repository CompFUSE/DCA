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

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation_gpu.hpp"

#include <cuda_runtime.h>

#include "dca/linalg/util/error_cuda.hpp"
#include "dca/util/cuda_blocks.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <typename Real>
__global__ void g0InterpolationTestKernel(Real tau, const int lindex,
                                          DeviceInterpolationData<Real> g0, Real* result) {
  *result = g0(tau, lindex);
}

template <typename Real>
Real G0Interpolation<linalg::GPU, Real>::operator()(Real tau, int lindex) const {
  Real* d_result;
  Real result;
  cudaMalloc((void**)&d_result, sizeof(Real));

  g0InterpolationTestKernel<<<1, 1>>>(tau, lindex, *this, d_result);

  assert(cudaSuccess == cudaPeekAtLastError());
  cudaMemcpy(&result, d_result, sizeof(Real), cudaMemcpyDeviceToHost);
  cudaFree(d_result);
  return result;
}

template class G0Interpolation<linalg::GPU, float>;
template class G0Interpolation<linalg::GPU, double>;

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
