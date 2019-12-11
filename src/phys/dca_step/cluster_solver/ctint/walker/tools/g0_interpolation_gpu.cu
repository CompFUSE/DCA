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

__device__ double DeviceInterpolationData::operator()(double tau, int lindex) const {
  assert(tau >= -beta_ && tau <= beta_);

  if (tau == 0)  // returns G0(tau = 0+)
    return g0_minus_[lindex];

  short int factor = 1;
  if (tau < 0) {
    tau += beta_;
    factor = -1;
  }

  // Scale tau in [0, n_time_slices). Assume even spacing in time.
  const double scaled_tau = tau * n_div_beta_;
  const int tau_index(scaled_tau);
  const double delta_tau = scaled_tau - tau_index;

  // Get the pointer to the first akima coeff.
  const double* coeff_ptr = &values_[tau_index * coeff_size_ + lindex * stride_];

  // Return akima interpolation.
  return factor *
         (coeff_ptr[0] +
          delta_tau * (coeff_ptr[1] + delta_tau * (coeff_ptr[2] + delta_tau * coeff_ptr[3])));
}

__global__ void g0InterpolationTestKernel(double tau, const int lindex,
                                          DeviceInterpolationData g0, double* result) {
  *result = g0(tau, lindex);
}

double G0Interpolation<linalg::GPU>::operator()(double tau, int lindex) const {
  double* d_result;
  double result;
  cudaMalloc((void**)&d_result, sizeof(double));

  g0InterpolationTestKernel<<<1, 1>>>(tau, lindex, *this, d_result);

  assert(cudaSuccess == cudaPeekAtLastError());
  cudaMemcpy(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(d_result);
  return result;
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
