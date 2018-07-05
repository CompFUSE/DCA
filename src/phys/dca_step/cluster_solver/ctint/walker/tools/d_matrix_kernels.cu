// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the device methods of DMatrixBuilder, G0Interpolation<GPU> and
// GlobalMemoryManager.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"

#include <climits>
#include <cuda.h>
#include <cuda_runtime.h>

#include "dca/linalg/util/error_cuda.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/device_memory/global_memory_manager.hpp"
#include "dca/util/cuda_blocks.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace global {
// dca::phys::solver::ctint::global
// ********** Global Memory *********
constexpr int static MAX_CLUSTER_SIZE = 64;
__constant__ uint g0_paremeters_step;
__constant__ ushort cluster_size;
__constant__ ushort cluster_site_diff[MAX_CLUSTER_SIZE * MAX_CLUSTER_SIZE];
__constant__ ushort sbdm_step[3];
}  // global
namespace details {
// dca::phys::solver::ctint::details::

// ********** D Matrix Builder *********
__global__ void buildG0MatrixKernel(MatrixView G0, const int n_init, const bool right_section,
                                    Configuration config, Interpolation g0_interp) {
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

  const int delta_r =
      global::cluster_site_diff[config.getLeftR(i) + config.getRightR(j) * global::cluster_size];
  const int label = b_i + b_j * global::sbdm_step[1] + delta_r * global::sbdm_step[2];

  G0(id_i, id_j) = g0_interp(tau_i - tau_j, label);
}

void buildG0Matrix(MatrixView G0, const int n_init, const bool right_section, Configuration config,
                   Interpolation g0_interp, cudaStream_t stream) {
  assert(GlobalMemoryManager::isInitialized());
  const auto blocks = dca::util::getBlockSize(G0.nrRows(), G0.nrCols());

  buildG0MatrixKernel<<<blocks[0], blocks[1], 0, stream>>>(G0, n_init, right_section, config,
                                                           g0_interp);
  checkErrorsCudaDebug();
}

// ************  G0 Interpolation.  **************
__device__ double DeviceInterpolationData::operator()(double tau, int lindex) const {
  assert(tau >= -beta_ and tau <= beta_);

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
  const double* coeff_ptr = &values_[tau_index * coeff_size_ + lindex * global::g0_paremeters_step];
  // Return akima interpolation.
  return factor *
         (coeff_ptr[0] +
          delta_tau * (coeff_ptr[1] + delta_tau * (coeff_ptr[2] + delta_tau * coeff_ptr[3])));
}

__global__ void g0InterpolationTestKernel(double tau, const int lindex, Interpolation g0,
                                          double* result) {
  *result = g0(tau, lindex);
}

double deviceInterpolationTest(Interpolation g0, double tau, int lindex) {
  double* d_result;
  double result;
  cudaMalloc((void**)&d_result, sizeof(double));

  g0InterpolationTestKernel<<<1, 1>>>(tau, lindex, g0, d_result);

  assert(cudaSuccess == cudaPeekAtLastError());
  cudaMemcpy(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(d_result);
  return result;
}
}
// dca::phys::solver::ctint::

// ************  Global Memory Manager  **************
bool GlobalMemoryManager::cluster_initialized_ = false;
bool GlobalMemoryManager::interpolation_initialized_ = false;

void GlobalMemoryManager::initialize(const linalg::MatrixView<int, linalg::CPU>& site_diff_matrix,
                                     const std::vector<int>& sbdm_step, const int parameters_step,
                                     bool override) {
  initializeCluster(site_diff_matrix, sbdm_step, override);
  initializeInterpolation(parameters_step, override);
}

void GlobalMemoryManager::initializeCluster(const linalg::MatrixView<int, linalg::CPU>& site_diff_matrix,
                                            const std::vector<int>& sbdm_step, bool override) {
  if (cluster_initialized_ and not override)
    return;
  constexpr int max_size = global::MAX_CLUSTER_SIZE;

  ushort compact_diff[max_size * max_size];
  const ushort cluster_size = site_diff_matrix.nrRows();

  if (cluster_size > max_size) {
    throw(std::bad_alloc());  // Not enough constant memory allocated on device.
  }

  // Copy matrix values discarding the padding.
  for (int j = 0; j < site_diff_matrix.nrCols(); j++)
    for (int i = 0; i < site_diff_matrix.nrRows(); i++)
      compact_diff[i + j * cluster_size] = site_diff_matrix(i, j);

  std::vector<ushort> sbdm_step_short(sbdm_step.size());
  for (size_t i = 0; i < sbdm_step_short.size(); i++)
    sbdm_step_short[i] = (ushort)sbdm_step[i];

  // copy values to GPU
  cudaMemcpyToSymbol(global::sbdm_step, sbdm_step_short.data(),
                     sbdm_step_short.size() * sizeof(ushort));
  cudaMemcpyToSymbol(global::cluster_site_diff, compact_diff,
                     cluster_size * cluster_size * sizeof(ushort));
  cudaMemcpyToSymbol(global::cluster_size, &cluster_size, sizeof(ushort));

  cluster_initialized_ = true;
}

void GlobalMemoryManager::initializeInterpolation(const int parameters_step, bool override) {
  if (interpolation_initialized_ and not override)
    return;
  cudaMemcpyToSymbol(global::g0_paremeters_step, &parameters_step, sizeof(ushort));
  interpolation_initialized_ = true;
}

}  // ctint
}  // solver
}  // phys
}  // dca
