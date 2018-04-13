// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the device methods of DMatrixBuilder, G0Interpolation<GPU> and
// GlobalMemoryManager.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"

#include <climits>
#include <cuda.h>
#include <cuda_runtime.h>

#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/device_memory/global_memory_manager.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#include "dca/util/integer_division.hpp"

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
using linalg::GPU;
using linalg::CPU;
using Configuration = DeviceConfiguration;
using Interpolation = details::DeviceInterpolationData;
using MatrixView = linalg::MatrixView<double, linalg::GPU>;

__global__ void computeDKernel(MatrixView Q, MatrixView R, MatrixView S, const int delta,
                               const double alpha1, const double alpha2, const double alpha3,
                               const Configuration config, const Interpolation g0) {
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;

  const int n = Q.nrRows();
  // Check boundaries.
  if ((id_i >= n + delta) or (id_j >= 2 * delta))
    return;

  // Determine matrix indices and write location.
  int i, j;
  double* write_location;
  if (id_j >= delta) {
    if (id_i >= n)
      return;
    else {  // Write to R
      i = id_j - delta + n;
      j = id_i;
      write_location = R.ptr(id_j - delta, id_i);
    }
  }
  else {
    if (id_i < n) {  // Write to Q.
      i = id_i;
      j = id_j + n;
      write_location = Q.ptr(id_i, id_j);
    }
    else {  // Write to S
      i = id_i;
      j = id_j + n;
      write_location = S.ptr(id_i - n, id_j);
    }
  }

  auto alphaField = [&](const int type) {
    switch (std::abs(type)) {
      case 1:
        return type < 0 ? (0.5 + alpha1) : (0.5 - alpha1);
      case 2:
        return type < 0 ? (0.5 + alpha2) : (0.5 - alpha2);
      case 3:
        return type < 0 ? alpha3 : -alpha3;
    }
    return 0.;
  };

  const ushort b1 = config.getLeftB(i);
  const ushort b2 = config.getRightB(j);

  using global::sbdm_step;
  const ushort delta_r =
      global::cluster_site_diff[config.getLeftR(i) + config.getRightR(j) * global::cluster_size];
  const ushort p_index = b1 + b2 * sbdm_step[1] + delta_r * sbdm_step[2];
  const double delta_tau = config.getTau(i) - config.getTau(j);
  *write_location = g0(delta_tau, p_index);
  if (i == j)
    *write_location -= alphaField(config.getAuxFieldType(i));
}

void computeD(linalg::MatrixView<double, linalg::GPU> Q, linalg::MatrixView<double, linalg::GPU> R,
              linalg::MatrixView<double, linalg::GPU> S, int n, int delta, double alpha_1,
              double alpha_2, double alpha_3, Configuration config, DeviceInterpolationData g0_data,
              cudaStream_t stream) {
  const int number_of_threads = std::min(128, n);  // maximum number of thread rows.
  dim3 threads, blocks;
  if (n) {
    threads = dim3(number_of_threads, delta);
    blocks = dim3(util::ceilDiv(n + delta, number_of_threads), 2);
  }
  else {
    threads = dim3(delta, delta);
    blocks = dim3(1, 1);
  }

  assert(GlobalMemoryManager::isInitialized());
  computeDKernel<<<blocks, threads, 0, stream>>>(Q, R, S, delta, alpha_1, alpha_2, alpha_3, config,
                                                 g0_data);
  assert(cudaDeviceSynchronize() == cudaSuccess);
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
