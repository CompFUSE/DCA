// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements G4Helper::set.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/g4_helper.cuh"

#include <algorithm>
#include <array>
#include <mutex>
#include <stdexcept>
#include <iostream>

#include "dca/platform/dca_gpu.h"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
// dca::phys::solver::accumulator::details::

__device__ __constant__ G4Helper g4_helper;

void G4Helper::set(int nb, int nk, int nw, const std::vector<int>& delta_k,
                   const std::vector<int>& delta_w, const int extension_offset, const int* add_k, int lda, const int* sub_k,
                   int lds) {
  // Initialize the reciprocal cluster if not done already.
  solver::details::ClusterHelper::set(nk, add_k, lda, sub_k, lds, true);

  G4Helper host_helper;
  host_helper.nb_ = nb;
  host_helper.nc_ = nk;
  host_helper.nw_ = nw;
  host_helper.n_k_ex_ = delta_k.size();
  host_helper.n_w_ex_ = delta_w.size();

  host_helper.ext_size_ = extension_offset;
  // for (const auto idx : delta_w)
  //   host_helper.ext_size_ = std::max(host_helper.ext_size_, static_cast<int>(std::abs(idx)));

  // compute strides
  const std::array<int, 10> sizes{nb,
                                  nb,
                                  nb,
                                  nb,
                                  nk,
                                  host_helper.nw_,
                                  nk,
                                  host_helper.nw_,
                                  static_cast<int>(delta_k.size()),
                                  static_cast<int>(delta_w.size())};

  host_helper.sbdm_steps_[0] = 1;
  for (std::size_t i = 1; i < sizes.size(); ++i)
    host_helper.sbdm_steps_[i] = host_helper.sbdm_steps_[i - 1] * sizes[i - 1];

  cudaMalloc(&host_helper.w_ex_indices_, sizeof(int) * delta_w.size());
  cudaMemcpy(const_cast<int*>(host_helper.w_ex_indices_), delta_w.data(),
             sizeof(int) * delta_w.size(), cudaMemcpyHostToDevice);

  cudaMalloc(&host_helper.k_ex_indices_, sizeof(int) * delta_k.size());
  cudaMemcpy(const_cast<int*>(host_helper.k_ex_indices_), delta_k.data(),
             sizeof(int) * delta_k.size(), cudaMemcpyHostToDevice);

#ifndef NDEBUG
  cudaMalloc(&host_helper.bad_indicies_, sizeof(int) * 1024);
#endif
  cudaMemcpyToSymbol(g4_helper, &host_helper, sizeof(G4Helper));
}

__device__ bool G4Helper::extendGIndices(int& k1, int& k2, int& w1, int& w2) const {
  const int extension_offset = ext_size_;
  w1 += extension_offset;
  w2 += extension_offset;
  const int n_w_ext = nw_ + ext_size_ + 1;
#ifndef NDEBUG
  if (w1 >= n_w_ext || w2 >= n_w_ext) {
    if (w1 >= n_w_ext)
      bad_indicies_[threadIdx.x + threadIdx.y * 32] = w1;
    else
      bad_indicies_[threadIdx.x + threadIdx.y * 32] = w2;
  }
#endif
  return false;
}

__device__ bool G4Helper::extendGIndicesMultiBand(int& k1 [[maybe_unused]],
                                                  int& k2 [[maybe_unused]], int& w1, int& w2) const {
  const int extension_offset = ext_size_;
  w1 += extension_offset;
  w2 += extension_offset;
  const int n_w_ext = nw_ + ext_size_ + 1;
#ifndef NDEBUG
  if (w1 >= n_w_ext || w2 >= n_w_ext) {
    if (w1 >= n_w_ext)
      bad_indicies_[threadIdx.x + threadIdx.y * 32] = w1;
    else
      bad_indicies_[threadIdx.x + threadIdx.y * 32] = w2;
  }
#endif
  return false;
}

}  // namespace details
}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca
