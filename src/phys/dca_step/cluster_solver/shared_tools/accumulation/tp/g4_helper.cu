// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
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

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
// dca::phys::solver::accumulator::details::

__device__ __constant__ G4Helper g4_helper;

void G4Helper::set(unsigned int nb, unsigned int nk, unsigned int nw_pos,
                   const std::vector<int>& delta_k, const std::vector<int>& delta_w,
                   const int* add_k, unsigned int lda, const int* sub_k, unsigned int lds) {
  static std::once_flag flag;

  std::call_once(flag, [=]() {
    // Initialize the reciprocal cluster if not done already.
    solver::details::ClusterHelper::set(nk, add_k, lda, sub_k, lds, true);

    G4Helper host_helper;
    host_helper.nb_ = nb;
    host_helper.nc_ = nk;
    host_helper.nw_pos_ = nw_pos;
    host_helper.nw_ = 2 * nw_pos;
    host_helper.n_k_ex_ = delta_k.size();
    host_helper.n_w_ex_ = delta_w.size();

    host_helper.ext_size_ = 0;
    for (const auto idx : delta_w)
      host_helper.ext_size_ = std::max(host_helper.ext_size_, static_cast<unsigned>(std::abs(idx)));

    // compute strides
    const std::array<unsigned, 10> sizes{nb,
                                         nb,
                                         nb,
                                         nb,
                                         nk,
                                         host_helper.nw_,
                                         nk,
                                         host_helper.nw_,
                                         static_cast<unsigned>(delta_k.size()),
                                         static_cast<unsigned>(delta_w.size())};

    host_helper.sbdm_steps_[0] = 1;
    for (std::size_t i = 1; i < sizes.size(); ++i)
      host_helper.sbdm_steps_[i] = host_helper.sbdm_steps_[i - 1] * sizes[i - 1];

    cudaMalloc(&host_helper.w_ex_indices_, sizeof(int) * delta_w.size());
    cudaMemcpy(const_cast<int*>(host_helper.w_ex_indices_), delta_w.data(),
               sizeof(int) * delta_w.size(), cudaMemcpyHostToDevice);

    cudaMalloc(&host_helper.k_ex_indices_, sizeof(int) * delta_k.size());
    cudaMemcpy(const_cast<int*>(host_helper.k_ex_indices_), delta_k.data(),
               sizeof(int) * delta_k.size(), cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(g4_helper, &host_helper, sizeof(G4Helper));
  });
}

}  // namespace details
}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca
