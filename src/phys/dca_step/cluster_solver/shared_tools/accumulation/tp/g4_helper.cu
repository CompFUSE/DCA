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

void G4Helper::set(int nb, int nk, int nw_pos, const std::vector<int>& delta_k,
                   const std::vector<int>& delta_w, const int* add_k, int lda, const int* sub_k,
                   int lds, int k0) {
  static std::once_flag flag;

  std::call_once(flag, [=]() {
    G4Helper host_helper;
    host_helper.lda_ = lda;
    host_helper.lds_ = lds;
    host_helper.nw_pos_ = nw_pos;
    host_helper.k0_ = k0;

    host_helper.ext_size_ = 0;
    for (const int idx : delta_w)
      host_helper.ext_size_ = std::max(host_helper.ext_size_, std::abs(idx));

    const int nb4 = nb * nb * nb * nb;
    const int nk3 = nk * nk * delta_k.size();
    const int nw = 2 * nw_pos;
    const std::array<int, 10> steps{1,
                                    nb,
                                    nb * nb,
                                    nb * nb * nb,
                                    nb4,
                                    nb4 * nk,
                                    nb4 * nk * nk,
                                    nb4 * nk3,
                                    nb4 * nk3 * nw,
                                    nb4 * nk3 * nw * nw};
    std::copy_n(steps.data(), steps.size(), host_helper.sbdm_steps_);

    cudaMalloc(&host_helper.add_matrix_, sizeof(int) * lda * nk);
    cudaMemcpy(host_helper.add_matrix_, add_k, sizeof(int) * lda * nk, cudaMemcpyHostToDevice);

    cudaMalloc(&host_helper.sub_matrix_, sizeof(int) * lds * nk);
    cudaMemcpy(host_helper.sub_matrix_, sub_k, sizeof(int) * lds * nk, cudaMemcpyHostToDevice);

    cudaMalloc(&host_helper.w_ex_indices_, sizeof(int) * delta_w.size());
    cudaMemcpy(host_helper.w_ex_indices_, delta_w.data(), sizeof(int) * delta_w.size(),
               cudaMemcpyHostToDevice);

    cudaMalloc(&host_helper.k_ex_indices_, sizeof(int) * delta_k.size());
    cudaMemcpy(host_helper.k_ex_indices_, delta_k.data(), sizeof(int) * delta_k.size(),
               cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(g4_helper, &host_helper, sizeof(G4Helper));
  });
}

}  // namespace details
}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca
