// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements ClusterHelper::set.

#include "dca/phys/dca_step/cluster_solver/shared_tools/cluster_helper.cuh"

#include <mutex>

namespace dca {
namespace phys {
namespace solver {
namespace details {
// dca::phys::solver::details::

__device__ __constant__ ClusterHelper cluster_real_helper;
__device__ __constant__ ClusterHelper cluster_momentum_helper;

void ClusterHelper::set(int nc, const int* add, int lda, const int* sub, int lds, int id_0,
                        bool momentum) {
  static std::array<std::once_flag, 2> flags;

  std::call_once(flags[momentum], [=]() {
    ClusterHelper host_helper;
    host_helper.lda_ = lda;
    host_helper.lds_ = lds;
    host_helper.id_0_ = id_0;

    cudaMalloc(&host_helper.add_matrix_, sizeof(int) * lda * nc);
    cudaMemcpy(const_cast<int*>(host_helper.add_matrix_), const_cast<int*>(add),
               sizeof(int) * lda * nc, cudaMemcpyHostToDevice);

    cudaMalloc(&host_helper.sub_matrix_, sizeof(int) * lds * nc);
    cudaMemcpy(const_cast<int*>(host_helper.sub_matrix_), const_cast<int*>(sub),
               sizeof(int) * lds * nc, cudaMemcpyHostToDevice);

    if (momentum) {
      cudaMemcpyToSymbol(cluster_momentum_helper, &host_helper, sizeof(ClusterHelper));
    }
    else {
      cudaMemcpyToSymbol(cluster_real_helper, &host_helper, sizeof(ClusterHelper));
    }
  });
}

}  // namespace details
}  // namespace solver
}  // namespace phys
}  // namespace dca
