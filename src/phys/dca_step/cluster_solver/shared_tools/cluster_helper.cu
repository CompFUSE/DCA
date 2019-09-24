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

#include "dca/linalg/util/allocators/vectors_typedefs.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace details {
// dca::phys::solver::details::

__device__ __constant__ ClusterHelper cluster_real_helper;
__device__ __constant__ ClusterHelper cluster_momentum_helper;

void ClusterHelper::set(int nc, const int* add, int lda, const int* sub, int lds, bool momentum) {
  static std::array<std::once_flag, 2> flags;

  std::call_once(flags[momentum], [=]() {
    ClusterHelper host_helper;
    host_helper.nc_ = nc;

    auto compact_transfer = [=](const int* matrix, int ldm, int** dest) {
      linalg::util::HostVector<int> compact(nc * nc);
      for (int j = 0; j < nc; ++j)
        for (int i = 0; i < nc; ++i)
          compact[i + nc * j] = matrix[i + ldm * j];

      cudaMalloc(dest, sizeof(int) * lds * nc);
      cudaMemcpy(*dest, compact.data(), sizeof(int) * nc * nc, cudaMemcpyHostToDevice);
    };

    compact_transfer(add, lda, const_cast<int**>(&host_helper.add_matrix_));
    compact_transfer(sub, lds, const_cast<int**>(&host_helper.sub_matrix_));

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
