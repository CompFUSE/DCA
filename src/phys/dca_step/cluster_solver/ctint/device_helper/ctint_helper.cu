// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the static methods of CtintHelper.

#include "dca/phys/dca_step/cluster_solver/ctint/device_helper/ctint_helper.cuh"

#include <stdexcept>
#include <mutex>

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

// Global helper instance.
__device__ __constant__ CtintHelper ctint_helper;

void CtintHelper::set(const int* add_r, int lda, const int* sub_r, int lds, const int nb,
                      const int nc, const int r0) {
  static std::once_flag flag;
  std::call_once(flag, [&] {
    // Initialize real space cluster.
    solver::details::ClusterHelper::set(nc, add_r, lda, sub_r, lds, 0);

    CtintHelper host_helper;
    host_helper.subdm_step_[0] = nb;
    host_helper.subdm_step_[1] = nb * nb;

    cudaMemcpyToSymbol(ctint_helper, &host_helper, sizeof(CtintHelper));
  });
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
