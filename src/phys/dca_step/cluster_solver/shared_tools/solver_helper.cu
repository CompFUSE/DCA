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

#include "dca/phys/dca_step/cluster_solver/shared_tools/solver_helper.cuh"

#include <stdexcept>
#include <mutex>

#include "dca/phys/dca_step/cluster_solver/shared_tools/cluster_helper.cuh"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

// Global helper instance.
__device__ __constant__ SolverHelper solver_helper;

bool SolverHelper::initialized_ = false;

void SolverHelper::set(const int* add_r, int lda, const int* sub_r, int lds, const int nb,
                       const int nc, const int r0) {
  static std::once_flag flag;
  std::call_once(flag, [&] {
    // Initialize real space cluster.
    solver::ClusterHelper::set(nc, add_r, lda, sub_r, lds, 0);

    SolverHelper host_helper;
    host_helper.subdm_step_[0] = nb;
    host_helper.subdm_step_[1] = nb * nb;

    cudaMemcpyToSymbol(solver_helper, &host_helper, sizeof(SolverHelper));
    initialized_ = true;
  });
}

}  // namespace solver
}  // namespace phys
}  // namespace dca
