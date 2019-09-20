// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Helper class used by the interpolation and DMatrixBuilder classes.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DEVICE_HELPER_CTINT_HELPER_CUH
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DEVICE_HELPER_CTINT_HELPER_CUH

#include <vector>

#include <cuda.h>

#include "dca/phys/dca_step/cluster_solver/shared_tools/cluster_helper.cuh"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

class CtintHelper {
public:
  static void set(const int* sum_r, int lda, const int* sub_r, int lds, int nb, int nc, int r0);

  // Return the index of a single particle function of b1, b2, r1 - r2.
  __device__ std::size_t index(int b1, int b2, int r1, int r2) const;

private:
  std::size_t subdm_step_[2];
};

// Global instance.
extern __device__ __constant__ CtintHelper ctint_helper;

__device__ inline std::size_t CtintHelper::index(int b1, int b2, int r1, int r2) const {
  const int delta_r = solver::details::cluster_real_helper.subtract(r2, r1);
  return b1 + b2 * subdm_step_[0] + delta_r * subdm_step_[1];
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DEVICE_HELPER_CTINT_HELPER_CUH
