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

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_SOLVER_HELPER_CUH
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_SOLVER_HELPER_CUH

#ifdef DCA_HAVE_CUDA
#include <cuda.h>

#include "dca/phys/dca_step/cluster_solver/shared_tools/cluster_helper.cuh"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#endif

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

#ifdef DCA_HAVE_CUDA
class SolverHelper {
public:
  static void set(const int* sum_r, int lda, const int* sub_r, int lds, int nb, int nc, int r0);

  template <class RDmn, class BDmn>
  static void set();

  // Return the index of a single particle function of b1, b2, r1 - r2.
  __device__ std::size_t index(int b1, int b2, int r1, int r2) const;

private:
  std::size_t subdm_step_[2];
};

// Global instance.
extern __device__ __constant__ SolverHelper solver_helper;

__device__ inline std::size_t SolverHelper::index(int b1, int b2, int r1, int r2) const {
  const int delta_r = solver::cluster_real_helper.subtract(r2, r1);
  return b1 + b2 * subdm_step_[0] + delta_r * subdm_step_[1];
}

template <class RDmn, class BDmn>
void SolverHelper::set() {
  using Cluster = typename RDmn::parameter_type;
  static_assert(Cluster::REPRESENTATION == domains::REAL_SPACE, "Domain mismatch.");
  const auto& add_matrix = Cluster::get_add_matrix();
  const auto& sub_matrix = Cluster::get_subtract_matrix();

  set(add_matrix.ptr(), add_matrix.leadingDimension(), sub_matrix.ptr(),
      sub_matrix.leadingDimension(), BDmn::dmn_size(), RDmn::dmn_size(), Cluster::origin_index());
}

#else  // !DCA_HAVE_CUDA
// No-op version.
class SolverHelper {
public:
  template <class RDmn, class BDmn>
  static void set() {}

  constexpr static bool initialized(){
      return false;
  }
};
#endif  // DCA_HAVE_CUDA

}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_SOLVER_HELPER_CUH
