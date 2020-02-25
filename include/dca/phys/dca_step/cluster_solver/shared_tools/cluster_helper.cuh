// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Helper class for adding and subtracting cluster elements on the device.

#ifndef DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_CLUSTER_HELPER_CUH
#define DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_CLUSTER_HELPER_CUH

#include <vector>

#include <cuda.h>

namespace dca {
namespace phys {
namespace solver {
namespace details {
// dca::phys::solver::details::

class ClusterHelper {
public:
  // Initialize real reciprocal cluster if mometnum == true.
  static void set(int nc, const int* add, int lda, const int* sub, int lds, bool momentum);

  // Returns the index of id_1 + id_2.
  __device__ inline int add(int id_1, int id_2) const;
  // Returns the index of id_2 - id_1.
  __device__ inline int subtract(int id_1, int id_2) const;
  // returns the index of -id
  __device__ inline int minus(int id) const;

private:
  int nc_;
  const int* add_matrix_;
  const int* sub_matrix_;
};

// Global instance for real space and momentum clusters.
extern __device__ __constant__ ClusterHelper cluster_real_helper;
extern __device__ __constant__ ClusterHelper cluster_momentum_helper;

inline __device__ int ClusterHelper::add(const int id_1, const int id_2) const {
  return add_matrix_[id_1 + nc_ * id_2];
}

inline __device__ int ClusterHelper::subtract(int id_1, int id_2) const {
  return sub_matrix_[id_1 + nc_ * id_2];
}

inline __device__ int ClusterHelper::minus(const int id) const {
  const int id_0 = sub_matrix_[0];
  return sub_matrix_[id + nc_ * id_0];
}

}  // namespace details
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_CLUSTER_HELPER_CUH
