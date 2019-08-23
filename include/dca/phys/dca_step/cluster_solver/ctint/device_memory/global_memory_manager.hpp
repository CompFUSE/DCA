// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Manges memory contained in ctint/device_memory/global_memory.cuh.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DEVICE_MEMORY_GLOBAL_MEMORY_MANAGER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DEVICE_MEMORY_GLOBAL_MEMORY_MANAGER_HPP

#include <vector>

#include "dca/linalg/matrix_view.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

class GlobalMemoryManager {
public:
  // Initialize const memory on the device.
  // If already initialized and not forced, of if DCA does not have CUDA, do nothing.
  static void initialize(const linalg::MatrixView<int, linalg::CPU>& site_diff_matrix,
                         const std::vector<std::size_t>& sbdm_step, const int parameters_step,
                         bool override = false);

  static void initializeInterpolation(const int parameters_step, bool override = false);
  static void initializeCluster(const linalg::MatrixView<int, linalg::CPU>& site_diff_matrix,
                                const std::vector<std::size_t>& sbdm_step, bool override = false);

  static bool isInitialized() {
    return cluster_initialized_ and interpolation_initialized_;
  }

private:
  static bool cluster_initialized_;
  static bool interpolation_initialized_;
};

#ifndef DCA_HAVE_CUDA
inline void GlobalMemoryManager::initialize(
    const linalg::MatrixView<int, linalg::CPU>& /*site_diff_matrix*/,
    const std::vector<int>& /*sbdm_step*/, const int /*parameters_step*/, bool /*force*/) {}
inline void GlobalMemoryManager::initializeCluster(
    const linalg::MatrixView<int, linalg::CPU>& /*site_diff_matrix*/,
    const std::vector<int>& /*sbdm_step*/, bool /*override*/) {}
inline void GlobalMemoryManager::initializeInterpolation(const int /*parameters_step*/,
                                                         bool /*override*/) {}
#endif  // DCA_HAVE_CUDA

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DEVICE_MEMORY_GLOBAL_MEMORY_MANAGER_HPP
