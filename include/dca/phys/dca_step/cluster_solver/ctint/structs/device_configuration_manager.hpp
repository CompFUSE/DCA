// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This class provides an interface to upload  to the GPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_DEVICE_CONFIGURATION_MANAGER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_DEVICE_CONFIGURATION_MANAGER_HPP

#ifndef DCA_HAVE_CUDA
#error "This file requires CUDA support."
#endif

#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/solver_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/sector_entry.hpp"
#include "dca/util/cuda_definitions.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

class DeviceConfigurationManager {
public:
  DeviceConfigurationManager() = default;

  inline void upload(const SolverConfiguration& config, int thread_id, int spin);
  inline void upload(const SolverConfiguration& config, int thread_id);

  details::DeviceConfiguration getDeviceData(int s) const {
    return device_pointers_[s];
  }

private:
  std::array<linalg::Vector<details::SectorEntry, linalg::GPU>, 2> device_entries_;
  std::array<details::DeviceConfiguration, 2> device_pointers_;
};

void DeviceConfigurationManager::upload(const SolverConfiguration& config, int thread_id, int spin) {
  assert(spin >= 0 and spin < 2);
  const auto& entries = config.getSector(spin).entries_;
  device_entries_[spin].setAsync(entries, linalg::util::getStream(thread_id, spin));
  device_pointers_[spin].data = device_entries_[spin].ptr();
}

void DeviceConfigurationManager::upload(const SolverConfiguration& config, int thread_id) {
  upload(config, thread_id, 0);
  upload(config, thread_id, 1);
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_DEVICE_CONFIGURATION_MANAGER_HPP
