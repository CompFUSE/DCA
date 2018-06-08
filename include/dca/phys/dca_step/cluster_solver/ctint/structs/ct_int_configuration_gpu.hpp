// Copyright_ (C) 2009-2016 ETH Zurich
// Copyright_ (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All right_s reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author:Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This class organizes the vertex configuration for ct-int.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_CONFIGURATION_GPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_CONFIGURATION_GPU_HPP
#ifdef DCA_HAVE_CUDA

#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/sector_entry.hpp"
#include "dca/util/cuda_definitions.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <>
struct SolverConfiguration<linalg::GPU> : public SolverConfiguration<linalg::CPU> {
  using BaseClass = SolverConfiguration<linalg::CPU>;

  SolverConfiguration(double beta, int n_bands, const InteractionVertices& H_int,
                      double double_move_prob = 0)
      : SolverConfiguration<linalg::CPU>(beta, n_bands, H_int, double_move_prob) {}

  inline SolverConfiguration() = default;

  inline SolverConfiguration(const SolverConfiguration<linalg::GPU>& other) = default;

  //  inline SolverConfiguration<linalg::GPU>& operator=(const SolverConfiguration<linalg::GPU>&
  //  rhs);

  inline void upload(int spin, int thread_id);
  inline void upload(int thread_id);

  using BaseClass::pop;
  using BaseClass::push_back;
  using BaseClass::back;

  details::DeviceConfiguration getDeviceData(int s) const {
    return device_pointers_[s];
  }

private:
  using BaseClass::vertices_;
  std::array<linalg::Vector<details::SectorEntry, linalg::GPU>, 2> device_entries_;
  std::array<details::DeviceConfiguration, 2> device_pointers_;
};

void SolverConfiguration<linalg::GPU>::upload(int spin, int thread_id) {
  assert(spin >= 0 and spin < 2);
  const auto& entries = BaseClass::BaseClass::getEntries(spin);
  device_entries_[spin].setAsync(entries, linalg::util::getStream(thread_id, spin));
  device_pointers_[spin].data = device_entries_[spin].ptr();
}

void SolverConfiguration<linalg::GPU>::upload(int thread_id){
    upload(0, thread_id);
    upload(1, thread_id);
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_CONFIGURATION_GPU_HPP
