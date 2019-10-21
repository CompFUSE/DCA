// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This structure organize the data of MatrixConfiguration for each sector in {up, down}

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_SECTOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_SECTOR_HPP

#include <array>
#include <cassert>
#include <vector>

#include "dca/linalg/util/allocators/vectors_typedefs.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/sector_entry.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

// TODO: replace by vector
class Sector {
public:
  Sector() : entries_(){};

  friend class SolverConfiguration;
  friend class MatrixConfiguration;
  friend class DeviceConfigurationManager;

  ushort size() const {
    return entries_.size();
  }

  const details::SectorEntry& operator[](const std::size_t idx) const {
    assert(idx < size());
    return entries_[idx];
  }

  inline ushort getLeftB(const int i) const {
    return entries_[i].b_left_;
  }

  inline ushort getRightB(const int i) const {
    return entries_[i].b_right_;
  }

  inline ushort getLeftR(const int i) const {
    return entries_[i].r_left_;
  }

  inline ushort getRightR(const int i) const {
    return entries_[i].r_right_;
  }

  inline double getTau(const int i) const {
    return entries_[i].tau_;
  }

  inline short getAuxFieldType(int i) const {
    return entries_[i].aux_field_type_;
  }

  auto getTag(int i) const {
    return tags_[i];
  }

  bool operator==(const Sector& rhs) const {
    return entries_ == rhs.entries_;
  }

private:
  linalg::util::HostVector<details::SectorEntry> entries_;
  std::vector<std::uint64_t> tags_;
};

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_SECTOR_HPP
