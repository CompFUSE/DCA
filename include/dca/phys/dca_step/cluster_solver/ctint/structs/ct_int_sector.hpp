// Copyright_ (C) 2009-2016 ETH Zurich
// Copyright_ (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All right_s reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author:Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This structure organize the data of MatrixConfiguration for each sector in {up, down}

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_SECTOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_SECTOR_HPP

#include <array>
#include <cassert>
#include <vector>

#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/sector_entry.hpp"
#include "dca/linalg/util/allocators.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

// TODO: replace by vector
class Sector {
public:
  Sector() : entries_(){};

  friend class MatrixConfiguration;

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

protected:
  linalg::util::HostVector<details::SectorEntry> entries_;
};

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_SECTOR_HPP
