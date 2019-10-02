// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This structure describes an element of a Sector of the CT_INT's configuration.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_SECTOR_ENTRY_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_SECTOR_ENTRY_HPP

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

struct SectorEntry {
  int get_left_band() const {
    return b_left_;
  }

  int get_right_band() const {
    return b_right_;
  }

  int get_left_site() const {
    return r_left_;
  }

  int get_right_site() const {
    return r_right_;
  }

  double get_tau() const {
    return tau_;
  }

  bool operator==(const SectorEntry& rhs) const {
    return b_left_ == rhs.b_left_ && r_left_ == rhs.r_left_ && b_right_ == rhs.b_right_ &&
           r_left_ == rhs.r_left_ && tau_ == rhs.tau_ && aux_field_type_ == rhs.aux_field_type_;
  }

  ushort b_left_;
  ushort r_left_;
  ushort b_right_;
  ushort r_right_;
  double tau_;
  short aux_field_type_;
};

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_SECTOR_ENTRY_HPP
