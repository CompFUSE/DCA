// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This class process the information contained in the SolverConfiguration class into a
// representation of the M matrix.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_MATRIX_CONFIGURATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_MATRIX_CONFIGURATION_HPP

#include <array>
#include <vector>

#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_sector.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

// Expansion term.
struct Vertex {
  bool aux_spin;
  ushort interaction_id;
  std::uint64_t tag;
  double tau;

  bool operator==(const Vertex& b) const {
    return aux_spin == b.aux_spin && interaction_id == b.interaction_id && tau == b.tau;
  }
};

class MatrixConfiguration {
public:
  MatrixConfiguration() = default;
  MatrixConfiguration(const MatrixConfiguration& rhs) = default;
  inline MatrixConfiguration& operator=(const MatrixConfiguration& rhs);
  inline MatrixConfiguration& operator=(MatrixConfiguration&& rhs);

  inline void swapSectorLabels(int a, int b, int s);

  const Sector& getSector(int s) const {
    return sectors_[s];
  }

  std::size_t size() const {
    return (size(0) + size(1)) / 2;
  }

  std::size_t size(int s) const {
    return sectors_[s].size();
  }

  const std::array<Sector, 2>& get_sectors() const {
    return sectors_;
  }

  bool operator==(const MatrixConfiguration& rhs) const {
    return sectors_ == rhs.sectors_;
  }

  // TODO try to return on the stack or preallocated mem.
  inline std::vector<int> findIndices(const uint64_t tag, const int s) const;

protected:
  inline MatrixConfiguration(const InteractionVertices* H_int, int bands);

  inline std::array<int, 2> addVertex(const Vertex& v);

  inline void pop(ushort idx_up, ushort idx_down);

  const auto& getEntries(const int s) const {
    return sectors_[s].entries_;
  }
  auto& getEntries(const int s) {
    return sectors_[s].entries_;
  }

protected:
  const InteractionVertices* H_int_ = nullptr;
  const int n_bands_ = -1;
  std::array<Sector, 2> sectors_;
};

MatrixConfiguration::MatrixConfiguration(const InteractionVertices* H_int, const int bands)
    : H_int_(H_int), n_bands_(bands), sectors_{Sector(), Sector()} {}

MatrixConfiguration& MatrixConfiguration::operator=(const MatrixConfiguration& rhs) {
  sectors_ = rhs.sectors_;
  return *this;
}

MatrixConfiguration& MatrixConfiguration::operator=(MatrixConfiguration&& rhs) {
  sectors_ = std::move(rhs.sectors_);
  return *this;
}

std::array<int, 2> MatrixConfiguration::addVertex(const Vertex& v) {
  auto spin = [=](const int nu) { return nu >= n_bands_; };
  auto band = [=](const int nu) -> ushort { return nu - n_bands_ * spin(nu); };

  auto field_type = [&](const Vertex& v, const int leg) -> short {
    const short sign = v.aux_spin ? 1 : -1;
    const InteractionElement& elem = (*H_int_)[v.interaction_id];
    if (elem.partners_id.size())
      return leg == 1 ? -3 * sign : 3 * sign;  // non density-density.
    else if (elem.w > 0)
      return leg == 1 ? -1 * sign : 1 * sign;  // positive dd interaction.
    else
      return 2 * sign;  // negative dd interaction.
  };

  std::array<int, 2> sizes{0, 0};
  const auto& nu = (*H_int_)[v.interaction_id].nu;
  const auto& r = (*H_int_)[v.interaction_id].r;

  for (ushort leg = 0; leg < 2; ++leg) {
    assert(spin(nu[0 + 2 * leg]) == spin(nu[1 + 2 * leg]));
    const short s = spin(nu[0 + 2 * leg]);
    Sector& sector = sectors_[s];
    ++sizes[s];
    sector.entries_.emplace_back(details::SectorEntry{band(nu[0 + 2 * leg]), r[0 + 2 * leg],
                                                      band(nu[1 + 2 * leg]), r[1 + 2 * leg], v.tau,
                                                      field_type(v, leg)});
    sector.tags_.push_back(v.tag);
  }
  assert(sectors_[0].entries_.size() == sectors_[0].tags_.size());
  assert(sectors_[1].entries_.size() == sectors_[1].tags_.size());
  return sizes;
}

void MatrixConfiguration::pop(ushort n_up, ushort n_down) {
  assert(n_up <= sectors_[0].size() and n_down <= sectors_[1].size());
  sectors_[0].entries_.erase(sectors_[0].entries_.end() - n_up, sectors_[0].entries_.end());
  sectors_[0].tags_.erase(sectors_[0].tags_.end() - n_up, sectors_[0].tags_.end());
  sectors_[1].entries_.erase(sectors_[1].entries_.end() - n_down, sectors_[1].entries_.end());
  sectors_[1].tags_.erase(sectors_[1].tags_.end() - n_down, sectors_[1].tags_.end());
}

std::vector<int> MatrixConfiguration::findIndices(const uint64_t tag, const int s) const {
  std::vector<int> indices;

  auto start = sectors_[s].tags_.begin();
  const auto end = sectors_[s].tags_.end();
  while (true) {
    start = std::find(start, end, tag);
    if (start != end) {
      indices.push_back(start - sectors_[s].tags_.begin());
      ++start;
    }
    else
      break;
  }

  return indices;
}

void MatrixConfiguration::swapSectorLabels(const int a, const int b, const int s) {
  std::swap(sectors_[s].entries_[a], sectors_[s].entries_[b]);
  std::swap(sectors_[s].tags_[a], sectors_[s].tags_[b]);
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_MATRIX_CONFIGURATION_HPP
