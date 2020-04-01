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
  Vertex() = default;
  Vertex(bool _aux_spin, unsigned short _interaction_id, std::uint64_t _tag, double _tau)
      : aux_spin(_aux_spin), interaction_id(_interaction_id), tag(_tag), tau(_tau) {
    spins.fill(-1);
    matrix_config_indices.fill(-1);
  }

  bool aux_spin;
  unsigned short interaction_id;
  std::uint64_t tag;
  double tau;

  // Reference to the matrix config
  std::array<std::uint8_t, 2> spins;
  std::array<unsigned, 2> matrix_config_indices;

  // Marks if this vertex is part of the accepted configuration.
  bool annihilatable = false;

  bool operator==(const Vertex& b) const {
    return aux_spin == b.aux_spin && interaction_id == b.interaction_id && tau == b.tau;
  }
};

struct ConfigRef {
  ConfigRef() = default;
  ConfigRef(unsigned _config_id, std::uint8_t _leg_id) : config_id(_config_id), leg_id(_leg_id) {}

  unsigned config_id;  // Index of the interaction in the SolverConfiguration.
  std::uint8_t leg_id;  // In {0, 1}. Stores if this is the first or second leg of an interaction vertex.
};

class MatrixConfiguration {
public:
  MatrixConfiguration() = default;
  MatrixConfiguration(const MatrixConfiguration& rhs) = default;
  MatrixConfiguration& operator=(const MatrixConfiguration& rhs);
  MatrixConfiguration& operator=(MatrixConfiguration&& rhs);

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

protected:
  MatrixConfiguration(const InteractionVertices* H_int, int bands);

  // In/Out: v. The vertex to be added, the matrix configuration spins and indices are updated.
  // In. config_id. Position of the vertex in the solver configuration.
  // Out: config_refs. References from matrix configuration to solver configuration to be updated.
  void addVertex(Vertex& v, unsigned config_id, std::array<std::vector<ConfigRef>, 2>& config_refs);

  //  inline void pop(unsigned short idx_up, unsigned short idx_down);

  const auto& getEntries(const int s) const {
    return sectors_[s].entries_;
  }
  auto& getEntries(const int s) {
    return sectors_[s].entries_;
  }

  const InteractionVertices* H_int_ = nullptr;
  int n_bands_ = -1;
  std::array<Sector, 2> sectors_;
};

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_MATRIX_CONFIGURATION_HPP
