// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This class organizes the vertex configuration for ct-int.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_SOLVER_CONFIGURATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_SOLVER_CONFIGURATION_HPP

#include <array>
#include <numeric>
#include <vector>
#include <vector>

#include "dca/io/buffer.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/utils.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/util/containers/random_access_map.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

class SolverConfiguration : public MatrixConfiguration {
public:
  using BaseClass = MatrixConfiguration;
  template <class T>
  using HostVector = dca::linalg::util::HostVector<T>;

  SolverConfiguration() {}
  SolverConfiguration(const SolverConfiguration& other) = default;
  SolverConfiguration(SolverConfiguration&& other) = default;
  SolverConfiguration(double beta, int n_bands, const InteractionVertices& H_int,
                      double double_move_prob = 0);

  SolverConfiguration& operator=(const SolverConfiguration& other) = default;
  SolverConfiguration& operator=(SolverConfiguration&& other) = default;

  bool operator==(const SolverConfiguration& rhs) const;

  template <class RngType>
  void insertRandom(RngType& rng);

  // Returns the indices of the removal candidates. -1 stands for a missing candidate.
  std::array<int, 2> randomRemovalCandidate(const std::array<double, 3>& rvals);
  // Similar to the method above.
  // Precondition: no double moves are proposed.
  std::array<int, 2> randomRemovalCandidate(double rval) {
    return randomRemovalCandidate({rval, -1, -1});
  }

  // Out: indices. Appends the result of the search to indices.
  template <class Alloc>
  void findIndices(std::vector<int, Alloc>& indices, unsigned config_index, int s) const;

  // Remove elements with index in 'remove' by copying elements from the end.
  // In: remove. List of vertices to be removed. Sorted on exit.
  // Out: sector remove, sector_from. Lists of matrix indices to be moved to maintain consistency.
  void moveAndShrink(std::array<HostVector<int>, 2>& sector_from,
                     std::array<HostVector<int>, 2>& sector_remove, std::vector<int>& remove);

  // In/Out: v. Matrix indices are updated.
  void push_back(Vertex& v);

  int nPartners(int vertex_index) const;

  void commitInsertion(int idx);
  void markForRemoval(int idx);

  std::size_t size() const {
    return vertices_.size();
  }
  const Vertex& operator[](const int i) const {
    assert(i < (int)size());
    return vertices_[i];
  }

  const Vertex& back() const {
    assert(vertices_.size());
    return vertices_.back();
  }

  auto getTag(int i) const {
    return vertices_[i].tag;
  }

  inline double getStrength(int vertex_index) const;
  inline short getSign(int vertex_index) const;

  double getDoubleUpdateProb() const {
    return double_insertion_prob_;
  }

  unsigned short lastInsertionSize() const {
    return last_insertion_size_;
  }

  std::array<int, 2> sizeIncrease() const;

  bool checkConsistency() const;
  // Return index corresponding to tag.
  int findTag(std::uint64_t tag) const;

  auto possiblePartners(unsigned idx) const {
    assert(idx < vertices_.size());
    return H_int_->possiblePartners(vertices_[idx].interaction_id);
  }

  friend io::Buffer& operator<<(io::Buffer& buff, const SolverConfiguration& config);
  friend io::Buffer& operator>>(io::Buffer& buff, SolverConfiguration& config);

private:
  const InteractionElement& coordinates(int v_idx) const;
  void addSectorSizes(int idx, std::array<int, 2>& sizes) const;

  template <class Rng>
  bool doDoubleUpdate(Rng& rng) const;

  const double double_insertion_prob_ = 0;

  std::vector<Vertex> vertices_;
  // Connection from Vertices to MatrixConfiguration elements.
  std::array<std::vector<ConfigRef>, 2> matrix_config_indices_;

  using BaseClass::H_int_;

  std::vector<details::VertexTypeList> existing_;

  // Temporary storage for the removal candidate method.
  std::vector<const details::VertexTypeList*> partners_lists_;

  unsigned short last_insertion_size_ = 1;
  const double max_tau_ = 0;
  const int n_bands_ = 0;

  dca::util::RandomAccessMap<std::size_t, unsigned> anhilatable_indices_;
  std::uint64_t current_tag_ = 0;
};

template <class Rng>
void SolverConfiguration::insertRandom(Rng& rng) {
  const std::pair<short, short> indices = H_int_->getInsertionIndices(rng, double_insertion_prob_);
  const double tau = rng() * max_tau_;
  const bool aux_spin = rng() > 0.5;

  Vertex v(aux_spin, indices.first, current_tag_++, tau);
  push_back(v);

  if (indices.second != -1) {
    const double tau2 = rng() * max_tau_;
    const bool aux_spin2 = rng() > 0.5;
    Vertex v2(aux_spin2, indices.second, current_tag_++, tau2);
    push_back(v2);
    last_insertion_size_ = 2;
  }
  else {
    last_insertion_size_ = 1;
  }
  assert(2 * size() == getSector(0).size() + getSector(1).size());
}

template <class Rng>
bool SolverConfiguration::doDoubleUpdate(Rng& rng) const {
  if (double_insertion_prob_ == 0)
    return false;
  else if (double_insertion_prob_ == 1)
    return true;
  else
    return rng() < double_insertion_prob_;
}

template <class Alloc>
void SolverConfiguration::findIndices(std::vector<int, Alloc>& indices, unsigned config_index,
                                      int s) const {
  const auto& v = vertices_[config_index];
  for (int leg = 0; leg < 2; ++leg) {
    if (v.spins[leg] == s) {
      indices.push_back(v.matrix_config_indices[leg]);
    }
  }
}

inline short SolverConfiguration::getSign(const int vertex_index) const {
  return getStrength(vertex_index) > 0 ? 1 : -1;
}

inline double SolverConfiguration::getStrength(int vertex_index) const {
  assert(vertex_index >= 0 && vertex_index < (int)size());
  return (*H_int_)[vertices_[vertex_index].interaction_id].w;
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_SOLVER_CONFIGURATION_HPP
