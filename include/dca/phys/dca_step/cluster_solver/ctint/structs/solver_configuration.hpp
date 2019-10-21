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

#include "dca/io/buffer.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/utils.hpp"
#include "dca/linalg/device_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

class SolverConfiguration : public MatrixConfiguration {
public:
  using BaseClass = MatrixConfiguration;

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

  // TODO: pass output as argument.
  template <class RngType>
  std::vector<int> randomRemovalCandidate(RngType& rng, double removal_rand);
  template <class RngType>
  std::vector<int> randomRemovalCandidate(RngType& rng);

  // Remove elements with id remove by copying elements from the end.
  // Postcondition: from and sector_from contains the indices that where moved into the place
  //                of the old elements.
  // Postcondition: all vectors are ordered, resized, and ready to be used for moving the rows of
  //                the M matrix
  // In/Out: sector_remove, remove.
  // Out: sector_from, from.
  template <class Alloc>
  void moveAndShrink(std::array<std::vector<int, Alloc>, 2>& sector_from,
                     std::array<std::vector<int, Alloc>, 2>& sector_remove, std::vector<int>& from,
                     std::vector<int>& remove);

  inline void swapVertices(short i, short j);

  inline void pop(int n = 1);
  inline void push_back(const Vertex& v);

  inline int nPartners(int vertex_index) const;

  void prepareForSubmatrixUpdate() {
    removable_.resize(vertices_.size());
    std::iota(removable_.begin(), removable_.end(), 0);
  }

  void commitInsertion(int idx) {
    assert(idx < size() && idx >= removable_.size());
    removable_.push_back(idx);
    if (double_insertion_prob_) {
      const auto tag = vertices_[idx].tag;
      auto& list = existing_[vertices_[idx].interaction_id];
      list.push_back(tag);
    }
  }

  void markForRemoval(int idx) {
    removable_.erase(std::find(removable_.begin(), removable_.end(), idx));
    if (double_insertion_prob_) {
      const auto tag = vertices_[idx].tag;
      auto& list = existing_[vertices_[idx].interaction_id];
      list.erase(std::find(list.begin(), list.end(), tag));
    }
  }

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

  int get_bands() const {
    return n_bands_;
  }

  auto getTag(int i) const {
    return vertices_[i].tag;
  }

  inline double getStrength(int vertex_index) const;
  inline short getSign(int vertex_index) const;

  ushort lastInsertionSize() const {
    return last_insertion_size_;
  }

  inline std::array<int, 2> sizeIncrease() const;

  bool checkConsistency() const;
  // Return index corresponding to tag.
  int findTag(std::uint64_t tag) const;

  auto possiblePartners() const {
    return H_int_->possiblePartners();
  }

  friend io::Buffer& operator<<(io::Buffer& buff, const SolverConfiguration& config);
  friend io::Buffer& operator>>(io::Buffer& buff, SolverConfiguration& config);

private:
  const InteractionElement& coordinates(int v_idx) const {
    assert(v_idx >= 0 and v_idx < (int)size());
    return (*H_int_)[vertices_[v_idx].interaction_id];
  }
  void addSectorSizes(int idx, std::array<int, 2>& sizes) const;

  template <class Rng>
  bool doDoubleUpdate(Rng& rng) const {
    if (double_insertion_prob_ == 0)
      return false;
    else if (double_insertion_prob_ == 1)
      return true;
    else
      return rng() < double_insertion_prob_;
  }

  // List of points entering into the first or second member of g0
  const double double_insertion_prob_ = 0;

  std::vector<Vertex> vertices_;

  const InteractionVertices* H_int_ = nullptr;
  std::vector<std::vector<std::uint64_t>> existing_;
  std::vector<const std::vector<std::size_t>*> partners_lists_;
  std::vector<int> removable_;
  ushort last_insertion_size_ = 1;
  const double max_tau_ = 0;
  const int n_bands_ = 0;

  std::uint64_t current_tag_ = 0;
};

inline SolverConfiguration::SolverConfiguration(const double beta, const int n_bands,
                                                const InteractionVertices& H_int,
                                                const double double_insertion)
    : MatrixConfiguration(&H_int, n_bands),
      double_insertion_prob_(double_insertion),

      H_int_(&H_int),
      existing_(double_insertion ? H_int_->size() : 0),
      max_tau_(beta),
      n_bands_(n_bands) {}

template <class RngType>
void SolverConfiguration::insertRandom(RngType& rng) {
  const std::pair<short, short> indices = H_int_->getInsertionIndices(rng, double_insertion_prob_);
  const double tau = rng() * max_tau_;
  const bool aux_spin = rng() > 0.5;

  push_back(Vertex{aux_spin, ushort(indices.first), current_tag_++, tau});

  if (double_insertion_prob_) {
    // TODO: generalize to multiband n_bands > 2
    if (indices.second != -1 && double_insertion_prob_) {
      const double tau2 = rng() * max_tau_;
      const bool aux_spin2 = rng() > 0.5;
      // TODO: decide if tag is same or not.
      push_back(Vertex{aux_spin2, ushort(indices.second), current_tag_++, tau2});
      last_insertion_size_ = 2;
    }
    else
      last_insertion_size_ = 1;
  }
  assert(2 * size() == getSector(0).size() + getSector(1).size());
}

template <class RngType>
std::vector<int> SolverConfiguration::randomRemovalCandidate(RngType& rng) {
  return randomRemovalCandidate(rng, rng());
}

template <class RngType>
std::vector<int> SolverConfiguration::randomRemovalCandidate(RngType& rng, double removal_rand) {
  std::vector<int> candidates;
  if (removable_.size() == 0)
    return candidates;

  candidates.push_back(removable_[removal_rand * removable_.size()]);

  if ((*H_int_)[vertices_[candidates[0]].interaction_id].partners_id.size() &&
      doDoubleUpdate(rng)) {  // Double removal.
    partners_lists_.clear();
    for (const auto& partner_id : (*H_int_)[vertices_[candidates[0]].interaction_id].partners_id)
      partners_lists_.push_back(&existing_[partner_id]);

    const auto tag = details::getRandomElement(partners_lists_, rng());
    candidates.push_back(findTag(tag));
    assert(candidates[1] < int(size()) && candidates[1] >= 0);
  }

  assert(candidates[0] < int(size()));
  return candidates;
}

void SolverConfiguration::push_back(const Vertex& v) {
  vertices_.push_back(v);
  BaseClass::addVertex(v);
}

void SolverConfiguration::pop(const int n) {
  assert(n <= (int)size());

  const int first_idx = size() - n;
  std::array<int, 2> removal_size{0, 0};
  for (std::size_t i = first_idx; i < size(); ++i)
    addSectorSizes(i, removal_size);
  vertices_.erase(vertices_.begin() + first_idx, vertices_.end());
  BaseClass::pop(removal_size[0], removal_size[1]);
}

std::array<int, 2> SolverConfiguration::sizeIncrease() const {
  assert(last_insertion_size_ > 0);
  std::array<int, 2> result{0, 0};
  for (std::size_t i = size() - last_insertion_size_; i < size(); ++i)
    addSectorSizes(i, result);

  return result;
}

inline void SolverConfiguration::addSectorSizes(int idx, std::array<int, 2>& sizes) const {
  auto spin = [=](ushort nu) { return nu >= n_bands_; };
  const auto& nu = coordinates(idx).nu;
  ++sizes[spin(nu[0])];
  ++sizes[spin(nu[2])];
}

inline short SolverConfiguration::getSign(const int vertex_index) const {
  return getStrength(vertex_index) > 0 ? 1 : -1;
}

inline double SolverConfiguration::getStrength(int vertex_index) const {
  assert(vertex_index >= 0 && vertex_index < (int)size());
  return (*H_int_)[vertices_[vertex_index].interaction_id].w;
}

int SolverConfiguration::nPartners(int vertex_index) const {
  assert(vertex_index < vertices_.size());
  const auto& partners_id = (*H_int_)[vertices_[vertex_index].interaction_id].partners_id;
  int n_partners = 0;
  for (const auto partner_id : partners_id)
    n_partners += existing_[partner_id].size();
  return n_partners;
}

template <class Alloc>
inline void SolverConfiguration::moveAndShrink(std::array<std::vector<int, Alloc>, 2>& sector_from,
                                               std::array<std::vector<int, Alloc>, 2>& sector_remove,
                                               std::vector<int>& from, std::vector<int>& remove) {
  for (int s = 0; s < 2; ++s) {
    // Sort and prepare source array.
    auto& sector = BaseClass::sectors_[s].entries_;
    std::sort(sector_remove[s].begin(), sector_remove[s].end());
    auto& tags = BaseClass::sectors_[s].tags_;
    sector_from[s].clear();
    int source_idx = sector.size() - sector_remove[s].size();
    for (; source_idx < sector.size(); ++source_idx) {
      while (std::binary_search(sector_remove[s].begin(), sector_remove[s].end(), source_idx))
        ++source_idx;
      if (source_idx < sector.size())
        sector_from[s].push_back(source_idx);
    }

    // Move configuration elements.
    for (int i = 0; i < sector_from[s].size(); ++i) {
      sector[sector_remove[s][i]] = sector[sector_from[s][i]];
      tags[sector_remove[s][i]] = tags[sector_from[s][i]];
    }
    // Shrink sector configuration.
    sector.erase(sector.end() - sector_remove[s].size(), sector.end());
    tags.erase(tags.end() - sector_remove[s].size(), tags.end());
  }

  std::sort(remove.begin(), remove.end());
  from.clear();
  int source_idx = size() - remove.size();
  for (; source_idx < size(); ++source_idx) {
    while (std::binary_search(remove.begin(), remove.end(), source_idx))
      ++source_idx;
    if (source_idx < size())
      from.push_back(source_idx);
  }

  // Move and shrink configuration.
  for (int i = 0; i < from.size(); ++i)
    vertices_[remove[i]] = vertices_[from[i]];
  vertices_.erase(vertices_.end() - remove.size(), vertices_.end());

  assert(checkConsistency());
}

inline bool SolverConfiguration::operator==(const SolverConfiguration& rhs) const {
  bool result = true;
  result &= vertices_ == rhs.vertices_;
  result &= existing_ == rhs.existing_;
  result &= max_tau_ == rhs.max_tau_;
  result &= n_bands_ == rhs.n_bands_;
  result &= double_insertion_prob_ == rhs.double_insertion_prob_;

  result &= static_cast<BaseClass>(*this) == static_cast<BaseClass>(rhs);

  return result;
}

void SolverConfiguration::swapVertices(const short i, const short j) {
  std::swap(vertices_[i], vertices_[j]);
}

inline bool SolverConfiguration::checkConsistency() const {
  for (auto v : vertices_) {
    int count = 0;
    for (int s = 0; s < 2; ++s) {
      auto indices = findIndices(v.tag, s);
      count += indices.size();
      for (auto idx : indices) {
        const auto& matrix_el = sectors_[s].entries_[idx];
        if (matrix_el.tau_ != v.tau)
          return false;
      }
    }
    if ((double_insertion_prob_ == 0 && count != 2) || !count)
      return false;
  }

  if (double_insertion_prob_) {
    for (const auto& v : vertices_) {
      const auto& list = existing_[v.interaction_id];
      if (std::find(list.begin(), list.end(), v.tag) == list.end())
        return false;
      // Check total size.
      int size_sum = 0;
      for (const auto& list : existing_)
        size_sum += list.size();
      if (size_sum != vertices_.size())
        return false;
    }
  }
  return true;
}

inline int SolverConfiguration::findTag(std::uint64_t tag) const {
  for (int i = 0; i < vertices_.size(); ++i)
    if (vertices_[i].tag == tag)
      return i;
  return -1;
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_SOLVER_CONFIGURATION_HPP
