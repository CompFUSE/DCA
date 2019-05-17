// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author:Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This class organizes the vertex configuration for ct-int.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_CONFIGURATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_CONFIGURATION_HPP

#include <array>
#include <vector>

#include "dca/io/buffer.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"
#include "dca/linalg/device_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

class SolverConfiguration : public MatrixConfiguration {
public:
  using BaseClass = MatrixConfiguration;

  inline SolverConfiguration() {}
  inline SolverConfiguration(const SolverConfiguration& other) = default;
  inline SolverConfiguration(double beta, int n_bands, const InteractionVertices& H_int,
                             double double_move_prob = 0);

  inline SolverConfiguration& operator=(const SolverConfiguration& other);
  inline SolverConfiguration& operator=(SolverConfiguration&& other);

  bool operator==(const SolverConfiguration& rhs) const;

  template <class RngType>
  void insertRandom(RngType& rng);

  template <class RngType>
  std::pair<short, short> randomRemovalCandidate(RngType& rng) const;

  // Move elements from 'from' to 'to' and shrink at the back by 'from' elements.
  template <class Alloc>
  void moveAndShrink(const std::array<std::vector<int, Alloc>, 2>& sector_from,
                     const std::array<std::vector<int, Alloc>, 2>& sector_to,
                     const std::vector<int>& from, const std::vector<int>& to);

  inline void swapVertices(short i, short j);

  inline void pop(int n = 1);
  inline void push_back(const Vertex& v);

  inline int occupationNumber(int vertex_index) const;
  inline int occupationNumber(const Vertex& vertex) const;

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

  const InteractionElement& coordinates(int v_idx) const {
    assert(v_idx >= 0 and v_idx < (int)size());
    return (*H_int_)[vertices_[v_idx].interaction_id];
  }

  int get_bands() const {
    return n_bands_;
  }

  const std::vector<ushort>& getExisting(const short id) const {
    return existing_[id];
  }

  inline ushort getLeftNu(const int matrix_index) const;
  inline ushort getRightNu(const int matrix_index) const;
  inline ushort getLeftR(const int matrix_index) const;
  inline ushort getRightR(const int matrix_index) const;
  inline double getTau(const int matrix_index) const;
  inline bool getAuxSpin(int matrix_index) const;
  auto getTag(int i) const {
    return vertices_[i].tag;
  }

  inline double getStrength(int vertex_index) const;
  inline short getSign(int vertex_index) const;

  ushort lastInsertionSize() const {
    return last_insertion_size_;
  }

  inline std::array<int, 2> sizeIncrease() const;

  friend io::Buffer& operator<<(io::Buffer& buff, const SolverConfiguration& config);
  friend io::Buffer& operator>>(io::Buffer& buff, SolverConfiguration& config);

private:
  inline void addSectorSizes(int idx, std::array<int, 2>& sizes) const;

protected:
  // List of points entering into the first or second member of g0
  const double double_insertion_prob_ = 0;

  std::vector<Vertex> vertices_;

  const InteractionVertices* H_int_ = nullptr;
  std::vector<std::vector<ushort>> existing_;
  ushort last_insertion_size_ = 1;
  const double max_tau_ = 0;
  const int n_bands_ = 0;

  std::uint64_t current_tag_ = 0;
};

SolverConfiguration& SolverConfiguration::operator=(const SolverConfiguration& other) {
  BaseClass::operator=(other);
  H_int_ = other.H_int_;
  vertices_ = other.vertices_;
  existing_ = other.existing_;
  return *this;
}

SolverConfiguration& SolverConfiguration::operator=(SolverConfiguration&& other) {
  H_int_ = other.H_int_;
  BaseClass::operator=(other);
  vertices_ = std::move(other.vertices_);
  return *this;
}

SolverConfiguration::SolverConfiguration(const double beta, const int n_bands,
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
  const std::pair<short, short> indices = H_int_->getInsertionIndices(rng());
  const double tau = rng() * max_tau_;
  const bool aux_spin = rng() > 0.5;

  push_back(Vertex{aux_spin, ushort(indices.first), ++current_tag_, tau});

  if (double_insertion_prob_) {
    existing_[indices.first].push_back(vertices_.size() - 1);
    // TODO: generalize to multiband n_bands > 2
    if (indices.second != -1 and rng() < double_insertion_prob_) {
      const double tau2 = rng() * max_tau_;
      const bool aux_spin2 = rng() > 0.5;
      // TODO: decide if tag is same or not.
      push_back(Vertex{aux_spin2, ushort(indices.second), ++current_tag_, tau2});
      existing_[indices.second].push_back(vertices_.size() - 1);
      last_insertion_size_ = 2;
    }
    else
      last_insertion_size_ = 1;
  }
  assert(2 * size() == getSector(0).size() + getSector(1).size());
}

template <class RngType>
std::pair<short, short> SolverConfiguration::randomRemovalCandidate(RngType& rng) const {
  std::pair<short, short> candidates(rng() * size(), -1);
  if (double_insertion_prob_) {
    const int partner_id = (*H_int_)[vertices_[candidates.first].interaction_id].partner_id;
    if (partner_id != -1 and rng() < double_insertion_prob_) {
      const std::size_t n_partners = existing_[partner_id].size();
      if (n_partners)
        candidates.second = existing_[partner_id][int(rng() * n_partners)];
    }
  }
  assert(candidates.first < size() and candidates.second < int(size()));
  return candidates;
}

void SolverConfiguration::push_back(const Vertex& v) {
  vertices_.push_back(v);
  BaseClass::addVertex(v);

  if (double_insertion_prob_)
    existing_[v.interaction_id].push_back(size() - 1);
}

void SolverConfiguration::pop(const int n) {
  assert(n <= (int)size());
  assert(n == 1 or (*H_int_)[vertices_[size() - 2].interaction_id].partner_id ==
                       vertices_[size() - 1].interaction_id);
  if (double_insertion_prob_) {
    for (int i = 1; i <= n; ++i) {
      // TODO improve.
      const int type = vertices_[size() - i].interaction_id;
      auto& list = existing_[type];
      const auto iter = std::find(list.begin(), list.end(), size() - i);
      list.erase(iter);
    }
  }
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

void SolverConfiguration::addSectorSizes(int idx, std::array<int, 2>& sizes) const {
  auto spin = [=](ushort nu) { return nu >= n_bands_; };
  const auto& nu = coordinates(idx).nu;
  ++sizes[spin(nu[0])];
  ++sizes[spin(nu[2])];
}

short SolverConfiguration::getSign(const int vertex_index) const {
  return (*H_int_)[vertices_[vertex_index].interaction_id].w >= 0 ? 1 : -1;
}

double SolverConfiguration::getStrength(int vertex_index) const {
  assert(vertex_index >= 0 && vertex_index < (int)size());
  return (*H_int_)[vertices_[vertex_index].interaction_id].w;
}

int SolverConfiguration::occupationNumber(int vertex_index) const {
  return existing_[vertices_[vertex_index].interaction_id].size();
}
int SolverConfiguration::occupationNumber(const Vertex& vertex) const {
  return existing_[vertex.interaction_id].size();
}

template <class Alloc>
inline void SolverConfiguration::moveAndShrink(
    const std::array<std::vector<int, Alloc>, 2>& sector_from,
    const std::array<std::vector<int, Alloc>, 2>& sector_to, const std::vector<int>& from,
    const std::vector<int>& to) {
  for (int s = 0; s < 2; ++s) {
    auto& sector = BaseClass::getEntries(s);
    auto& tags = BaseClass::sectors_[s].tags_;
    for (int i = 0; i < sector_from[s].size(); ++i) {
      sector[sector_to[s][i]] = sector[sector_from[s][i]];
      tags[sector_to[s][i]] = tags[sector_from[s][i]];
    }
    sector.erase(sector.end() - sector_to[s].size(), sector.end());
    tags.erase(tags.end() - sector_to[s].size(), tags.end());
  }

  if (double_insertion_prob_) {
    for (int i = 0; i < to.size(); ++i) {
      const int type = vertices_[size() - i].interaction_id;
      auto& list = existing_[type];
      const auto iter = std::find(list.begin(), list.end(), size() - i);
      list.erase(iter);
    }
  }

  for (int i = 0; i < from.size(); ++i)
    vertices_[to[i]] = vertices_[from[i]];
  vertices_.erase(vertices_.end() - to.size(), vertices_.end());
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
  // Swap occupation list entries.
  if (double_insertion_prob_) {
    auto& list1 = existing_[vertices_[i].interaction_id];
    auto& list2 = existing_[vertices_[j].interaction_id];

    auto it1 = std::find(list1.begin(), list1.end(), i);
    assert(it1 != list1.end());
    *it1 = j;
    auto it2 = std::find(list2.begin(), list2.end(), j);
    assert(it2 != list1.end());
    *it2 = i;
  }
  // Swap the vertices.
  std::swap(vertices_[i], vertices_[j]);
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_CT_INT_CONFIGURATION_HPP
