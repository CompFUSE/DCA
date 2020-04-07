// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This file implements solver_configuration.hpp.

#include "dca/phys/dca_step/cluster_solver/ctint/structs/solver_configuration.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

SolverConfiguration::SolverConfiguration(const double beta, const int n_bands,
                                         const InteractionVertices& H_int,
                                         const double double_insertion)
    : MatrixConfiguration(&H_int, n_bands),
      double_insertion_prob_(double_insertion),

      existing_(double_insertion ? H_int_->size() : 0),
      max_tau_(beta),
      n_bands_(n_bands) {}

std::array<int, 2> SolverConfiguration::randomRemovalCandidate(const std::array<double, 3>& rvals) {
  std::array<int, 2> candidates{-1, -1};
  if (anhilatable_indices_.size() == 0)
    return candidates;

  const std::size_t index = rvals[0] * anhilatable_indices_.size();
  candidates[0] = anhilatable_indices_[index];

  if (rvals[1] < double_insertion_prob_ &&
      (*H_int_)[vertices_[candidates[0]].interaction_id].partners_id.size()) {  // Double removal.
    partners_lists_.clear();
    for (const auto& partner_id : (*H_int_)[vertices_[candidates[0]].interaction_id].partners_id)
      partners_lists_.push_back(&existing_[partner_id]);

    const auto tag = details::getRandomElement(partners_lists_, rvals[2]);
    if (tag != -1) {
      candidates[1] = findTag(tag);
      assert(candidates[1] < int(size()) && candidates[1] >= 0);
      assert(vertices_[candidates[1]].annihilatable);
    }
  }

  assert(candidates[0] < int(size()));
  return candidates;
}

void SolverConfiguration::push_back(Vertex& v) {
  BaseClass::addVertex(v, size(), matrix_config_indices_);
  vertices_.push_back(v);
}

std::array<int, 2> SolverConfiguration::sizeIncrease() const {
  assert(last_insertion_size_ > 0);
  std::array<int, 2> result{0, 0};
  for (std::size_t i = size() - last_insertion_size_; i < size(); ++i)
    addSectorSizes(i, result);

  return result;
}

void SolverConfiguration::addSectorSizes(int idx, std::array<int, 2>& sizes) const {
  auto spin = [=](unsigned short nu) { return nu >= n_bands_; };
  const auto& nu = coordinates(idx).nu;
  ++sizes[spin(nu[0])];
  ++sizes[spin(nu[2])];
}

const InteractionElement& SolverConfiguration::coordinates(int v_idx) const {
  assert(v_idx >= 0 and v_idx < (int)size());
  return (*H_int_)[vertices_[v_idx].interaction_id];
}

int SolverConfiguration::nPartners(int vertex_index) const {
  assert(vertex_index < vertices_.size());
  const auto& partners_id = (*H_int_)[vertices_[vertex_index].interaction_id].partners_id;
  int n_partners = 0;
  for (const auto partner_id : partners_id)
    n_partners += existing_[partner_id].size();
  return n_partners;
}

void SolverConfiguration::commitInsertion(int idx) {
  assert(vertices_[idx].annihilatable == false);
  vertices_[idx].annihilatable = true;

  anhilatable_indices_.insert(vertices_[idx].tag, idx);

  if (double_insertion_prob_) {
    const auto tag = vertices_[idx].tag;
    auto& list = existing_[vertices_[idx].interaction_id];
    // TODO: use set instead of map
    list.insert(tag, tag);
  }
}

void SolverConfiguration::markForRemoval(int idx) {
  assert(vertices_[idx].annihilatable == true);
  anhilatable_indices_.erase(vertices_[idx].tag);

  vertices_[idx].annihilatable = false;

  if (double_insertion_prob_) {
    const auto tag = vertices_[idx].tag;
    auto& list = existing_[vertices_[idx].interaction_id];
    list.erase(tag);
  }
}

void SolverConfiguration::moveAndShrink(std::array<HostVector<int>, 2>& sector_from,
                                        std::array<HostVector<int>, 2>& sector_remove,
                                        std::vector<int>& remove) {
  for (int s = 0; s < 2; ++s) {
    // Prepare sorted matrix removal indices.
    sector_remove[s].clear();
    sector_from[s].clear();
    for (int idx : remove) {
      findIndices(sector_remove[s], idx, s);
    }
    std::sort(sector_remove[s].begin(), sector_remove[s].end());

    auto& sector = BaseClass::sectors_[s].entries_;
    auto& matrix_config_indices = matrix_config_indices_[s];
    const unsigned expected_new_size = sector.size() - sector_remove[s].size();

    int right_index = sector_remove[s].size() - 1;
    int living_index = sector.size() - 1;
    for (int left_index = 0; left_index < sector_remove[s].size(); ++left_index) {
      const int dead_index = sector_remove[s][left_index];
      while (right_index >= left_index &&
             living_index == sector_remove[s].back()) {  // Living not found: remove from back.
        sector_remove[s].pop_back();
        sector.pop_back();
        matrix_config_indices.pop_back();
        --living_index;
        --right_index;
      }
      if (dead_index >= living_index)
        break;

      // Update matrix config index
      const auto [config_id, leg] = matrix_config_indices[living_index];
      vertices_[config_id].spins[leg] = s;
      vertices_[config_id].matrix_config_indices[leg] = dead_index;

      // Move living spin into dead spin.
      sector[dead_index] = sector[living_index];
      matrix_config_indices[dead_index] = matrix_config_indices[living_index];
      sector.pop_back();
      matrix_config_indices.pop_back();

      sector_from[s].push_back(living_index);
      --living_index;
    }

    assert(sector_from[s].size() == sector_remove[s].size());
    if (sector.size() != expected_new_size)
      throw(std::logic_error("Size mismatch."));
  }

  std::sort(remove.begin(), remove.end());
  int right_index = remove.size() - 1;
  int living_index = size() - 1;

  for (int left_index = 0; left_index <= right_index; ++left_index) {
    const int dead_index = remove[left_index];
    assert(vertices_[dead_index].annihilatable == false);
    while (right_index >= left_index &&
           living_index == remove[right_index]) {  // Living not found: remove from back.
      vertices_.pop_back();
      --right_index;
      --living_index;
    }
    if (dead_index >= living_index)
      break;

    // Move living vertex in place of the dead one.
    vertices_[dead_index] = vertices_[living_index];
    vertices_.pop_back();

    // Update config index
    const auto& v = vertices_[dead_index];
    for (int leg = 0; leg < 2; ++leg) {
      matrix_config_indices_[v.spins[leg]].at(v.matrix_config_indices[leg]).config_id = dead_index;
    }
    // Update tag index
    if (v.annihilatable)
      anhilatable_indices_.find(v.tag) = dead_index;

    --living_index;
  }

  assert(checkConsistency());
}

bool SolverConfiguration::operator==(const SolverConfiguration& rhs) const {
  bool result = true;
  result &= vertices_ == rhs.vertices_;
  result &= max_tau_ == rhs.max_tau_;
  result &= n_bands_ == rhs.n_bands_;

  result &= static_cast<BaseClass>(*this) == static_cast<BaseClass>(rhs);

  return result;
}

bool SolverConfiguration::checkConsistency() const {
  if (2 * size() != sectors_[0].size() + sectors_[1].size())
    return false;

  for (int i = 0; i < vertices_.size(); ++i) {
    const auto& v = vertices_[i];
    for (int leg = 0; leg < 2; ++leg) {
      const auto& config_id = matrix_config_indices_[v.spins[leg]][v.matrix_config_indices[leg]];
      if (config_id.leg_id != leg)
        return false;
      if (config_id.config_id != i)
        return false;

      const auto& matrix_elem = sectors_[v.spins[leg]][v.matrix_config_indices[leg]];
      if (v.tau != matrix_elem.get_tau())
        return false;
    }
  }

  // Check annihilatable.
  bool annhilatable_consitency = true;
  unsigned idx = 0;
  for (const auto& v : vertices_) {
    if (v.annihilatable) {
      if (idx != anhilatable_indices_.find(v.tag))
        annhilatable_consitency = false;
    }
    else {  // !v.annihilatable
      try {
        anhilatable_indices_.find(v.tag);  // must throw.
        annhilatable_consitency = false;
      }
      catch (...) {
      }
    }
    ++idx;
  }

  if (!annhilatable_consitency) {
    std::cerr << "Non consistant annihilatable tags." << std::endl;
    return false;
  }

  if (double_insertion_prob_) {
    for (const auto& v : vertices_) {
      // check tags.
      if (v.annihilatable) {
        const auto& list = existing_[v.interaction_id];
        try {
          list.find(v.tag);
        }
        catch (...) {
          return false;
        }
      }
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

int SolverConfiguration::findTag(std::uint64_t tag) const {
  return anhilatable_indices_.find(tag);
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
