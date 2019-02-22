// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a base class containing the shared methods between the CPU and GPU
// implementations of the class CachedNdft.

#ifndef DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_BASE_HPP
#define DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_BASE_HPP

#include <array>
#include <cassert>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/linalg/util/allocators/vectors_typedefs.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/triple.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
class CachedNdftBase {
protected:
  CachedNdftBase();

  template <class Configuration>
  void sortConfiguration(const Configuration& configuration);

protected:
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using BRDmn = func::dmn_variadic<BDmn, RDmn>;
  using Triple = details::Triple<ScalarType>;
  template <class T>
  using HostVector = linalg::util::HostVector<T>;

  HostVector<ScalarType> w_;

  std::array<HostVector<Triple>, 2> indexed_config_;

  std::array<std::vector<int>, 2> start_index_;
  std::array<std::vector<int>, 2> end_index_;

  const HostVector<Triple>& config_left_;
  const HostVector<Triple>& config_right_;
  std::vector<int>& start_index_left_;
  std::vector<int>& start_index_right_;
  std::vector<int>& end_index_left_;
  std::vector<int>& end_index_right_;

  const int n_orbitals_;
};

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
CachedNdftBase<ScalarType, RDmn, WDmn, WPosDmn, non_density_density>::CachedNdftBase()
    : w_(),
      config_left_(indexed_config_[0]),
      config_right_(non_density_density ? indexed_config_[1] : indexed_config_[0]),
      start_index_left_(start_index_[0]),
      start_index_right_(non_density_density ? start_index_[1] : start_index_[0]),
      end_index_left_(end_index_[0]),
      end_index_right_(non_density_density ? end_index_[1] : end_index_[0]),
      n_orbitals_(BDmn::dmn_size() * RDmn::dmn_size()) {
  for (const auto elem : WDmn::parameter_type::get_elements())
    w_.push_back(static_cast<ScalarType>(elem));

  start_index_left_.resize(n_orbitals_);
  start_index_right_.resize(n_orbitals_);
  end_index_right_.resize(n_orbitals_);
  end_index_left_.resize(n_orbitals_);
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
template <class Configuration>
void CachedNdftBase<ScalarType, RDmn, WDmn, WPosDmn, non_density_density>::sortConfiguration(
    const Configuration& configuration) {
  const int n_b = BDmn::dmn_size();
  const int n_v = configuration.size();
  constexpr int n_sides = non_density_density ? 2 : 1;

  auto orbital = [&](const int i, const int side) {
    const auto& vertex = configuration[i];
    if (side)  // Switch left and right band as M is an inverse matrix.
      return vertex.get_left_band() + n_b * vertex.get_left_site();
    else
      return vertex.get_right_band() + n_b * vertex.get_right_site();
  };

  for (int side = 0; side < n_sides; ++side) {
    auto& config_side = indexed_config_[side];
    config_side.resize(n_v);

    for (int l = 0; l < n_v; ++l) {
      config_side[l].orbital = orbital(l, side);
      config_side[l].tau = configuration[l].get_tau();
      config_side[l].idx = l;
    }

    sort(config_side.begin(), config_side.end());

    for (int orb = 0; orb < n_orbitals_; ++orb) {
      details::Triple<ScalarType> trp{orb, 0, 0};

      start_index_[side][orb] =
          lower_bound(config_side.begin(), config_side.end(), trp) - config_side.begin();

      trp.idx = n_v;

      end_index_[side][orb] =
          upper_bound(config_side.begin(), config_side.end(), trp) - config_side.begin();
    }
  }
}

}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_BASE_HPP
