// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
//
// This function creates a set of disorder configurations based on the
// sites in the cluster and requested density

#ifndef DCA_PHYS_DCA_LOOP_MAKE_DISORDER_HPP
#define DCA_PHYS_DCA_LOOP_MAKE_DISORDER_HPP

#include <algorithm>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>
#include "dca/phys/types/dca_shared_types.hpp"

namespace dca::phys {

namespace detail {
struct VectorHash {
  size_t operator()(const std::vector<int>& v) const {
    size_t seed = 0;
    for (int x : v) {
      // XOR current hash with hash of element + constant + bit rotations
      // This helps avoid collisions for similar vectors (e.g.,
      // [1,-1] vs [-1,1])
      // Based on the boost hash_combine function.
      constexpr size_t golden_ratio = sizeof(int) == 8 ? 0x9e3779b97f4a7c15 : 0x9e3779b9;
      seed ^= std::hash<int>{}(x) + golden_ratio + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

}  // namespace detail

template <typename Parameters>
class MakeDisorderConfigurations {
  using DST = DcaSharedTypes<Parameters>;
  using DisorderConfiguration = typename DST::DisorderConfiguration;
  using DDA = typename DST::DDA;
  using RDmn = typename DDA::RClusterDmn;
  using BDmn = typename DDA::BDmn;

public:
  void operator()(Parameters& parameters, std::vector<DisorderConfiguration>& disorder_configurations,
                  std::vector<double>& disorder_weights) {
    const int m_configs = parameters.get_disorder_num_configurations();
    assert(m_configs > 0);
    disorder_configurations.resize(m_configs);
    disorder_weights.resize(m_configs);

    const int n_bands = BDmn::dmn_size();
    const int n_rsites = RDmn::dmn_size();

    const std::string& distribution = parameters.get_disorder_distribution();
    const double potential = parameters.get_disorder_potential();
    const double density = parameters.get_disorder_density();
    const bool enforce_unique = parameters.get_disorder_unique_configs();

    const bool is_box = (distribution == "box");

    // Seed from the Monte Carlo seed, shifted by a fixed offset so the disorder realizations
    // use a stream distinct from the QMC walkers (which derive their seeds from get_seed()).
    // The seed is read on the root rank and broadcast as part of Parameters, so it is identical
    // on every rank (every rank generates the same disorder ensemble) and fixed for the whole
    // run (the ensemble is the same in every DCA iteration, i.e. quenched disorder). No separate
    // seed broadcast is needed here.
    constexpr unsigned disorder_seed_offset = 0x9e3779b9u;  // golden ratio; arbitrary but fixed
    std::mt19937 gen(static_cast<unsigned>(parameters.get_seed()) + disorder_seed_offset);
    // box:    site potential ~ Uniform[-W/2, W/2], with W = potential.
    std::uniform_real_distribution<double> box_dist(-0.5 * potential, 0.5 * potential);
    std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

    // Draw one disorder value for a single (band, site):
    //  - box:    uniform in [-W/2, W/2].
    //  - binary: +V/2 with probability `density`, else -V/2 (V = potential).
    auto draw_site_value = [&]() -> double {
      if (is_box)
        return box_dist(gen);
      return (unit_dist(gen) < density ? 1.0 : -1.0) * 0.5 * potential;
    };

    // Uniqueness only makes sense for the discrete binary distribution; a continuous
    // box draw is distinct with probability one.
    const bool tracking_unique = enforce_unique && !is_box;
    std::unordered_set<std::vector<int>, detail::VectorHash> seen;

    int generated = 0;
    int attempts = 0;
    const int max_attempts = 1000;
    while (generated < m_configs && attempts < max_attempts) {
      auto& config = disorder_configurations[generated];

      std::vector<int> pattern;  // sign pattern, only populated for the uniqueness check
      if (tracking_unique)
        pattern.reserve(static_cast<std::size_t>(n_bands) * n_rsites);

      // Fill the per-(band, site) random potential, spin-symmetrically.
      for (int isite = 0; isite < n_rsites; ++isite)
        for (int iband = 0; iband < n_bands; ++iband) {
          const double value = draw_site_value();
          config(iband, 0, isite) = value;
          config(iband, 1, isite) = value;
          if (tracking_unique)
            pattern.push_back(value > 0.0 ? 1 : -1);
        }

      if (tracking_unique && !seen.insert(pattern).second) {
        ++attempts;  // duplicate binary configuration, redraw
        continue;
      }
      ++generated;
      attempts = 0;
    }
    // If uniqueness was requested but the cluster is too small to supply m_configs
    // distinct configurations, keep only the ones actually generated.
    if (generated < m_configs) {
      disorder_configurations.resize(generated);
      disorder_weights.resize(generated);
    }

    // Equal weight for every configuration in the disorder average.
    const double weight = 1.0 / static_cast<double>(disorder_configurations.size());
    std::fill(disorder_weights.begin(), disorder_weights.end(), weight);
  }
};
}  // namespace dca::phys

#endif
