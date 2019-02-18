// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the methods of momentum_exchange_domain.hpp.

#include "dca/phys/domains/cluster/momentum_exchange_domain.hpp"

#include <iostream>
#include <numeric>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

// Static members initialization.
std::vector<int> MomentumExchangeDomain::elements_;
bool MomentumExchangeDomain::initialized_ = false;

void MomentumExchangeDomain::initialize(bool compute_all_transfers, int transfer_index,
                                        const std::vector<std::vector<double>>& cluster_elements,
                                        const std::vector<std::array<int, 2>>& symmetries,
                                        bool verbose) {
  initialized_ = true;

  if (!compute_all_transfers) {
    elements_ = std::vector<int>{transfer_index};
    return;
  }

  elements_.resize(cluster_elements.size());
  std::iota(elements_.begin(), elements_.end(), 0);

  for (const auto& symmetry : symmetries) {
    if (symmetry[0] != symmetry[1]) {
      const int to_remove = std::max(symmetry[0], symmetry[1]);
      elements_.erase(std::remove(elements_.begin(), elements_.end(), to_remove), elements_.end());
    }
  }

  if (verbose) {
    std::cout << "Independent reciprocal cluster sites:\n";
    for (const auto id : elements_) {
      std::cout << id << ":\t";
      for (auto k : cluster_elements[id])
        std::cout << k << "\t";
      std::cout << "\n";
    }
  }
}

}  // domains
}  // phys
}  // dca
