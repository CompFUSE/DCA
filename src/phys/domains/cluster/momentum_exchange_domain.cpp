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

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

// Static members initialization.
std::vector<int> MomentumExchangeDomain::elements_;
bool MomentumExchangeDomain::initialized_ = false;

void MomentumExchangeDomain::initialize(const bool compute_all_transfers, const int transfer_index,
                                        const int cluster_size) {
  if (compute_all_transfers) {
    ;
    elements_.resize(cluster_size);
    int idx_value = 0;
    for (int& elem : elements_)
      elem = idx_value++;
  }

  else {
    elements_ = std::vector<int>{transfer_index};
  }

  initialized_ = true;
}

}  // domains
}  // phys
}  // dca
