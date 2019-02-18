// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This files implements the methods of frequency_exchange_domain.hpp.

#include <stdexcept>

#include "dca/phys/domains/time_and_frequency/frequency_exchange_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

// Static members initialization.
std::vector<int> FrequencyExchangeDomain::elements_;
bool FrequencyExchangeDomain::initialized_ = false;
int FrequencyExchangeDomain::extension_size_ = -1;

void FrequencyExchangeDomain::initialize(const bool compute_all_transfers,
                                         const int frequency_transfer) {
  if (compute_all_transfers) {
    if (frequency_transfer < 0)
      throw(
          std::logic_error("The frequency transfer for multiple transfers must be non-negative."));
    elements_.resize(frequency_transfer + 1);
    int idx_value = 0;
    for (int& elem : elements_)
      elem = idx_value++;
  }

  else {
    elements_ = std::vector<int>{frequency_transfer};
  }

  // Compute the extension size.
  extension_size_ = 0;
  for (auto el : elements_)
    extension_size_ = std::max(extension_size_, std::abs(el));

  initialized_ = true;
}

}  // domains
}  // phys
}  // dca
