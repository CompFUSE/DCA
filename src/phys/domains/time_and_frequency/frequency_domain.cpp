// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements frequency_domain.hpp.

#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

#include <cmath>  // M_PI
#include <stdexcept>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

bool frequency_domain::initialized_ = false;
const std::string frequency_domain::name_ = "frequency-domain";
std::vector<frequency_domain::element_type> frequency_domain::elements_;
std::vector<int> frequency_domain::indices_;

void frequency_domain::initialize(const ScalarType beta, const int num_freqs) {
  if (initialized_)
    throw std::logic_error("frequency_domain has already been initialzed.");

  const int size = 2 * num_freqs;

  indices_.resize(size);
  elements_.resize(size);

  for (int l = 0; l < num_freqs; ++l) {
    const int index = 2 * l + 1;

    indices_[size / 2 + 0 + l] = index;
    indices_[size / 2 - 1 - l] = -index;

    elements_[size / 2 + 0 + l] = index * M_PI / beta;
    elements_[size / 2 - 1 - l] = -index * M_PI / beta;
  }

  initialized_ = true;
}

}  // domains
}  // phys
}  // dca
