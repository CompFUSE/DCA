// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
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

void frequency_domain::initialize(const double beta, const int num_freqs) {
  if (initialized_)
    throw std::logic_error("frequency_domain has already been initialzed.");

  get_basis()[0] = (2. * M_PI) / beta;
  get_inverse_basis()[0] = beta / (2. * M_PI);

  get_size() = 2 * num_freqs;

  get_elements().resize(get_size());
  get_integer_wave_vectors().resize(get_size());

  for (int l = 0; l < num_freqs; ++l) {
    get_elements()[get_size() / 2 + 0 + l] = M_PI / beta * (1 + 2 * l);
    get_elements()[get_size() / 2 - 1 - l] = -M_PI / beta * (1 + 2 * l);

    get_integer_wave_vectors()[get_size() / 2 + 0 + l] = (1 + 2 * l);
    get_integer_wave_vectors()[get_size() / 2 - 1 - l] = -(1 + 2 * l);
  }

  initialized_ = true;
}

}  // domains
}  // phys
}  // dca
