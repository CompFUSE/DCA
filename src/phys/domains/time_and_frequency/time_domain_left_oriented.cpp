// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements time_domain_left_oriented.hpp.

#include "dca/phys/domains/time_and_frequency/time_domain_left_oriented.hpp"
#include <stdexcept>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

bool time_domain_left_oriented::initialized_ = false;
const std::string time_domain_left_oriented::name_ = "time-domain-left-oriented";
std::vector<time_domain_left_oriented::element_type> time_domain_left_oriented::elements_;

void time_domain_left_oriented::initialize() {
  if (initialized_)
    throw std::logic_error("time_domain_left_oriented has already been initialized.");

  if (!time_domain::is_initialized())
    throw std::logic_error("time_domain must be initialized first.");

  elements_.resize(time_domain::get_size() - 2);

  for (std::size_t t_ind = 0; t_ind < time_domain::get_size() / 2 - 1; t_ind++)
    elements_[t_ind] = time_domain::get_elements()[t_ind];

  for (std::size_t t_ind = time_domain::get_size() / 2; t_ind < time_domain::get_size() - 1; t_ind++)
    // The conditional replaces the element that has value 'eps' with zero.
    elements_[t_ind - 1] =
        (t_ind == time_domain::get_size() / 2) ? 0 : time_domain::get_elements()[t_ind];

  initialized_ = true;
}

}  // domains
}  // phys
}  // dca
