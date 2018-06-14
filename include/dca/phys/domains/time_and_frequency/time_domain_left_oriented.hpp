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
// This class parametrizes the imaginary time domain with elements defined in the interval
// [-beta, beta), i.e. it leaves out the right edge.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_LEFT_ORIENTED_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_LEFT_ORIENTED_HPP

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>

#include "dca/phys/domains/time_and_frequency/time_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class time_domain_left_oriented {
public:
  using element_type = time_domain::element_type;

  static bool is_initialized() {
    return initialized_;
  }

  static const std::string& get_name() {
    return name_;
  }

  static std::size_t get_size() {
    assert(initialized_);
    return elements_.size();
  }

  // TODO: Add const qualifier when rest of the code is fixed.
  static /*const*/ std::vector<element_type>& get_elements() {
    assert(initialized_);
    return elements_;
  }

  // Initializes the elements of the imaginary time domain (without the right edge) as follows,
  // [-beta+eps, -beta+step, ..., -step, 0 step, ..., beta-step],
  // where step = beta/time_slices is defined in time_domain::initialize.
  // Reuses the elements of time_domain.
  // Preconditions: time_domain is initialized.
  // INTERNAL: Do we want the shift at -beta?
  static void initialize();

  // This method is provided for a consistent initialization interface for all domains.
  template <class parameters_t>
  static void initialize(const parameters_t& /*parameters*/) {
    initialize();
  }

private:
  static bool initialized_;
  const static std::string name_;
  static std::vector<element_type> elements_;
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_LEFT_ORIENTED_HPP
