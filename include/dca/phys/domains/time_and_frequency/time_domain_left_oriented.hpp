// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class parametrizes the time domain with N intervals, but it leaves out the right edge.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_LEFT_ORIENTED_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_LEFT_ORIENTED_HPP

#include <vector>
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class time_domain_left_oriented {
public:
  typedef double element_type;
  typedef time_domain_left_oriented this_type;

  static int& get_size() {
    static int SIZE = time_domain::get_size() - 2;
    return SIZE;
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> elements(get_size(), 0);
    return elements;
  }

  template <class parameters_t>
  static void initialize(parameters_t& parameters);
};

template <class parameters_t>
void time_domain_left_oriented::initialize(parameters_t& /*parameters*/) {
  for (int t_ind = 0; t_ind < time_domain::get_size() / 2 - 1; t_ind++)
    time_domain_left_oriented::get_elements()[t_ind] = time_domain::get_elements()[t_ind];

  for (int t_ind = time_domain::get_size() / 2; t_ind < time_domain::get_size() - 1; t_ind++)
    time_domain_left_oriented::get_elements()[t_ind - 1] =
        t_ind == time_domain::get_size() / 2 ? 0 : time_domain::get_elements()[t_ind];
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_LEFT_ORIENTED_HPP
