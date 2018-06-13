// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// HS spin domain.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_HS_SPIN_DOMAIN_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_HS_SPIN_DOMAIN_HPP

#include <vector>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

enum HS_spin_states { HS_DN = -1, HS_ZERO = 0, HS_UP = 1 };
typedef HS_spin_states HS_spin_states_type;

class HS_spin_domain {
public:
  typedef HS_spin_states_type element_type;

  static int get_size() {
    return 3;
  }

  static std::vector<HS_spin_states_type>& get_elements() {
    static std::vector<HS_spin_states_type> v = initialize_elements();
    return v;
  }

  static int to_coordinate(element_type spin);

private:
  static std::vector<HS_spin_states_type> initialize_elements();
};

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_HS_SPIN_DOMAIN_HPP
