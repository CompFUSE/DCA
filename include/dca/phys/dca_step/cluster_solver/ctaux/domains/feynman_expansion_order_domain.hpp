// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Feynman expansion order domain.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_FEYNMAN_EXPANSION_ORDER_DOMAIN_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_FEYNMAN_EXPANSION_ORDER_DOMAIN_HPP

#include <vector>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

class Feynman_expansion_order_domain {
  const static int MAX_ORDER_SQUARED = 10000;

public:
  typedef int element_type;

  static int get_size() {
    const static int size = MAX_ORDER_SQUARED;
    return size;
  }

  static std::vector<int>& get_elements() {
    static std::vector<int>& v = initialize_elements();
    return v;
  }

private:
  static std::vector<int>& initialize_elements();
};

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_FEYNMAN_EXPANSION_ORDER_DOMAIN_HPP
