// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// HS field sign domain.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_HS_FIELD_SIGN_DOMAIN_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_HS_FIELD_SIGN_DOMAIN_HPP

#include <vector>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

enum HS_field_sign { HS_FIELD_DN = -1, HS_FIELD_UP = 1 };
typedef HS_field_sign HS_field_sign_type;

class HS_field_sign_domain {
public:
  typedef HS_field_sign_type element_type;

  static int get_size() {
    return 2;
  }

  static std::vector<HS_field_sign_type>& get_elements() {
    static std::vector<HS_field_sign_type> v = initialize_elements();
    return v;
  }

  static int to_coordinate(HS_field_sign sign);

private:
  static std::vector<HS_field_sign_type> initialize_elements();
};

}  // ctaux
}  // phys
}  // solver
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_HS_FIELD_SIGN_DOMAIN_HPP
