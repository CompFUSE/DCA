// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements hs_field_sign_domain.hpp.

#include "dca/phys/dca_step/cluster_solver/ctaux/domains/hs_field_sign_domain.hpp"

#include <stdexcept>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

std::vector<HS_field_sign_type> HS_field_sign_domain::initialize_elements() {
  static std::vector<HS_field_sign_type> v(0);

  v.push_back(HS_FIELD_DN);
  v.push_back(HS_FIELD_UP);

  return v;
}

int HS_field_sign_domain::to_coordinate(HS_field_sign sign) {
  switch (sign) {
    case HS_FIELD_DN:
      return 0;
      break;

    case HS_FIELD_UP:
      return 1;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

}  // ctaux
}  // phys
}  // solver
}  // dca
