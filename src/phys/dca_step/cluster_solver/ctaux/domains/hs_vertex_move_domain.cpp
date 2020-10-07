// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements hs_vertex_move_domain.hpp.

#include "dca/phys/dca_step/cluster_solver/ctaux/domains/hs_vertex_move_domain.hpp"

#include <stdexcept>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

std::vector<HS_vertex_move_type> HS_vertex_move_domain::initialize_elements() {
  return std::vector<HS_vertex_move_type>{ANNIHILATION, STATIC, CREATION};
}

int HS_vertex_move_domain::to_coordinate(element_type vertex_move) {
  switch (vertex_move) {
    case ANNIHILATION:
      return 0;
      break;

    case STATIC:
      return 1;
      break;

    case CREATION:
      return 2;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca
