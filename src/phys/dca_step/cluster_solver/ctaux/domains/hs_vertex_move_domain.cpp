// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
  static std::vector<HS_vertex_move_type> v(0);

  v.push_back(ANNIHILATION);
  v.push_back(STATIC);
  v.push_back(CREATION);

  return v;
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

}  // ctaux
}  // solver
}  // phys
}  // dca
