// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_DOMAINS_HS_VERTEX_MOVE_DOMAIN_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_DOMAINS_HS_VERTEX_MOVE_DOMAIN_H

#include <stdexcept>
#include <vector>

namespace DCA {
namespace QMCI {
// DCA::QMCI::

enum HS_vertex_move { ANNIHILATION = -1, STATIC = 0, CREATION = 1 };
typedef HS_vertex_move HS_vertex_move_type;

class HS_vertex_move_domain {
public:
  typedef HS_vertex_move_type element_type;

public:
  static int get_size();
  static std::vector<HS_vertex_move_type>& get_elements();

  static int to_coordinate(element_type vertex_move);

private:
  static std::vector<HS_vertex_move_type> initialize_elements();
};

int HS_vertex_move_domain::get_size() {
  return 3;
}

std::vector<HS_vertex_move_type>& HS_vertex_move_domain::get_elements() {
  static std::vector<HS_vertex_move_type> v = initialize_elements();
  return v;
}

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

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_DOMAINS_HS_VERTEX_MOVE_DOMAIN_H
