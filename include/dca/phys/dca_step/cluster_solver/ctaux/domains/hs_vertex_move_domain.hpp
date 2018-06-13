// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// HS vertex move domain.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_HS_VERTEX_MOVE_DOMAIN_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_HS_VERTEX_MOVE_DOMAIN_HPP

#include <vector>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

enum HS_vertex_move { ANNIHILATION = -1, STATIC = 0, CREATION = 1 };
typedef HS_vertex_move HS_vertex_move_type;

class HS_vertex_move_domain {
public:
  typedef HS_vertex_move_type element_type;

  static int get_size() {
    return 3;
  }

  static std::vector<HS_vertex_move_type>& get_elements() {
    static std::vector<HS_vertex_move_type> v = initialize_elements();
    return v;
  }

  static int to_coordinate(element_type vertex_move);

private:
  static std::vector<HS_vertex_move_type> initialize_elements();
};

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_DOMAINS_HS_VERTEX_MOVE_DOMAIN_HPP
