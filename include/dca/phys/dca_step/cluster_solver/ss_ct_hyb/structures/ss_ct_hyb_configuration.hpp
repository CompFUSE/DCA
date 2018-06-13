// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Bart Ydens
//
// This class organizes the configuration space in the single-site hybridization QMC.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_STRUCTURES_SS_CT_HYB_CONFIGURATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_STRUCTURES_SS_CT_HYB_CONFIGURATION_HPP

#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/structures/hybridization_vertex.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

class SS_CT_HYB_configuration {
public:
  using this_type = SS_CT_HYB_configuration;
  using orbital_configuration_type = std::vector<Hybridization_vertex>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  SS_CT_HYB_configuration();

  int size();

  void initialize();

  orbital_configuration_type& get_vertices(int i) {
    return vertices(i);
  }

  bool& get_full_line(int i) {
    return has_full_line(i);
  }

  void copy_from(this_type& other_configuration);

  void print();

private:
  func::function<orbital_configuration_type, nu> vertices;
  func::function<bool, nu> has_full_line;

  int N_spin_orbitals;
};

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_STRUCTURES_SS_CT_HYB_CONFIGURATION_HPP
