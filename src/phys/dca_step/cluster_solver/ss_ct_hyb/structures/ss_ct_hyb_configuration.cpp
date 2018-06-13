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
// This file implements ss_ct_hyb_configuration.hpp.

#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/structures/ss_ct_hyb_configuration.hpp"

#include <iostream>

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

SS_CT_HYB_configuration::SS_CT_HYB_configuration()
    : vertices("SS_CT_HYB_vertices"), has_full_line("SS_CT_HYB_lines") {
  N_spin_orbitals = s::dmn_size() * b::dmn_size();

  for (int i = 0; i < N_spin_orbitals; i++)
    has_full_line(i) = false;
}

int SS_CT_HYB_configuration::size() {
  int size = 0;

  for (int l = 0; l < N_spin_orbitals; l++)
    size += get_vertices(l).size();

  return size;
}

void SS_CT_HYB_configuration::initialize() {
  for (int i = 0; i < has_full_line.size(); i++)
    has_full_line(i) = false;

  for (int i = 0; i < vertices.size(); i++)
    vertices(i).resize(0);
}

void SS_CT_HYB_configuration::copy_from(this_type& other_configuration) {
  for (int l = 0; l < nu::dmn_size(); l++) {
    orbital_configuration_type& other_vertices = other_configuration.get_vertices(l);

    vertices(l).resize(other_vertices.size());

    for (std::size_t i = 0; i < other_vertices.size(); i++) {
      vertices(l)[i].set_t_end(other_vertices[i].t_end());
      vertices(l)[i].set_t_start(other_vertices[i].t_start());
    }
  }

  for (int l = 0; l < nu::dmn_size(); l++)
    has_full_line(l) = other_configuration.get_full_line(l);
}

void SS_CT_HYB_configuration::print() {
  for (int i = 0; i < N_spin_orbitals; i++) {
    std::cout << i;
    if (vertices(i).size() > 0)
      for (size_t j = 0; j < vertices(i).size(); j++)
        std::cout << "\t {" << vertices(i)[j].t_start() << " , " << vertices(i)[j].t_end() << " }";
    else if (has_full_line(i))
      std::cout << "\t full line";
    std::cout << "\n";
  }
  std::cout << "\n";
}

}  // cthyb
}  // solver
}  // phys
}  // dca
