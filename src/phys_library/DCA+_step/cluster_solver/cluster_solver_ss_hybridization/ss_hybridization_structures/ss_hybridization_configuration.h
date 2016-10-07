// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Bart Ydens
//
// This class organizes the configuration space in the single-site hybridization QMC.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_STRUCTURES_SS_HYBRIDIZATION_CONFIGURATION_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_STRUCTURES_SS_HYBRIDIZATION_CONFIGURATION_H

#include <iostream>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_structures/ss_hybridization_vertex.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

namespace DCA {

class SS_CT_HYB_configuration {
public:
  using this_type = SS_CT_HYB_configuration;
  using orbital_configuration_type = std::vector<Hybridization_vertex>;

  using b = func::dmn_0<electron_band_domain>;
  using s = func::dmn_0<electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

public:
  SS_CT_HYB_configuration();

  int size();

  void initialize();

  orbital_configuration_type& get_vertices(int i);
  bool& get_full_line(int i);

  void copy_from(this_type& other_configuration);

  void print();

private:
  func::function<orbital_configuration_type, nu> vertices;
  func::function<bool, nu> has_full_line;

  int N_spin_orbitals;
};

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

std::vector<Hybridization_vertex>& SS_CT_HYB_configuration::get_vertices(int i) {
  return vertices(i);
}

bool& SS_CT_HYB_configuration::get_full_line(int i) {
  return has_full_line(i);
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

}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_STRUCTURES_SS_HYBRIDIZATION_CONFIGURATION_H
