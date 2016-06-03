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

#ifndef PHYS_LIBRARY_PARAMETERS_MODELS_MATERIAL_HAMILTONIANS_MATERIAL_LATTICE_TEMPLATE_H
#define PHYS_LIBRARY_PARAMETERS_MODELS_MATERIAL_HAMILTONIANS_MATERIAL_LATTICE_TEMPLATE_H

#include <complex>
#include <utility>
#include <vector>

#include "comp_library/function_library/include_function_library.h"
#include "phys_library/parameters/models/material_hamiltonians/material_names.hpp"
#include "phys_library/domains/cluster/symmetries/point_groups/No_symmetry.h"

template <material_name_type name, typename point_group_type>
class material_lattice {
public:
  const static int DIMENSION = -1;
  const static int BANDS = -1;

  typedef no_symmetry<DIMENSION> LDA_point_group;
  typedef point_group_type DCA_point_group;

  static double* initialize_r_DCA_basis();
  static double* initialize_r_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> get_orbital_permutations();

  template <class domain, class parameters_type>
  static void initialize_H_interaction(FUNC_LIB::function<double, domain>& H_interaction,
                                       parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(FUNC_LIB::function<int, domain>& H_symmetry);

  template <class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters, std::vector<double> k,
                                                   int b1, int s1, int b2, int s2);
};

#endif  // PHYS_LIBRARY_PARAMETERS_MODELS_MATERIAL_HAMILTONIANS_MATERIAL_LATTICE_TEMPLATE_H
