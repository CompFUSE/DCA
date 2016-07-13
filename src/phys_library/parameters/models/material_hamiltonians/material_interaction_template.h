// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef PHYS_LIBRARY_PARAMETERS_MODELS_MATERIAL_HAMILTONIANS_MATERIAL_INTERACTION_TEMPLATE_H
#define PHYS_LIBRARY_PARAMETERS_MODELS_MATERIAL_HAMILTONIANS_MATERIAL_INTERACTION_TEMPLATE_H

#include "phys_library/parameters/models/material_hamiltonians/material_names.hpp"
#include "comp_library/function_library/include_function_library.h"

template <material_name_type name, typename point_group_type>
class material_interaction {
public:
  template <class domain, class parameters_type>
  static void initialize_H_interaction(FUNC_LIB::function<double, domain>& H_interaction,
                                       parameters_type& parameters);
};

#endif  // PHYS_LIBRARY_PARAMETERS_MODELS_MATERIAL_HAMILTONIANS_MATERIAL_INTERACTION_TEMPLATE_H
