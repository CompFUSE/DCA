// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class defines common types for the Single-Site Hybridization Monte Carlo Integrator.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_TYPE_DEFINITIONS_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_TYPE_DEFINITIONS_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/mc_type_definitions.hpp"
#include "comp_library/linalg/linalg.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_structures/ss_hybridization_configuration.h"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <class parameters_type, class MOMS_type>
class MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type> {
public:
  // Types that define the profiling.
  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_type;

  // Types that define the scalar type and matrix type.
  typedef double scalartype;
  // typedef resizeable_square_matrix<scalartype> vertex_vertex_matrix_type;
  typedef LIN_ALG::matrix<scalartype, LIN_ALG::CPU> vertex_vertex_matrix_type;

  // Types that define the vertex and configuration type.
  typedef SS_CT_HYB_configuration configuration_type;
  typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;
};

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_TYPE_DEFINITIONS_H
