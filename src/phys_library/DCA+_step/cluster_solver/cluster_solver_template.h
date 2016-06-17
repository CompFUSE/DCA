// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This file provides an empty template class for cluster solvers.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_TEMPLATE_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_TEMPLATE_H

#include <string>

#include "comp_library/IO_library/IO.hpp"
#include "comp_library/linalg/linalg_device_types.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_name.hpp"

namespace DCA {

template <ClusterSolverName solver_name, LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
class cluster_solver {
public:
  cluster_solver(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  void initialize(int dca_iteration);

  void execute();

  template <typename dca_info_struct_t>
  void finalize(dca_info_struct_t& dca_info_struct);

  void read(std::string filename);

  void write(std::string filename);

  template <IO::FORMAT DATA_FORMAT>
  void read(IO::reader<DATA_FORMAT>& reader);

  template <IO::FORMAT DATA_FORMAT>
  void write(IO::writer<DATA_FORMAT>& reader);
};

}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_TEMPLATE_H
