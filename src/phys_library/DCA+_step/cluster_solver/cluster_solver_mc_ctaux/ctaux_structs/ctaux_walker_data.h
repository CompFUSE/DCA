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

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_STRUCTS_CTAUX_WALKER_DATA_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_STRUCTS_CTAUX_WALKER_DATA_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/mc_walker_data.hpp"
#include <utility>
#include "comp_library/linalg/linalg.hpp"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <dca::linalg::DeviceType device_t, typename parameters_type>
class MC_walker_data<CT_AUX_SOLVER, device_t, parameters_type> {
protected:
  const static int MAX_VERTEX_SINGLETS = 4;

public:
  MC_walker_data(parameters_type& parameters, int id);

public:
  int thread_id;

  dca::linalg::Matrix<double, device_t> N_up;
  dca::linalg::Matrix<double, device_t> N_dn;

  dca::linalg::Matrix<double, device_t> G0_up;
  dca::linalg::Matrix<double, device_t> G0_dn;

  dca::linalg::Matrix<double, device_t> Gamma_up;
  dca::linalg::Matrix<double, device_t> Gamma_dn;

  dca::linalg::Matrix<double, device_t> stored_Gamma_up;
  dca::linalg::Matrix<double, device_t> stored_Gamma_dn;

  dca::linalg::Matrix<double, device_t> G_up;
  dca::linalg::Matrix<double, device_t> G_dn;
};

template <dca::linalg::DeviceType device_t, typename parameters_type>
MC_walker_data<CT_AUX_SOLVER, device_t, parameters_type>::MC_walker_data(parameters_type& parameters,
                                                                         int id)
    : thread_id(id),

      N_up("N_up", 0, parameters.get_initial_matrix_size()),
      N_dn("N_dn", 0, parameters.get_initial_matrix_size()),

      G0_up("G0_up", 0, parameters.get_initial_matrix_size()),
      G0_dn("G0_dn", 0, parameters.get_initial_matrix_size()),

      Gamma_up("Gamma_up", 0, MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),
      Gamma_dn("Gamma_dn", 0, MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),

      stored_Gamma_up("stored_Gamma_up", 0, MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),
      stored_Gamma_dn("stored_Gamma_dn", 0, MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),

      G_up("G_up", std::pair<int, int>(0, 0),
           std::pair<int, int>(parameters.get_initial_matrix_size(),
                               MAX_VERTEX_SINGLETS * parameters.get_K_PHANI())),
      G_dn("G_dn", std::pair<int, int>(0, 0),
           std::pair<int, int>(parameters.get_initial_matrix_size(),
                               MAX_VERTEX_SINGLETS * parameters.get_K_PHANI())) {
  N_up.setThreadAndStreamId(thread_id, 0);
  N_dn.setThreadAndStreamId(thread_id, 0);

  G0_up.setThreadAndStreamId(thread_id, 0);
  G0_dn.setThreadAndStreamId(thread_id, 0);

  Gamma_up.setThreadAndStreamId(thread_id, 0);
  Gamma_dn.setThreadAndStreamId(thread_id, 0);

  stored_Gamma_up.setThreadAndStreamId(thread_id, 0);
  stored_Gamma_dn.setThreadAndStreamId(thread_id, 0);

  G_up.setThreadAndStreamId(thread_id, 0);
  G_dn.setThreadAndStreamId(thread_id, 0);
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_STRUCTS_CTAUX_WALKER_DATA_H
