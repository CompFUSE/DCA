// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class stores the data of a CT-AUX walker.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_CTAUX_WALKER_DATA_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_CTAUX_WALKER_DATA_HPP

#include <utility>

#include "dca/linalg/matrix.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <dca::linalg::DeviceType device_t, typename parameters_type>
class CtauxWalkerData {
protected:
  const static int MAX_VERTEX_SINGLETS = 2;

public:
  CtauxWalkerData(parameters_type& parameters, int id);

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
CtauxWalkerData<device_t, parameters_type>::CtauxWalkerData(parameters_type& parameters, int id)
    : thread_id(id),

      N_up("N_up", 0, parameters.get_initial_matrix_size()),
      N_dn("N_dn", 0, parameters.get_initial_matrix_size()),

      G0_up("G0_up", 0, parameters.get_initial_matrix_size()),
      G0_dn("G0_dn", 0, parameters.get_initial_matrix_size()),

      Gamma_up("Gamma_up", 0, MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),
      Gamma_dn("Gamma_dn", 0, MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),

      stored_Gamma_up("stored_Gamma_up", 0,
                      MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),
      stored_Gamma_dn("stored_Gamma_dn", 0,
                      MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),

      G_up("G_up", std::pair<int, int>(0, 0),
           std::pair<int, int>(parameters.get_initial_matrix_size(),
                               MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size())),
      G_dn("G_dn", std::pair<int, int>(0, 0),
           std::pair<int, int>(parameters.get_initial_matrix_size(),
                               MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size())) {}

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_CTAUX_WALKER_DATA_HPP
