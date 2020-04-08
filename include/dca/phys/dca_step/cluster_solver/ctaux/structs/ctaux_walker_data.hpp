// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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

template <dca::linalg::DeviceType device_t, typename Parameters, typename Real>
class CtauxWalkerData {
protected:
  const static int MAX_VERTEX_SINGLETS = 2;

public:
  CtauxWalkerData(Parameters& parameters, int id);

  std::size_t deviceFingerprint() const {
    return N_up.deviceFingerprint() + N_dn.deviceFingerprint() + G0_up.deviceFingerprint() +
           G0_dn.deviceFingerprint() + Gamma_up.deviceFingerprint() + G0_dn.deviceFingerprint() +
           stored_Gamma_up.deviceFingerprint() + stored_Gamma_dn.deviceFingerprint() +
           G_up.deviceFingerprint() + G_dn.deviceFingerprint();
  }

  int thread_id;

  dca::linalg::Matrix<Real, device_t> N_up;
  dca::linalg::Matrix<Real, device_t> N_dn;

  dca::linalg::Matrix<Real, device_t> G0_up;
  dca::linalg::Matrix<Real, device_t> G0_dn;

  dca::linalg::Matrix<Real, device_t> Gamma_up;
  dca::linalg::Matrix<Real, device_t> Gamma_dn;

  dca::linalg::Matrix<Real, device_t> stored_Gamma_up;
  dca::linalg::Matrix<Real, device_t> stored_Gamma_dn;

  dca::linalg::Matrix<Real, device_t> G_up;
  dca::linalg::Matrix<Real, device_t> G_dn;
};

template <dca::linalg::DeviceType device_t, typename Parameters, typename Real>
CtauxWalkerData<device_t, Parameters, Real>::CtauxWalkerData(Parameters& parameters, int id)
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
