// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Jérémie Bouquet (bouquetj@gmail.com)
//
//

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_HPP

#include <array>
#include <cassert>
#include <type_traits>

#include "dca/linalg/linalg.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/device_memory/global_memory_manager.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <linalg::DeviceType device_t>
class DMatrixBuilder;

template <linalg::DeviceType device_t>
using Matrix = linalg::Matrix<double, device_t>;

template <>
class DMatrixBuilder<linalg::CPU> {
private:
  using Matrix = linalg::Matrix<double, linalg::CPU>;
  using MatrixPair = std::array<linalg::Matrix<double, linalg::CPU>, 2>;

public:
  DMatrixBuilder(const G0Interpolation<linalg::CPU>& g0,
                 const linalg::Matrix<int, linalg::CPU>& site_diff,
                 const std::vector<int>& sbdm_step, const std::array<double, 3>& alphas);

  void buildSQR(MatrixPair& S, MatrixPair& Q, MatrixPair& R,
                const SolverConfiguration& config) const;

  const G0Interpolation<linalg::CPU>& getG0() const {
    return g0_ref_;
  }

  double computeD(const int i, const int j, const Sector& config) const;
  double computeAlpha(const int aux_spin_type) const;
  double computeDSubmatrix(const int i, const int j, const Sector& configuration) const;
  double computeF(const double alpha) const;
  double computeF(const int i, const Sector& configuration) const;
  double computeF(const int aux_spin_type) const;
  double computeG(const int i, const int j, const Sector& configuration, const Matrix& M) const;
  double computeGFast(const int i, const int j, const int aux_spin_type, const double M_ij) const;
  double computeG0(const int i, const int j, const Sector& configuration) const;
  double computeGamma(const int aux_spin_type, const int new_aux_spin_type) const;
  void computeG0Init(Matrix& G0, const Sector& configuration, const int n_init, const int n_max) const;
  void computeG0(linalg::Matrix<double, linalg::CPU>& G0, const Sector& configuration,
                 const int n_init, const int n_max, const int which_section) const;

private:
  int label(const int nu1, const int nu2, const int r) const;

protected:
  const G0Interpolation<linalg::CPU>& g0_ref_;
  const double alpha_1_ = 0;
  const double alpha_2_ = 0;
  const double alpha_3_ = 0;
  const int n_bands_ = -1;
  const std::vector<int>& sbdm_step_;
  const linalg::Matrix<int, linalg::CPU>& site_diff_;
};

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_HPP