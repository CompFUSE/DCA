// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//          Jérémie Bouquet (bouquetj@gmail.com)
//
// Class to compute new entries of the D matrix.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_HPP

#include <array>
#include <cassert>
#include <type_traits>

#include "dca/linalg/linalg.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/solver_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"

#ifdef DCA_HAVE_CUDA
#include <cuda.h>
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <linalg::DeviceType device_t, typename Real>
class DMatrixBuilder;

template <typename Real>
class DMatrixBuilder<linalg::CPU, Real> {
private:
  using Matrix = linalg::Matrix<Real, linalg::CPU>;
  using MatrixPair = std::array<Matrix, 2>;

public:
  template <class RDmn>
  DMatrixBuilder(const G0Interpolation<linalg::CPU, Real>& g0, int nb, const RDmn& /*r_dmn*/);

  DMatrixBuilder(const G0Interpolation<linalg::CPU, Real>& g0,
                 const linalg::Matrix<int, linalg::CPU>& site_diff, int nb);

  virtual ~DMatrixBuilder() {}

  // Set the parameters of the alpha field (see next function).
  void setAlphas(const std::array<double, 3>& alphas_base, bool adjust_dd);

  // Computes the right (Q), bottom(R), and right-bottom (S) blocks of the D matrix (M^-1) with elements
  // D_i,j = G0(v_i, v_j) + delta(i, j) alpha_field(v_i), where v_i is an interaction vertex.
  void buildSQR(MatrixPair& S, MatrixPair& Q, MatrixPair& R, const SolverConfiguration& config) const;

  const G0Interpolation<linalg::CPU, Real>& getG0() const {
    return g0_ref_;
  }

  // Computes the i, j entry of the D matrix.
  Real computeD(const int i, const int j, const Sector& config) const;

  // Computes the alpha field given auxiliary spin and band.
  Real computeAlpha(const int aux_spin_type, const int b) const;

  // Computes the quantity f = alpha_field / (1 - alpha_field). Used by the submatrix solver.
  Real computeF(const Real alpha) const;
  Real computeF(const int aux_spin_type, int b) const;

  // computes the gamma factor: gamma = (f(new_spin) - f(old_spin)) / f(old_spin)
  Real computeGamma(int aux_spin_type, int new_aux_spin_type, int b) const;

  // Computes a sector of the G0 matrix.
  // If which_section == 0, computes the bottom sector. If which_sector == 1 computes the right sector.
  void computeG0(Matrix& G0, const Sector& configuration, const int n_init, const int n_max,
                 const int which_section) const;

#ifdef DCA_HAVE_CUDA
  virtual void computeG0(linalg::Matrix<Real, linalg::GPU>& /*G0*/,
                         const details::DeviceConfiguration& /*configuration*/, int /*n_init*/,
                         bool /*right_section*/, cudaStream_t /*stream*/) const {
    throw(std::runtime_error("Not implemented."));
  }
#endif  // DCA_HAVE_CUDA

private:
  int label(const int nu1, const int nu2, const int r) const;

protected:
  const G0Interpolation<linalg::CPU, Real>& g0_ref_;

  std::vector<Real> alpha_dd_;  // alpha parameters for density-density interactions. One for each band.
  Real alpha_dd_neg_ = 0;  // parameter for density-density interactions with negative coupling.
  Real alpha_ndd_ = 0;     // parameter for non-density-density interactions.

  const int n_bands_ = -1;
  std::array<int, 2> sbdm_step_;
  // Note: site_diff is a matrix where site_diff(i,j) = r_j - r_i.
  const linalg::Matrix<int, linalg::CPU>& site_diff_;
};

template <typename Real>
template <class RDmn>
DMatrixBuilder<linalg::CPU, Real>::DMatrixBuilder(const G0Interpolation<linalg::CPU, Real>& g0,
                                                  int nb, const RDmn& /*r_dmn*/)
    : DMatrixBuilder(g0, RDmn::parameter_type::get_subtract_matrix(), nb) {}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_HPP
