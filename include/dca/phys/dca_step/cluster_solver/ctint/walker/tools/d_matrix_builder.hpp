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
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"

#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu.h"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#include "dca/linalg/util/gpu_stream.hpp"

#endif  // DCA_HAVE_GPU

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <linalg::DeviceType device_t, typename Scalar>
class DMatrixBuilder;

template <typename Scalar>
class DMatrixBuilder<linalg::CPU, Scalar> {
private:
  using Matrix = linalg::Matrix<Scalar, linalg::CPU>;
  using MatrixPair = std::array<Matrix, 2>;
  using Real = dca::util::RealAlias<Scalar>;
  using GpuStream = dca::linalg::util::GpuStream;
public:
  template <class RDmn>
  DMatrixBuilder(const G0Interpolation<linalg::CPU, Scalar>& g0, int nb, const RDmn& /*r_dmn*/);

  DMatrixBuilder(const G0Interpolation<linalg::CPU, Scalar>& g0,
                 const linalg::Matrix<int, linalg::CPU>& site_diff, int nb);

  virtual ~DMatrixBuilder() {}

  // Set the parameters of the alpha field (see next function).
  void setAlphas(const std::array<double, 3>& alphas_base, bool adjust_dd);

  // Computes the right (Q), bottom(R), and right-bottom (S) blocks of the D matrix (M^-1) with elements
  // D_i,j = G0(v_i, v_j) + delta(i, j) alpha_field(v_i), where v_i is an interaction vertex.
  void buildSQR(MatrixPair& S, MatrixPair& Q, MatrixPair& R, const SolverConfiguration& config) const;

  const G0Interpolation<linalg::CPU, Scalar>& getG0() const {
    return g0_ref_;
  }

  // Computes the i, j entry of the D matrix.
  Scalar computeD(const int i, const int j, const Sector& config) const;

  // Computes the alpha field given auxiliary spin and band.
  auto computeAlpha(const int aux_spin_type, const int b) const -> Real;

  // Computes the quantity f = alpha_field / (1 - alpha_field). Used by the submatrix solver.
  auto computeF(const Real alpha) const -> Real;
  auto computeF(const int aux_spin_type, int b) const -> Real;

  // computes the gamma factor: gamma = (f(new_spin) - f(old_spin)) / f(old_spin)
  auto computeGamma(int aux_spin_type, int new_aux_spin_type, int b) const -> Real;

  // Computes a sector of the G0 matrix.
  // If which_section == 0, computes the bottom sector. If which_sector == 1 computes the right sector.
  void computeG0(Matrix& G0, const Sector& configuration, const int n_init, const int n_max,
                 const int which_section) const;

#ifdef DCA_HAVE_GPU
  virtual void computeG0(linalg::Matrix<Scalar, linalg::GPU>& /*G0*/,
                         const details::DeviceConfiguration& /*configuration*/, int /*n_init*/,
                         bool /*right_section*/, const GpuStream& /*stream*/) const {
    throw(std::runtime_error("Not implemented."));
  }
#endif  // DCA_HAVE_GPU

private:
  int label(const int nu1, const int nu2, const int r) const;

protected:
  const G0Interpolation<linalg::CPU, Scalar>& g0_ref_;

  std::vector<Real> alpha_dd_;  // alpha parameters for density-density interactions. One for each band.
  Real alpha_dd_neg_ = 0;  // parameter for density-density interactions with negative coupling.
  Real alpha_ndd_ = 0;     // parameter for non-density-density interactions.

  const int n_bands_ = -1;
  std::array<int, 2> sbdm_step_;
  // Note: site_diff is a matrix where site_diff(i,j) = r_j - r_i.
  const linalg::Matrix<int, linalg::CPU>& site_diff_;
};

template <typename Scalar>
template <class RDmn>
DMatrixBuilder<linalg::CPU, Scalar>::DMatrixBuilder(const G0Interpolation<linalg::CPU, Scalar>& g0,
                                                  int nb, const RDmn& /*r_dmn*/)
    : DMatrixBuilder(g0, RDmn::parameter_type::get_subtract_matrix(), nb) {}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_HPP
