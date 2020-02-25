// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//          Jérémie Bouquet (bouquetj@gmail.com).
//
// CPU implementation of d_matrix_builder.hpp

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

using linalg::CPU;
using linalg::GPU;

DMatrixBuilder<CPU>::DMatrixBuilder(const G0Interpolation<CPU>& g0,
                                    const linalg::Matrix<int, linalg::CPU>& site_diff, const int nb)
    : g0_ref_(g0), n_bands_(nb), sbdm_step_{nb, nb * nb}, site_diff_(site_diff) {}

void DMatrixBuilder<CPU>::setAlphas(const std::array<double, 3>& alphas_base, bool adjust_dd) {
  alpha_dd_neg_ = alphas_base[1];
  alpha_ndd_ = alphas_base[2];

  const int r0 = site_diff_(0, 0);
  alpha_dd_.resize(n_bands_);
  for (int b = 0; b < n_bands_; ++b) {
    alpha_dd_[b] = alphas_base[0];
    // TODO: check if shifting by g0(0+) or g0(0-).
    if (adjust_dd)  // shift by g0(0).
      alpha_dd_[b] += std::abs(g0_ref_(0., label(b, b, r0)));
  }
}

void DMatrixBuilder<CPU>::buildSQR(MatrixPair& S, MatrixPair& Q, MatrixPair& R,
                                   const SolverConfiguration& config) const {
  std::array<int, 2> size_increase = config.sizeIncrease();

  for (int s = 0; s < 2; ++s) {
    const int delta = size_increase[s];
    const Sector& sector = config.getSector(s);
    const int n = sector.size() - delta;

    Q[s].resizeNoCopy(std::make_pair(n, delta));
    R[s].resizeNoCopy(std::make_pair(delta, n));
    S[s].resizeNoCopy(std::make_pair(delta, delta));

    for (int j = 0; j < delta; j++)
      for (int i = 0; i < n; i++)
        Q[s](i, j) = computeD(i, n + j, sector);
    for (int j = 0; j < n; j++)
      for (int i = 0; i < delta; i++)
        R[s](i, j) = computeD(i + n, j, sector);
    for (int j = 0; j < delta; j++)
      for (int i = 0; i < delta; i++)
        S[s](i, j) = computeD(n + i, n + j, sector);
  }
}

double DMatrixBuilder<CPU>::computeD(const int i, const int j, const Sector& configuration) const {
  assert(configuration.size() > i and configuration.size() > j);

  const int b1 = configuration.getLeftB(i);
  const int b2 = configuration.getRightB(j);
  const int delta_r = site_diff_(configuration.getRightR(j), configuration.getLeftR(i));
  const int p_index = label(b1, b2, delta_r);
  const double delta_tau = configuration.getTau(i) - configuration.getTau(j);
  const double g0_val = g0_ref_(delta_tau, p_index);
  if (i == j)
    return g0_val - computeAlpha(configuration.getAuxFieldType(i), b1);
  else
    return g0_val;
}

int DMatrixBuilder<CPU>::label(const int b1, const int b2, const int r) const {
  return b1 + b2 * sbdm_step_[0] + r * sbdm_step_[1];
}

double DMatrixBuilder<CPU>::computeAlpha(const int aux_spin_type, const int b) const {
  assert(alpha_dd_.size());
  switch (std::abs(aux_spin_type)) {
    case 1:
      return aux_spin_type < 0 ? (0.5 + alpha_dd_[b]) : (0.5 - alpha_dd_[b]);
    case 2:
      return aux_spin_type < 0 ? (0.5 + alpha_dd_neg_) : (0.5 - alpha_dd_neg_);
    case 3:
      return aux_spin_type < 0 ? alpha_ndd_ : -alpha_ndd_;
    default:
      throw(std::logic_error("type not recognized."));
  }
}

double DMatrixBuilder<CPU>::computeF(const double alpha) const {
  return alpha / (alpha - 1);
}

double DMatrixBuilder<CPU>::computeF(const int aux_spin_type, int b) const {
  if (aux_spin_type == 0)
    return 1;
  else
    return computeF(computeAlpha(aux_spin_type, b));
}

double DMatrixBuilder<CPU>::computeGamma(const int aux_spin_type, const int new_aux_spin_type,
                                         int b) const {
  return (computeF(new_aux_spin_type, b) - computeF(aux_spin_type, b)) / computeF(aux_spin_type, b);
}

// Compute only the parts of G0 required at a given moment. (Re)Computing every element is not needed in most situations.
void DMatrixBuilder<CPU>::computeG0(Matrix& G0, const Sector& configuration, const int n_init,
                                    const int n_max, const int which_section) const {
  int b_i, b_j, r_i, r_j;
  double tau_i, tau_j;

  if (which_section == 0) {
    for (int i = n_init; i < n_max; ++i) {
      b_i = configuration.getLeftB(i);
      tau_i = configuration.getTau(i);
      r_i = configuration.getLeftR(i);

      for (int j = 0; j < n_init; ++j) {
        b_j = configuration.getRightB(j);
        tau_j = configuration.getTau(j);
        r_j = configuration.getRightR(j);

        G0(i - n_init, j) = g0_ref_(tau_i - tau_j, label(b_i, b_j, site_diff_(r_j, r_i)));
      }
    }
  }

  else {
    for (int i = 0; i < n_max; ++i) {
      b_i = configuration.getLeftB(i);
      tau_i = configuration.getTau(i);
      r_i = configuration.getLeftR(i);

      for (int j = n_init; j < n_max; ++j) {
        b_j = configuration.getRightB(j);
        tau_j = configuration.getTau(j);
        r_j = configuration.getRightR(j);

        G0(i, j - n_init) = g0_ref_(tau_i - tau_j, label(b_i, b_j, site_diff_(r_j, r_i)));
      }
    }
  }
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
