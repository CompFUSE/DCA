// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"

#include "dca/linalg/make_constant_view.hpp"
#include "dca/linalg/matrix_view.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/device_memory/global_memory_manager.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

using linalg::CPU;
using linalg::GPU;

DMatrixBuilder<CPU>::DMatrixBuilder(const G0Interpolation<CPU>& g0,
                                    const linalg::Matrix<int, linalg::CPU>& site_diff,
                                    const std::vector<int>& sbdm_step,
                                    const std::array<double, 3>& alphas)
    : g0_ref_(g0),
      alpha_1_(alphas[0]),
      alpha_2_(alphas[1]),
      alpha_3_(alphas[2]),
      n_bands_(sbdm_step[1]),
      sbdm_step_(sbdm_step),
      site_diff_(site_diff) {
  assert(sbdm_step.size() == 3);
}

void DMatrixBuilder<CPU>::buildSQR(MatrixPair<linalg::CPU>& S, MatrixPair<linalg::CPU>& Q,
                                   MatrixPair<linalg::CPU>& R,
                                   const SolverConfiguration<linalg::CPU>& config) const {
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

  auto alphaField = [&](const int type) {
    switch (std::abs(type)) {
      case 1:
        return type < 0 ? (0.5 + alpha_1_) : (0.5 - alpha_1_);
      case 2:
        return type < 0 ? (0.5 + alpha_2_) : (0.5 - alpha_2_);
      case 3:
        return type < 0 ? alpha_3_ : -alpha_3_;
      default:
        throw(std::logic_error("type not recognized."));
    }
  };

  const int b1 = configuration.getLeftB(i);
  const int b2 = configuration.getRightB(j);
  const int delta_r = site_diff_(configuration.getLeftR(i), configuration.getRightR(j));
  const int p_index = label(b1, b2, delta_r);
  const double delta_tau = configuration.getTau(i) - configuration.getTau(j);
  const double g0_val = g0_ref_(delta_tau, p_index);
  if (i == j)
    return g0_val - alphaField(configuration.getAuxFieldType(i));
  else
    return g0_val;
}

int DMatrixBuilder<CPU>::label(const int b1, const int b2, const int r) const {
  return b1 + b2 * sbdm_step_[1] + r * sbdm_step_[2];
}

#ifdef DCA_HAVE_CUDA

DMatrixBuilder<linalg::GPU>::DMatrixBuilder(const G0Interpolation<GPU>& g0,
                                            const linalg::Matrix<int, linalg::CPU>& site_diff,
                                            const std::vector<int>& sbdm_step,
                                            const std::array<double, 3>& alpha)
    : g0_ref_(g0), alpha_1_(alpha[0]), alpha_2_(alpha[1]), alpha_3_(alpha[2]), n_bands_(sbdm_step[1]) {
  assert(sbdm_step.size() == 3);

  ctint::GlobalMemoryManager::initializeCluster(*linalg::makeConstantView(site_diff),
                                                sbdm_step);
}

void DMatrixBuilder<linalg::GPU>::buildSQR(MatrixPair<linalg::CPU>& S, MatrixPair<linalg::CPU>& Q,
                                           MatrixPair<linalg::CPU>& R, DeviceWorkspace& devspace,
                                           SolverConfiguration<linalg::GPU>& config, int thread_id,
                                           int stream_id) const {
  std::array<int, 2> size_increase = config.sizeIncrease();

  for (int s = 0; s < 2; ++s) {
    const Sector& sector = config.getSector(s);
    const int delta = size_increase[s];
    const int n = sector.size() - delta;
    const int spin_stream_id = 2 * stream_id + s;
    const auto stream = linalg::util::getStream(thread_id, spin_stream_id);

    config.upload(s, thread_id, spin_stream_id);
    assert(cudaDeviceSynchronize() == cudaSuccess);

    auto& Q_dev = devspace.Q[s];
    auto& R_dev = devspace.R[s];
    auto& S_dev = devspace.S[s];
    Q_dev.resizeNoCopy(std::make_pair(n, delta));
    R_dev.resizeNoCopy(std::make_pair(delta, n));
    S_dev.resizeNoCopy(std::make_pair(delta, delta));

    using MatrixView = linalg::MatrixView<double, linalg::GPU>;
    assert(GlobalMemoryManager::isInitialized());
    details::computeD(MatrixView(Q_dev), MatrixView(R_dev), MatrixView(S_dev), n, delta, alpha_1_,
                      alpha_2_, alpha_3_, config.getDeviceData(s), g0_ref_, stream);

    assert(cudaDeviceSynchronize() == cudaSuccess);

    // Download to host.
    Q[s].setAsync(Q_dev, stream);
    R[s].setAsync(R_dev, stream);
    S[s].setAsync(S_dev, stream);
  }
}

#endif  // DCA_HAVE_CUDA

}  // ctint
}  // solver
}  // phys
}  // dca
