// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch).

// This class organizes the MC walker in the CT-INT QMC using a GPU accelerator.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_SUBMATRIX_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_SUBMATRIX_HPP

#include "dca/config/haves_defines.hpp"

#ifndef DCA_HAVE_GPU
#error "This file requires GPU support."
#endif

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"

#include <cassert>
#include <stdexcept>
#include <vector>

#include "dca/linalg/util/gpu_event.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_submatrix_base.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration_manager.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/kernels_interface.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {

template <class Parameters, DistType DIST>
class CtintWalkerSubmatrixGpu : public CtintWalkerSubmatrixBase<Parameters, DIST> {
public:
  using this_type = CtintWalkerSubmatrixGpu<Parameters, DIST>;
  using SubmatrixBase = CtintWalkerSubmatrixBase<Parameters, DIST>;
  using BaseClass = CtintWalkerBase<Parameters, DIST>;

  using typename BaseClass::Data;
  using typename BaseClass::Profiler;
  using typename BaseClass::Rng;
  using typename BaseClass::GpuStream;
  using typename BaseClass::Real;
  using typename BaseClass::Scalar;

  using Resource = DMatrixBuilder<linalg::GPU, Scalar>;

  template <linalg::DeviceType dev>
  using MatrixView = linalg::MatrixView<Scalar, dev>;
  template <linalg::DeviceType dev>
  using Matrix = linalg::Matrix<Scalar, dev>;
  template <linalg::DeviceType device_t>
  using MatrixPair = std::array<linalg::Matrix<Scalar, device_t>, 2>;

  constexpr static linalg::DeviceType device = linalg::GPU;

  CtintWalkerSubmatrixGpu(const Parameters& pars_ref, const Data& /*data_ref*/, Rng& rng_ref,
                          DMatrixBuilder<linalg::GPU, Scalar>& d_matrix_builder, int id = 0);

  void computeM(MatrixPair<linalg::GPU>& m_accum);

  MatrixPair<linalg::CPU> getM();

  void doSweep() override;

  void synchronize() const override;

  using BaseClass::order;
  using BaseClass::get_stream;

  std::size_t deviceFingerprint() const override;

  void setMFromConfig() override;
  void uploadConfiguration();

protected:
  // For testing purposes:
  void doStep(int n_moves_to_delay);

private:
  void doStep() override;

  void computeMInit();
  void computeGInit();
  void updateM() override;

protected:
  using BaseClass::configuration_;
  using BaseClass::M_;
  using BaseClass::sweeps_per_meas_;
  using BaseClass::partial_order_avg_;
  using BaseClass::thermalization_steps_;
  using BaseClass::order_avg_;
  using BaseClass::sign_avg_;
  using BaseClass::n_steps_;
  using BaseClass::mc_log_weight_;
  using BaseClass::n_accepted_;

private:
  using BaseClass::concurrency_;
  using BaseClass::thread_id_;
  using SubmatrixBase::removal_list_;
  using SubmatrixBase::source_list_;
  using SubmatrixBase::conf_removal_list_;
  using SubmatrixBase::parameters_;
  using SubmatrixBase::rng_;
  using SubmatrixBase::Gamma_inv_;
  using SubmatrixBase::gamma_values_;
  using SubmatrixBase::f_;
  // Initial and current sector sizes.
  using SubmatrixBase::n_init_;
  // Maximal sector size after submatrix update.
  using SubmatrixBase::n_max_;
  using SubmatrixBase::n_bands_;
  using SubmatrixBase::gamma_;
  using SubmatrixBase::G_;
  using SubmatrixBase::move_indices_;
  using SubmatrixBase::flop_;
  using SubmatrixBase::prob_const_;
  DeviceConfigurationManager device_config_;

  std::array<linalg::Vector<int, linalg::GPU>, 2> removal_list_dev_;
  std::array<linalg::Vector<int, linalg::GPU>, 2> source_list_dev_;

  MatrixPair<linalg::GPU> M_dev_;
  MatrixPair<linalg::GPU> Gamma_inv_dev_;
  MatrixPair<linalg::GPU> D_dev_;
  MatrixPair<linalg::GPU> G_dev_;
  MatrixPair<linalg::GPU> G0_dev_;

  std::array<linalg::Vector<Real, linalg::GPU>, 2> f_dev_;
  std::array<linalg::Vector<Real, linalg::CPU>, 2> f_values_;

  std::array<linalg::Vector<int, linalg::GPU>, 2> move_indices_dev_;
  std::array<linalg::util::HostVector<std::pair<int, Real>>, 2> gamma_index_;
  std::array<linalg::Vector<std::pair<int, Real>, linalg::GPU>, 2> gamma_index_dev_;

  std::array<linalg::util::GpuEvent, 2> config_copied_;

  DMatrixBuilder<linalg::GPU, Scalar>& d_matrix_builder_;
};

template <class Parameters, DistType DIST>
CtintWalkerSubmatrixGpu<Parameters, DIST>::CtintWalkerSubmatrixGpu(
    const Parameters& pars_ref, const Data& data, Rng& rng_ref,
    DMatrixBuilder<linalg::GPU, Scalar>& d_matrix_builder, int id)
    : SubmatrixBase(pars_ref, data, rng_ref, id), d_matrix_builder_(d_matrix_builder) {
  for (int b = 0; b < n_bands_; ++b) {
    for (int i = 1; i <= 3; ++i) {
      f_[i][b] = d_matrix_builder_.computeF(i, b);
      f_[-i][b] = d_matrix_builder_.computeF(-i, b);

      gamma_values_[std::make_pair(0, i)][b] = d_matrix_builder_.computeGamma(0, i, b);
      gamma_values_[std::make_pair(0, -i)][b] = d_matrix_builder_.computeGamma(0, -i, b);
      gamma_values_[std::make_pair(i, 0)][b] = d_matrix_builder_.computeGamma(i, 0, b);
      gamma_values_[std::make_pair(-i, 0)][b] = d_matrix_builder_.computeGamma(-i, 0, b);

      prob_const_[i][b] = prob_const_[-i][b] = -1. / (f_[i][b] - 1) / (f_[-i][b] - 1);
    }
    f_[0][b] = 1;
  }

  if (concurrency_.id() == concurrency_.first() && thread_id_ == 0)
    std::cout << "\nCT-INT submatrix walker extended to GPU." << std::endl;
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixGpu<Parameters, DIST>::setMFromConfig() {
  BaseClass::setMFromConfigImpl(d_matrix_builder_);
  for (int s = 0; s < 2; ++s) {
    M_dev_[s].setAsync(M_[s], get_stream(s));
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixGpu<Parameters, DIST>::synchronize() const {
  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  checkRC(cudaStreamSynchronize(get_stream(0)));
  checkRC(cudaStreamSynchronize(get_stream(1)));
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixGpu<Parameters, DIST>::doSweep() {
  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  SubmatrixBase::doSteps();
  uploadConfiguration();
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixGpu<Parameters, DIST>::doStep(const int n_moves_to_delay) {
  SubmatrixBase::nbr_of_moves_to_delay_ = n_moves_to_delay;
  doStep();
  uploadConfiguration();
  doStep();
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixGpu<Parameters, DIST>::doStep() {
  SubmatrixBase::generateDelayedMoves(SubmatrixBase::nbr_of_moves_to_delay_);
  uploadConfiguration();

  computeMInit();
  computeGInit();
  synchronize();
  SubmatrixBase::mainSubmatrixProcess();
  updateM();
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixGpu<Parameters, DIST>::uploadConfiguration() {
  for (int s = 0; s < 2; ++s)
    config_copied_[s].block();

  // Upload configuration and f values.
  device_config_.upload(configuration_, thread_id_);
  for (int s = 0; s < 2; ++s) {
    auto& values = f_values_[s];
    const auto& sector = configuration_.getSector(s);
    values.resizeNoCopy(sector.size());
    for (int i = 0; i < sector.size(); ++i) {
      const auto field_type = sector.getAuxFieldType(i);
      const auto b = sector.getLeftB(i);
      values[i] = f_[field_type][b];
    }
    f_dev_[s].setAsync(values, get_stream(s));
  }

  for (int s = 0; s < 2; ++s)
    config_copied_[s].record(get_stream(s));
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixGpu<Parameters, DIST>::computeMInit() {
  //  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s)
    M_dev_[s].resize(n_max_[s]);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];
    if (delta > 0) {
      D_dev_[s].resizeNoCopy(std::make_pair(delta, n_init_[s]));

      if (delta == 0 || n_init_[s] == 0)
        throw std::runtime_error(
            "expansion factor dropped to 0 or below use a higher beta or larger interaction!");

      d_matrix_builder_.computeG0(D_dev_[s], device_config_.getDeviceData(s), n_init_[s], false,
                                  get_stream(s));
      flop_ += D_dev_[s].nrCols() * D_dev_[s].nrRows() * 10;

      MatrixView<linalg::GPU> D_view(D_dev_[s]);
      details::multiplyByFColFactor(D_view, f_dev_[s].ptr(), get_stream(s));

      MatrixView<linalg::GPU> M(M_dev_[s], 0, 0, n_init_[s], n_init_[s]);
      MatrixView<linalg::GPU> D_M(M_dev_[s], n_init_[s], 0, delta, n_init_[s]);

      linalg::matrixop::gemm(D_dev_[s], M, D_M, thread_id_, s);
      flop_ += 2 * D_dev_[s].nrRows() * D_dev_[s].nrCols() * M.nrCols();

      details::setRightSectorToId(M_dev_[s].ptr(), M_dev_[s].leadingDimension(), n_init_[s],
                                  n_max_[s], get_stream(s));
    }
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixGpu<Parameters, DIST>::computeGInit() {
  //  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];
    auto& f_dev = f_dev_[s];

    G_dev_[s].resizeNoCopy(n_max_[s]);

    MatrixView<linalg::GPU> G(G_dev_[s]);
    const MatrixView<linalg::GPU> M(M_dev_[s]);
    details::computeGLeft(G, M, f_dev.ptr(), n_init_[s], get_stream(s));
    flop_ += n_init_[s] * n_max_[s] * 2;

    if (delta > 0) {
      G0_dev_[s].resizeNoCopy(std::make_pair(n_max_[s], delta));
      // This doesn't do flops but the g0 interp data it uses does somewhere.
      d_matrix_builder_.computeG0(G0_dev_[s], device_config_.getDeviceData(s), n_init_[s], true,
                                  get_stream(s));
      flop_ += G0_dev_[s].nrCols() * G0_dev_[s].nrRows() * 10;
      MatrixView<linalg::GPU> G(G_dev_[s], 0, n_init_[s], n_max_[s], delta);
      // compute G right.
      linalg::matrixop::gemm(M_dev_[s], G0_dev_[s], G, thread_id_, s);
      flop_ += 2 * M_dev_[s].nrRows() * M_dev_[s].nrCols() * G0_dev_[s].nrCols();
    }
    G_[s].setAsync(G_dev_[s], get_stream(s));
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixGpu<Parameters, DIST>::updateM() {
  //  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s)
    Gamma_inv_dev_[s].setAsync(Gamma_inv_[s], get_stream(s));

  // Copy gamma factors.
  for (int s = 0; s < 2; ++s) {
    auto& gamma_index = gamma_index_[s];
    const int n = gamma_[s].size();
    gamma_index.resize(n);
    for (int i = 0; i < n; ++i)
      gamma_index[i] = std::make_pair(move_indices_[s][i], gamma_[s][i]);
    gamma_index_dev_[s].setAsync(gamma_index, get_stream(s));
  }

  for (int s = 0; s < 2; ++s) {
    if (gamma_[s].size() == 0)
      continue;

    // Reuse previously allocated memory as workspace.
    auto& old_M = D_dev_[s];
    auto& old_G = G0_dev_[s];

    const int gamma_size = gamma_[s].size();
    old_G.resizeNoCopy(std::make_pair(n_max_[s], gamma_size));
    old_M.resizeNoCopy(std::make_pair(gamma_size, n_max_[s]));

    assert(dca::linalg::util::getStream(thread_id_, s) == get_stream(s));
    move_indices_dev_[s].setAsync(move_indices_[s], get_stream(s));
    // Note: an event synchronization might be necessary if the order of operation is changed.
    linalg::matrixop::copyCols(G_dev_[s], move_indices_dev_[s], old_G, thread_id_, s);
    linalg::matrixop::copyRows(M_dev_[s], move_indices_dev_[s], old_M, thread_id_, s);

    auto& tmp = G_dev_[s];
    // Note: the following resize is safe as it does not deallocate.
    tmp.resizeNoCopy(std::make_pair(gamma_size, n_max_[s]));
    linalg::matrixop::gemm(Gamma_inv_dev_[s], old_M, tmp, thread_id_, s);
    linalg::matrixop::gemm(Scalar(-1.), old_G, tmp, Scalar(1.), M_dev_[s], thread_id_, s);
    flop_ += 2 * Gamma_inv_dev_[s].nrRows() * Gamma_inv_dev_[s].nrCols() * old_M.nrCols();
    flop_ += 2 * old_G.nrRows() * old_G.nrCols() * tmp.nrCols();

    details::divideByGammaFactor(MatrixView<linalg::GPU>(M_dev_[s]), gamma_index_dev_[s].ptr(),
                                 gamma_size, get_stream(s));
    flop_ += gamma_size * M_dev_[s].nrCols() * 2;
  }

  // Remove non-interacting rows and columns.
  configuration_.moveAndShrink(source_list_, removal_list_, conf_removal_list_);
  for (int s = 0; s < 2; ++s) {
    removal_list_dev_[s].setAsync(removal_list_[s], get_stream(s));
    source_list_dev_[s].setAsync(source_list_[s], get_stream(s));
    config_copied_[s].record(get_stream(s));
    linalg::matrixop::copyRows(M_dev_[s], source_list_dev_[s], M_dev_[s], removal_list_dev_[s],
                               thread_id_, s);
    linalg::matrixop::copyCols(M_dev_[s], source_list_dev_[s], M_dev_[s], removal_list_dev_[s],
                               thread_id_, s);
    M_dev_[s].resize(configuration_.getSector(s).size());
  }

  assert(configuration_.getSector(0).size() == M_dev_[0].nrRows());
  assert(configuration_.getSector(1).size() == M_dev_[1].nrRows());
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixGpu<Parameters, DIST>::computeM(
    std::array<dca::linalg::Matrix<Scalar, linalg::GPU>, 2>& m_accum) {
  for (int s = 0; s < 2; ++s)
    m_accum[s].resizeNoCopy(M_dev_[s].size());

  for (int s = 0; s < 2; ++s) {
    MatrixView<linalg::GPU> m_in(M_dev_[s]);
    MatrixView<linalg::GPU> m_out(m_accum[s]);
    details::multiplyByInverseFFactor(m_in, m_out, f_dev_[s].ptr(), get_stream(s));
    flop_ += m_in.nrRows() * m_out.nrCols();
  }

  // TODO: understand why this is necessary.
  config_copied_[0].block();
  config_copied_[1].block();
}

template <class Parameters, DistType DIST>
std::size_t CtintWalkerSubmatrixGpu<Parameters, DIST>::deviceFingerprint() const {
  std::size_t res = 0;
  for (int s = 0; s < 2; ++s) {
    res += M_dev_[s].deviceFingerprint();
    res += Gamma_inv_dev_[s].deviceFingerprint();
    res += D_dev_[s].deviceFingerprint();
    res += G_dev_[s].deviceFingerprint();
    res += G0_dev_[s].deviceFingerprint();
  }
  return res;
}

template <class Parameters, DistType DIST>
CtintWalkerSubmatrixGpu<Parameters, DIST>::MatrixPair<linalg::CPU> CtintWalkerSubmatrixGpu<
    Parameters, DIST>::getM() {
  std::array<dca::linalg::Matrix<Scalar, device>, 2> M;

  computeM(M);
  checkRC(cudaDeviceSynchronize());

  std::array<dca::linalg::Matrix<Scalar, linalg::CPU>, 2> M_copy{M[0].size(), M[1].size()};
  linalg::util::memoryCopyD2H(M_copy[0].ptr(), M_copy[0].leadingDimension(), M[0].ptr(), M[0].leadingDimension(),
                M[0].size());
  linalg::util::memoryCopyD2H(M_copy[1].ptr(), M_copy[1].leadingDimension(), M[1].ptr(), M[1].leadingDimension(),
                M[1].size());
  return M_copy;
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_SUBMATRIX_HPP
