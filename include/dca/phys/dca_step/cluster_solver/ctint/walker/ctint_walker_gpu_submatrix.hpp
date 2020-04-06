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

#ifndef DCA_HAVE_CUDA
#error "This file requires CUDA support."
#endif

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"

#include <cassert>
#include <stdexcept>
#include <vector>

#include "dca/linalg/util/cuda_event.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration_manager.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/kernels_interface.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {

template <class Parameters, typename Real, bool fix_rng_order>
class CtintWalkerSubmatrixGpu : public CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order> {
public:
  using this_type = CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>;
  using BaseClass = CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>;
  using RootClass = CtintWalkerBase<Parameters, Real>;

  using typename BaseClass::Data;
  using typename BaseClass::Profiler;
  using typename BaseClass::Rng;
  using typename BaseClass::CudaStream;

  template <linalg::DeviceType dev>
  using MatrixView = linalg::MatrixView<Real, dev>;
  template <linalg::DeviceType dev>
  using Matrix = linalg::Matrix<Real, dev>;
  template <linalg::DeviceType device_t>
  using MatrixPair = std::array<linalg::Matrix<Real, device_t>, 2>;

  constexpr static linalg::DeviceType device = linalg::GPU;

  CtintWalkerSubmatrixGpu(const Parameters& pars_ref, const Data& /*data_ref*/, Rng& rng_ref,
                          int id = 0);

  void computeM(MatrixPair<linalg::GPU>& m_accum);

  void doSweep() override;

  void synchronize();

  using BaseClass::order;
  using RootClass::get_stream;

  std::size_t deviceFingerprint() const;

protected:
  // For testing purposes:
  void doStep(int n_moves_to_delay);

  void setMFromConfig() override;

private:
  void doStep() override;

  void computeMInit();
  void computeGInit();
  void updateM();

  void uploadConfiguration();

protected:
  using BaseClass::configuration_;
  using BaseClass::M_;

private:
  using BaseClass::concurrency_;
  using BaseClass::thread_id_;
  using BaseClass::removal_list_;
  using BaseClass::source_list_;
  using BaseClass::conf_removal_list_;
  using RootClass::d_builder_ptr_;
  using BaseClass::parameters_;
  using BaseClass::rng_;
  using BaseClass::Gamma_inv_;
  using BaseClass::f_;
  // Initial and current sector sizes.
  using BaseClass::n_init_;
  // Maximal sector size after submatrix update.
  using BaseClass::n_max_;
  using BaseClass::gamma_;
  using BaseClass::G_;
  using BaseClass::move_indices_;
  using BaseClass::flop_;

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

  std::array<linalg::util::CudaEvent, 2> config_copied_;
};

template <class Parameters, typename Real, bool fix_rng_order>
CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::CtintWalkerSubmatrixGpu(
    const Parameters& pars_ref, const Data& data, Rng& rng_ref, int id)
    : BaseClass(pars_ref, data, rng_ref, id) {
  if (concurrency_.id() == concurrency_.first() && thread_id_ == 0)
    std::cout << "\nCT-INT submatrix walker extended to GPU." << std::endl;
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::setMFromConfig() {
  BaseClass::setMFromConfig();
  for (int s = 0; s < 2; ++s) {
    M_dev_[s].setAsync(M_[s], get_stream(s));
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::synchronize() {
  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  cudaStreamSynchronize(get_stream(0));
  cudaStreamSynchronize(get_stream(1));
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::doSweep() {
  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  BaseClass::doSteps();
  uploadConfiguration();
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::doStep(const int n_moves_to_delay) {
  BaseClass::nbr_of_moves_to_delay_ = n_moves_to_delay;
  doStep();
  uploadConfiguration();
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::doStep() {
  BaseClass::generateDelayedMoves(BaseClass::nbr_of_moves_to_delay_);
  uploadConfiguration();

  computeMInit();
  computeGInit();
  synchronize();
  BaseClass::mainSubmatrixProcess();
  updateM();
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::uploadConfiguration() {
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

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::computeMInit() {
  //  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s)
    M_dev_[s].resize(n_max_[s]);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];
    if (delta > 0) {
      D_dev_[s].resizeNoCopy(std::make_pair(delta, n_init_[s]));
      d_builder_ptr_->computeG0(D_dev_[s], device_config_.getDeviceData(s), n_init_[s], false,
                                get_stream(s));

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

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::computeGInit() {
  //  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];
    auto& f_dev = f_dev_[s];

    G_dev_[s].resizeNoCopy(n_max_[s]);

    MatrixView<linalg::GPU> G(G_dev_[s]);
    const MatrixView<linalg::GPU> M(M_dev_[s]);
    details::computeGLeft(G, M, f_dev.ptr(), n_init_[s], get_stream(s));

    if (delta > 0) {
      G0_dev_[s].resizeNoCopy(std::make_pair(n_max_[s], delta));
      d_builder_ptr_->computeG0(G0_dev_[s], device_config_.getDeviceData(s), n_init_[s], true,
                                get_stream(s));

      MatrixView<linalg::GPU> G(G_dev_[s], 0, n_init_[s], n_max_[s], delta);
      // compute G right.
      linalg::matrixop::gemm(M_dev_[s], G0_dev_[s], G, thread_id_, s);
      flop_ += 2 * M_dev_[s].nrRows() * M_dev_[s].nrCols() * G0_dev_[s].nrCols();
    }
    G_[s].setAsync(G_dev_[s], get_stream(s));
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::updateM() {
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

    move_indices_dev_[s].setAsync(move_indices_[s], get_stream(s));
    // Note: an event synchronization might be necessary if the order of operation is changed.
    linalg::matrixop::copyCols(G_dev_[s], move_indices_dev_[s], old_G, thread_id_, s);
    linalg::matrixop::copyRows(M_dev_[s], move_indices_dev_[s], old_M, thread_id_, s);

    auto& tmp = G_dev_[s];
    // Note: the following resize is safe as it does not deallocate.
    tmp.resizeNoCopy(std::make_pair(gamma_size, n_max_[s]));
    linalg::matrixop::gemm(Gamma_inv_dev_[s], old_M, tmp, thread_id_, s);
    linalg::matrixop::gemm(Real(-1.), old_G, tmp, Real(1.), M_dev_[s], thread_id_, s);
    flop_ += 2 * Gamma_inv_dev_[s].nrRows() * Gamma_inv_dev_[s].nrCols() * old_M.nrCols();
    flop_ += 2 * old_G.nrRows() * old_G.nrCols() * tmp.nrCols();

    details::divideByGammaFactor(MatrixView<linalg::GPU>(M_dev_[s]), gamma_index_dev_[s].ptr(),
                                 gamma_size, get_stream(s));
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

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::computeM(
    std::array<dca::linalg::Matrix<Real, linalg::GPU>, 2>& m_accum) {
  for (int s = 0; s < 2; ++s)
    m_accum[s].resizeNoCopy(M_dev_[s].size());

  for (int s = 0; s < 2; ++s) {
    MatrixView<linalg::GPU> m_in(M_dev_[s]);
    MatrixView<linalg::GPU> m_out(m_accum[s]);
    details::multiplyByInverseFFactor(m_in, m_out, f_dev_[s].ptr(), get_stream(s));
  }

  // TODO: understand why this is necessary.
  config_copied_[0].block();
  config_copied_[1].block();
}

template <class Parameters, typename Real, bool fix_rng_order>
std::size_t CtintWalkerSubmatrixGpu<Parameters, Real, fix_rng_order>::deviceFingerprint() const {
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

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_SUBMATRIX_HPP
