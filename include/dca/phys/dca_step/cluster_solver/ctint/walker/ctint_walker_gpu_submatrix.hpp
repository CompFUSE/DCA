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

template <class Parameters>
class CtintWalkerSubmatrix<linalg::GPU, Parameters>
    : public CtintWalkerSubmatrix<linalg::CPU, Parameters> {
public:
  using this_type = CtintWalkerSubmatrix<linalg::GPU, Parameters>;
  using BaseClass = CtintWalkerSubmatrix<linalg::CPU, Parameters>;

  using typename BaseClass::Profiler;
  using typename BaseClass::Rng;

  CtintWalkerSubmatrix(Parameters& pars_ref, Rng& rng_ref, const InteractionVertices& vertices,
                       const DMatrixBuilder<linalg::GPU>& builder_ref, int id = 0);

public:
  void doSweep() override;

  void synchronize() const;

  using BaseClass::order;

private:
  void doStep() override;

  void computeMInit();
  void computeGInit();
  void updateM();
  void pushToEnd();

  void uploadConfiguration();

protected:
  // For testing purposes:
  void doStep(int n_moves_to_delay);

protected:
  using BaseClass::configuration_;
  using BaseClass::M_;

private:
  template <linalg::DeviceType dev = linalg::GPU>
  using MatrixView = linalg::MatrixView<double, dev>;
  using Matrix = linalg::Matrix<double, linalg::CPU>;

  DeviceConfigurationManager device_config_;

  using BaseClass::removal_list_;
  using BaseClass::conf_removal_list_;

  const DMatrixBuilder<linalg::GPU>& d_builder_;

  using BaseClass::parameters_;
  using BaseClass::rng_;

private:
  template <linalg::DeviceType device_t>
  using MatrixPair = std::array<linalg::Matrix<double, device_t>, 2>;

  using BaseClass::G_;
  MatrixPair<linalg::GPU> M_dev_;
  MatrixPair<linalg::GPU> Gamma_inv_dev_;
  MatrixPair<linalg::GPU> D_dev_;
  MatrixPair<linalg::GPU> G_dev_;
  MatrixPair<linalg::GPU> G0_dev_;

  using BaseClass::Gamma_inv_;
  using BaseClass::f_;

  std::array<linalg::Vector<double, linalg::GPU>, 2> f_dev_;
  std::array<linalg::Vector<double, linalg::CPU>, 2> f_values_;

  using BaseClass::gamma_;
  using BaseClass::move_indices_;
  std::array<linalg::util::HostVector<std::pair<int, double>>, 2> gamma_index_;
  std::array<linalg::Vector<std::pair<int, double>, linalg::GPU>, 2> gamma_index_dev_;

  std::array<linalg::util::CudaEvent, 2> config_copied_;

  using BaseClass::concurrency_;
  using BaseClass::thread_id_;
  std::array<cudaStream_t, 2> stream_;

  // Initial and current sector sizes.
  using BaseClass::n_init_;

  // Maximal sector size after submatrix update.
  using BaseClass::n_max_;

  // If true the m matrix is multiplied by a f-factor.
  bool m_transformed_ = false;
};

template <class Parameters>
CtintWalkerSubmatrix<linalg::GPU, Parameters>::CtintWalkerSubmatrix(
    Parameters& parameters_ref, Rng& rng_ref, const InteractionVertices& vertices,
    const DMatrixBuilder<linalg::GPU>& builder_ref, int id)
    : BaseClass(parameters_ref, rng_ref, vertices, builder_ref, id), d_builder_(builder_ref) {
  for (int s = 0; s < 2; ++s) {
    stream_[s] = linalg::util::getStream(thread_id_, s);
    M_dev_[s].setAsync(M_[s], stream_[s]);
  }
  uploadConfiguration();

  if (concurrency_.id() == concurrency_.first() && thread_id_ == 0)
    std::cout << "\nCT-INT submatrix walker extended to GPU." << std::endl;
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::synchronize() const {
  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  cudaStreamSynchronize(stream_[0]);
  cudaStreamSynchronize(stream_[1]);
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::doSweep() {
  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    MatrixView<linalg::GPU> M(M_dev_[s]);
    details::multiplyByFFactor(M, f_dev_[s].ptr(), true, true, stream_[s]);
  }

  BaseClass::doSteps();
  uploadConfiguration();

  for (int s = 0; s < 2; ++s) {
    MatrixView<linalg::GPU> M(M_dev_[s]);
    details::multiplyByFFactor(M, f_dev_[s].ptr(), false, true, stream_[s]);

    if (BaseClass::thermalized_)  // TODO: do not download if accumulator accepts GPU matrices.
      M_[s].setAsync(M_dev_[s], stream_[s]);
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::doStep(const int n_moves_to_delay) {
  for (int s = 0; s < 2; ++s) {
    MatrixView<linalg::GPU> M(M_dev_[s]);
    details::multiplyByFFactor(M, f_dev_[s].ptr(), true, true, stream_[s]);
  }

  BaseClass::nbr_of_moves_to_delay_ = n_moves_to_delay;
  doStep();
  uploadConfiguration();

  for (int s = 0; s < 2; ++s) {
    MatrixView<linalg::GPU> M(M_dev_[s]);
    details::multiplyByFFactor(M, f_dev_[s].ptr(), false, true, stream_[s]);
    M_[s].setAsync(M_dev_[s], stream_[s]);
  }
  synchronize();
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::doStep() {
  BaseClass::generateDelayedMoves(BaseClass::nbr_of_moves_to_delay_);
  uploadConfiguration();

  computeMInit();
  computeGInit();
  synchronize();
  BaseClass::mainSubmatrixProcess();
  updateM();
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::uploadConfiguration() {
  for (int s = 0; s < 2; ++s)
    config_copied_[s].block();

  // Upload configuration and f values.
  device_config_.upload(configuration_, thread_id_);
  for (int s = 0; s < 2; ++s) {
    auto& values = f_values_[s];
    const auto& sector = configuration_.getSector(s);
    values.resizeNoCopy(sector.size());
    for (int i = 0; i < sector.size(); ++i)
      values[i] = f_[sector.getAuxFieldType(i)];
    f_dev_[s].setAsync(values, stream_[s]);
  }

  for (int s = 0; s < 2; ++s)
    config_copied_[s].record(stream_[s]);
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::computeMInit() {
  //  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s)
    M_dev_[s].resize(n_max_[s]);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];
    if (delta > 0) {
      D_dev_[s].resizeNoCopy(std::make_pair(delta, n_init_[s]));
      d_builder_.computeG0(D_dev_[s], device_config_.getDeviceData(s), n_init_[s], false, stream_[s]);

      MatrixView<linalg::GPU> D_view(D_dev_[s]);
      details::multiplyByFFactor(D_view, f_dev_[s].ptr(), false, false, stream_[s]);

      MatrixView<linalg::GPU> M(M_dev_[s], 0, 0, n_init_[s], n_init_[s]);
      MatrixView<linalg::GPU> D_M(M_dev_[s], n_init_[s], 0, delta, n_init_[s]);

      linalg::matrixop::gemm(D_dev_[s], M, D_M, thread_id_, s);

      details::setRightSectorToId(M_dev_[s].ptr(), M_dev_[s].leadingDimension(), n_init_[s],
                                  n_max_[s], stream_[s]);
    }
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::computeGInit() {
  //  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];
    auto& f_dev = f_dev_[s];

    G_dev_[s].resizeNoCopy(n_max_[s]);

    MatrixView<> G(G_dev_[s]);
    const MatrixView<> M(M_dev_[s]);
    details::computeGLeft(G, M, f_dev.ptr(), n_init_[s], stream_[s]);

    if (delta > 0) {
      G0_dev_[s].resizeNoCopy(std::make_pair(n_max_[s], delta));
      d_builder_.computeG0(G0_dev_[s], device_config_.getDeviceData(s), n_init_[s], true, stream_[s]);

      MatrixView<linalg::GPU> G(G_dev_[s], 0, n_init_[s], n_max_[s], delta);
      // compute G right.
      linalg::matrixop::gemm(M_dev_[s], G0_dev_[s], G, thread_id_, s);
    }
    G_[s].setAsync(G_dev_[s], stream_[s]);
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::updateM() {
  //  Profiler profiler(__FUNCTION__, "CT-INT GPU walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s)
    Gamma_inv_dev_[s].setAsync(Gamma_inv_[s], stream_[s]);

  // Copy gamma factors.
  for (int s = 0; s < 2; ++s) {
    auto& gamma_index = gamma_index_[s];
    const int n = gamma_[s].size();
    gamma_index.resize(n);
    for (int i = 0; i < n; ++i)
      gamma_index[i] = std::make_pair(move_indices_[s][i], gamma_[s][i]);
    gamma_index_dev_[s].setAsync(gamma_index, stream_[s]);
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

    // TODO: this can be done in a single kernel with the copied memory.
    int p;
    for (int j = 0; j < gamma_size; ++j) {
      p = move_indices_[s][j];

      linalg::matrixop::copyCol(G_dev_[s], p, old_G, j, thread_id_, s);
      linalg::matrixop::copyRow(M_dev_[s], p, old_M, j, thread_id_, s);
    }

    auto& tmp = G_dev_[s];
    // Note: the following resize is safe as it does not deallocate.
    tmp.resizeNoCopy(std::make_pair(gamma_size, n_max_[s]));
    linalg::matrixop::gemm(Gamma_inv_dev_[s], old_M, tmp, thread_id_, s);
    linalg::matrixop::gemm(-1., old_G, tmp, 1., M_dev_[s], thread_id_, s);

    details::divideByGammaFactor(MatrixView<linalg::GPU>(M_dev_[s]), gamma_index_dev_[s].ptr(),
                                 gamma_size, stream_[s]);
  }

  // Remove "non-interacting" rows and columns.
  pushToEnd();

  for (int s = 0; s < 2; ++s)
    M_dev_[s].resize(n_max_[s] - BaseClass::removal_list_[s].size());

  const int n_new_vertices = (M_dev_[0].size().first + M_dev_[1].size().first) / 2;
  while (configuration_.size() > n_new_vertices)
    configuration_.pop();

  assert(configuration_.getSector(0).size() == M_dev_[0].nrRows());
  assert(configuration_.getSector(1).size() == M_dev_[1].nrRows());
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::pushToEnd() {
  int source;
  int destination;

  for (int s = 0; s < 2; ++s) {
    // Sort in reverse order.
    std::sort(removal_list_[s].rbegin(), removal_list_[s].rend());

    destination = n_max_[s] - 1;

    for (int i = 0; i < removal_list_[s].size(); ++i) {
      source = removal_list_[s][i];

      removal_list_[s][i] = destination;

      linalg::matrixop::swapRowAndCol(M_dev_[s], source, destination);
      configuration_.swapSectorLabels(source, destination, s);

      --destination;
    }
  }

  std::sort(conf_removal_list_.rbegin(), conf_removal_list_.rend());

  destination = configuration_.size() - 1;

  for (int i = 0; i < conf_removal_list_.size(); ++i) {
    source = conf_removal_list_[i];

    conf_removal_list_[i] = destination;

    configuration_.swapVertices(source, destination);

    --destination;
  }
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_SUBMATRIX_HPP
