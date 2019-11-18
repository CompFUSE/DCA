// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements on the GPU a 2D NDFT from imaginary time to Matsubara frequency, applied
// independently to each pair of orbitals, where an orbital is a combination of cluster site and
// band.

#ifndef DCA_HAVE_CUDA
#pragma error "GPU algorithm requested but DCA_HAVE_CUDA is not defined."
#endif  // DCA_HAVE_CUDA

#ifndef DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_GPU_HPP
#define DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_GPU_HPP

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft_base.hpp"

#include <complex>
#include <memory>

#include "dca/config/mc_options.hpp"
#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/reshapable_matrix.hpp"
#include "dca/linalg/vector.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/util/magma_vbatched_gemm.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft_template.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/kernels_interface.hpp"
#include "dca/util/ignore.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
class CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>
    : public CachedNdftBase<Real, RDmn, WDmn, WPosDmn, non_density_density> {
private:
  using BaseClass = CachedNdftBase<Real, RDmn, WDmn, WPosDmn, non_density_density>;

  using typename BaseClass::BDmn;

  using Complex = std::complex<Real>;
  using Matrix = linalg::Matrix<Complex, dca::linalg::GPU>;
  using RMatrix = linalg::ReshapableMatrix<Complex, dca::linalg::GPU,
                                           config::McOptions::TpAllocator<Complex>>;
  using MatrixHost = linalg::Matrix<Complex, dca::linalg::CPU>;

public:
  CachedNdft(magma_queue_t queue);

  // For each pair of orbitals, performs the non-uniform 2D Fourier Transform from time to frequency
  // defined as M(w1, w2) = \sum_{t1, t2} exp(i (w1 t1 - w2 t2)) M(t1, t2).
  // Out: M_r_r_w_w.
  template <class Configuration, typename RealIn>
  float execute(const Configuration& configuration, const linalg::Matrix<RealIn, linalg::GPU>& M,
                RMatrix& M_r_r_w_w);

  void setWorkspace(const std::shared_ptr<RMatrix>& workspace) {
    workspace_ = workspace;
  }

  cudaStream_t get_stream() const {
    return stream_;
  }

  void synchronizeCopy() {
    copy_event_.block();
  }

  std::size_t deviceFingerprint() const;

private:
  template <typename RealIn>
  void sortM(const linalg::Matrix<RealIn, linalg::GPU>& M, RMatrix& M_sorted) const;
  void computeT();
  double performFT(RMatrix& work);
  void rearrangeOutput(RMatrix& output);

private:
  using BaseClass::config_left_;
  using BaseClass::config_right_;
  using BaseClass::end_index_;
  using BaseClass::end_index_left_;
  using BaseClass::end_index_right_;
  using BaseClass::indexed_config_;
  using BaseClass::n_orbitals_;
  using BaseClass::start_index_;
  using BaseClass::start_index_left_;
  using BaseClass::start_index_right_;
  using BaseClass::w_;

  linalg::Vector<Real, linalg::GPU> w_dev_;
  magma_queue_t magma_queue_;
  cudaStream_t stream_;
  linalg::util::CudaEvent copy_event_;
  std::array<linalg::Vector<details::Triple<Real>, linalg::GPU>, 2> config_dev_;

  std::shared_ptr<RMatrix> workspace_;

  RMatrix T_l_dev_;
  RMatrix T_r_dev_;

  linalg::util::MagmaVBatchedGemm<Complex> magma_plan1_;
  linalg::util::MagmaVBatchedGemm<Complex> magma_plan2_;
};

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::CachedNdft(
    const magma_queue_t queue)
    : BaseClass(),
      magma_queue_(queue),
      stream_(magma_queue_get_cuda_stream(magma_queue_)),
      magma_plan1_(n_orbitals_, magma_queue_),
      magma_plan2_(n_orbitals_, magma_queue_) {
  workspace_ = std::make_shared<RMatrix>();
  w_dev_.setAsync(w_, stream_);
}

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
template <class Configuration, typename RealIn>
float CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::execute(
    const Configuration& configuration, const linalg::Matrix<RealIn, linalg::GPU>& M, RMatrix& M_out) {
  float flop = 0.;

  if (configuration.size() == 0) {  // The result is zero
    M_out.resizeNoCopy(std::make_pair(w_.size() / 2 * n_orbitals_, w_.size() * n_orbitals_));
    M_out.setToZero(stream_);
    return flop;
  }

  BaseClass::sortConfiguration(configuration);
  config_dev_[0].setAsync(config_left_, stream_);
  config_dev_[1].setAsync(config_right_, stream_);
  assert(cudaPeekAtLastError() == cudaSuccess);

  // Allocate enough memory.
  const int nw = w_.size();
  const int order = std::max(config_dev_[0].size(), config_dev_[1].size());
  const int g_first_size = nw / 2 * n_orbitals_;
  M_out.reserveNoCopy(g_first_size * std::max(order, 2 * g_first_size));
  workspace_->reserveNoCopy(std::max(order * order, 2 * g_first_size * g_first_size));

  sortM(M, *workspace_);
  copy_event_.record(stream_);

  computeT();
  flop += performFT(M_out);
  rearrangeOutput(M_out);

  return flop;
}

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
template <typename RealIn>
void CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::sortM(
    const linalg::Matrix<RealIn, linalg::GPU>& M, RMatrix& M_sorted) const {
  M_sorted.resizeNoCopy(M.size());
  details::sortM(M.nrCols(), M.ptr(), M.leadingDimension(), M_sorted.ptr(),
                 M_sorted.leadingDimension(), config_dev_[0].ptr(), config_dev_[1].ptr(), stream_);

  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
void CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::computeT() {
  const int nw = w_.size();
  const int nw_pos = nw / 2;
  const int k = indexed_config_[0].size();
  T_l_dev_.resizeNoCopy(std::make_pair(nw_pos, k));
  details::computeT(nw_pos, k, T_l_dev_.ptr(), T_l_dev_.leadingDimension(), config_dev_[0].ptr(),
                    w_dev_.ptr() + nw_pos, false, stream_);

  T_r_dev_.resizeNoCopy(std::make_pair(k, nw));
  details::computeT(k, nw, T_r_dev_.ptr(), T_r_dev_.leadingDimension(), config_dev_[1].ptr(),
                    w_dev_.ptr(), true, stream_);

  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
double CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::performFT(
    RMatrix& M_out) {
  double flop = 0.;
  const auto& M_t_t = *workspace_;
  const int nw = w_.size();
  const int order = indexed_config_[0].size();

  auto& T_times_M = M_out;
  bool realloc = T_times_M.resizeNoCopy(std::make_pair(nw / 2 * n_orbitals_, order));
  util::ignoreUnused(realloc);
  assert(!realloc);
  T_times_M.setToZero(stream_);

  {
    // Performs T_times_M = T_l * M_t_t.
    const int lda = T_l_dev_.leadingDimension();
    const int ldb = M_t_t.leadingDimension();
    const int ldc = T_times_M.leadingDimension();

    magma_plan1_.synchronizeCopy();
    for (int i = 0; i < n_orbitals_; ++i) {
      const int n_i = end_index_left_[i] - start_index_left_[i];
      if (!n_i)
        continue;
      magma_plan1_.addGemm(nw / 2, order, n_i, T_l_dev_.ptr(0, start_index_left_[i]), lda,
                           M_t_t.ptr(start_index_left_[i], 0), ldb, T_times_M.ptr(i * nw / 2, 0),
                           ldc);

      flop += 8. * nw / 2 * order * n_i;
    }
    magma_plan1_.execute('N', 'N');
  }

  auto& T_times_M_times_T = *workspace_;
  realloc = T_times_M_times_T.resizeNoCopy(std::make_pair(nw / 2 * n_orbitals_, nw * n_orbitals_));
  assert(!realloc);
  T_times_M_times_T.setToZero(stream_);

  {
    // Performs T_times_M_times_T = T_times_M * T_r.
    const int lda = T_times_M.leadingDimension();
    const int ldb = T_r_dev_.leadingDimension();
    const int ldc = T_times_M_times_T.leadingDimension();
    magma_plan2_.synchronizeCopy();
    for (int j = 0; j < n_orbitals_; ++j) {
      const int n_j = end_index_right_[j] - start_index_right_[j];
      if (!n_j)
        continue;
      magma_plan2_.addGemm(T_times_M.nrRows(), nw, n_j, T_times_M.ptr(0, start_index_right_[j]),
                           lda, T_r_dev_.ptr(start_index_right_[j], 0), ldb,
                           T_times_M_times_T.ptr(0, nw * j), ldc);

      flop += 8. * T_times_M.nrRows() * nw * n_j;
    }
    magma_plan2_.execute('N', 'N');
  }

  return flop;
}

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
void CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::rearrangeOutput(
    RMatrix& M_out) {
  const auto& M_w_w = *workspace_;

  const int nw = w_.size();
  M_out.resizeNoCopy(M_w_w.size());
  const int n_bands = BDmn::dmn_size();

  // Rearranges the index order, from fast to slow, from {frequency, band, site} to { site, band,
  // frequency}.
  details::rearrangeOutput(nw, n_orbitals_, n_bands, M_w_w.ptr(), M_w_w.leadingDimension(),
                           M_out.ptr(), M_out.leadingDimension(), stream_);

  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
std::size_t CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::deviceFingerprint()
    const {
  std::size_t res(0);
  if (workspace_.unique())
    res += workspace_->deviceFingerprint();
  res += T_l_dev_.deviceFingerprint();
  res += T_r_dev_.deviceFingerprint();
  return res;
}

}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_GPU_HPP
