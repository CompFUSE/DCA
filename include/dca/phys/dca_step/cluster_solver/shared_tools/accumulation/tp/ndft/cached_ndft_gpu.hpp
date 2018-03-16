// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements a 2d NDFT while tracking the band indices of each measurement.

#ifndef DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_GPU_HPP
#define DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_GPU_HPP

#ifdef DCA_HAVE_CUDA

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft.hpp"

#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/vector.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/util/magma_vbatched_gemm.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/kernels_interface.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
class CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>
    : private CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density> {
private:
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using ClusterDmn = typename RDmn::parameter_type;
  using BRDmn = func::dmn_variadic<BDmn, RDmn>;
  using Complex = std::complex<Real>;
  using Matrix = linalg::Matrix<Complex, dca::linalg::GPU>;
  using MatrixHost = linalg::Matrix<Complex, dca::linalg::CPU>;

public:
  using BaseClass = CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>;

  CachedNdft(magma_queue_t queue);

  // For each pair of orbitals, performs the non-uniform 2D Fourier Transform from time to frequency
  // defined as M(w1, w2) = \sum_{t1, t2} exp(i (w1 t1 - w2 t2)) M(t1, t2).
  // In: M: input matrix provided by the walker.
  // Out: M_r_r_w_w: stores the result of the computation on the device.
  template <class Configuration, typename ScalarInp>
  void execute(const Configuration& configuration, const linalg::Matrix<ScalarInp, linalg::CPU>& M,
               Matrix& M_r_r_w_w);

  cudaStream_t get_stream() const {
    return stream_;
  }

private:
  template <class Configuration, class OutDmn>
  void executeImpl(const Configuration& configuration, const MatrixHost& M,
                   func::function<std::complex<Real>, OutDmn>& M_r_r_w_w, const int spin);

  using BaseClass::sortConfiguration;
  void sortM(const Matrix& M, Matrix& M_sorted);
  void computeT();
  void performFT(const Matrix& M_t_t, Matrix& M_w_w);
  void rearrangeOutput(const Matrix& M_w_w, Matrix& output);

private:
  using BaseClass::w_;
  linalg::Vector<Real, linalg::GPU> w_dev_;

  using BaseClass::indexed_config_;

  using BaseClass::start_index_;
  using BaseClass::end_index_;

  using BaseClass::n_orbitals_;
  using BaseClass::config_left_;
  using BaseClass::config_right_;
  using BaseClass::start_index_left_;
  using BaseClass::start_index_right_;
  using BaseClass::end_index_left_;
  using BaseClass::end_index_right_;

  magma_queue_t magma_queue_;
  cudaStream_t stream_;
  linalg::util::CudaEvent copy_event_;
  std::array<linalg::Vector<details::Triple<Real>, linalg::GPU>, 2> config_dev_;

  MatrixHost M_;
  Matrix workspace_;
  Matrix work1_;
  Matrix work2_;
  Matrix T_times_M_;

  Matrix T_l_dev_;
  Matrix T_r_dev_;

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
  w_dev_.setAsync(w_, stream_);
}

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
template <class Configuration, typename ScalarInp>
void CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::execute(
    const Configuration& configuration, const linalg::Matrix<ScalarInp, linalg::CPU>& M,
    Matrix& M_out) {
  if (configuration.size() == 0) {  // The result is zero
    M_out.resizeNoCopy(std::make_pair(w_.size() / 2 * n_orbitals_, w_.size() * n_orbitals_));
    M_out.setToZero(stream_);
    return;
  }

  copy_event_.block();
  M_.resizeNoCopy(M.size());
  for (int j = 0; j < M.nrCols(); ++j)
    for (int i = 0; i < M.nrRows(); ++i)
      M_(i, j) = Complex(M(i, j));

  M_out.setAsync(M_, stream_);

  sortConfiguration(configuration);
  // upload configuration.
  config_dev_[0].setAsync(config_left_, stream_);
  config_dev_[1].setAsync(config_right_, stream_);
  assert(cudaPeekAtLastError() == cudaSuccess);

  copy_event_.record(stream_);

  sortM(M_out, work1_);
  computeT();
  performFT(work1_, work2_);
  rearrangeOutput(work2_, M_out);
}

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
void CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::sortM(const Matrix& M,
                                                                                    Matrix& M_sorted) {
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
void CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::performFT(
    const Matrix& M_in, Matrix& M_out) {
  const auto& M_t_t = M_in;
  const int nw = w_.size();
  const int order = indexed_config_[0].size();
  T_times_M_.resizeNoCopy(std::make_pair(nw / 2 * n_orbitals_, order));
  cudaMemsetAsync(T_times_M_.ptr(), 0,
                  T_times_M_.leadingDimension() * T_times_M_.nrCols() * sizeof(Complex), stream_);
  auto& T_times_M_times_T = M_out;
  T_times_M_times_T.resizeNoCopy(std::make_pair(nw / 2 * n_orbitals_, nw * n_orbitals_));
  cudaMemsetAsync(
      T_times_M_times_T.ptr(), 0,
      T_times_M_times_T.leadingDimension() * T_times_M_times_T.nrCols() * sizeof(Complex), stream_);

  {
    // Performs T_l * M_t_t
    const int lda = T_l_dev_.leadingDimension();
    const int ldb = M_t_t.leadingDimension();
    const int ldc = T_times_M_.leadingDimension();

    magma_plan1_.synchronizeCopy();
    for (int i = 0; i < n_orbitals_; ++i) {
      const int n_i = end_index_left_[i] - start_index_left_[i];
      if (!n_i)
        continue;
      magma_plan1_.addGemm(nw / 2, order, n_i, T_l_dev_.ptr(0, start_index_left_[i]), lda,
                           M_t_t.ptr(start_index_left_[i], 0), ldb, T_times_M_.ptr(i * nw / 2, 0),
                           ldc);
    }
    magma_plan1_.execute('N', 'N');
  }
  {
    // Performs T_times_M_ * T_r
    const int lda = T_times_M_.leadingDimension();
    const int ldb = T_r_dev_.leadingDimension();
    const int ldc = T_times_M_times_T.leadingDimension();
    magma_plan2_.synchronizeCopy();
    for (int j = 0; j < n_orbitals_; ++j) {
      const int n_j = end_index_right_[j] - start_index_right_[j];
      if (!n_j)
        continue;
      magma_plan2_.addGemm(T_times_M_.nrRows(), nw, n_j, T_times_M_.ptr(0, start_index_right_[j]),
                           lda, T_r_dev_.ptr(start_index_right_[j], 0), ldb,
                           T_times_M_times_T.ptr(0, nw * j), ldc);
    }
    magma_plan2_.execute('N', 'N');
  }
}

template <typename Real, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
void CachedNdft<Real, RDmn, WDmn, WPosDmn, linalg::GPU, non_density_density>::rearrangeOutput(
    const Matrix& M_w_w, Matrix& M_out) {
  const int nw = w_.size();
  M_out.resizeNoCopy(M_w_w.size());
  const int n_bands = BDmn::dmn_size();
  details::rearrangeOutput(nw, n_orbitals_, n_bands, M_w_w.ptr(), M_w_w.leadingDimension(),
                           M_out.ptr(), M_out.leadingDimension(), stream_);

  assert(cudaPeekAtLastError() == cudaSuccess);
}

}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_GPU_HPP
