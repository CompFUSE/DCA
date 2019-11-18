// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class performs on the GPU the 2D transform from space to momentum used by the tp
// accumulation.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TRANSFORM_2D_GPU
#define DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TRANSFORM_2D_GPU

#ifndef DCA_HAVE_CUDA
#pragma error "This file requires CUDA support."
#endif

#include "dca/math/function_transform/special_transforms/space_transform_2D_gpu.hpp"

#include <array>
#include <memory>

#include "dca/config/mc_options.hpp"
#include "dca/linalg/reshapable_matrix.hpp"
#include "dca/linalg/util/magma_batched_gemm.hpp"
#include "dca/math/function_transform/special_transforms/kernels_interface.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

template <class RDmn, class KDmn, typename Real = double>
class SpaceTransform2DGpu : private SpaceTransform2D<RDmn, KDmn, Real> {
private:
  using BaseClass = SpaceTransform2D<RDmn, KDmn, Real>;

  using Complex = std::complex<Real>;
  using MatrixDev = linalg::Matrix<Complex, linalg::GPU>;
  using VectorDev = linalg::Vector<Complex, linalg::GPU>;
  using RMatrix =
      linalg::ReshapableMatrix<Complex, linalg::GPU, config::McOptions::TpAllocator<Complex>>;

public:
  // Constructor
  // In: nw_pos: number of extended positive frequencies.
  // In: queue: the magma queue on which execute will run.
  SpaceTransform2DGpu(int nw_pos, magma_queue_t queue);

  // Performs the 2D fourier transform from real to momentum space in place and rearranges the
  // order of M's labels from (r, b, w) to (b, r, w).
  // See space_transform_2D.hpp for a definition of the transform.
  // In/Out: M
  // Returns: number of flop.
  float execute(RMatrix& M);

  void setWorkspace(const std::shared_ptr<RMatrix>& workspace) {
    workspace_ = workspace;
  }

  // Returns: the stream associated with the magma queue.
  cudaStream_t get_stream() const {
    return stream_;
  }

  std::size_t deviceFingerprint() const {
    std::size_t res(0);

    if (workspace_.unique())
      res += workspace_->deviceFingerprint();
    return res;
  }

private:
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;

  static const MatrixDev& get_T_matrix();
  const auto& getPhaseFactors();

  void phaseFactorsAndRearrange(const RMatrix& in, RMatrix& out);

  const int n_bands_;
  const int nw_;
  const int nc_;

  magma_queue_t queue_;
  cudaStream_t stream_;

  std::shared_ptr<RMatrix> workspace_;

  linalg::util::MagmaBatchedGemm<Complex> plan1_;
  linalg::util::MagmaBatchedGemm<Complex> plan2_;
};

template <class RDmn, class KDmn, typename Real>
SpaceTransform2DGpu<RDmn, KDmn, Real>::SpaceTransform2DGpu(const int nw_pos, magma_queue_t queue)
    : n_bands_(BDmn::dmn_size()),
      nw_(2 * nw_pos),
      nc_(RDmn::dmn_size()),
      queue_(queue),
      stream_(magma_queue_get_cuda_stream(queue_)),
      plan1_(queue_),
      plan2_(queue_) {
  workspace_ = std::make_shared<RMatrix>();
}

template <class RDmn, class KDmn, typename Real>
float SpaceTransform2DGpu<RDmn, KDmn, Real>::execute(RMatrix& M) {
  float flop = 0.;

  auto& T_times_M = *(workspace_);
  auto& T_times_M_times_T = M;

  T_times_M.resizeNoCopy(M.size());
  const auto& T = get_T_matrix();

  {
    const int n_trafo = M.nrRows() / nc_;

    plan1_.synchronizeCopy();
    for (int i = 0; i < n_trafo; ++i)
      plan1_.addGemm(T.ptr(), M.ptr(i * nc_, 0), T_times_M.ptr(i * nc_, 0));

    const int lda = T.leadingDimension();
    const int ldb = M.leadingDimension();
    const int ldc = T_times_M.leadingDimension();
    plan1_.execute('N', 'N', nc_, M.nrCols(), nc_, Complex(1), Complex(0), lda, ldb, ldc);
    flop += n_trafo * 8. * nc_ * M.nrCols() * nc_;
  }

  {
    const int n_trafo = M.nrCols() / nc_;
    plan2_.synchronizeCopy();
    for (int j = 0; j < n_trafo; ++j)
      plan2_.addGemm(T_times_M.ptr(0, j * nc_), T.ptr(), T_times_M_times_T.ptr(0, j * nc_));

    const int lda = T_times_M.leadingDimension();
    const int ldb = T.leadingDimension();
    const int ldc = T_times_M_times_T.leadingDimension();
    const Complex norm(1. / nc_);
    plan2_.execute('N', 'C', M.nrRows(), nc_, nc_, norm, Complex(0), lda, ldb, ldc);
    flop += n_trafo * 8. * M.nrRows() * nc_ * nc_;
  }

  phaseFactorsAndRearrange(T_times_M_times_T, *workspace_);
  M.swap(*workspace_);

  return flop;
}

template <class RDmn, class KDmn, typename Real>
void SpaceTransform2DGpu<RDmn, KDmn, Real>::phaseFactorsAndRearrange(const RMatrix& in, RMatrix& out) {
  out.resizeNoCopy(in.size());
  const Complex* const phase_factors_ptr =
      BaseClass::hasPhaseFactors() ? getPhaseFactors().ptr() : nullptr;
  details::phaseFactorsAndRearrange(in.ptr(), in.leadingDimension(), out.ptr(),
                                    out.leadingDimension(), n_bands_, nc_, nw_, phase_factors_ptr,
                                    stream_);
}

template <class RDmn, class KDmn, typename Real>
const linalg::Matrix<std::complex<Real>, linalg::GPU>& SpaceTransform2DGpu<RDmn, KDmn,
                                                                           Real>::get_T_matrix() {
  auto initialize_T_matrix = []() {
    const auto T_host = BaseClass::get_T_matrix();
    return MatrixDev(T_host);
  };

  static const MatrixDev T(initialize_T_matrix(), "T_space_2D_device");
  return T;
}

template <class RDmn, class KDmn, typename Real>
const auto& SpaceTransform2DGpu<RDmn, KDmn, Real>::getPhaseFactors() {
  auto initialize = []() {
    const auto& phase_factors = BaseClass::getPhaseFactors();
    linalg::Vector<std::complex<Real>, linalg::CPU> host_vector(phase_factors.size());
    std::copy_n(phase_factors.values(), phase_factors.size(), host_vector.ptr());
    return VectorDev(host_vector);
  };

  static const VectorDev phase_factors_dev(initialize());

  return phase_factors_dev;
}

}  // namespace transform
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TRANSFORM_2D_GPU
