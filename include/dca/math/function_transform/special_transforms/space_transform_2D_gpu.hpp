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
  using Complex = std::complex<Real>;
  using MatrixDev = linalg::Matrix<Complex, linalg::GPU>;

public:
  // Constructor
  // In: nw_pos: number of extended positive frequencies.
  // In: queue: the magma queue on which execute will run.
  SpaceTransform2DGpu(int nw_pos, magma_queue_t queue);

  // Performs the 2D fourier transform from real to momentum space in place and rearranges the
  // order of M's labels from (r, b, w) to (b, r, w).
  // The transform is equivalent to M(k1, k2) = \sum_{r1, r2} exp(i(k1 * r1 - k2 * r2)) M(r1, r2)
  // In/Out: M
  void execute(MatrixDev& M);

  // Returns: the stream associated with the magma queue.
  cudaStream_t get_stream() const {
    return stream_;
  }

  std::size_t deviceFingerprint() const {
    return T_times_M_.deviceFingerprint() + T_times_M_times_T_.deviceFingerprint();
  }

private:
  using BaseClass = SpaceTransform2D<RDmn, KDmn, Real>;
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;

  static const MatrixDev& get_T_matrix();
  void rearrangeResult(const MatrixDev& in, MatrixDev& out);

  const int n_bands_;
  const int nw_;
  const int nc_;

  magma_queue_t queue_;
  cudaStream_t stream_;

  MatrixDev T_times_M_;
  MatrixDev T_times_M_times_T_;

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
      plan2_(queue_) {}

template <class RDmn, class KDmn, typename Real>
void SpaceTransform2DGpu<RDmn, KDmn, Real>::execute(MatrixDev& M) {
  T_times_M_.resizeNoCopy(M.size());
  const auto& T = get_T_matrix();

  {
    const int n_trafo = M.nrRows() / nc_;

    plan1_.synchronizeCopy();
    for (int i = 0; i < n_trafo; ++i)
      plan1_.addGemm(T.ptr(), M.ptr(i * nc_, 0), T_times_M_.ptr(i * nc_, 0));

    const int lda = T.leadingDimension();
    const int ldb = M.leadingDimension();
    const int ldc = T_times_M_.leadingDimension();
    plan1_.execute('N', 'N', nc_, M.nrCols(), nc_, Complex(1), Complex(0), lda, ldb, ldc);
  }

  T_times_M_times_T_.resizeNoCopy(M.size());
  {
    const int n_trafo = M.nrCols() / nc_;
    plan2_.synchronizeCopy();
    for (int j = 0; j < n_trafo; ++j)
      plan2_.addGemm(T_times_M_.ptr(0, j * nc_), T.ptr(), T_times_M_times_T_.ptr(0, j * nc_));

    const int lda = T_times_M_.leadingDimension();
    const int ldb = T.leadingDimension();
    const int ldc = T_times_M_times_T_.leadingDimension();
    const Complex norm(1. / nc_);
    plan2_.execute('N', 'C', M.nrRows(), nc_, nc_, norm, Complex(0), lda, ldb, ldc);
  }

  rearrangeResult(T_times_M_times_T_, M);
}

template <class RDmn, class KDmn, typename Real>
void SpaceTransform2DGpu<RDmn, KDmn, Real>::rearrangeResult(const MatrixDev& in, MatrixDev& out) {
  out.resizeNoCopy(in.size());
  details::rearrangeResult(in.ptr(), in.leadingDimension(), out.ptr(), out.leadingDimension(),
                           n_bands_, nc_, nw_, stream_);
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

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TRANSFORM_2D_GPU
